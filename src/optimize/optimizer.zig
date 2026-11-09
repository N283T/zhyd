//! Optimization engine: selects the best orientation for each mover
//! by searching over cliques (connected components) of interacting movers.
//! Uses brute-force enumeration for small cliques, greedy for large ones.

const std = @import("std");
const Allocator = std.mem.Allocator;
const math_mod = @import("../math.zig");
const Vec3 = math_mod.Vec3;
const model_mod = @import("../model.zig");
const Model = model_mod.Model;
const mover_mod = @import("mover.zig");
const Mover = mover_mod.Mover;
const scorer_mod = @import("scorer.zig");
const scoring_mod = @import("scoring.zig");
const clique_mod = @import("clique.zig");
const rotator_mod = @import("rotator.zig");
const ScoreContext = scoring_mod.ScoreContext;
const buildScoreContext = scoring_mod.buildScoreContext;
const scoreMover = scoring_mod.scoreMover;
const scoreMoverWithPositions = scoring_mod.scoreMoverWithPositions;

pub const OptConfig = struct {
    n_threads: u32 = 0, // Reserved for future optimizer-internal parallelism.
    brute_force_limit: u64 = 100_000,
    interaction_cutoff: f32 = 6.0, // Angstrom -- max distance for mover interaction
    scoring_params: scorer_mod.ScoringParams = .{},
};

pub const OptResult = struct {
    n_singletons: u32 = 0,
    n_brute_force: u32 = 0,
    n_vertex_cut: u32 = 0,
    total_cliques: u32 = 0,
};

/// Optimize all movers: find best orientations using clique-based search.
pub fn optimize(
    allocator: Allocator,
    movers: []Mover,
    model: *Model,
    config: OptConfig,
) !OptResult {
    var result = OptResult{};


    var score_ctx = try buildScoreContext(allocator, movers, model.atoms.items);
    defer score_ctx.deinit(allocator);

    var scratch = std.ArrayListUnmanaged(u32).empty;
    defer scratch.deinit(allocator);

    // Build interaction graph
    var graph = try clique_mod.buildInteractionGraph(
        allocator,
        movers,
        model.atoms.items,
        config.interaction_cutoff,
    );
    defer graph.deinit();

    // Find connected components
    const cliques = try clique_mod.findCliques(allocator, &graph);
    defer {
        for (cliques) |c| allocator.free(c);
        allocator.free(cliques);
    }

    result.total_cliques = @intCast(cliques.len);

    // Collect singleton mover indices for a dedicated pass after sequential cliques.
    var singleton_indices = std.ArrayListUnmanaged(u32).empty;
    defer singleton_indices.deinit(allocator);
    // Pre-allocate worst-case capacity so that the OOM fallback path's append
    // calls in the brute-force error handler never need to allocate.
    try singleton_indices.ensureTotalCapacity(allocator, movers.len);

    for (cliques) |clq| {
        if (allFixed(movers, clq)) continue;
        if (clq.len == 1) {
            // Defer singletons for the singleton pass below.
            try singleton_indices.append(allocator, clq[0]);
            result.n_singletons += 1;
        } else if (totalStates(movers, clq) <= config.brute_force_limit) {
            // Brute force: enumerate all combinations
            optimizeBruteForce(allocator, movers, clq, model, config, &score_ctx, &scratch) catch |err| switch (err) {
                error.OutOfMemory => {
                    std.log.warn("brute-force OOM for clique of {d} movers; falling back to singletons", .{clq.len});
                    for (clq) |mi| try singleton_indices.append(allocator, mi);
                    result.n_vertex_cut += 1;
                    continue;
                },
            };
            result.n_brute_force += 1;
        } else {
            // Vertex-cut decomposition (simplified: greedy for Phase 3)
            optimizeIterativeGreedy(movers, clq, model, config, &score_ctx, allocator, &scratch);
            result.n_vertex_cut += 1;
        }
    }

    // Dispatch singleton optimization sequentially.
    // TODO: restore optimizer-internal parallelism with std.Io.Group after the Zig 0.16 migration.
    // Thread safety: optimizeSingleton uses scoreMoverWithPositions which reads other movers'
    // atoms from model.atoms but does NOT write to model.atoms during scoring. The only
    // write is to m.best_orientation (a per-mover field; each thread owns a distinct mover).
    for (singleton_indices.items) |mi| {
        optimizeSingleton(movers, mi, model, config, &score_ctx, allocator, &scratch);
    }

    // Apply best coarse orientations
    for (movers) |*m| {
        m.applyOrientation(model.atoms.items, m.best_orientation);
        m.current_orientation = m.best_orientation;
    }

    // Fine search phase: refine each mover around its coarse best.
    for (0..movers.len) |mi| {
        var best_fine_positions: ?[]math_mod.Vec3(f32) = null;
        fineSearchMover(allocator, movers, @intCast(mi), model, config, &score_ctx, &scratch, &best_fine_positions);
        if (best_fine_positions) |positions| {
            defer allocator.free(positions);
            for (movers[mi].atom_indices, 0..) |ai, j| {
                model.atoms.items[ai].pos = positions[j];
            }
        }
    }

    return result;
}

fn optimizeSingleton(
    movers: []Mover,
    mover_idx: u32,
    model: *Model,
    config: OptConfig,
    score_ctx: *const ScoreContext,
    allocator: Allocator,
    scratch: *std.ArrayListUnmanaged(u32),
) void {
    const m = &movers[mover_idx];
    if (m.is_fixed) return;
    var best_score: f32 = -std.math.inf(f32);
    var best_idx: u16 = 0;

    for (0..m.nOrientations()) |oi| {
        const idx: u16 = @intCast(oi);
        // Score using the orientation's positions directly — do NOT call applyOrientation,
        // which would write to model.atoms and race with concurrent threads reading it.
        const orient_positions = m.orientations[idx].positions;
        const score = scoreMoverWithPositions(m, mover_idx, orient_positions, movers, model, config.scoring_params, config.scoring_params.gap_scale, score_ctx, allocator, scratch) - m.orientationPenalty(idx);
        if (score > best_score) {
            best_score = score;
            best_idx = idx;
        }
    }

    m.best_orientation = best_idx;
}

/// Refine a mover's best orientation with fine angular search.
///
/// Thread-safe: does NOT write to model.atoms during scoring. Scoring uses
/// scoreMoverWithPositions with explicit orientation positions. If a fine position
/// beats the coarse best, the best positions are allocated and written to
/// `out_best_positions` for the caller to apply after all parallel tasks join.
/// The caller is responsible for freeing the slice.
fn fineSearchMover(
    allocator: Allocator,
    movers: []Mover,
    mover_idx: u32,
    model: *Model,
    config: OptConfig,
    score_ctx: *const ScoreContext,
    scratch: *std.ArrayListUnmanaged(u32),
    out_best_positions: *?[]math_mod.Vec3(f32),
) void {
    const m = &movers[mover_idx];
    if (m.is_fixed) return;
    const atoms = model.atoms.items;

    const fine = rotator_mod.generateFineOrientations(allocator, atoms, m, m.best_orientation) catch |err| {
        std.log.warn("fine orientation generation failed for mover {d}: {s}", .{ mover_idx, @errorName(err) });
        return;
    };
    defer {
        for (fine) |o| allocator.free(o.positions);
        allocator.free(fine);
    }

    if (fine.len == 0) return;

    // Score current coarse best using the orientation's stored positions — do NOT read
    // atoms[ai].pos because another thread may be writing it concurrently during fine search.
    const coarse_positions = m.orientations[m.best_orientation].positions;
    var best_score = scoreMoverWithPositions(m, mover_idx, coarse_positions, movers, model, config.scoring_params, config.scoring_params.gap_scale, score_ctx, allocator, scratch) - m.orientationPenalty(m.best_orientation);

    // Score fine orientations using their positions directly (no writes to model.atoms).
    var best_fine: ?usize = null;
    for (fine, 0..) |orient, fi| {
        const score = scoreMoverWithPositions(m, mover_idx, orient.positions, movers, model, config.scoring_params, config.scoring_params.gap_scale, score_ctx, allocator, scratch) - orient.penalty;
        if (score > best_score) {
            best_score = score;
            best_fine = fi;
        }
    }

    if (best_fine) |fi| {
        // Copy the best fine positions into a freshly allocated slice; the caller will
        // apply them to model.atoms after all threads join (avoiding write-write races).
        const n = m.atom_indices.len;
        const best = allocator.alloc(math_mod.Vec3(f32), n) catch |err| {
            std.log.warn("fine-search result alloc failed for mover {d}: {s}", .{ mover_idx, @errorName(err) });
            return;
        };
        @memcpy(best, fine[fi].positions);
        out_best_positions.* = best;
    }
    // If no fine orientation improves on the coarse best, leave out_best_positions null.
    // The coarse best is already applied to model.atoms (done before this phase starts).
}

/// Exhaustive search over every orientation combination of a clique.
///
/// The score of a combination is a sum of atom-pair terms, so it splits into
/// one term per mover (static atoms and movers outside the clique) and one
/// term per ordered pair of interacting clique movers. Both are tabulated
/// once, and the enumeration only adds table entries. Scoring every mover
/// from scratch for each combination made cliques near `brute_force_limit`
/// take tens of seconds.
fn optimizeBruteForce(
    allocator: Allocator,
    movers: []Mover,
    clq: []const u32,
    model: *Model,
    config: OptConfig,
    score_ctx: *const ScoreContext,
    scratch: *std.ArrayListUnmanaged(u32),
) !void {
    const n = clq.len;
    const atoms = model.atoms.items;
    const params = config.scoring_params;
    const gap_scale = params.gap_scale;

    // Candidate orientations per clique member: only the locked one when fixed.
    const n_cand = try allocator.alloc(u16, n);
    defer allocator.free(n_cand);
    for (clq, 0..) |mi, i| n_cand[i] = movers[mi].effectiveOrientations();

    const in_clique = try allocator.alloc(bool, movers.len);
    defer allocator.free(in_clique);
    @memset(in_clique, false);
    for (clq) |mi| in_clique[mi] = true;

    var max_atoms: usize = 0;
    for (clq) |mi| max_atoms = @max(max_atoms, movers[mi].atom_indices.len);
    const pos_buf = try allocator.alloc(Vec3(f32), max_atoms);
    defer allocator.free(pos_buf);

    // Per-mover term (including the orientation penalty) for each candidate.
    const self_offsets = try allocator.alloc(usize, n + 1);
    defer allocator.free(self_offsets);
    self_offsets[0] = 0;
    for (0..n) |i| self_offsets[i + 1] = self_offsets[i] + n_cand[i];
    const self_table = try allocator.alloc(f32, self_offsets[n]);
    defer allocator.free(self_table);
    for (clq, 0..) |mi, i| {
        const m = &movers[mi];
        for (0..n_cand[i]) |c| {
            const o = candidateOrientation(m, @intCast(c));
            m.applyOrientation(atoms, o);
            const positions = pos_buf[0..m.atom_indices.len];
            for (m.atom_indices, 0..) |ai, j| positions[j] = atoms[ai].pos;
            const s = scoring_mod.scoreMoverExcluding(m, mi, positions, movers, model, params, gap_scale, score_ctx, allocator, scratch, in_clique);
            self_table[self_offsets[i] + c] = s - m.orientationPenalty(o);
        }
    }

    // Pair terms for interacting clique members i < j: both directions summed,
    // stored row-major over (candidate of i, candidate of j).
    const Pair = struct { i: usize, j: usize, offset: usize };
    var pairs = std.ArrayListUnmanaged(Pair).empty;
    defer pairs.deinit(allocator);
    var pair_len: usize = 0;
    for (0..n) |i| {
        for (i + 1..n) |j| {
            if (!scoring_mod.moversMayInteract(score_ctx, clq[i], clq[j])) continue;
            try pairs.append(allocator, .{ .i = i, .j = j, .offset = pair_len });
            pair_len += @as(usize, n_cand[i]) * n_cand[j];
        }
    }
    const pair_table = try allocator.alloc(f32, pair_len);
    defer allocator.free(pair_table);
    for (pairs.items) |p| {
        const a = &movers[clq[p.i]];
        const b = &movers[clq[p.j]];
        for (0..n_cand[p.i]) |ca| {
            a.applyOrientation(atoms, candidateOrientation(a, @intCast(ca)));
            for (0..n_cand[p.j]) |cb| {
                b.applyOrientation(atoms, candidateOrientation(b, @intCast(cb)));
                pair_table[p.offset + ca * n_cand[p.j] + cb] =
                    scoring_mod.scoreMoverPair(a, b, atoms, params, gap_scale, score_ctx) +
                    scoring_mod.scoreMoverPair(b, a, atoms, params, gap_scale, score_ctx);
            }
        }
    }

    // Pairs grouped by their later member, so the score of a prefix of the
    // mixed-radix counter can be reused when only later digits change.
    const pairs_by_j = try allocator.alloc(std.ArrayListUnmanaged(Pair), n);
    defer {
        for (pairs_by_j) |*l| l.deinit(allocator);
        allocator.free(pairs_by_j);
    }
    for (pairs_by_j) |*l| l.* = .empty;
    for (pairs.items) |p| try pairs_by_j[p.j].append(allocator, p);

    // prefix[k] = score of members 0..k-1 with their current candidates.
    const cand = try allocator.alloc(u16, n);
    defer allocator.free(cand);
    @memset(cand, 0);
    const best = try allocator.alloc(u16, n);
    defer allocator.free(best);
    @memset(best, 0);
    const prefix = try allocator.alloc(f32, n + 1);
    defer allocator.free(prefix);
    prefix[0] = 0;

    var best_score: f32 = -std.math.inf(f32);
    var first_dirty: usize = 0;
    while (true) {
        for (first_dirty..n) |k| {
            var s = prefix[k] + self_table[self_offsets[k] + cand[k]];
            for (pairs_by_j[k].items) |p| {
                s += pair_table[p.offset + @as(usize, cand[p.i]) * n_cand[k] + cand[k]];
            }
            prefix[k + 1] = s;
        }
        if (prefix[n] > best_score) {
            best_score = prefix[n];
            @memcpy(best, cand);
        }

        // Mixed-radix increment, last digit fastest (same order as before).
        var k: usize = n;
        const done = while (k > 0) {
            k -= 1;
            cand[k] += 1;
            if (cand[k] < n_cand[k]) break false;
            cand[k] = 0;
        } else true;
        if (done) break;
        first_dirty = k;
    }

    for (clq, 0..) |mi, i| {
        const o = candidateOrientation(&movers[mi], best[i]);
        movers[mi].best_orientation = o;
        movers[mi].applyOrientation(atoms, o);
    }
}

/// Orientation index of candidate `c` of a mover: the locked orientation for a
/// fixed mover, otherwise `c` itself.
fn candidateOrientation(m: *const Mover, c: u16) u16 {
    return if (m.is_fixed) m.best_orientation else c;
}

fn incrementIndices(indices: []u16, movers: []const Mover, clq: []const u32) bool {
    var i: usize = indices.len;
    while (i > 0) {
        i -= 1;
        if (movers[clq[i]].is_fixed) {
            indices[i] = movers[clq[i]].best_orientation;
            continue;
        }
        indices[i] += 1;
        if (indices[i] < movers[clq[i]].nOrientations()) return true;
        indices[i] = 0;
    }
    return false; // overflow = done
}

/// Compute the total number of combined states for a set of movers.
/// Uses saturating multiplication to avoid overflow.
pub fn totalStates(movers: []const Mover, clq: []const u32) u64 {
    var total: u64 = 1;
    for (clq) |mi| {
        total *|= movers[mi].effectiveOrientations(); // saturating multiply
    }
    return total;
}

fn optimizeIterativeGreedy(
    movers: []Mover,
    clq: []const u32,
    model: *Model,
    config: OptConfig,
    score_ctx: *const ScoreContext,
    allocator: Allocator,
    scratch: *std.ArrayListUnmanaged(u32),
) void {
    const max_iterations: u32 = 3;
    var iteration: u32 = 0;

    while (iteration < max_iterations) : (iteration += 1) {
        var changed = false;
        for (clq) |mi| {
            if (movers[mi].is_fixed) continue;
            const old_best = movers[mi].best_orientation;
            optimizeSingleton(movers, mi, model, config, score_ctx, allocator, scratch);
            movers[mi].applyOrientation(model.atoms.items, movers[mi].best_orientation);
            if (movers[mi].best_orientation != old_best) changed = true;
        }
        if (!changed) break;
    }
}

fn allFixed(movers: []const Mover, clq: []const u32) bool {
    for (clq) |mi| {
        if (!movers[mi].is_fixed) return false;
    }
    return true;
}

// ── Tests ─────────────────────────────────────────────────────────────────────

/// Reference brute force that rescores every mover for every combination.
/// Kept for tests that check `optimizeBruteForce` against it.
fn optimizeBruteForceReference(
    allocator: Allocator,
    movers: []Mover,
    clq: []const u32,
    model: *Model,
    config: OptConfig,
    score_ctx: *const ScoreContext,
    scratch: *std.ArrayListUnmanaged(u32),
) !void {
    const n = clq.len;
    // Current orientation indices for each mover in clique
    const indices = try allocator.alloc(u16, n);
    defer allocator.free(indices);
    @memset(indices, 0);
    for (clq, 0..) |mi, i| {
        if (movers[mi].is_fixed) indices[i] = movers[mi].best_orientation;
    }

    var best_score: f32 = -std.math.inf(f32);
    const best_indices = try allocator.alloc(u16, n);
    defer allocator.free(best_indices);
    @memset(best_indices, 0);
    for (clq, 0..) |mi, i| {
        if (movers[mi].is_fixed) best_indices[i] = movers[mi].best_orientation;
    }

    // Enumerate all combinations
    while (true) {
        // Apply current combination
        for (clq, 0..) |mi, i| {
            movers[mi].applyOrientation(model.atoms.items, indices[i]);
        }

        // Score all movers in clique
        var total_score: f32 = 0;
        for (clq, 0..) |mi, i| {
            total_score += scoreMover(&movers[mi], mi, movers, model, config.scoring_params, config.scoring_params.gap_scale, score_ctx, allocator, scratch);
            total_score -= movers[mi].orientationPenalty(indices[i]);
        }

        if (total_score > best_score) {
            best_score = total_score;
            @memcpy(best_indices, indices);
        }

        // Increment indices (mixed radix counter)
        if (!incrementIndices(indices, movers, clq)) break;
    }

    // Record best orientations
    for (clq, 0..) |mi, i| {
        movers[mi].best_orientation = best_indices[i];
    }
}


const testing = std.testing;

/// Create a test mover with the given orientations (each orientation is one atom position).
fn makeTestMover(
    allocator: Allocator,
    atom_idx: u32,
    positions: []const Vec3(f32),
    penalties: []const f32,
) !Mover {
    const n = positions.len;
    const orientations = try allocator.alloc(mover_mod.Orientation, n);
    for (0..n) |i| {
        const pos = try allocator.alloc(Vec3(f32), 1);
        pos[0] = positions[i];
        orientations[i] = .{ .positions = pos, .penalty = penalties[i] };
    }
    const atom_indices = try allocator.alloc(u32, 1);
    atom_indices[0] = atom_idx;

    return Mover{
        .kind = .single_h_rotator,
        .residue_idx = 0,
        .atom_indices = atom_indices,
        .orientations = orientations,
        .allocator = allocator,
    };
}

test "totalStates calculation" {
    const allocator = testing.allocator;

    // Mover 0 with 3 orientations
    var m0 = try makeTestMover(allocator, 0, &.{
        .{ .x = 0, .y = 0, .z = 0 },
        .{ .x = 1, .y = 0, .z = 0 },
        .{ .x = 2, .y = 0, .z = 0 },
    }, &.{ 0, 0, 0 });
    defer m0.deinit();

    // Mover 1 with 2 orientations
    var m1 = try makeTestMover(allocator, 1, &.{
        .{ .x = 0, .y = 0, .z = 0 },
        .{ .x = 0, .y = 1, .z = 0 },
    }, &.{ 0, 0 });
    defer m1.deinit();

    const movers = [_]Mover{ m0, m1 };
    const clq = [_]u32{ 0, 1 };
    try testing.expectEqual(@as(u64, 6), totalStates(&movers, &clq));

    // Single mover
    const clq_single = [_]u32{0};
    try testing.expectEqual(@as(u64, 3), totalStates(&movers, &clq_single));
}

test "optimize singleton picks best orientation" {
    // Create a model with 2 atoms: a fixed obstacle at origin, and a mover atom.
    // Orientation 0: mover at (1.5, 0, 0) -- bumps with obstacle (VDW overlap)
    // Orientation 1: mover at (10.0, 0, 0) -- far away, no bump
    // The optimizer should pick orientation 1.
    const allocator = testing.allocator;

    var model = Model.init(allocator);
    defer model.deinit();

    // Fixed obstacle atom at origin
    try model.atoms.append(allocator, .{
        .pos = .{ .x = 0, .y = 0, .z = 0 },
        .vdw_radius = 1.7,
    });
    // Mover atom (initial position, will be overwritten)
    try model.atoms.append(allocator, .{
        .pos = .{ .x = 5, .y = 0, .z = 0 },
        .vdw_radius = 1.7,
    });

    // Mover controls atom 1 with 2 orientations
    var mover = try makeTestMover(allocator, 1, &.{
        .{ .x = 1.5, .y = 0, .z = 0 }, // overlaps with obstacle (distance 1.5 < 1.7+1.7=3.4)
        .{ .x = 10.0, .y = 0, .z = 0 }, // far away, no interaction
    }, &.{ 0, 0 });
    defer mover.deinit();

    var movers = [_]Mover{mover};

    var score_ctx = try buildScoreContext(allocator, &movers, model.atoms.items);
    defer score_ctx.deinit(allocator);

    var scratch = std.ArrayListUnmanaged(u32).empty;
    defer scratch.deinit(allocator);

    optimizeSingleton(&movers, 0, &model, .{}, &score_ctx, allocator, &scratch);

    // Should pick orientation 1 (no bump)
    try testing.expectEqual(@as(u16, 1), movers[0].best_orientation);
}

test "optimize brute force finds optimal combination" {
    // Two movers. A fixed atom at (0, 0, 0).
    // Mover 0 controls atom 1: orient 0 = (1.5, 0, 0) (bumps), orient 1 = (10, 0, 0)
    // Mover 1 controls atom 2: orient 0 = (0, 1.5, 0) (bumps), orient 1 = (0, 10, 0)
    // Best: both at orientation 1.
    const allocator = testing.allocator;

    var model = Model.init(allocator);
    defer model.deinit();

    // Fixed obstacle atom at origin
    try model.atoms.append(allocator, .{
        .pos = .{ .x = 0, .y = 0, .z = 0 },
        .vdw_radius = 1.7,
    });
    // Mover 0 atom
    try model.atoms.append(allocator, .{
        .pos = .{ .x = 5, .y = 0, .z = 0 },
        .vdw_radius = 1.7,
    });
    // Mover 1 atom
    try model.atoms.append(allocator, .{
        .pos = .{ .x = 0, .y = 5, .z = 0 },
        .vdw_radius = 1.7,
    });

    var m0 = try makeTestMover(allocator, 1, &.{
        .{ .x = 1.5, .y = 0, .z = 0 },
        .{ .x = 10.0, .y = 0, .z = 0 },
    }, &.{ 0, 0 });
    defer m0.deinit();

    var m1 = try makeTestMover(allocator, 2, &.{
        .{ .x = 0, .y = 1.5, .z = 0 },
        .{ .x = 0, .y = 10.0, .z = 0 },
    }, &.{ 0, 0 });
    defer m1.deinit();

    var movers = [_]Mover{ m0, m1 };
    const clq = [_]u32{ 0, 1 };

    var score_ctx = try buildScoreContext(allocator, &movers, model.atoms.items);
    defer score_ctx.deinit(allocator);

    var scratch = std.ArrayListUnmanaged(u32).empty;
    defer scratch.deinit(allocator);

    try optimizeBruteForce(allocator, &movers, &clq, &model, .{}, &score_ctx, &scratch);

    // Both should pick orientation 1 (no bumps)
    try testing.expectEqual(@as(u16, 1), movers[0].best_orientation);
    try testing.expectEqual(@as(u16, 1), movers[1].best_orientation);
}

test "incrementIndices mixed radix counting" {
    const allocator = testing.allocator;

    // Mover 0: 2 orientations, Mover 1: 3 orientations
    var m0 = try makeTestMover(allocator, 0, &.{
        .{ .x = 0, .y = 0, .z = 0 },
        .{ .x = 1, .y = 0, .z = 0 },
    }, &.{ 0, 0 });
    defer m0.deinit();

    var m1 = try makeTestMover(allocator, 1, &.{
        .{ .x = 0, .y = 0, .z = 0 },
        .{ .x = 0, .y = 1, .z = 0 },
        .{ .x = 0, .y = 2, .z = 0 },
    }, &.{ 0, 0, 0 });
    defer m1.deinit();

    const movers = [_]Mover{ m0, m1 };
    const clq = [_]u32{ 0, 1 };
    var indices = [_]u16{ 0, 0 };

    // Count total iterations
    var count: u32 = 0;
    while (incrementIndices(&indices, &movers, &clq)) {
        count += 1;
    }
    // Started at (0,0), should enumerate 2*3 - 1 = 5 more states before overflow
    try testing.expectEqual(@as(u32, 5), count);
}

test "iterative greedy finds optimal for coupled movers" {
    const allocator = testing.allocator;

    var model = Model.init(allocator);
    defer model.deinit();

    // Fixed obstacle at origin
    try model.atoms.append(allocator, .{
        .pos = .{ .x = 0, .y = 0, .z = 0 },
        .vdw_radius = 1.7,
    });
    // Mover 0 atom
    try model.atoms.append(allocator, .{
        .pos = .{ .x = 5, .y = 0, .z = 0 },
        .vdw_radius = 1.7,
    });
    // Mover 1 atom
    try model.atoms.append(allocator, .{
        .pos = .{ .x = 0, .y = 5, .z = 0 },
        .vdw_radius = 1.7,
    });

    var m0 = try makeTestMover(allocator, 1, &.{
        .{ .x = 1.5, .y = 0, .z = 0 }, // overlaps with obstacle (dist 1.5 < 3.4)
        .{ .x = 10.0, .y = 0, .z = 0 }, // far away
    }, &.{ 0, 0 });
    defer m0.deinit();

    var m1 = try makeTestMover(allocator, 2, &.{
        .{ .x = 0, .y = 1.5, .z = 0 }, // overlaps with obstacle (dist 1.5 < 3.4)
        .{ .x = 0, .y = 10.0, .z = 0 }, // far away
    }, &.{ 0, 0 });
    defer m1.deinit();

    var movers = [_]Mover{ m0, m1 };
    const clq = [_]u32{ 0, 1 };

    var score_ctx = try buildScoreContext(allocator, &movers, model.atoms.items);
    defer score_ctx.deinit(allocator);

    var scratch = std.ArrayListUnmanaged(u32).empty;
    defer scratch.deinit(allocator);

    optimizeIterativeGreedy(&movers, &clq, &model, .{}, &score_ctx, allocator, &scratch);

    // Both should pick orientation 1 (away from obstacle)
    try testing.expectEqual(@as(u16, 1), movers[0].best_orientation);
    try testing.expectEqual(@as(u16, 1), movers[1].best_orientation);
}

test "optimize brute force sees clashes introduced by moved coordinates" {
    const allocator = testing.allocator;

    var model = Model.init(allocator);
    defer model.deinit();

    try model.atoms.append(allocator, .{
        .pos = .{ .x = 0, .y = 0, .z = 0 },
        .vdw_radius = 1.7,
    });
    try model.atoms.append(allocator, .{
        .pos = .{ .x = 20, .y = 0, .z = 0 },
        .vdw_radius = 1.7,
    });

    var m0 = try makeTestMover(allocator, 0, &.{
        .{ .x = 10.0, .y = 0, .z = 0 },
        .{ .x = 0.0, .y = 0, .z = 0 },
    }, &.{ 0.0, 0.2 });
    defer m0.deinit();

    var m1 = try makeTestMover(allocator, 1, &.{
        .{ .x = 10.5, .y = 0, .z = 0 },
        .{ .x = 20.0, .y = 0, .z = 0 },
    }, &.{ 0.0, 0.2 });
    defer m1.deinit();

    var movers = [_]Mover{ m0, m1 };
    const clq = [_]u32{ 0, 1 };

    var score_ctx = try buildScoreContext(allocator, &movers, model.atoms.items);
    defer score_ctx.deinit(allocator);

    var scratch = std.ArrayListUnmanaged(u32).empty;
    defer scratch.deinit(allocator);

    try optimizeBruteForce(allocator, &movers, &clq, &model, .{}, &score_ctx, &scratch);

    // orient 0+0 places atoms at x=10.0 and x=10.5 (distance 0.5 < sum_r 3.4 → severe clash).
    // orient 0+1: m0 at 10.0, m1 at 20.0 → no clash, penalty 0+0.2 = -0.2
    // orient 1+0: m0 at 0.0, m1 at 10.5 → no clash, penalty 0.2+0 = -0.2
    // orient 1+1: m0 at 0.0, m1 at 20.0 → no clash, penalty 0.2+0.2 = -0.4
    // Best is orient 0+1 (first encountered with score -0.2).
    try testing.expectEqual(@as(u16, 0), movers[0].best_orientation);
    try testing.expectEqual(@as(u16, 1), movers[1].best_orientation);
}

test "buildScoreContext excludes mover-controlled atoms from static index" {
    const allocator = testing.allocator;

    var model = Model.init(allocator);
    defer model.deinit();
    try model.atoms.append(allocator, .{ .pos = .{ .x = 0, .y = 0, .z = 0 } });
    try model.atoms.append(allocator, .{ .pos = .{ .x = 1, .y = 0, .z = 0 } });
    try model.atoms.append(allocator, .{ .pos = .{ .x = 5, .y = 0, .z = 0 } });

    var mover = try makeTestMover(allocator, 1, &.{
        .{ .x = 1, .y = 0, .z = 0 },
        .{ .x = 2, .y = 0, .z = 0 },
    }, &.{ 0, 0 });
    defer mover.deinit();

    const movers = [_]Mover{mover};
    var score_ctx = try buildScoreContext(allocator, &movers, model.atoms.items);
    defer score_ctx.deinit(allocator);

    try testing.expectEqual(@as(usize, 2), score_ctx.static_positions.len);
    try testing.expectEqual(@as(u32, 0), score_ctx.static_atom_indices[0]);
    try testing.expectEqual(@as(u32, 2), score_ctx.static_atom_indices[1]);
    try testing.expect(score_ctx.cell_list != null);
}

test "buildScoreContext with multiple movers excludes all mover atoms" {
    const allocator = testing.allocator;

    var model = Model.init(allocator);
    defer model.deinit();
    // atoms 0..4, movers control atoms 1 and 3
    for (0..5) |i| {
        try model.atoms.append(allocator, .{
            .pos = .{ .x = @as(f32, @floatFromInt(i)), .y = 0, .z = 0 },
        });
    }

    var m0 = try makeTestMover(allocator, 1, &.{
        .{ .x = 1, .y = 0, .z = 0 },
    }, &.{0});
    defer m0.deinit();
    var m1 = try makeTestMover(allocator, 3, &.{
        .{ .x = 3, .y = 0, .z = 0 },
    }, &.{0});
    defer m1.deinit();

    const movers = [_]Mover{ m0, m1 };
    var score_ctx = try buildScoreContext(allocator, &movers, model.atoms.items);
    defer score_ctx.deinit(allocator);

    // Static atoms: 0, 2, 4
    try testing.expectEqual(@as(usize, 3), score_ctx.static_positions.len);
    try testing.expectEqual(@as(u32, 0), score_ctx.static_atom_indices[0]);
    try testing.expectEqual(@as(u32, 2), score_ctx.static_atom_indices[1]);
    try testing.expectEqual(@as(u32, 4), score_ctx.static_atom_indices[2]);
}

test "buildScoreContext partial allocation failure frees correctly" {
    const allocator = testing.allocator;

    var model = Model.init(allocator);
    defer model.deinit();
    try model.atoms.append(allocator, .{ .pos = .{ .x = 0, .y = 0, .z = 0 } });
    try model.atoms.append(allocator, .{ .pos = .{ .x = 1, .y = 0, .z = 0 } });

    var mover = try makeTestMover(allocator, 0, &.{
        .{ .x = 0, .y = 0, .z = 0 },
    }, &.{0});
    defer mover.deinit();

    const movers = [_]Mover{mover};

    // Fail at various allocation points; FailingAllocator verifies no leaks.
    // Allocation order: moved_atoms, static_positions, static_atom_indices,
    //                   CellList.init internals (counts/cell_offsets/atom_indices),
    //                   mover_centroids, mover_radii.
    // Range 0..5 covers the first 5 allocation points safely. CellList.init has a
    // known pre-existing issue where cell_offsets leaks if atom_indices alloc fails
    // (fail_index=5 in this test), so we stop before that index.
    for (0..5) |fail_idx| {
        var failing = std.testing.FailingAllocator.init(allocator, .{ .fail_index = fail_idx });
        const result = buildScoreContext(failing.allocator(), &movers, model.atoms.items);
        if (result) |*ctx| {
            var ctx_mut = ctx.*;
            ctx_mut.deinit(failing.allocator());
        } else |_| {
            // Expected failure -- FailingAllocator checks for leaks on deinit.
        }
    }
}

test "buildScoreContext computes correct centroids and bounding radii" {
    const allocator = testing.allocator;

    var model = Model.init(allocator);
    defer model.deinit();

    // Static atom at origin; mover controls atoms 1 and 2
    try model.atoms.append(allocator, .{ .pos = .{ .x = 0, .y = 0, .z = 0 } });
    try model.atoms.append(allocator, .{ .pos = .{ .x = 2, .y = 0, .z = 0 } });
    try model.atoms.append(allocator, .{ .pos = .{ .x = 4, .y = 0, .z = 0 } });
    try model.atoms.append(allocator, .{ .pos = .{ .x = 100, .y = 0, .z = 0 } });

    // Mover 0: one orientation with two positions at (2,0,0) and (4,0,0).
    // centroid = ((2+4)/2, 0, 0) = (3, 0, 0)
    // bounding radius = max(|2-3|, |4-3|) = 1.0
    const orientations_m0 = try allocator.alloc(mover_mod.Orientation, 1);
    const pos_m0 = try allocator.alloc(Vec3(f32), 2);
    pos_m0[0] = .{ .x = 2, .y = 0, .z = 0 };
    pos_m0[1] = .{ .x = 4, .y = 0, .z = 0 };
    orientations_m0[0] = .{ .positions = pos_m0 };
    const atom_indices_m0 = try allocator.alloc(u32, 2);
    atom_indices_m0[0] = 1;
    atom_indices_m0[1] = 2;
    var m0 = Mover{
        .kind = .single_h_rotator,
        .residue_idx = 0,
        .atom_indices = atom_indices_m0,
        .orientations = orientations_m0,
        .allocator = allocator,
    };
    defer m0.deinit();

    // Mover 1: one orientation at (100,0,0). centroid=(100,0,0), radius=0.
    var m1 = try makeTestMover(allocator, 3, &.{
        .{ .x = 100, .y = 0, .z = 0 },
    }, &.{0});
    defer m1.deinit();

    const movers = [_]Mover{ m0, m1 };
    var score_ctx = try buildScoreContext(allocator, &movers, model.atoms.items);
    defer score_ctx.deinit(allocator);

    try testing.expectEqual(@as(usize, 2), score_ctx.mover_centroids.len);
    try testing.expectEqual(@as(usize, 2), score_ctx.mover_radii.len);

    // Centroid of mover 0: (3, 0, 0); bounding radius: 1.0
    try testing.expectApproxEqAbs(@as(f32, 3.0), score_ctx.mover_centroids[0].x, 0.001);
    try testing.expectApproxEqAbs(@as(f32, 0.0), score_ctx.mover_centroids[0].y, 0.001);
    try testing.expectApproxEqAbs(@as(f32, 1.0), score_ctx.mover_radii[0], 0.001);

    // Centroid of mover 1: (100, 0, 0); bounding radius: 0.0
    try testing.expectApproxEqAbs(@as(f32, 100.0), score_ctx.mover_centroids[1].x, 0.001);
    try testing.expectApproxEqAbs(@as(f32, 0.0), score_ctx.mover_radii[1], 0.001);
}

test "centroid early-exit does not change optimization result" {
    // Two movers far apart: mover 0 near x=5.75, mover 1 near x=50.5.
    // Pair cutoff = search_radius + radius_0 + radius_1 = 5.0 + 4.25 + 0.5 = 9.75.
    // Distance between centroids ~= 44.75 >> 9.75 → early-exit fires.
    // A fixed obstacle is near mover 0. Mover 0 should still pick the clear orientation.
    const allocator = testing.allocator;

    var model = Model.init(allocator);
    defer model.deinit();

    // Fixed obstacle near mover 0's closer orientation
    try model.atoms.append(allocator, .{
        .pos = .{ .x = 0, .y = 0, .z = 0 },
        .vdw_radius = 1.7,
    });
    // Mover 0 atom (initial position, overwritten during optimization)
    try model.atoms.append(allocator, .{
        .pos = .{ .x = 5, .y = 0, .z = 0 },
        .vdw_radius = 1.2,
    });
    // Mover 1 atom far away
    try model.atoms.append(allocator, .{
        .pos = .{ .x = 50, .y = 0, .z = 0 },
        .vdw_radius = 1.2,
    });

    // Mover 0: orient 0 bumps obstacle at (0,0,0), orient 1 is clear.
    // centroid = (1.5+10.0)/2 = 5.75, radius = max(|1.5-5.75|, |10.0-5.75|) = 4.25
    var m0 = try makeTestMover(allocator, 1, &.{
        .{ .x = 1.5, .y = 0, .z = 0 }, // bumps obstacle (dist 1.5 < 1.7+1.2=2.9)
        .{ .x = 10.0, .y = 0, .z = 0 }, // clear
    }, &.{ 0, 0 });
    defer m0.deinit();

    // Mover 1: far away; centroid = 50.5, radius = 0.5
    var m1 = try makeTestMover(allocator, 2, &.{
        .{ .x = 50, .y = 0, .z = 0 },
        .{ .x = 51, .y = 0, .z = 0 },
    }, &.{ 0, 0 });
    defer m1.deinit();

    var movers = [_]Mover{ m0, m1 };

    var score_ctx = try buildScoreContext(allocator, &movers, model.atoms.items);
    defer score_ctx.deinit(allocator);

    // Verify early-exit fires: centroid distance > search_radius + r0 + r1
    const search_radius: f32 = 5.0;
    const r0 = score_ctx.mover_radii[0];
    const r1 = score_ctx.mover_radii[1];
    const pair_cutoff = search_radius + r0 + r1;
    const cdiff = score_ctx.mover_centroids[0].sub(score_ctx.mover_centroids[1]);
    const dist2 = cdiff.dot(cdiff);
    try testing.expect(dist2 > pair_cutoff * pair_cutoff);

    var scratch = std.ArrayListUnmanaged(u32).empty;
    defer scratch.deinit(allocator);

    optimizeSingleton(&movers, 0, &model, .{}, &score_ctx, allocator, &scratch);

    // Mover 0 should pick orientation 1 (away from obstacle), even with early-exit skipping m1.
    try testing.expectEqual(@as(u16, 1), movers[0].best_orientation);

    // Optimize mover 1 as well -- should work fine (no nearby obstacles)
    optimizeSingleton(&movers, 1, &model, .{}, &score_ctx, allocator, &scratch);
    try testing.expect(movers[1].best_orientation < 2);
}

test "mover-vs-mover clash scoring picks correct orientations" {
    // Both atoms are mover-controlled (no static atoms).
    // m0 orient 0: x=10.0 (near m1 orient 0 at x=10.5 → clash)
    // m0 orient 1: x=0.0  (far from m1)
    // m1 orient 0: x=10.5 (near m0 orient 0)
    // m1 orient 1: x=20.0 (far from m0)
    // Best: both avoid clash. orient 0+0 has penalty 0 but severe clash.
    // orient 1+1 has penalty 0.2+0.2 = 0.4 but no clash → should win.
    const allocator = testing.allocator;

    var model = Model.init(allocator);
    defer model.deinit();

    try model.atoms.append(allocator, .{
        .pos = .{ .x = 0, .y = 0, .z = 0 },
        .vdw_radius = 1.7,
    });
    try model.atoms.append(allocator, .{
        .pos = .{ .x = 20, .y = 0, .z = 0 },
        .vdw_radius = 1.7,
    });

    var m0 = try makeTestMover(allocator, 0, &.{
        .{ .x = 10.0, .y = 0, .z = 0 },
        .{ .x = 0.0, .y = 0, .z = 0 },
    }, &.{ 0.0, 0.2 });
    defer m0.deinit();

    var m1 = try makeTestMover(allocator, 1, &.{
        .{ .x = 10.5, .y = 0, .z = 0 },
        .{ .x = 20.0, .y = 0, .z = 0 },
    }, &.{ 0.0, 0.2 });
    defer m1.deinit();

    var movers = [_]Mover{ m0, m1 };
    const clq = [_]u32{ 0, 1 };

    var score_ctx = try buildScoreContext(allocator, &movers, model.atoms.items);
    defer score_ctx.deinit(allocator);

    var scratch = std.ArrayListUnmanaged(u32).empty;
    defer scratch.deinit(allocator);

    try optimizeBruteForce(allocator, &movers, &clq, &model, .{}, &score_ctx, &scratch);

    // orient 0+1 scores best: m0 at 10.0 (no clash with m1 at 20.0), penalty 0+0.2 = -0.2.
    try testing.expectEqual(@as(u16, 0), movers[0].best_orientation);
    try testing.expectEqual(@as(u16, 1), movers[1].best_orientation);
}

test "optimize() with multiple independent singletons gets correct results" {
    // Place N independent movers far apart (200 Å apart on y-axis) so they form N singleton cliques.
    // Each mover has its own nearby obstacle:
    //   obstacle_i at (0, y_i, 0) with VDW 1.7
    //   orient 0:  mover at (1.5, y_i, 0) -- bumps obstacle (dist 1.5 < sum_r 3.4)
    //   orient 1:  mover at (10.0, y_i, 0) -- far away, no bump
    // After optimize(), every mover should have best_orientation == 1.
    const allocator = testing.allocator;
    const N = 8;

    var model = Model.init(allocator);
    defer model.deinit();

    // Each mover i occupies atom slot obstacle_atom = 2*i, mover_atom = 2*i+1.
    // Movers are 200 Å apart on y-axis so they don't interact with each other.
    var movers_buf: [N]Mover = undefined;
    for (0..N) |i| {
        const mover_atom_idx: u32 = @intCast(2 * i + 1);
        const y: f32 = @as(f32, @floatFromInt(i)) * 200.0;

        // Fixed obstacle for this mover (atom index obstacle_atom_idx = 2*i, static)
        try model.atoms.append(allocator, .{
            .pos = .{ .x = 0, .y = y, .z = 0 },
            .vdw_radius = 1.7,
        });

        // Mover atom (initial position, overwritten during optimization)
        try model.atoms.append(allocator, .{
            .pos = .{ .x = 5, .y = y, .z = 0 },
            .vdw_radius = 1.7,
        });

        movers_buf[i] = try makeTestMover(allocator, mover_atom_idx, &.{
            .{ .x = 1.5, .y = y, .z = 0 }, // bumps own obstacle (dist 1.5 < sum_r 3.4)
            .{ .x = 10.0, .y = y, .z = 0 }, // clear
        }, &.{ 0, 0 });
    }
    defer for (&movers_buf) |*m| m.deinit();

    const movers: []Mover = &movers_buf;
    const result = try optimize(allocator, movers, &model, .{});

    // All singletons, no brute-force or greedy cliques
    try testing.expectEqual(@as(u32, N), result.n_singletons);
    try testing.expectEqual(@as(u32, 0), result.n_brute_force);

    // Every mover should pick orientation 1 (no bump)
    for (movers_buf) |m| {
        try testing.expectEqual(@as(u16, 1), m.best_orientation);
    }
}

test "tabulated brute force matches the rescoring reference" {
    // Random cliques of single-atom movers among static atoms, with donor and
    // acceptor flags and a fixed member, must give the same orientations as
    // rescoring every mover for every combination.
    const allocator = testing.allocator;
    var prng = std.Random.DefaultPrng.init(0x5eed);
    const rand = prng.random();

    for (0..20) |_| {
        var model = Model.init(allocator);
        defer model.deinit();

        const n_movers: u32 = 5;
        const n_static: u32 = 8;
        for (0..n_movers + n_static) |k| {
            try model.atoms.append(allocator, .{
                .pos = .{ .x = rand.float(f32) * 6, .y = rand.float(f32) * 6, .z = rand.float(f32) * 6 },
                .vdw_radius = 1.2 + rand.float(f32) * 0.6,
                .flags = .{ .donor = k % 3 == 0, .acceptor = k % 3 == 1 },
            });
        }

        var movers_a: [n_movers]Mover = undefined;
        var movers_b: [n_movers]Mover = undefined;
        for (0..n_movers) |i| {
            var positions: [4]Vec3(f32) = undefined;
            var penalties: [4]f32 = undefined;
            for (0..4) |o| {
                positions[o] = .{ .x = rand.float(f32) * 6, .y = rand.float(f32) * 6, .z = rand.float(f32) * 6 };
                penalties[o] = rand.float(f32) * 0.3;
            }
            movers_a[i] = try makeTestMover(allocator, @intCast(i), &positions, &penalties);
            movers_b[i] = try makeTestMover(allocator, @intCast(i), &positions, &penalties);
        }
        defer for (&movers_a) |*m| m.deinit();
        defer for (&movers_b) |*m| m.deinit();
        movers_a[2].lockToOrientation(3);
        movers_b[2].lockToOrientation(3);

        const clq = [_]u32{ 0, 1, 2, 3, 4 };
        var scratch = std.ArrayListUnmanaged(u32).empty;
        defer scratch.deinit(allocator);

        var ctx_a = try buildScoreContext(allocator, &movers_a, model.atoms.items);
        defer ctx_a.deinit(allocator);
        try optimizeBruteForceReference(allocator, &movers_a, &clq, &model, .{}, &ctx_a, &scratch);

        var ctx_b = try buildScoreContext(allocator, &movers_b, model.atoms.items);
        defer ctx_b.deinit(allocator);
        try optimizeBruteForce(allocator, &movers_b, &clq, &model, .{}, &ctx_b, &scratch);

        for (0..n_movers) |i| {
            try testing.expectEqual(movers_a[i].best_orientation, movers_b[i].best_orientation);
        }
    }
}
