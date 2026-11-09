//! Scoring infrastructure for the optimization engine.
//! Contains ScoreContext, pair-scoring functions, and mover-level scoring.
//! Separated from optimizer.zig to keep each file focused and within size limits.

const std = @import("std");
const Allocator = std.mem.Allocator;
const math_mod = @import("../math.zig");
const Vec3 = math_mod.Vec3;
const model_mod = @import("../model.zig");
const Model = model_mod.Model;
const Atom = model_mod.Atom;
const CellList = model_mod.CellList;
const mover_mod = @import("mover.zig");
const Mover = mover_mod.Mover;
const scorer_mod = @import("scorer.zig");
const element = @import("../element.zig");

pub const ScoreContext = struct {
    cell_list: ?CellList,
    static_positions: []math_mod.Vec3(f32),
    static_atom_indices: []u32,
    /// VDW radii for static atoms (parallel to static_positions).
    static_radii: []f32,
    /// AtomFlags for static atoms (parallel to static_positions).
    static_flags: []element.AtomFlags,
    /// Mean position of each mover's atoms across all its orientations.
    mover_centroids: []math_mod.Vec3(f32),
    /// Bounding radius of each mover: max distance from centroid to any orientation position.
    mover_radii: []f32,
    /// VDW radii for mover-controlled atoms, indexed directly by atom index.
    /// Only entries for atoms belonging to movers are populated; others are undefined.
    /// Length equals the total number of atoms in the model.
    /// VDW radii are immutable after placement, so this cache is always valid.
    mover_atom_radii: []f32,
    /// Movers that pass the centroid early-exit test against each mover, in
    /// ascending index order: `mover_neighbors[mover_neighbor_offsets[i]..mover_neighbor_offsets[i + 1]]`.
    /// Precomputed so that scoring a mover does not scan every mover in the model.
    mover_neighbor_offsets: []u32,
    mover_neighbors: []u32,

    pub fn neighborsOf(self: *const ScoreContext, mover_idx: u32) []const u32 {
        return self.mover_neighbors[self.mover_neighbor_offsets[mover_idx]..self.mover_neighbor_offsets[mover_idx + 1]];
    }

    pub fn deinit(self: *ScoreContext, allocator: Allocator) void {
        if (self.cell_list) |*cl| cl.deinit();
        allocator.free(self.static_positions);
        allocator.free(self.static_atom_indices);
        allocator.free(self.static_radii);
        allocator.free(self.static_flags);
        allocator.free(self.mover_centroids);
        allocator.free(self.mover_radii);
        allocator.free(self.mover_atom_radii);
        allocator.free(self.mover_neighbor_offsets);
        allocator.free(self.mover_neighbors);
    }
};

/// Build a ScoreContext from the current model state.
/// Partitions atoms into static (non-mover) and mover-controlled sets,
/// builds a CellList spatial index over static atoms, and computes
/// per-mover centroids and bounding radii for early-exit distance tests.
pub fn buildScoreContext(
    allocator: Allocator,
    movers: []const Mover,
    atoms: []const Atom,
) !ScoreContext {
    // Mark mover-controlled atoms (temporary, freed before return).
    const moved_atoms = try allocator.alloc(bool, atoms.len);
    defer allocator.free(moved_atoms);
    @memset(moved_atoms, false);
    for (movers) |m| {
        for (m.atom_indices) |ai| {
            // Debug-only: verify mover atom indices are disjoint across movers.
            if (std.debug.runtime_safety and moved_atoms[ai]) {
                std.debug.panic("mover atom index {d} appears in multiple movers", .{ai});
            }
            moved_atoms[ai] = true;
        }
    }

    var static_count: usize = 0;
    for (moved_atoms) |is_moved| {
        if (!is_moved) static_count += 1;
    }

    const static_positions = try allocator.alloc(math_mod.Vec3(f32), static_count);
    errdefer allocator.free(static_positions);
    const static_atom_indices = try allocator.alloc(u32, static_count);
    errdefer allocator.free(static_atom_indices);
    const static_radii = try allocator.alloc(f32, static_count);
    errdefer allocator.free(static_radii);
    const static_flags = try allocator.alloc(element.AtomFlags, static_count);
    errdefer allocator.free(static_flags);

    var out_i: usize = 0;
    for (atoms, 0..) |a, i| {
        if (moved_atoms[i]) continue;
        static_positions[out_i] = a.pos;
        static_atom_indices[out_i] = @intCast(i);
        static_radii[out_i] = a.vdw_radius;
        static_flags[out_i] = a.flags;
        out_i += 1;
    }

    // Build SoA VDW radii cache for mover-controlled atoms.
    // Positions are not cached here because they change as orientations are applied.
    // Flags are not cached because flip movers (His) update them per orientation.
    // VDW radii are immutable after placement, so this cache is always valid and
    // reduces cache pressure in the mover-vs-mover scoring loop.
    const mover_atom_radii = try allocator.alloc(f32, atoms.len);
    errdefer allocator.free(mover_atom_radii);
    for (movers) |m| {
        for (m.atom_indices) |ai| {
            mover_atom_radii[ai] = atoms[ai].vdw_radius;
        }
    }

    // CellList may fail with GridTooLarge for very spread-out structures;
    // fall back to null (pairwise scoring only) rather than aborting.
    var cell_list: ?CellList = CellList.init(allocator, static_positions, 5.0) catch |err| blk: {
        switch (err) {
            error.GridTooLarge, error.NonFiniteCoordinate => {
                std.log.warn("CellList init failed ({s}) for static atoms; falling back to pairwise scoring", .{@errorName(err)});
                break :blk null;
            },
            error.OutOfMemory => return err,
        }
    };
    errdefer if (cell_list) |*cl| cl.deinit();

    // Compute centroid (mean position) and bounding radius for each mover.
    // Centroid is averaged over all orientation positions; bounding radius is the
    // max distance from the centroid to any orientation position.
    // Used by scoreMover for early-exit of distant mover-vs-mover pairs:
    //   skip if dist(centroid_a, centroid_b) > search_radius + radius_a + radius_b
    const mover_centroids = try allocator.alloc(math_mod.Vec3(f32), movers.len);
    errdefer allocator.free(mover_centroids);
    const mover_radii = try allocator.alloc(f32, movers.len);
    errdefer allocator.free(mover_radii);

    for (movers, 0..) |m, mi| {
        var sum = math_mod.Vec3(f32).zero;
        var count: usize = 0;
        for (m.orientations) |orient| {
            for (orient.positions) |pos| {
                if (mover_mod.isAbsentPos(pos)) continue;
                sum = sum.add(pos);
                count += 1;
            }
        }
        const centroid = if (count > 0) sum.scale(1.0 / @as(f32, @floatFromInt(count))) else math_mod.Vec3(f32).zero;
        mover_centroids[mi] = centroid;

        var max_r2: f32 = 0.0;
        for (m.orientations) |orient| {
            for (orient.positions) |pos| {
                if (mover_mod.isAbsentPos(pos)) continue;
                const d = pos.sub(centroid);
                const r2 = d.dot(d);
                if (r2 > max_r2) max_r2 = r2;
            }
        }
        mover_radii[mi] = @sqrt(max_r2);
    }

    const neighbors = try buildMoverNeighbors(allocator, mover_centroids, mover_radii);

    return .{
        .mover_neighbor_offsets = neighbors.offsets,
        .mover_neighbors = neighbors.list,
        .cell_list = cell_list,
        .static_positions = static_positions,
        .static_atom_indices = static_atom_indices,
        .static_radii = static_radii,
        .static_flags = static_flags,
        .mover_centroids = mover_centroids,
        .mover_radii = mover_radii,
        .mover_atom_radii = mover_atom_radii,
    };
}

/// Search radius around each mover atom used by mover scoring, in Angstrom.
const score_search_radius: f32 = 5.0;

/// Movers whose bounding radius is above this are paired against every mover
/// instead of through the sweep. Absent H positions are left out of the radius,
/// so this only guards unusual movers.
const sweep_radius_limit: f32 = 10.0;

fn mayInteract(centroids: []const math_mod.Vec3(f32), radii: []const f32, a: usize, b: usize) bool {
    const cutoff = score_search_radius + radii[a] + radii[b];
    const d = centroids[a].sub(centroids[b]);
    return d.dot(d) <= cutoff * cutoff;
}

/// Build the neighbor lists of `ScoreContext.mover_neighbors`. Pairs of
/// ordinary movers are found with a sweep along x; movers with a large
/// bounding radius are tested against all others. Every candidate pair goes
/// through the same test as the scoring loop used before, so the lists hold
/// exactly the movers that loop would have scored.
fn buildMoverNeighbors(
    allocator: Allocator,
    centroids: []const math_mod.Vec3(f32),
    radii: []const f32,
) !struct { offsets: []u32, list: []u32 } {
    const n = centroids.len;
    var edges = std.ArrayListUnmanaged([2]u32).empty;
    defer edges.deinit(allocator);

    const order = try allocator.alloc(u32, n);
    defer allocator.free(order);
    var n_small: usize = 0;
    for (0..n) |i| {
        if (radii[i] <= sweep_radius_limit) {
            order[n_small] = @intCast(i);
            n_small += 1;
        }
    }
    const small = order[0..n_small];
    std.mem.sort(u32, small, centroids, struct {
        fn lessThan(c: []const math_mod.Vec3(f32), a: u32, b: u32) bool {
            return c[a].x < c[b].x;
        }
    }.lessThan);
    const window = score_search_radius + 2 * sweep_radius_limit;
    for (small, 0..) |a, k| {
        for (small[k + 1 ..]) |b| {
            if (centroids[b].x - centroids[a].x > window) break;
            if (mayInteract(centroids, radii, a, b)) try edges.append(allocator, .{ a, b });
        }
    }
    for (0..n) |a| {
        if (radii[a] <= sweep_radius_limit) continue;
        for (0..n) |b| {
            if (b == a) continue;
            // A pair of two large movers is visited twice; keep it once.
            if (radii[b] > sweep_radius_limit and b < a) continue;
            if (mayInteract(centroids, radii, a, b)) try edges.append(allocator, .{ @intCast(a), @intCast(b) });
        }
    }

    const offsets = try allocator.alloc(u32, n + 1);
    errdefer allocator.free(offsets);
    @memset(offsets, 0);
    for (edges.items) |e| {
        offsets[e[0] + 1] += 1;
        offsets[e[1] + 1] += 1;
    }
    for (0..n) |i| offsets[i + 1] += offsets[i];
    const list = try allocator.alloc(u32, offsets[n]);
    errdefer allocator.free(list);
    const fill = try allocator.alloc(u32, n);
    defer allocator.free(fill);
    @memcpy(fill, offsets[0..n]);
    for (edges.items) |e| {
        list[fill[e[0]]] = e[1];
        fill[e[0]] += 1;
        list[fill[e[1]]] = e[0];
        fill[e[1]] += 1;
    }
    for (0..n) |i| std.mem.sort(u32, list[offsets[i]..offsets[i + 1]], {}, std.sort.asc(u32));
    return .{ .offsets = offsets, .list = list };
}

/// Check whether a given atom index belongs to a specific mover.
pub fn isMoverAtom(m: *const Mover, atom_idx: u32) bool {
    for (m.atom_indices) |ai| {
        if (ai == atom_idx) return true;
    }
    return false;
}

/// Score an atom pair from two Atom structs.
/// Returns a positive score for contacts/H-bonds and negative for bumps.
pub fn scorePair(a: Atom, other: Atom, scoring_params: scorer_mod.ScoringParams, gap_scale: f32) f32 {
    const diff = a.pos.sub(other.pos);
    const dist2 = diff.dot(diff);
    const sum_r = a.vdw_radius + other.vdw_radius;
    const sum_r2 = sum_r * sum_r;

    if (dist2 < sum_r2) {
        const dist = @sqrt(dist2);
        const gap = dist - sum_r;
        if (scorer_mod.isHBond(a.flags, other.flags, gap, scoring_params)) {
            return scoring_params.hb_weight * (-0.5 * gap);
        }
        return -scoring_params.bump_weight * (-0.5 * gap);
    }

    const threshold = sum_r + 0.5;
    if (dist2 < threshold * threshold) {
        const dist = @sqrt(dist2);
        const gap = dist - sum_r;
        const ratio = gap / gap_scale;
        return math_mod.fastExp(-ratio * ratio);
    }
    return 0.0;
}

/// Like scorePair but takes individual fields instead of Atom structs.
/// Used in the static-atom scoring path to read from compact SoA arrays,
/// reducing cache pressure compared to loading full Atom structs (~64B each).
pub fn scorePairSoA(
    pos_a: math_mod.Vec3(f32),
    radius_a: f32,
    flags_a: element.AtomFlags,
    pos_b: math_mod.Vec3(f32),
    radius_b: f32,
    flags_b: element.AtomFlags,
    scoring_params: scorer_mod.ScoringParams,
    gap_scale: f32,
) f32 {
    const diff = pos_a.sub(pos_b);
    const dist2 = diff.dot(diff);
    const sum_r = radius_a + radius_b;
    const sum_r2 = sum_r * sum_r;

    if (dist2 < sum_r2) {
        const dist = @sqrt(dist2);
        const gap = dist - sum_r;
        if (scorer_mod.isHBond(flags_a, flags_b, gap, scoring_params)) {
            return scoring_params.hb_weight * (-0.5 * gap);
        }
        return -scoring_params.bump_weight * (-0.5 * gap);
    }

    const threshold = sum_r + 0.5;
    if (dist2 < threshold * threshold) {
        const dist = @sqrt(dist2);
        const gap = dist - sum_r;
        const ratio = gap / gap_scale;
        return math_mod.fastExp(-ratio * ratio);
    }
    return 0.0;
}

/// Score a mover against the model using explicit positions for its atoms.
///
/// `positions_override` provides the positions for the mover's atoms in the same
/// order as `m.atom_indices`. This avoids writing to `model.atoms` during scoring,
/// which is required for thread-safe parallel optimization: concurrent threads can
/// each call this function with different orientation positions without data races.
///
/// Other movers' atoms are read from `model.atoms` at their current positions
/// (which is their initial or last-applied coarse orientation — acceptable for the
/// parallel coarse/fine search phases).
///
/// Static (non-mover) atoms are queried via CellList spatial index (O(nearby)).
/// Mover-controlled atoms are scored by direct iteration.
pub fn scoreMoverWithPositions(
    m: *const Mover,
    mover_idx: u32,
    positions_override: []const math_mod.Vec3(f32),
    movers: []const Mover,
    model: *const Model,
    scoring_params: scorer_mod.ScoringParams,
    gap_scale: f32,
    score_ctx: *const ScoreContext,
    allocator: Allocator,
    scratch: *std.ArrayListUnmanaged(u32),
) f32 {
    return scoreMoverExcluding(m, mover_idx, positions_override, movers, model, scoring_params, gap_scale, score_ctx, allocator, scratch, null);
}

/// Like `scoreMoverWithPositions`, but skips every other mover whose entry in
/// `skip_movers` is true. Brute-force search uses it for the part of a mover's
/// score that does not depend on the other movers of its clique.
pub fn scoreMoverExcluding(
    m: *const Mover,
    mover_idx: u32,
    positions_override: []const math_mod.Vec3(f32),
    movers: []const Mover,
    model: *const Model,
    scoring_params: scorer_mod.ScoringParams,
    gap_scale: f32,
    score_ctx: *const ScoreContext,
    allocator: Allocator,
    scratch: *std.ArrayListUnmanaged(u32),
    skip_movers: ?[]const bool,
) f32 {
    var total: f32 = 0;
    const atoms = model.atoms.items;
    const search_radius = score_search_radius;

    for (m.atom_indices, 0..) |ai, j| {
        // Use the override position for this mover's atom instead of atoms[ai].pos.
        // Read radius from the SoA cache (immutable). Read flags from atoms[ai] because
        // applyOrientation may have updated them (e.g. His flip changes donor/acceptor state).
        const a_pos = positions_override[j];
        if (mover_mod.isAbsentPos(a_pos)) continue;
        const a_radius = score_ctx.mover_atom_radii[ai];
        const a_flags = atoms[ai].flags;

        // Static atoms: use spatial index when available, else brute-force scan.
        // Use SoA arrays (static_radii/static_flags) to avoid loading full Atom structs.
        if (score_ctx.cell_list) |cl| {
            scratch.clearRetainingCapacity();
            cl.neighborsInRadius(a_pos, search_radius, scratch, allocator, score_ctx.static_positions) catch |err| {
                std.log.warn("neighbor query failed for atom {d}: {s}; returning -inf", .{ ai, @errorName(err) });
                return -std.math.inf(f32);
            };
            for (scratch.items) |static_idx| {
                total += scorePairSoA(
                    a_pos, a_radius, a_flags,
                    score_ctx.static_positions[static_idx],
                    score_ctx.static_radii[static_idx],
                    score_ctx.static_flags[static_idx],
                    scoring_params,
                    gap_scale,
                );
            }
        } else {
            // Fallback: pairwise scan of all static atoms (no CellList available).
            for (0..score_ctx.static_atom_indices.len) |static_idx| {
                total += scorePairSoA(
                    a_pos, a_radius, a_flags,
                    score_ctx.static_positions[static_idx],
                    score_ctx.static_radii[static_idx],
                    score_ctx.static_flags[static_idx],
                    scoring_params,
                    gap_scale,
                );
            }
        }

        // Mover-controlled atoms: score directly with current coordinates.
        // Only movers that pass the centroid early-exit test (their bounding
        // spheres can reach this mover's scoring range) are in the neighbor list.
        // Radii come from the SoA cache (immutable). Positions and flags are read
        // from model.atoms because flip movers (His) may update flags per orientation.
        for (score_ctx.neighborsOf(mover_idx)) |other_idx| {
            if (skip_movers) |skip| if (skip[other_idx]) continue;
            for (movers[other_idx].atom_indices) |oi| {
                if (oi == ai) continue;
                if (mover_mod.isAbsentPos(atoms[oi].pos)) continue;
                if (isMoverAtom(m, oi)) continue;
                total += scorePairSoA(
                    a_pos, a_radius, a_flags,
                    atoms[oi].pos,
                    score_ctx.mover_atom_radii[oi],
                    atoms[oi].flags,
                    scoring_params,
                    gap_scale,
                );
            }
        }
    }
    return total;
}

/// Whether two movers can interact at all, using the same centroid early-exit
/// test as `scoreMoverExcluding`.
pub fn moversMayInteract(score_ctx: *const ScoreContext, a: u32, b: u32) bool {
    return mayInteract(score_ctx.mover_centroids, score_ctx.mover_radii, a, b);
}

/// Score mover `a`'s atoms against mover `b`'s atoms at their current
/// positions and flags in `atoms`. This is the term that `scoreMoverExcluding`
/// adds for `b` when scoring `a`; it is not symmetric, because H-bond
/// detection reads the flags of the scored atom first.
pub fn scoreMoverPair(
    a: *const Mover,
    b: *const Mover,
    atoms: []const Atom,
    scoring_params: scorer_mod.ScoringParams,
    gap_scale: f32,
    score_ctx: *const ScoreContext,
) f32 {
    var total: f32 = 0;
    for (a.atom_indices) |ai| {
        const a_pos = atoms[ai].pos;
        if (mover_mod.isAbsentPos(a_pos)) continue;
        const a_radius = score_ctx.mover_atom_radii[ai];
        const a_flags = atoms[ai].flags;
        for (b.atom_indices) |oi| {
            if (oi == ai) continue;
            if (mover_mod.isAbsentPos(atoms[oi].pos)) continue;
            if (isMoverAtom(a, oi)) continue;
            total += scorePairSoA(
                a_pos, a_radius, a_flags,
                atoms[oi].pos,
                score_ctx.mover_atom_radii[oi],
                atoms[oi].flags,
                scoring_params,
                gap_scale,
            );
        }
    }
    return total;
}

/// Score a mover's current orientation against the model.
/// Reads atom positions from model.atoms — only safe to call sequentially
/// (e.g. brute-force and greedy clique search where orientations are applied first).
/// For parallel singleton/fine-search use scoreMoverWithPositions instead.
pub fn scoreMover(
    m: *const Mover,
    mover_idx: u32,
    movers: []const Mover,
    model: *const Model,
    scoring_params: scorer_mod.ScoringParams,
    gap_scale: f32,
    score_ctx: *const ScoreContext,
    allocator: Allocator,
    scratch: *std.ArrayListUnmanaged(u32),
) f32 {
    return scoreMoverWithPositions(
        m,
        mover_idx,
        blk: {
            // Build a slice of current positions from model.atoms for this mover.
            // This is a small, bounded allocation on the stack (up to ~16 atoms per mover).
            // NOTE: pos_buf lives on this call frame; the slice is only valid for this
            // synchronous call — do NOT pass it to a spawned thread.
            var pos_buf: [64]math_mod.Vec3(f32) = undefined;
            const n = m.atom_indices.len;
            if (n > pos_buf.len) {
                std.log.err("mover has {d} atoms, exceeding stack buffer of {d}", .{ n, pos_buf.len });
                return -std.math.inf(f32);
            }
            for (m.atom_indices, 0..) |ai, i| {
                pos_buf[i] = model.atoms.items[ai].pos;
            }
            break :blk pos_buf[0..n];
        },
        movers,
        model,
        scoring_params,
        gap_scale,
        score_ctx,
        allocator,
        scratch,
    );
}

// ── Tests ─────────────────────────────────────────────────────────────────────

const testing = std.testing;

test "scorePair returns correct values for all branches" {
    const params = scorer_mod.ScoringParams{};
    const gap_scale = params.gap_scale;

    // Branch 1: No interaction (atoms far apart) -- returns 0.0
    {
        const a = Atom{ .pos = .{ .x = 0, .y = 0, .z = 0 }, .vdw_radius = 1.7 };
        const b = Atom{ .pos = .{ .x = 10, .y = 0, .z = 0 }, .vdw_radius = 1.7 };
        try testing.expectEqual(@as(f32, 0.0), scorePair(a, b, params, gap_scale));
    }

    // Branch 2: Contact within threshold (dist between sum_r and sum_r + 0.5)
    // sum_r = 1.7 + 1.7 = 3.4, threshold = 3.9
    // place at distance 3.6: gap = 0.2
    {
        const a = Atom{ .pos = .{ .x = 0, .y = 0, .z = 0 }, .vdw_radius = 1.7 };
        const b = Atom{ .pos = .{ .x = 3.6, .y = 0, .z = 0 }, .vdw_radius = 1.7 };
        const score = scorePair(a, b, params, gap_scale);
        // gap = 0.2, ratio = 0.2 / 0.25 = 0.8, fastExp(-0.64) ≈ 0.527
        try testing.expect(score > 0.0);
        try testing.expect(@abs(score - math_mod.fastExp(@as(f32, -0.64))) < 0.001);
    }

    // Branch 3: Overlap without H-bond (bump) -- no donor/acceptor flags
    // distance 3.0 < sum_r 3.4, gap = 3.0 - 3.4 = -0.4
    {
        const a = Atom{ .pos = .{ .x = 0, .y = 0, .z = 0 }, .vdw_radius = 1.7 };
        const b = Atom{ .pos = .{ .x = 3.0, .y = 0, .z = 0 }, .vdw_radius = 1.7 };
        const score = scorePair(a, b, params, gap_scale);
        // bump: -bump_weight * (-0.5 * gap) = -10.0 * (-0.5 * -0.4) = -10.0 * 0.2 = -2.0
        try testing.expect(score < 0.0);
        try testing.expect(@abs(score - (-2.0)) < 0.001);
    }

    // Branch 4: Overlap with H-bond (donor + acceptor flags)
    // distance 3.0 < sum_r 3.4, gap = -0.4 → -gap = 0.4 <= min_reg_hb_gap(0.6) → H-bond
    {
        const a = Atom{
            .pos = .{ .x = 0, .y = 0, .z = 0 },
            .vdw_radius = 1.7,
            .flags = .{ .donor = true },
        };
        const b = Atom{
            .pos = .{ .x = 3.0, .y = 0, .z = 0 },
            .vdw_radius = 1.7,
            .flags = .{ .acceptor = true },
        };
        const score = scorePair(a, b, params, gap_scale);
        // hb: hb_weight * (-0.5 * gap) = 4.0 * (-0.5 * -0.4) = 4.0 * 0.2 = 0.8
        try testing.expect(score > 0.0);
        try testing.expect(@abs(score - 0.8) < 0.001);
    }
}

test "scorePairSoA produces identical results to scorePair" {
    const params = scorer_mod.ScoringParams{};
    const gap_scale = params.gap_scale;

    // Test all 4 branches with the same atom configurations
    const TestCase = struct {
        pos_a: math_mod.Vec3(f32), r_a: f32, flags_a: element.AtomFlags,
        pos_b: math_mod.Vec3(f32), r_b: f32, flags_b: element.AtomFlags,
    };
    const test_cases = [_]TestCase{
        // Branch 1: No interaction (far apart)
        .{ .pos_a = .{ .x = 0, .y = 0, .z = 0 }, .r_a = 1.7, .flags_a = .{},
           .pos_b = .{ .x = 10, .y = 0, .z = 0 }, .r_b = 1.7, .flags_b = .{} },
        // Branch 2: Contact (within threshold, sum_r=3.4, threshold=3.9)
        .{ .pos_a = .{ .x = 0, .y = 0, .z = 0 }, .r_a = 1.7, .flags_a = .{},
           .pos_b = .{ .x = 3.6, .y = 0, .z = 0 }, .r_b = 1.7, .flags_b = .{} },
        // Branch 3: Bump (overlap, no H-bond)
        .{ .pos_a = .{ .x = 0, .y = 0, .z = 0 }, .r_a = 1.7, .flags_a = .{},
           .pos_b = .{ .x = 3.0, .y = 0, .z = 0 }, .r_b = 1.7, .flags_b = .{} },
        // Branch 4: H-bond (overlap with donor + acceptor)
        .{ .pos_a = .{ .x = 0, .y = 0, .z = 0 }, .r_a = 1.7, .flags_a = .{ .donor = true },
           .pos_b = .{ .x = 3.0, .y = 0, .z = 0 }, .r_b = 1.7, .flags_b = .{ .acceptor = true } },
    };

    for (test_cases) |tc| {
        const atom_a = Atom{ .pos = tc.pos_a, .vdw_radius = tc.r_a, .flags = tc.flags_a };
        const atom_b = Atom{ .pos = tc.pos_b, .vdw_radius = tc.r_b, .flags = tc.flags_b };
        const score_pair_result = scorePair(atom_a, atom_b, params, gap_scale);
        const score_soa = scorePairSoA(tc.pos_a, tc.r_a, tc.flags_a, tc.pos_b, tc.r_b, tc.flags_b, params, gap_scale);
        try testing.expectEqual(score_pair_result, score_soa);
    }
}

test "buildMoverNeighbors matches the all-pairs centroid test" {
    const allocator = testing.allocator;
    var prng = std.Random.DefaultPrng.init(42);
    const rand = prng.random();
    const n = 300;
    var centroids: [n]math_mod.Vec3(f32) = undefined;
    var radii: [n]f32 = undefined;
    for (0..n) |i| {
        centroids[i] = .{ .x = rand.float(f32) * 80, .y = rand.float(f32) * 40, .z = rand.float(f32) * 40 };
        // A few movers with an absent H far away get a very large radius.
        radii[i] = if (i % 37 == 0) 900 + rand.float(f32) * 100 else rand.float(f32) * 4;
    }
    const nb = try buildMoverNeighbors(allocator, &centroids, &radii);
    defer allocator.free(nb.offsets);
    defer allocator.free(nb.list);

    var expected = std.ArrayListUnmanaged(u32).empty;
    defer expected.deinit(allocator);
    for (0..n) |a| {
        expected.clearRetainingCapacity();
        for (0..n) |b| {
            if (b != a and mayInteract(&centroids, &radii, a, b)) try expected.append(allocator, @intCast(b));
        }
        try testing.expectEqualSlices(u32, expected.items, nb.list[nb.offsets[a]..nb.offsets[a + 1]]);
    }
}
