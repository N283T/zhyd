//! Interaction graph and clique (connected component) detection for movers.
//! Two movers interact if any of their controlled atoms are within a cutoff distance.
//! Absent H atoms (at `Mover.ABSENT_H_POS`) are ignored.

const std = @import("std");
const Allocator = std.mem.Allocator;
const math_mod = @import("../math.zig");
const Vec3 = math_mod.Vec3;
const model_mod = @import("../model.zig");
const Atom = model_mod.Atom;
const mover_mod = @import("mover.zig");
const Mover = mover_mod.Mover;
const CellList = @import("../model/neighbor.zig").CellList;

pub const InteractionGraph = struct {
    /// Adjacency list: adj[i] = list of mover indices that interact with mover i
    adjacency: [][]u32,
    n_movers: u32,
    allocator: Allocator,

    pub fn deinit(self: *InteractionGraph) void {
        for (self.adjacency) |neighbors| {
            self.allocator.free(neighbors);
        }
        self.allocator.free(self.adjacency);
    }
};

/// Build interaction graph: edge between movers whose atoms can interact.
/// cutoff: maximum distance between any pair of atoms for two movers to be considered interacting.
/// Typically VDW_radius_sum + probe_radius + overlap_margin (~6-8 Å).
pub fn buildInteractionGraph(
    allocator: Allocator,
    movers: []const Mover,
    atoms: []const Atom,
    cutoff: f32,
) !InteractionGraph {
    const n = movers.len;
    var adj_lists = try allocator.alloc(std.ArrayListUnmanaged(u32), n);
    for (adj_lists) |*al| al.* = .empty;
    errdefer {
        for (adj_lists) |*al| al.deinit(allocator);
        allocator.free(adj_lists);
    }

    // Use CellList spatial index to avoid checking all O(M^2) mover pairs.
    // Compute a centroid for each mover and find the max distance from any
    // centroid to its atoms (max_radius). Then build a CellList over centroids
    // with cell_size = cutoff + 2*max_radius so that querying within that radius
    // gives exactly the candidate pairs that could possibly interact.
    // Fall back to brute force if the grid would be too large.
    buildWithCellList(allocator, movers, atoms, cutoff, adj_lists) catch |err| switch (err) {
        error.GridTooLarge => {
            // Brute-force fallback
            for (0..n) |i| {
                for (i + 1..n) |j| {
                    if (moversInteract(movers[i], movers[j], atoms, cutoff)) {
                        try adj_lists[i].append(allocator, @intCast(j));
                        try adj_lists[j].append(allocator, @intCast(i));
                    }
                }
            }
        },
        else => return err,
    };

    // Convert to owned slices
    const adjacency = try allocator.alloc([]u32, n);
    for (0..n) |i| {
        adjacency[i] = try adj_lists[i].toOwnedSlice(allocator);
    }
    allocator.free(adj_lists);

    return InteractionGraph{
        .adjacency = adjacency,
        .n_movers = @intCast(n),
        .allocator = allocator,
    };
}

/// CellList-accelerated inner loop for buildInteractionGraph.
/// Returns error.GridTooLarge when CellList.init fails with that error, so the
/// caller can fall back to brute force. Other errors are propagated normally.
fn buildWithCellList(
    allocator: Allocator,
    movers: []const Mover,
    atoms: []const Atom,
    cutoff: f32,
    adj_lists: []std.ArrayListUnmanaged(u32),
) !void {
    const n = movers.len;
    if (n == 0) return;

    // 1. Compute centroid for each mover.
    const centroids = try allocator.alloc(Vec3(f32), n);
    defer allocator.free(centroids);

    // Absent H atoms sit at ABSENT_H_POS and are skipped throughout; a mover
    // with no present atom gets the origin, and moversInteract never pairs it.
    for (movers, 0..) |mover, i| {
        var sum = Vec3(f32){ .x = 0, .y = 0, .z = 0 };
        var n_present: usize = 0;
        for (mover.atom_indices) |ai| {
            if (mover_mod.isAbsentPos(atoms[ai].pos)) continue;
            sum = sum.add(atoms[ai].pos);
            n_present += 1;
        }
        const count: f32 = @floatFromInt(@max(n_present, 1));
        centroids[i] = Vec3(f32){
            .x = sum.x / count,
            .y = sum.y / count,
            .z = sum.z / count,
        };
    }

    // 2. Compute max_radius: the furthest any atom is from its mover's centroid.
    var max_radius: f32 = 0.0;
    for (movers, 0..) |mover, i| {
        for (mover.atom_indices) |ai| {
            if (mover_mod.isAbsentPos(atoms[ai].pos)) continue;
            const diff = atoms[ai].pos.sub(centroids[i]);
            const d = @sqrt(diff.dot(diff));
            if (d > max_radius) max_radius = d;
        }
    }

    // 3. Build CellList over mover centroids.
    //    Query radius = cutoff + 2*max_radius covers all pairs that might interact.
    const query_radius = cutoff + 2.0 * max_radius;
    const cell_size = @max(query_radius, 1.0); // one cell spans the whole query radius

    var cell_list = try CellList.init(allocator, centroids, cell_size);
    defer cell_list.deinit();

    // 4. For each mover i, query nearby centroids and run the detailed check.
    var candidates = std.ArrayListUnmanaged(u32).empty;
    defer candidates.deinit(allocator);

    for (0..n) |i| {
        candidates.clearRetainingCapacity();
        try cell_list.neighborsInRadius(centroids[i], query_radius, &candidates, allocator, centroids);

        for (candidates.items) |j_u32| {
            const j: usize = j_u32;
            // Only check each pair once (i < j).
            if (j <= i) continue;
            if (moversInteract(movers[i], movers[j], atoms, cutoff)) {
                try adj_lists[i].append(allocator, @intCast(j));
                try adj_lists[j].append(allocator, @intCast(i));
            }
        }
    }
}

fn moversInteract(a: Mover, b: Mover, atoms: []const Atom, cutoff: f32) bool {
    const cutoff2 = cutoff * cutoff;
    for (a.atom_indices) |ai| {
        if (mover_mod.isAbsentPos(atoms[ai].pos)) continue;
        for (b.atom_indices) |bi| {
            if (mover_mod.isAbsentPos(atoms[bi].pos)) continue;
            const diff = atoms[ai].pos.sub(atoms[bi].pos);
            if (diff.dot(diff) <= cutoff2) return true;
        }
    }
    return false;
}

/// Find connected components in the interaction graph.
/// Returns a list of cliques, where each clique is a list of mover indices.
pub fn findCliques(allocator: Allocator, graph: *const InteractionGraph) ![][]u32 {
    const n = graph.n_movers;
    const visited = try allocator.alloc(bool, n);
    defer allocator.free(visited);
    @memset(visited, false);

    var cliques = std.ArrayListUnmanaged([]u32).empty;
    errdefer {
        for (cliques.items) |c| allocator.free(c);
        cliques.deinit(allocator);
    }

    for (0..n) |i| {
        if (visited[i]) continue;
        // BFS from i
        var component = std.ArrayListUnmanaged(u32).empty;
        errdefer component.deinit(allocator);
        var queue = std.ArrayListUnmanaged(u32).empty;
        defer queue.deinit(allocator);

        try queue.append(allocator, @intCast(i));
        visited[i] = true;

        var read_idx: usize = 0;
        while (read_idx < queue.items.len) {
            const node = queue.items[read_idx];
            read_idx += 1;
            try component.append(allocator, node);
            for (graph.adjacency[node]) |neighbor| {
                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    try queue.append(allocator, neighbor);
                }
            }
        }

        const slice = try component.toOwnedSlice(allocator);
        errdefer allocator.free(slice);
        try cliques.append(allocator, slice);
    }

    return try cliques.toOwnedSlice(allocator);
}

// ── Tests ─────────────────────────────────────────────────────────────────────

const testing = std.testing;

/// Create a minimal Mover for testing. Only atom_indices matters for interaction checks.
fn makeMover(allocator: Allocator, atom_indices: []const u32) !Mover {
    const indices = try allocator.dupe(u32, atom_indices);
    const orientations = try allocator.alloc(mover_mod.Orientation, 0);
    return Mover{
        .kind = .single_h_rotator,
        .residue_idx = 0,
        .atom_indices = indices,
        .orientations = orientations,
        .allocator = allocator,
    };
}

test "single isolated mover forms own clique" {
    const allocator = testing.allocator;

    var atoms = [_]Atom{
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } },
    };

    var m0 = try makeMover(allocator, &.{0});
    defer m0.deinit();

    const movers = [_]Mover{m0};
    var graph = try buildInteractionGraph(allocator, &movers, &atoms, 6.0);
    defer graph.deinit();

    const cliques = try findCliques(allocator, &graph);
    defer {
        for (cliques) |c| allocator.free(c);
        allocator.free(cliques);
    }

    try testing.expectEqual(@as(usize, 1), cliques.len);
    try testing.expectEqual(@as(usize, 1), cliques[0].len);
    try testing.expectEqual(@as(u32, 0), cliques[0][0]);
}

test "two interacting movers form one clique" {
    const allocator = testing.allocator;

    // Two atoms close together (distance = 2.0, well within 6.0 cutoff)
    var atoms = [_]Atom{
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } },
        .{ .pos = .{ .x = 2.0, .y = 0.0, .z = 0.0 } },
    };

    var m0 = try makeMover(allocator, &.{0});
    defer m0.deinit();
    var m1 = try makeMover(allocator, &.{1});
    defer m1.deinit();

    const movers = [_]Mover{ m0, m1 };
    var graph = try buildInteractionGraph(allocator, &movers, &atoms, 6.0);
    defer graph.deinit();

    const cliques = try findCliques(allocator, &graph);
    defer {
        for (cliques) |c| allocator.free(c);
        allocator.free(cliques);
    }

    try testing.expectEqual(@as(usize, 1), cliques.len);
    try testing.expectEqual(@as(usize, 2), cliques[0].len);
}

test "two distant movers form two cliques" {
    const allocator = testing.allocator;

    // Two atoms far apart (distance = 100.0, beyond 6.0 cutoff)
    var atoms = [_]Atom{
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } },
        .{ .pos = .{ .x = 100.0, .y = 0.0, .z = 0.0 } },
    };

    var m0 = try makeMover(allocator, &.{0});
    defer m0.deinit();
    var m1 = try makeMover(allocator, &.{1});
    defer m1.deinit();

    const movers = [_]Mover{ m0, m1 };
    var graph = try buildInteractionGraph(allocator, &movers, &atoms, 6.0);
    defer graph.deinit();

    const cliques = try findCliques(allocator, &graph);
    defer {
        for (cliques) |c| allocator.free(c);
        allocator.free(cliques);
    }

    try testing.expectEqual(@as(usize, 2), cliques.len);
    try testing.expectEqual(@as(usize, 1), cliques[0].len);
    try testing.expectEqual(@as(usize, 1), cliques[1].len);
}

test "chain of interactions" {
    const allocator = testing.allocator;

    // 3 movers: A at x=0, B at x=4 (within 6 of A), C at x=8 (within 6 of B, not of A directly)
    // A-B distance = 4, B-C distance = 4, A-C distance = 8
    var atoms = [_]Atom{
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } }, // mover A
        .{ .pos = .{ .x = 4.0, .y = 0.0, .z = 0.0 } }, // mover B
        .{ .pos = .{ .x = 8.0, .y = 0.0, .z = 0.0 } }, // mover C
    };

    var m0 = try makeMover(allocator, &.{0});
    defer m0.deinit();
    var m1 = try makeMover(allocator, &.{1});
    defer m1.deinit();
    var m2 = try makeMover(allocator, &.{2});
    defer m2.deinit();

    const movers = [_]Mover{ m0, m1, m2 };
    var graph = try buildInteractionGraph(allocator, &movers, &atoms, 6.0);
    defer graph.deinit();

    // Verify A-C do not directly interact (distance 8 > 6 cutoff)
    try testing.expectEqual(@as(usize, 1), graph.adjacency[0].len); // A only connects to B
    try testing.expectEqual(@as(usize, 2), graph.adjacency[1].len); // B connects to A and C

    const cliques = try findCliques(allocator, &graph);
    defer {
        for (cliques) |c| allocator.free(c);
        allocator.free(cliques);
    }

    // All three form one clique through transitivity
    try testing.expectEqual(@as(usize, 1), cliques.len);
    try testing.expectEqual(@as(usize, 3), cliques[0].len);
}

test "absent H atoms do not connect distant movers" {
    const allocator = testing.allocator;

    // Two His-like movers 100 A apart, each with one H parked at ABSENT_H_POS.
    // The sentinel atoms coincide, but the movers must stay separate.
    var atoms = [_]Atom{
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } },
        .{ .pos = Mover.ABSENT_H_POS, .is_added = true },
        .{ .pos = .{ .x = 100.0, .y = 0.0, .z = 0.0 } },
        .{ .pos = Mover.ABSENT_H_POS, .is_added = true },
    };

    var m0 = try makeMover(allocator, &.{ 0, 1 });
    defer m0.deinit();
    var m1 = try makeMover(allocator, &.{ 2, 3 });
    defer m1.deinit();

    const movers = [_]Mover{ m0, m1 };
    var graph = try buildInteractionGraph(allocator, &movers, &atoms, 6.0);
    defer graph.deinit();

    const cliques = try findCliques(allocator, &graph);
    defer {
        for (cliques) |c| allocator.free(c);
        allocator.free(cliques);
    }

    try testing.expectEqual(@as(usize, 2), cliques.len);
}
