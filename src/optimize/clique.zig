//! Interaction graph and clique (connected component) detection for movers.
//! Two movers interact if some atom of one, in any orientation, is within
//! scoring range of some atom of the other. Absent H atoms are ignored.

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

/// Build the interaction graph: an edge joins two movers when some atom of one,
/// in some orientation, comes within `r_a + r_b + contact_margin` of some atom
/// of the other in some orientation. Beyond that distance `scorePair` is zero,
/// so movers in different connected components never change each other's
/// score and each component can be optimized on its own without loss.
///
/// Positions considered for each mover atom are its current position and its
/// position in every orientation. Absent H atoms (at `Mover.ABSENT_H_POS`)
/// are ignored. Radii are the atoms' `vdw_radius`.
pub fn buildInteractionGraph(
    allocator: Allocator,
    movers: []const Mover,
    atoms: []const Atom,
    contact_margin: f32,
) !InteractionGraph {
    const n = movers.len;
    var adj_lists = try allocator.alloc(std.ArrayListUnmanaged(u32), n);
    for (adj_lists) |*al| al.* = .empty;
    errdefer {
        for (adj_lists) |*al| al.deinit(allocator);
        allocator.free(adj_lists);
    }

    // Every position a mover atom can take, with its owner and radius.
    var points = std.ArrayListUnmanaged(Vec3(f32)).empty;
    defer points.deinit(allocator);
    var owner = std.ArrayListUnmanaged(u32).empty;
    defer owner.deinit(allocator);
    var radius = std.ArrayListUnmanaged(f32).empty;
    defer radius.deinit(allocator);
    var max_radius: f32 = 0;
    for (movers, 0..) |m, mi| {
        for (m.atom_indices, 0..) |ai, j| {
            const r = atoms[ai].vdw_radius;
            max_radius = @max(max_radius, r);
            try appendPoint(allocator, &points, &owner, &radius, atoms[ai].pos, @intCast(mi), r);
            for (m.orientations) |o| {
                try appendPoint(allocator, &points, &owner, &radius, o.positions[j], @intCast(mi), r);
            }
        }
    }

    // stamp[b] == a + 1 once the edge a-b has been recorded while visiting a.
    const stamp = try allocator.alloc(u32, n);
    defer allocator.free(stamp);
    @memset(stamp, 0);

    const max_reach = 2 * max_radius + contact_margin;
    var cell_list: ?CellList = CellList.init(allocator, points.items, @max(max_reach, 1.0)) catch |err| switch (err) {
        error.GridTooLarge, error.NonFiniteCoordinate => null,
        error.OutOfMemory => return err,
    };
    defer if (cell_list) |*cl| cl.deinit();

    var candidates = std.ArrayListUnmanaged(u32).empty;
    defer candidates.deinit(allocator);

    // Points of one mover are contiguous, so `stamp` only needs the owner of
    // the point being visited.
    for (points.items, 0..) |p, pi| {
        const a = owner.items[pi];
        candidates.clearRetainingCapacity();
        if (cell_list) |cl| {
            try cl.neighborsInRadius(p, radius.items[pi] + max_radius + contact_margin, &candidates, allocator, points.items);
        } else {
            // Fallback when the grid would be too large: scan every point.
            for (0..points.items.len) |qi| try candidates.append(allocator, @intCast(qi));
        }
        for (candidates.items) |qi| {
            const b = owner.items[qi];
            if (b <= a or stamp[b] == a + 1) continue;
            const reach = radius.items[pi] + radius.items[qi] + contact_margin;
            const d = p.sub(points.items[qi]);
            if (d.dot(d) >= reach * reach) continue;
            stamp[b] = a + 1;
            try adj_lists[a].append(allocator, b);
            try adj_lists[b].append(allocator, a);
        }
    }

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

fn appendPoint(
    allocator: Allocator,
    points: *std.ArrayListUnmanaged(Vec3(f32)),
    owner: *std.ArrayListUnmanaged(u32),
    radius: *std.ArrayListUnmanaged(f32),
    pos: Vec3(f32),
    mover_idx: u32,
    r: f32,
) !void {
    if (mover_mod.isAbsentPos(pos)) return;
    try points.append(allocator, pos);
    try owner.append(allocator, mover_idx);
    try radius.append(allocator, r);
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
/// Contact margin used by the tests (the scorer's value).
const test_margin: f32 = 0.5;

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
    var graph = try buildInteractionGraph(allocator, &movers, &atoms, test_margin);
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
    var graph = try buildInteractionGraph(allocator, &movers, &atoms, test_margin);
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
    var graph = try buildInteractionGraph(allocator, &movers, &atoms, test_margin);
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

    // 3 movers 3 A apart on a line; radii 1.7, so the reach is 3.9 A.
    // A-B and B-C are in range, A-C (6 A) is not.
    var atoms = [_]Atom{
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } }, // mover A
        .{ .pos = .{ .x = 3.0, .y = 0.0, .z = 0.0 } }, // mover B
        .{ .pos = .{ .x = 6.0, .y = 0.0, .z = 0.0 } }, // mover C
    };

    var m0 = try makeMover(allocator, &.{0});
    defer m0.deinit();
    var m1 = try makeMover(allocator, &.{1});
    defer m1.deinit();
    var m2 = try makeMover(allocator, &.{2});
    defer m2.deinit();

    const movers = [_]Mover{ m0, m1, m2 };
    var graph = try buildInteractionGraph(allocator, &movers, &atoms, test_margin);
    defer graph.deinit();

    // Verify A-C do not directly interact (distance 6 > 3.9 reach)
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
    var graph = try buildInteractionGraph(allocator, &movers, &atoms, test_margin);
    defer graph.deinit();

    const cliques = try findCliques(allocator, &graph);
    defer {
        for (cliques) |c| allocator.free(c);
        allocator.free(cliques);
    }

    try testing.expectEqual(@as(usize, 2), cliques.len);
}

test "movers just beyond scoring range do not interact" {
    const allocator = testing.allocator;

    // Radii 1.7 + 1.7 + margin 0.5 = 3.9 A; 4.0 A apart is out of range.
    var atoms = [_]Atom{
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } },
        .{ .pos = .{ .x = 4.0, .y = 0.0, .z = 0.0 } },
    };
    var m0 = try makeMover(allocator, &.{0});
    defer m0.deinit();
    var m1 = try makeMover(allocator, &.{1});
    defer m1.deinit();

    const movers = [_]Mover{ m0, m1 };
    var graph = try buildInteractionGraph(allocator, &movers, &atoms, test_margin);
    defer graph.deinit();
    try testing.expectEqual(@as(usize, 0), graph.adjacency[0].len);
}

test "movers interact through a non-current orientation" {
    const allocator = testing.allocator;

    // Mover 0 sits 10 A from mover 1, but one of its orientations reaches
    // within 3 A of it.
    var atoms = [_]Atom{
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } },
        .{ .pos = .{ .x = 10.0, .y = 0.0, .z = 0.0 } },
    };
    var positions = [_]Vec3(f32){ .{ .x = 0.0, .y = 0.0, .z = 0.0 }, .{ .x = 7.0, .y = 0.0, .z = 0.0 } };
    var orientations = [_]mover_mod.Orientation{ .{ .positions = positions[0..1] }, .{ .positions = positions[1..2] } };
    const m0 = Mover{
        .kind = .single_h_rotator,
        .residue_idx = 0,
        .atom_indices = @constCast(&[_]u32{0}),
        .orientations = &orientations,
        .allocator = allocator,
    };
    var m1 = try makeMover(allocator, &.{1});
    defer m1.deinit();

    const movers = [_]Mover{ m0, m1 };
    var graph = try buildInteractionGraph(allocator, &movers, &atoms, test_margin);
    defer graph.deinit();
    try testing.expectEqual(@as(usize, 1), graph.adjacency[0].len);
}
