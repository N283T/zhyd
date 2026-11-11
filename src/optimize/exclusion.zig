//! Bonded exclusions for mover scoring.
//!
//! Atoms within three covalent bonds of a mover atom (1-2, 1-3 and 1-4
//! partners) are left out of its score, as reduce and probe do. Without this a
//! rotating H is scored against its own parent and the atoms around it: the
//! 1-2 and 1-3 terms are constant, but the 1-4 terms change with the torsion
//! and push OH and NH3+ groups away from orientations that eclipse a heavy
//! atom, whatever H-bonds those orientations make.
//!
//! The model has no full bond table, so bonds are inferred from distance:
//! two atoms are bonded when they are closer than the sum of their covalent
//! radii plus `bond_tolerance`. H-H pairs and metals are never bonded.

const std = @import("std");
const Allocator = std.mem.Allocator;
const math_mod = @import("../math.zig");
const Vec3 = math_mod.Vec3;
const model_mod = @import("../model.zig");
const Atom = model_mod.Atom;
const CellList = model_mod.CellList;
const mover_mod = @import("mover.zig");
const Mover = mover_mod.Mover;

/// Slack added to the covalent radius sum when inferring a bond, in Angstrom.
pub const bond_tolerance: f32 = 0.4;
/// Pairs closer than this are overlapping atoms (e.g. altloc copies), not bonds.
const min_bond_length: f32 = 0.4;
/// Bond paths up to this length are excluded.
const max_bond_path: u32 = 3;
/// Upper bound on any inferred bond length, used as the search radius.
const max_bond_length: f32 = 2.6;

pub const Exclusions = struct {
    /// CSR over atom indices: atoms excluded from scoring atom `i` are
    /// `list[offsets[i]..offsets[i + 1]]`, sorted. Empty for non-mover atoms.
    offsets: []u32,
    list: []u32,

    pub fn of(self: *const Exclusions, atom_idx: u32) []const u32 {
        return self.list[self.offsets[atom_idx]..self.offsets[atom_idx + 1]];
    }

    pub fn excludes(self: *const Exclusions, atom_idx: u32, other_idx: u32) bool {
        return std.sort.binarySearch(u32, self.of(atom_idx), other_idx, struct {
            fn order(key: u32, item: u32) std.math.Order {
                return std.math.order(key, item);
            }
        }.order) != null;
    }

    pub fn deinit(self: *Exclusions, allocator: Allocator) void {
        allocator.free(self.offsets);
        allocator.free(self.list);
    }
};

fn canBond(a: Atom) bool {
    if (mover_mod.isAbsentPos(a.pos)) return false;
    return !a.element_type.info().flags.metallic;
}

fn bonded(a: Atom, b: Atom) bool {
    if (a.is_hydrogen and b.is_hydrogen) return false;
    const d = a.pos.distance(b.pos);
    if (d < min_bond_length) return false;
    const cut = a.element_type.info().covalent_radius + b.element_type.info().covalent_radius + bond_tolerance;
    return d < cut;
}

/// Build the exclusion lists for every atom controlled by a mover, from the
/// current atom positions.
pub fn build(allocator: Allocator, movers: []const Mover, atoms: []const Atom) !Exclusions {
    const n = atoms.len;

    // Atoms that take part in the bond graph, and their positions for the grid.
    var graph_atoms = std.ArrayListUnmanaged(u32).empty;
    defer graph_atoms.deinit(allocator);
    var graph_pos = std.ArrayListUnmanaged(Vec3(f32)).empty;
    defer graph_pos.deinit(allocator);
    for (atoms, 0..) |a, i| {
        if (!canBond(a)) continue;
        try graph_atoms.append(allocator, @intCast(i));
        try graph_pos.append(allocator, a.pos);
    }

    // Bonded neighbors (CSR over atom indices).
    var pairs = std.ArrayListUnmanaged([2]u32).empty;
    defer pairs.deinit(allocator);
    var cell_list: ?CellList = CellList.init(allocator, graph_pos.items, max_bond_length) catch |err| switch (err) {
        error.GridTooLarge, error.NonFiniteCoordinate => null,
        error.OutOfMemory => return err,
    };
    defer if (cell_list) |*cl| cl.deinit();
    var found = std.ArrayListUnmanaged(u32).empty;
    defer found.deinit(allocator);
    for (graph_atoms.items, 0..) |ai, gi| {
        found.clearRetainingCapacity();
        if (cell_list) |cl| {
            try cl.neighborsInRadius(graph_pos.items[gi], max_bond_length, &found, allocator, graph_pos.items);
        } else {
            for (0..graph_atoms.items.len) |gj| try found.append(allocator, @intCast(gj));
        }
        for (found.items) |gj| {
            const bj = graph_atoms.items[gj];
            if (bj <= ai) continue;
            if (bonded(atoms[ai], atoms[bj])) try pairs.append(allocator, .{ ai, bj });
        }
    }
    const nb_offsets = try allocator.alloc(u32, n + 1);
    defer allocator.free(nb_offsets);
    @memset(nb_offsets, 0);
    for (pairs.items) |p| {
        nb_offsets[p[0] + 1] += 1;
        nb_offsets[p[1] + 1] += 1;
    }
    for (0..n) |i| nb_offsets[i + 1] += nb_offsets[i];
    const nb = try allocator.alloc(u32, nb_offsets[n]);
    defer allocator.free(nb);
    {
        const fill = try allocator.alloc(u32, n);
        defer allocator.free(fill);
        @memcpy(fill, nb_offsets[0..n]);
        for (pairs.items) |p| {
            nb[fill[p[0]]] = p[1];
            fill[p[0]] += 1;
            nb[fill[p[1]]] = p[0];
            fill[p[1]] += 1;
        }
    }

    // Breadth-first search up to `max_bond_path` bonds from each mover atom.
    const is_mover_atom = try allocator.alloc(bool, n);
    defer allocator.free(is_mover_atom);
    @memset(is_mover_atom, false);
    for (movers) |m| {
        for (m.atom_indices) |ai| is_mover_atom[ai] = true;
    }

    const offsets = try allocator.alloc(u32, n + 1);
    errdefer allocator.free(offsets);
    var list = std.ArrayListUnmanaged(u32).empty;
    errdefer list.deinit(allocator);

    var frontier = std.ArrayListUnmanaged(u32).empty;
    defer frontier.deinit(allocator);
    var next = std.ArrayListUnmanaged(u32).empty;
    defer next.deinit(allocator);

    offsets[0] = 0;
    for (0..n) |i| {
        if (is_mover_atom[i]) {
            const start = list.items.len;
            frontier.clearRetainingCapacity();
            try frontier.append(allocator, @intCast(i));
            for (0..max_bond_path) |_| {
                next.clearRetainingCapacity();
                for (frontier.items) |a| {
                    for (nb[nb_offsets[a]..nb_offsets[a + 1]]) |b| {
                        if (b == i) continue;
                        if (std.mem.indexOfScalar(u32, list.items[start..], b) != null) continue;
                        try list.append(allocator, b);
                        try next.append(allocator, b);
                    }
                }
                std.mem.swap(std.ArrayListUnmanaged(u32), &frontier, &next);
            }
            std.mem.sort(u32, list.items[start..], {}, std.sort.asc(u32));
        }
        offsets[i + 1] = @intCast(list.items.len);
    }

    return .{ .offsets = offsets, .list = try list.toOwnedSlice(allocator) };
}

// ── Tests ─────────────────────────────────────────────────────────────────────

const testing = std.testing;

test "excludes atoms within three bonds of a mover atom" {
    // Chain CA-CB-OG-HG plus an unbonded O 2.5 A from HG. HG is the mover atom.
    const atoms = [_]Atom{
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 }, .element_type = .C }, // 0 CA
        .{ .pos = .{ .x = 1.53, .y = 0.0, .z = 0.0 }, .element_type = .C }, // 1 CB
        .{ .pos = .{ .x = 2.0, .y = 1.36, .z = 0.0 }, .element_type = .O }, // 2 OG
        .{ .pos = .{ .x = 2.95, .y = 1.5, .z = 0.0 }, .element_type = .Hpol, .is_hydrogen = true }, // 3 HG
        .{ .pos = .{ .x = 5.4, .y = 1.5, .z = 0.0 }, .element_type = .O }, // 4 acceptor
        .{ .pos = .{ .x = -0.5, .y = -1.4, .z = 0.0 }, .element_type = .N }, // 5 N, 4 bonds from HG
    };
    var positions = [_]Vec3(f32){atoms[3].pos};
    var orientations = [_]mover_mod.Orientation{.{ .positions = &positions }};
    const movers = [_]Mover{.{
        .kind = .single_h_rotator,
        .residue_idx = 0,
        .atom_indices = @constCast(&[_]u32{3}),
        .orientations = &orientations,
        .allocator = testing.allocator,
    }};

    var ex = try build(testing.allocator, &movers, &atoms);
    defer ex.deinit(testing.allocator);

    try testing.expectEqualSlices(u32, &.{ 0, 1, 2 }, ex.of(3));
    try testing.expect(ex.excludes(3, 0));
    try testing.expect(!ex.excludes(3, 4));
    try testing.expect(!ex.excludes(3, 5));
    try testing.expectEqual(@as(usize, 0), ex.of(0).len);
}
