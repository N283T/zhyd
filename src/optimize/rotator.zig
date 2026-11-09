//! Rotator movers: OH/SH (single H), NH3+ and NH2 rotators.
//! Each creates a Mover with discrete orientations by rotating hydrogens
//! around a bond axis.

const std = @import("std");
const math_mod = @import("../math.zig");
const model_mod = @import("../model.zig");
const mover_mod = @import("mover.zig");

const Vec3 = math_mod.Vec3;
const Atom = model_mod.Atom;
const Mover = mover_mod.Mover;
const Orientation = mover_mod.Orientation;

/// Create a SingleH rotator mover (OH, SH, SeH).
/// h_idx: index of the hydrogen atom to rotate
/// center_idx: index of the atom the H is bonded to (e.g., O in O-H)
/// axis_idx: index of the atom defining the rotation axis (e.g., C in C-O-H)
/// Produces 12 orientations at 30° intervals. No penalty for any orientation.
pub fn createSingleHRotator(
    allocator: std.mem.Allocator,
    atoms: []const Atom,
    h_idx: u32,
    center_idx: u32,
    axis_idx: u32,
    residue_idx: u32,
) !Mover {
    const n_orientations = 12;
    const n_atoms_per = 1;

    const center = atoms[center_idx].pos;
    const axis_atom = atoms[axis_idx].pos;
    const axis_dir = center.sub(axis_atom);
    const h_pos = atoms[h_idx].pos;

    // Allocate one contiguous block for all orientation positions.
    const all_positions = try allocator.alloc(Vec3(f32), n_orientations * n_atoms_per);
    errdefer allocator.free(all_positions);

    const orientations = try allocator.alloc(Orientation, n_orientations);
    errdefer allocator.free(orientations);

    for (0..n_orientations) |i| {
        const angle: f32 = @as(f32, @floatFromInt(i)) * 30.0;
        const rotated = math_mod.rotateAroundAxis(f32, h_pos, center, axis_dir, angle);
        const slice = all_positions[i * n_atoms_per .. (i + 1) * n_atoms_per];
        slice[0] = rotated;
        orientations[i] = .{ .positions = slice, .penalty = 0.0 };
    }

    const atom_indices = try allocator.alloc(u32, 1);
    errdefer allocator.free(atom_indices);
    atom_indices[0] = h_idx;

    return Mover{
        .kind = .single_h_rotator,
        .residue_idx = residue_idx,
        .atom_indices = atom_indices,
        .orientations = orientations,
        .allocator = allocator,
        .center_idx = center_idx,
        .axis_idx = axis_idx,
        .positions_backing = all_positions,
    };
}

/// Coarse NH3+ rotation offsets in degrees; see `createNH3Rotator`.
const nh3_angle_offsets = [4]f32{ 0.0, 30.0, 60.0, 90.0 };

/// Create an NH3+ rotator mover.
/// h_indices: indices of the 3 hydrogen atoms
/// center_idx: the N atom index
/// axis_idx: the atom defining rotation axis (CA usually)
/// NH3+ is 3-fold symmetric, so rotations repeat every 120°. Produces 4
/// orientations at 0°, 30°, 60° and 90° offsets, the distinct positions of a
/// 30° grid. (An earlier 0°/+60°/-60° set had two identical orientations.)
/// Non-default orientations have a small penalty of 0.05.
pub fn createNH3Rotator(
    allocator: std.mem.Allocator,
    atoms: []const Atom,
    h_indices: [3]u32,
    center_idx: u32,
    axis_idx: u32,
    residue_idx: u32,
) !Mover {
    const n_orientations = nh3_angle_offsets.len;
    const n_atoms_per = 3;

    const center = atoms[center_idx].pos;
    const axis_atom = atoms[axis_idx].pos;
    const axis_dir = center.sub(axis_atom);

    const angle_offsets = nh3_angle_offsets;
    const penalties = [nh3_angle_offsets.len]f32{ 0.0, 0.05, 0.05, 0.05 };

    // Allocate one contiguous block for all orientation positions.
    const all_positions = try allocator.alloc(Vec3(f32), n_orientations * n_atoms_per);
    errdefer allocator.free(all_positions);

    const orientations = try allocator.alloc(Orientation, n_orientations);
    errdefer allocator.free(orientations);

    for (0..n_orientations) |i| {
        const slice = all_positions[i * n_atoms_per .. (i + 1) * n_atoms_per];
        for (0..n_atoms_per) |j| {
            const h_pos = atoms[h_indices[j]].pos;
            slice[j] = math_mod.rotateAroundAxis(f32, h_pos, center, axis_dir, angle_offsets[i]);
        }
        orientations[i] = .{ .positions = slice, .penalty = penalties[i] };
    }

    const atom_indices = try allocator.alloc(u32, 3);
    errdefer allocator.free(atom_indices);
    for (0..3) |i| {
        atom_indices[i] = h_indices[i];
    }

    return Mover{
        .kind = .nh3_rotator,
        .residue_idx = residue_idx,
        .atom_indices = atom_indices,
        .orientations = orientations,
        .allocator = allocator,
        .center_idx = center_idx,
        .axis_idx = axis_idx,
        .positions_backing = all_positions,
    };
}

/// Create an NH2 rotator mover for a neutral primary amine or a PRO NH2+
/// secondary ammonium. Rotates a rigid pair of H atoms together around the
/// N-CA bond so the intrinsic HNH angle (set by the placer) is preserved.
/// h_indices: indices of the 2 hydrogen atoms
/// center_idx: the N atom index
/// axis_idx: the atom defining the rotation axis (CA)
/// Produces 12 orientations at 30° intervals. No penalty for any orientation.
/// Unlike NH3+, there is no 3-fold symmetry to exploit — every angle yields
/// a distinct configuration of the (H, H, lone pair) tetrahedron.
pub fn createNH2Rotator(
    allocator: std.mem.Allocator,
    atoms: []const Atom,
    h_indices: [2]u32,
    center_idx: u32,
    axis_idx: u32,
    residue_idx: u32,
) !Mover {
    const n_orientations = 12;
    const n_atoms_per = 2;

    const center = atoms[center_idx].pos;
    const axis_atom = atoms[axis_idx].pos;
    const axis_dir = center.sub(axis_atom);

    // Allocate one contiguous block for all orientation positions.
    const all_positions = try allocator.alloc(Vec3(f32), n_orientations * n_atoms_per);
    errdefer allocator.free(all_positions);

    const orientations = try allocator.alloc(Orientation, n_orientations);
    errdefer allocator.free(orientations);

    for (0..n_orientations) |i| {
        const angle: f32 = @as(f32, @floatFromInt(i)) * 30.0;
        const slice = all_positions[i * n_atoms_per .. (i + 1) * n_atoms_per];
        for (0..n_atoms_per) |j| {
            const h_pos = atoms[h_indices[j]].pos;
            slice[j] = math_mod.rotateAroundAxis(f32, h_pos, center, axis_dir, angle);
        }
        orientations[i] = .{ .positions = slice, .penalty = 0.0 };
    }

    const atom_indices = try allocator.alloc(u32, 2);
    errdefer allocator.free(atom_indices);
    for (0..2) |i| {
        atom_indices[i] = h_indices[i];
    }

    return Mover{
        .kind = .nh2_rotator,
        .residue_idx = residue_idx,
        .atom_indices = atom_indices,
        .orientations = orientations,
        .allocator = allocator,
        .center_idx = center_idx,
        .axis_idx = axis_idx,
        .positions_backing = all_positions,
    };
}

/// Generate fine angular orientations around a coarse best orientation.
/// Returns empty slice if mover has no center_idx/axis_idx (non-rotator).
pub fn generateFineOrientations(
    allocator: std.mem.Allocator,
    atoms: []const Atom,
    m: *const Mover,
    coarse_best_idx: u16,
) ![]Orientation {
    const ci = m.center_idx orelse return try allocator.alloc(Orientation, 0);
    const ai = m.axis_idx orelse return try allocator.alloc(Orientation, 0);

    const center_pos = atoms[ci].pos;
    const axis_dir = center_pos.sub(atoms[ai].pos);

    // Determine base angle from coarse_best_idx
    const base_angle: f32 = switch (m.kind) {
        .single_h_rotator, .nh2_rotator => @as(f32, @floatFromInt(coarse_best_idx)) * 30.0,
        .nh3_rotator => nh3_angle_offsets[coarse_best_idx],
        else => return try allocator.alloc(Orientation, 0),
    };

    // All rotators use a 30° coarse grid, so +-15° in 5° steps covers the gap
    // between neighbouring coarse orientations.
    const half_range: f32 = 15.0;
    const step: f32 = 5.0;

    // Collect fine offsets, skipping 0 (the coarse position)
    var offsets_buf: [12]f32 = undefined;
    var n_offsets: usize = 0;
    var off: f32 = -half_range;
    while (off <= half_range + 0.001) : (off += step) {
        if (@abs(off) < 0.001) continue; // skip 0 offset
        offsets_buf[n_offsets] = off;
        n_offsets += 1;
    }

    const n_atoms = m.orientations[0].positions.len;
    const original_positions = m.orientations[0].positions;

    const result = try allocator.alloc(Orientation, n_offsets);
    var allocated: usize = 0;
    errdefer {
        for (result[0..allocated]) |o| allocator.free(o.positions);
        allocator.free(result);
    }

    for (0..n_offsets) |i| {
        const angle = base_angle + offsets_buf[i];
        const positions = try allocator.alloc(Vec3(f32), n_atoms);
        for (0..n_atoms) |j| {
            positions[j] = math_mod.rotateAroundAxis(f32, original_positions[j], center_pos, axis_dir, angle);
        }
        result[i] = .{ .positions = positions, .penalty = 0.0 };
        allocated += 1;
    }

    return result;
}

// ── Tests ─────────────────────────────────────────────────────────────────────

const testing = std.testing;

test "single H rotator has 12 orientations" {
    const allocator = testing.allocator;

    // Geometry: axis at (0,0,0), center (O) at (1,0,0), H at (2,0,0)
    const center_pos = Vec3(f32){ .x = 1.0, .y = 0.0, .z = 0.0 };
    const atoms = [_]Atom{
        .{ .pos = .{ .x = 2.0, .y = 0.0, .z = 0.0 } }, // H at index 0
        .{ .pos = center_pos }, // center (O) at index 1
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } }, // axis atom at index 2
    };

    var mover = try createSingleHRotator(allocator, &atoms, 0, 1, 2, 0);
    defer mover.deinit();

    try testing.expectEqual(@as(usize, 12), mover.orientations.len);
    for (mover.orientations) |o| {
        const dist = o.positions[0].distance(center_pos);
        try testing.expectApproxEqAbs(@as(f32, 1.0), dist, 0.01);
    }
}

test "single H rotator has no penalties" {
    const allocator = testing.allocator;

    const atoms = [_]Atom{
        .{ .pos = .{ .x = 2.0, .y = 0.0, .z = 0.0 } },
        .{ .pos = .{ .x = 1.0, .y = 0.0, .z = 0.0 } },
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } },
    };

    var mover = try createSingleHRotator(allocator, &atoms, 0, 1, 2, 0);
    defer mover.deinit();

    for (mover.orientations) |o| {
        try testing.expectApproxEqAbs(@as(f32, 0.0), o.penalty, 0.001);
    }
}

test "NH3 rotator has 4 orientations with 3 atoms each" {
    const allocator = testing.allocator;

    // N at origin, CA at (-1,0,0), 3 H atoms placed around N
    const atoms = [_]Atom{
        .{ .pos = .{ .x = 0.0, .y = 1.0, .z = 0.0 } }, // H0 at index 0
        .{ .pos = .{ .x = 0.866, .y = -0.5, .z = 0.0 } }, // H1 at index 1
        .{ .pos = .{ .x = -0.866, .y = -0.5, .z = 0.0 } }, // H2 at index 2
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } }, // N at index 3
        .{ .pos = .{ .x = -1.0, .y = 0.0, .z = 0.0 } }, // CA at index 4
    };

    var mover = try createNH3Rotator(allocator, &atoms, .{ 0, 1, 2 }, 3, 4, 0);
    defer mover.deinit();

    try testing.expectEqual(@as(usize, 4), mover.orientations.len);
    for (mover.orientations) |o| try testing.expectEqual(@as(usize, 3), o.positions.len);
}

test "NH3 rotator orientations are distinct under 3-fold symmetry" {
    const allocator = testing.allocator;

    // Axis along x, H in the yz plane at 120° spacing.
    const atoms = [_]Atom{
        .{ .pos = .{ .x = 0.0, .y = 1.0, .z = 0.0 } },
        .{ .pos = .{ .x = 0.0, .y = -0.5, .z = 0.866 } },
        .{ .pos = .{ .x = 0.0, .y = -0.5, .z = -0.866 } },
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } },
        .{ .pos = .{ .x = -1.0, .y = 0.0, .z = 0.0 } },
    };

    var mover = try createNH3Rotator(allocator, &atoms, .{ 0, 1, 2 }, 3, 4, 0);
    defer mover.deinit();

    // Two orientations are the same configuration when every H of one sits on
    // some H of the other.
    for (mover.orientations, 0..) |oa, a| {
        for (mover.orientations[a + 1 ..]) |ob| {
            var all_matched = true;
            for (oa.positions) |pa| {
                var matched = false;
                for (ob.positions) |pb| {
                    if (pa.distance(pb) < 0.01) matched = true;
                }
                if (!matched) all_matched = false;
            }
            try testing.expect(!all_matched);
        }
    }
}

test "NH3 rotator first orientation has zero penalty, others have 0.05" {
    const allocator = testing.allocator;

    const atoms = [_]Atom{
        .{ .pos = .{ .x = 0.0, .y = 1.0, .z = 0.0 } },
        .{ .pos = .{ .x = 0.866, .y = -0.5, .z = 0.0 } },
        .{ .pos = .{ .x = -0.866, .y = -0.5, .z = 0.0 } },
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } },
        .{ .pos = .{ .x = -1.0, .y = 0.0, .z = 0.0 } },
    };

    var mover = try createNH3Rotator(allocator, &atoms, .{ 0, 1, 2 }, 3, 4, 0);
    defer mover.deinit();

    try testing.expectApproxEqAbs(@as(f32, 0.0), mover.orientations[0].penalty, 0.001);
    for (mover.orientations[1..]) |o| try testing.expectApproxEqAbs(@as(f32, 0.05), o.penalty, 0.001);
}

test "rotator atom distances preserved after rotation" {
    const allocator = testing.allocator;

    // H at (2,0,0), center at (1,0,0), axis at (0,0,0)
    // H is 1.0 away from center; this should be preserved for all orientations
    const center_pos = Vec3(f32){ .x = 1.0, .y = 0.0, .z = 0.0 };
    const atoms = [_]Atom{
        .{ .pos = .{ .x = 2.0, .y = 0.0, .z = 0.0 } },
        .{ .pos = center_pos },
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } },
    };

    var mover = try createSingleHRotator(allocator, &atoms, 0, 1, 2, 5);
    defer mover.deinit();

    try testing.expectEqual(@as(u32, 5), mover.residue_idx);
    for (mover.orientations) |o| {
        const d = o.positions[0].distance(center_pos);
        try testing.expectApproxEqAbs(@as(f32, 1.0), d, 0.01);
    }
}

test "generateFineOrientations produces 6 fine positions for single H" {
    const allocator = testing.allocator;
    const atoms = [_]Atom{
        .{ .pos = .{ .x = 2.0, .y = 0.0, .z = 0.0 } },
        .{ .pos = .{ .x = 1.0, .y = 0.0, .z = 0.0 } },
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } },
    };

    var mover = try createSingleHRotator(allocator, &atoms, 0, 1, 2, 0);
    defer mover.deinit();

    const fine = try generateFineOrientations(allocator, &atoms, &mover, 0);
    defer {
        for (fine) |o| allocator.free(o.positions);
        allocator.free(fine);
    }

    try testing.expectEqual(@as(usize, 6), fine.len);
    // Distance to center preserved
    const center_pos = atoms[1].pos;
    for (fine) |o| {
        const dist = o.positions[0].distance(center_pos);
        try testing.expectApproxEqAbs(@as(f32, 1.0), dist, 0.01);
    }
}

test "generateFineOrientations produces 6 fine positions for NH3" {
    const allocator = testing.allocator;
    const atoms = [_]Atom{
        .{ .pos = .{ .x = 0.0, .y = 1.0, .z = 0.0 } },
        .{ .pos = .{ .x = 0.866, .y = -0.5, .z = 0.0 } },
        .{ .pos = .{ .x = -0.866, .y = -0.5, .z = 0.0 } },
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } },
        .{ .pos = .{ .x = -1.0, .y = 0.0, .z = 0.0 } },
    };

    var mover = try createNH3Rotator(allocator, &atoms, .{ 0, 1, 2 }, 3, 4, 0);
    defer mover.deinit();

    const fine = try generateFineOrientations(allocator, &atoms, &mover, 0);
    defer {
        for (fine) |o| allocator.free(o.positions);
        allocator.free(fine);
    }

    try testing.expectEqual(@as(usize, 6), fine.len);
    try testing.expectEqual(@as(usize, 3), fine[0].positions.len);
}
