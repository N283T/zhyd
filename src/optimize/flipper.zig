//! Flip movers for Asn/Gln (amide flip) and His (histidine ring flip).

const std = @import("std");
const math_mod = @import("../math.zig");
const model_mod = @import("../model.zig");
const mover_mod = @import("mover.zig");

const Vec3 = math_mod.Vec3;
const Atom = model_mod.Atom;
const Mover = mover_mod.Mover;
const Orientation = mover_mod.Orientation;
const element = @import("../element.zig");
const AtomFlags = element.AtomFlags;
const bond_policy = @import("../place/bond_policy.zig");
const geometry = @import("../place/geometry.zig");

/// Compute planar NH2 hydrogen positions given:
///   n_pos: the nitrogen position
///   c_pos: the carbon bonded to N (e.g. CG or CD)
///   o_pos: the oxygen bonded to C (for reference plane)
/// Returns [2]Vec3(f32): the two H positions on the amide nitrogen.
fn computeAmideNH2(n_pos: Vec3(f32), c_pos: Vec3(f32), o_pos: Vec3(f32), mode: bond_policy.BondLengthMode) [2]Vec3(f32) {
    const bond_len = bond_policy.adjustedBondLength(mode, 1.00, .N, .Hpol);
    const half_angle_deg: f32 = 60.0; // half of 120°

    // Vector from C to N
    const cn = n_pos.sub(c_pos).normalize();
    // Vector from C to O (to define the plane)
    const co = o_pos.sub(c_pos).normalize();

    // Normal to the C-N-O plane
    const normal = cn.cross(co).normalize();

    // The two H atoms are placed symmetrically around the C-N direction
    // in the plane of the amide group, at bond_len from N.
    // We rotate cn by +half_angle and -half_angle around the normal.
    const h1_dir = math_mod.rotateAroundAxis(f32, n_pos.add(cn.scale(bond_len)), n_pos, normal, half_angle_deg);
    const h2_dir = math_mod.rotateAroundAxis(f32, n_pos.add(cn.scale(bond_len)), n_pos, normal, -half_angle_deg);

    return .{ h1_dir, h2_dir };
}

/// Create an amide flip mover for Asn or Gln.
/// o_idx: index of the O atom (OD1 or OE1)
/// n_idx: index of the N atom (ND2 or NE2)
/// h1_idx, h2_idx: indices of the two H atoms on N
/// c_idx: index of the C atom (CG or CD) for H recalculation
pub fn createAmideFlipper(
    allocator: std.mem.Allocator,
    atoms: []const Atom,
    o_idx: u32,
    n_idx: u32,
    h1_idx: u32,
    h2_idx: u32,
    c_idx: u32,
    residue_idx: u32,
    mode: bond_policy.BondLengthMode,
) !Mover {
    const n_orientations = 2;
    const n_atoms_per = 4;

    const o_pos = atoms[o_idx].pos;
    const n_pos = atoms[n_idx].pos;
    const h1_pos = atoms[h1_idx].pos;
    const h2_pos = atoms[h2_idx].pos;
    const c_pos = atoms[c_idx].pos;

    // Orientation 1: swapped O<->N, recompute H positions
    // After swap: new N is at old O position, new O is at old N position
    const new_n_pos = o_pos; // N moves to where O was
    const new_o_pos = n_pos; // O moves to where N was
    const new_hs = computeAmideNH2(new_n_pos, c_pos, new_o_pos, mode);

    // Allocate one contiguous block for all orientation positions.
    const all_positions = try allocator.alloc(Vec3(f32), n_orientations * n_atoms_per);
    errdefer allocator.free(all_positions);

    const orientations = try allocator.alloc(Orientation, n_orientations);
    errdefer allocator.free(orientations);

    // Orientation 0: original positions (O, N, H1, H2)
    const pos0 = all_positions[0..n_atoms_per];
    pos0[0] = o_pos;
    pos0[1] = n_pos;
    pos0[2] = h1_pos;
    pos0[3] = h2_pos;
    orientations[0] = .{ .positions = pos0, .penalty = 0.0 };

    // Orientation 1: swapped positions
    const pos1 = all_positions[n_atoms_per .. 2 * n_atoms_per];
    pos1[0] = new_o_pos;
    pos1[1] = new_n_pos;
    pos1[2] = new_hs[0];
    pos1[3] = new_hs[1];
    orientations[1] = .{ .positions = pos1, .penalty = 0.5 };

    const atom_indices = try allocator.alloc(u32, 4);
    errdefer allocator.free(atom_indices);
    atom_indices[0] = o_idx;
    atom_indices[1] = n_idx;
    atom_indices[2] = h1_idx;
    atom_indices[3] = h2_idx;

    return Mover{
        .kind = .amide_flip,
        .residue_idx = residue_idx,
        .atom_indices = atom_indices,
        .orientations = orientations,
        .allocator = allocator,
        .positions_backing = all_positions,
    };
}

const ABSENT_H_POS = Mover.ABSENT_H_POS;

/// Compute position of a ring H (HD1 on ND1 or HE2 on NE2) given:
///   n_pos: the nitrogen bearing the H
///   bonded1, bonded2: the two atoms bonded to N in the ring
fn computeRingNH(n_pos: Vec3(f32), bonded1: Vec3(f32), bonded2: Vec3(f32), mode: bond_policy.BondLengthMode) Vec3(f32) {
    const bond_len = bond_policy.adjustedBondLength(mode, 1.02, .N, .Hpol);
    // Place H in the plane opposite to the bisector of the two ring bonds
    const v1 = bonded1.sub(n_pos).normalize();
    const v2 = bonded2.sub(n_pos).normalize();
    const bisector = v1.add(v2).normalize();
    return n_pos.add(bisector.scale(-bond_len));
}

/// Create a His flip mover.
/// Neutral His: 4 orientations (2 flip states × HID/HIE). `charged` (HIP,
/// from a pKa-derived state): 2 orientations (2 flip states, both H).
///
/// Atom order in atom_indices / positions: [ND1, CD2, CE1, NE2, HD1, HE2]
/// hd1_idx / he2_idx may be null if not yet present in the model;
/// if null, a placeholder index (same as nd1_idx) is used and the
/// position is set to ABSENT_H_POS when H is absent.
pub fn createHisFlipper(
    allocator: std.mem.Allocator,
    atoms: []const Atom,
    nd1_idx: u32,
    cd2_idx: u32,
    ce1_idx: u32,
    ne2_idx: u32,
    hd1_idx: ?u32,
    he2_idx: ?u32,
    residue_idx: u32,
    mode: bond_policy.BondLengthMode,
    charged: bool,
) !Mover {
    const n_orientations: usize = if (charged) 2 else 4;
    const n_atoms_per = 6;

    // --- Original ring positions ---
    const nd1 = atoms[nd1_idx].pos;
    const cd2 = atoms[cd2_idx].pos;
    const ce1 = atoms[ce1_idx].pos;
    const ne2 = atoms[ne2_idx].pos;

    // HD1 and HE2 original positions (if atoms exist)
    const hd1_orig: Vec3(f32) = if (hd1_idx) |idx| atoms[idx].pos else ABSENT_H_POS;
    const he2_orig: Vec3(f32) = if (he2_idx) |idx| atoms[idx].pos else ABSENT_H_POS;

    // Compute H positions based on ring geometry (for orientations that need them)
    // HD1 sits on ND1, bonded to CE1 and CG (we don't have CG here; use CE1 and CD2 as the two ring neighbours of ND1)
    const hd1_computed = computeRingNH(nd1, ce1, cd2, mode);
    // HE2 sits on NE2, bonded to CE1 and CD2
    const he2_computed = computeRingNH(ne2, ce1, cd2, mode);

    // --- Flipped ring positions ---
    // Swap ND1<->CD2 and CE1<->NE2
    const nd1_f = cd2; // ND1 goes to CD2 position
    const cd2_f = nd1; // CD2 goes to ND1 position
    const ce1_f = ne2; // CE1 goes to NE2 position
    const ne2_f = ce1; // NE2 goes to CE1 position

    // Flipped H positions
    const hd1_f_computed = computeRingNH(nd1_f, ce1_f, cd2_f, mode);
    const he2_f_computed = computeRingNH(ne2_f, ce1_f, cd2_f, mode);

    // --- Build 4 orientations ---
    // At neutral pH, HIS is predominantly neutral (HID or HIE, pKa ~6.0).
    // HIP (doubly protonated, +1 charge) is excluded: penalty-based approaches
    // fail because the dot-sphere scoring rewards additional H-bond donors
    // without accounting for desolvation/charge costs, causing nearly all HIS
    // to be assigned HIP regardless of penalty magnitude.
    //
    // Orient 0: no flip, HIE (NE2-HE2 protonated, ND1 acceptor)
    // Orient 1: no flip, HID (ND1-HD1 protonated, NE2 acceptor)
    // Orient 2: flip, HIE
    // Orient 3: flip, HID
    const neutral_penalties = [4]f32{ 0.00, 0.00, 0.50, 0.50 };

    // For each orientation: positions = [ND1, CD2, CE1, NE2, HD1, HE2]
    const OrientSpec = struct {
        nd1_p: Vec3(f32),
        cd2_p: Vec3(f32),
        ce1_p: Vec3(f32),
        ne2_p: Vec3(f32),
        hd1_p: Vec3(f32),
        he2_p: Vec3(f32),
    };

    // Prefer original positions for H if they exist in the original model
    const hd1_no_flip = if (hd1_idx != null) hd1_orig else hd1_computed;
    const he2_no_flip = if (he2_idx != null) he2_orig else he2_computed;

    const neutral_specs = [4]OrientSpec{
        // 0: no flip, HIE (HE2 only)
        .{ .nd1_p = nd1, .cd2_p = cd2, .ce1_p = ce1, .ne2_p = ne2, .hd1_p = ABSENT_H_POS, .he2_p = he2_no_flip },
        // 1: no flip, HID (HD1 only)
        .{ .nd1_p = nd1, .cd2_p = cd2, .ce1_p = ce1, .ne2_p = ne2, .hd1_p = hd1_no_flip, .he2_p = ABSENT_H_POS },
        // 2: flip, HIE
        .{ .nd1_p = nd1_f, .cd2_p = cd2_f, .ce1_p = ce1_f, .ne2_p = ne2_f, .hd1_p = ABSENT_H_POS, .he2_p = he2_f_computed },
        // 3: flip, HID
        .{ .nd1_p = nd1_f, .cd2_p = cd2_f, .ce1_p = ce1_f, .ne2_p = ne2_f, .hd1_p = hd1_f_computed, .he2_p = ABSENT_H_POS },
    };

    // Flags per orientation: [ND1, CD2, CE1, NE2, HD1, HE2]
    // Aromatic ring carbons keep ARA (aromatic+acceptor) in all states.
    // Nitrogens switch donor/acceptor based on protonation.
    // H atoms are always donor when present.
    const ar_acc = AtomFlags{ .aromatic = true, .acceptor = true }; // ring C
    const ar_don = AtomFlags{ .aromatic = true, .donor = true }; // protonated N
    const ar_acc_n = AtomFlags{ .aromatic = true, .acceptor = true }; // unprotonated N
    const h_don = AtomFlags{ .donor = true }; // H atom present
    const h_absent = AtomFlags{}; // absent H (no flags)

    // [ND1, CD2, CE1, NE2, HD1, HE2]
    const neutral_flags = [4][6]AtomFlags{
        // 0: HIE (HE2 only) → ND1=acceptor, NE2=donor
        .{ ar_acc_n, ar_acc, ar_acc, ar_don, h_absent, h_don },
        // 1: HID (HD1 only) → ND1=donor, NE2=acceptor
        .{ ar_don, ar_acc, ar_acc, ar_acc_n, h_don, h_absent },
        // 2: flip, HIE → ND1=acceptor, NE2=donor (flipped positions)
        .{ ar_acc_n, ar_acc, ar_acc, ar_don, h_absent, h_don },
        // 3: flip, HID → ND1=donor, NE2=acceptor
        .{ ar_don, ar_acc, ar_acc, ar_acc_n, h_don, h_absent },
    };

    // HIP: both ring N protonated and positive, both H present.
    const ar_pos = AtomFlags{ .aromatic = true, .donor = true, .positive = true };
    const charged_specs = [2]OrientSpec{
        .{ .nd1_p = nd1, .cd2_p = cd2, .ce1_p = ce1, .ne2_p = ne2, .hd1_p = hd1_no_flip, .he2_p = he2_no_flip },
        .{ .nd1_p = nd1_f, .cd2_p = cd2_f, .ce1_p = ce1_f, .ne2_p = ne2_f, .hd1_p = hd1_f_computed, .he2_p = he2_f_computed },
    };
    const charged_flags = [2][6]AtomFlags{
        .{ ar_pos, ar_acc, ar_acc, ar_pos, h_don, h_don },
        .{ ar_pos, ar_acc, ar_acc, ar_pos, h_don, h_don },
    };
    const charged_penalties = [2]f32{ 0.00, 0.50 };
    const specs: []const OrientSpec = if (charged) &charged_specs else &neutral_specs;
    const orient_flags: []const [6]AtomFlags = if (charged) &charged_flags else &neutral_flags;
    const penalties: []const f32 = if (charged) &charged_penalties else &neutral_penalties;

    // Allocate one contiguous block for all orientation positions.
    const all_positions = try allocator.alloc(Vec3(f32), n_orientations * n_atoms_per);
    errdefer allocator.free(all_positions);

    const orientations = try allocator.alloc(Orientation, n_orientations);
    var allocated: usize = 0;
    errdefer {
        for (orientations[0..allocated]) |o| {
            if (o.flags) |f| allocator.free(f);
        }
        allocator.free(orientations);
    }
    for (specs, 0..) |spec, i| {
        const positions = all_positions[i * n_atoms_per .. (i + 1) * n_atoms_per];
        positions[0] = spec.nd1_p;
        positions[1] = spec.cd2_p;
        positions[2] = spec.ce1_p;
        positions[3] = spec.ne2_p;
        positions[4] = spec.hd1_p;
        positions[5] = spec.he2_p;

        const flags = try allocator.alloc(AtomFlags, 6);
        @memcpy(flags, &orient_flags[i]);

        orientations[i] = .{ .positions = positions, .flags = flags, .penalty = penalties[i] };
        allocated += 1;
    }

    // atom_indices: [ND1, CD2, CE1, NE2, HD1_or_placeholder, HE2_or_placeholder]
    const actual_hd1_idx = hd1_idx orelse nd1_idx;
    const actual_he2_idx = he2_idx orelse nd1_idx;

    const atom_indices = try allocator.alloc(u32, 6);
    errdefer allocator.free(atom_indices);
    atom_indices[0] = nd1_idx;
    atom_indices[1] = cd2_idx;
    atom_indices[2] = ce1_idx;
    atom_indices[3] = ne2_idx;
    atom_indices[4] = actual_hd1_idx;
    atom_indices[5] = actual_he2_idx;

    return Mover{
        .kind = .his_flip,
        .residue_idx = residue_idx,
        .atom_indices = atom_indices,
        .orientations = orientations,
        .allocator = allocator,
        .positions_backing = all_positions,
    };
}

/// Create a carboxyl flip mover for a protonated Asp/Glu whose proton's
/// oxygen is left open. 2 orientations: the H syn on the first oxygen, or
/// syn on the second; the unused H sits at the absent sentinel. Anti
/// carboxylic acid conformers are not modelled.
///
/// Atom order in atom_indices / positions: [O1, O2, H1, H2]
pub fn createCarboxylFlipper(
    allocator: std.mem.Allocator,
    atoms: []const Atom,
    o1_idx: u32,
    o2_idx: u32,
    c_idx: u32,
    h1_idx: u32,
    h2_idx: u32,
    residue_idx: u32,
    mode: bond_policy.BondLengthMode,
) !Mover {
    const n_orientations = 2;
    const n_atoms_per = 4;
    const o1 = atoms[o1_idx].pos;
    const o2 = atoms[o2_idx].pos;
    const c = atoms[c_idx].pos;
    const bond_len = bond_policy.adjustedBondLength(mode, 0.97, .O, .Hpol);
    const h1 = geometry.placeHXR2Planar(o1.cast(f64), c.cast(f64), o2.cast(f64), bond_len, 0.0).cast(f32);
    const h2 = geometry.placeHXR2Planar(o2.cast(f64), c.cast(f64), o1.cast(f64), bond_len, 0.0).cast(f32);

    const o_don = AtomFlags{ .donor = true };
    const o_acc = AtomFlags{ .acceptor = true };
    const h_don = AtomFlags{ .donor = true };
    const h_absent = AtomFlags{};
    const specs = [n_orientations][n_atoms_per]Vec3(f32){
        .{ o1, o2, h1, ABSENT_H_POS },
        .{ o1, o2, ABSENT_H_POS, h2 },
    };
    const orient_flags = [n_orientations][n_atoms_per]AtomFlags{
        .{ o_don, o_acc, h_don, h_absent },
        .{ o_acc, o_don, h_absent, h_don },
    };

    const all_positions = try allocator.alloc(Vec3(f32), n_orientations * n_atoms_per);
    errdefer allocator.free(all_positions);
    const orientations = try allocator.alloc(Orientation, n_orientations);
    var allocated: usize = 0;
    errdefer {
        for (orientations[0..allocated]) |o| {
            if (o.flags) |f| allocator.free(f);
        }
        allocator.free(orientations);
    }
    for (specs, orient_flags, 0..) |spec, flag_row, i| {
        const positions = all_positions[i * n_atoms_per .. (i + 1) * n_atoms_per];
        @memcpy(positions, &spec);
        const flags = try allocator.alloc(AtomFlags, n_atoms_per);
        @memcpy(flags, &flag_row);
        orientations[i] = .{ .positions = positions, .flags = flags };
        allocated += 1;
    }

    const atom_indices = try allocator.alloc(u32, n_atoms_per);
    errdefer allocator.free(atom_indices);
    @memcpy(atom_indices, &[n_atoms_per]u32{ o1_idx, o2_idx, h1_idx, h2_idx });

    return Mover{
        .kind = .carboxyl_flip,
        .residue_idx = residue_idx,
        .atom_indices = atom_indices,
        .orientations = orientations,
        .allocator = allocator,
        .positions_backing = all_positions,
    };
}

// ── Tests ─────────────────────────────────────────────────────────────────────

const testing = std.testing;

test "amide flipper has 2 orientations" {
    const allocator = testing.allocator;

    // Simple geometry: C at origin, O along +x, N along +y, H1/H2 near N
    const atoms = [_]Atom{
        .{ .pos = .{ .x = 1.5, .y = 0.0, .z = 0.0 } }, // 0: O (OD1)
        .{ .pos = .{ .x = 0.0, .y = 1.5, .z = 0.0 } }, // 1: N (ND2)
        .{ .pos = .{ .x = 0.2, .y = 2.3, .z = 0.5 } }, // 2: H1
        .{ .pos = .{ .x = -0.2, .y = 2.3, .z = -0.5 } }, // 3: H2
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } }, // 4: C (CG)
        .{ .pos = .{ .x = -1.0, .y = -1.0, .z = 0.0 } }, // 5: CB
    };

    var mover = try createAmideFlipper(allocator, &atoms, 0, 1, 2, 3, 4, 7, .neutron);
    defer mover.deinit();

    try testing.expectEqual(@as(usize, 2), mover.orientations.len);
    try testing.expectApproxEqAbs(@as(f32, 0.0), mover.orientations[0].penalty, 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 0.5), mover.orientations[1].penalty, 1e-6);
}

test "amide flipper swaps O and N positions" {
    const allocator = testing.allocator;

    const atoms = [_]Atom{
        .{ .pos = .{ .x = 1.5, .y = 0.0, .z = 0.0 } }, // 0: O
        .{ .pos = .{ .x = 0.0, .y = 1.5, .z = 0.0 } }, // 1: N
        .{ .pos = .{ .x = 0.2, .y = 2.3, .z = 0.5 } }, // 2: H1
        .{ .pos = .{ .x = -0.2, .y = 2.3, .z = -0.5 } }, // 3: H2
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } }, // 4: C
        .{ .pos = .{ .x = -1.0, .y = -1.0, .z = 0.0 } }, // 5: CB
    };

    var mover = try createAmideFlipper(allocator, &atoms, 0, 1, 2, 3, 4, 0, .neutron);
    defer mover.deinit();

    const orient0 = mover.orientations[0];
    const orient1 = mover.orientations[1];

    // O position in orient 0 == N position in orient 1
    try testing.expectApproxEqAbs(orient0.positions[0].x, orient1.positions[1].x, 1e-6);
    try testing.expectApproxEqAbs(orient0.positions[0].y, orient1.positions[1].y, 1e-6);
    try testing.expectApproxEqAbs(orient0.positions[0].z, orient1.positions[1].z, 1e-6);

    // N position in orient 0 == O position in orient 1
    try testing.expectApproxEqAbs(orient0.positions[1].x, orient1.positions[0].x, 1e-6);
    try testing.expectApproxEqAbs(orient0.positions[1].y, orient1.positions[0].y, 1e-6);
    try testing.expectApproxEqAbs(orient0.positions[1].z, orient1.positions[0].z, 1e-6);
}

test "his flipper has 4 orientations" {
    const allocator = testing.allocator;

    // Simple imidazole ring geometry
    const atoms = [_]Atom{
        .{ .pos = .{ .x = 0.0, .y = 1.38, .z = 0.0 } }, // 0: ND1
        .{ .pos = .{ .x = 1.19, .y = -0.38, .z = 0.0 } }, // 1: CD2
        .{ .pos = .{ .x = 1.21, .y = 0.97, .z = 0.0 } }, // 2: CE1
        .{ .pos = .{ .x = 0.0, .y = -1.38, .z = 0.0 } }, // 3: NE2
        .{ .pos = .{ .x = -0.9, .y = 2.0, .z = 0.0 } }, // 4: HD1
        .{ .pos = .{ .x = -0.9, .y = -2.0, .z = 0.0 } }, // 5: HE2
    };

    var mover = try createHisFlipper(allocator, &atoms, 0, 1, 2, 3, 4, 5, 3, .neutron, false);
    defer mover.deinit();

    try testing.expectEqual(@as(usize, 4), mover.orientations.len);
}

test "his flipper penalties match spec" {
    const allocator = testing.allocator;

    const atoms = [_]Atom{
        .{ .pos = .{ .x = 0.0, .y = 1.38, .z = 0.0 } }, // 0: ND1
        .{ .pos = .{ .x = 1.19, .y = -0.38, .z = 0.0 } }, // 1: CD2
        .{ .pos = .{ .x = 1.21, .y = 0.97, .z = 0.0 } }, // 2: CE1
        .{ .pos = .{ .x = 0.0, .y = -1.38, .z = 0.0 } }, // 3: NE2
        .{ .pos = .{ .x = -0.9, .y = 2.0, .z = 0.0 } }, // 4: HD1
        .{ .pos = .{ .x = -0.9, .y = -2.0, .z = 0.0 } }, // 5: HE2
    };

    var mover = try createHisFlipper(allocator, &atoms, 0, 1, 2, 3, 4, 5, 3, .neutron, false);
    defer mover.deinit();

    const expected_penalties = [4]f32{ 0.00, 0.00, 0.50, 0.50 };
    for (mover.orientations, 0..) |o, i| {
        try testing.expectApproxEqAbs(expected_penalties[i], o.penalty, 1e-6);
    }
}

test "his flipper orientations have chemistry flags" {
    const allocator = testing.allocator;

    const atoms = [_]Atom{
        .{ .pos = .{ .x = 0.0, .y = 1.38, .z = 0.0 } }, // 0: ND1
        .{ .pos = .{ .x = 1.19, .y = -0.38, .z = 0.0 } }, // 1: CD2
        .{ .pos = .{ .x = 1.21, .y = 0.97, .z = 0.0 } }, // 2: CE1
        .{ .pos = .{ .x = 0.0, .y = -1.38, .z = 0.0 } }, // 3: NE2
        .{ .pos = .{ .x = -0.9, .y = 2.0, .z = 0.0 } }, // 4: HD1
        .{ .pos = .{ .x = -0.9, .y = -2.0, .z = 0.0 } }, // 5: HE2
    };

    var mover = try createHisFlipper(allocator, &atoms, 0, 1, 2, 3, 4, 5, 3, .neutron, false);
    defer mover.deinit();

    // All 4 orientations should have flags
    for (mover.orientations) |o| {
        try testing.expect(o.flags != null);
        try testing.expectEqual(@as(usize, 6), o.flags.?.len);
    }

    // Orient 0: HIE (HE2 only) → ND1=acceptor (not donor), NE2=donor
    const flags0 = mover.orientations[0].flags.?;
    try testing.expect(flags0[0].acceptor); // ND1 acceptor
    try testing.expect(!flags0[0].donor); // ND1 not donor
    try testing.expect(flags0[3].donor); // NE2 donor
    try testing.expect(!flags0[3].acceptor); // NE2 not acceptor

    // Orient 1: HID (HD1 only) → ND1=donor, NE2=acceptor
    const flags1 = mover.orientations[1].flags.?;
    try testing.expect(flags1[0].donor); // ND1 donor
    try testing.expect(!flags1[0].acceptor); // ND1 not acceptor
    try testing.expect(flags1[3].acceptor); // NE2 acceptor
    try testing.expect(!flags1[3].donor); // NE2 not donor

    // No HIP orientations — all orientations have exactly one absent H
    for (mover.orientations) |o| {
        const hd1_absent = mover_mod.isAbsentH(.{ .pos = o.positions[4], .is_added = true });
        const he2_absent = mover_mod.isAbsentH(.{ .pos = o.positions[5], .is_added = true });
        // Exactly one H should be absent (HID xor HIE, never HIP)
        try testing.expect(hd1_absent != he2_absent);
    }
}

test "carboxyl flipper puts the proton on one oxygen at a time" {
    const allocator = testing.allocator;
    const atoms = [_]Atom{
        .{ .pos = .{ .x = 1.2, .y = 0.7, .z = 0.0 } }, // 0: OD1
        .{ .pos = .{ .x = -1.2, .y = 0.7, .z = 0.0 } }, // 1: OD2
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } }, // 2: CG
        .{ .pos = .{ .x = 1.5, .y = 1.5, .z = 0.0 } }, // 3: HD1
        .{ .pos = Mover.ABSENT_H_POS }, // 4: HD2
    };
    var mover = try createCarboxylFlipper(allocator, &atoms, 0, 1, 2, 3, 4, 7, .neutron);
    defer mover.deinit();
    try testing.expectEqual(mover_mod.MoverKind.carboxyl_flip, mover.kind);
    try testing.expectEqual(@as(u16, 2), mover.nOrientations());
    for (mover.orientations, 0..) |o, i| {
        const h1_absent = mover_mod.isAbsentPos(o.positions[2]);
        const h2_absent = mover_mod.isAbsentPos(o.positions[3]);
        try testing.expect(h1_absent != h2_absent);
        try testing.expectEqual(i == 1, h1_absent);
        // The protonated oxygen is the donor.
        try testing.expectEqual(!h1_absent, o.flags.?[0].donor);
        try testing.expectEqual(!h2_absent, o.flags.?[1].donor);
        // O-H bond length.
        const h = if (h1_absent) o.positions[3] else o.positions[2];
        const parent = if (h1_absent) o.positions[1] else o.positions[0];
        try testing.expectApproxEqAbs(@as(f32, 0.97), h.sub(parent).length(), 0.05);
    }
}

test "charged His flipper keeps both ring H in both flip states" {
    const allocator = testing.allocator;
    const atoms = [_]Atom{
        .{ .pos = .{ .x = 1.0, .y = 0.0, .z = 0.0 } }, // 0: ND1
        .{ .pos = .{ .x = -0.5, .y = 1.2, .z = 0.0 } }, // 1: CD2
        .{ .pos = .{ .x = 0.8, .y = 1.3, .z = 0.0 } }, // 2: CE1
        .{ .pos = .{ .x = 0.0, .y = 2.0, .z = 0.0 } }, // 3: NE2
        .{ .pos = .{ .x = 2.0, .y = -0.3, .z = 0.0 } }, // 4: HD1
        .{ .pos = .{ .x = -0.3, .y = 3.0, .z = 0.0 } }, // 5: HE2
    };
    var mover = try createHisFlipper(allocator, &atoms, 0, 1, 2, 3, 4, 5, 3, .neutron, true);
    defer mover.deinit();
    try testing.expectEqual(@as(u16, 2), mover.nOrientations());
    for (mover.orientations) |o| {
        try testing.expect(!mover_mod.isAbsentPos(o.positions[4]));
        try testing.expect(!mover_mod.isAbsentPos(o.positions[5]));
        try testing.expect(o.flags.?[0].positive and o.flags.?[3].positive);
    }
    try testing.expect(mover.orientations[1].penalty > mover.orientations[0].penalty);
}
