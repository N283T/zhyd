//! Integration tests for the hydrogen placement pipeline.
//!
//! Extracted from placer.zig to keep the production module focused.
//! Tests exercise the public API: addHydrogens, addHydrogensWithConfig,
//! applyChemistry, applyChemistryWithConfig.

const std = @import("std");
const testing = std.testing;
const mmcif = @import("../mmcif.zig");
const element = @import("../element.zig");
const math_mod = @import("../math.zig");
const Vec3f32 = math_mod.Vec3(f32);
const Model = @import("../model.zig").Model;
const placer = @import("placer.zig");
const addHydrogens = placer.addHydrogens;
const addHydrogensWithConfig = placer.addHydrogensWithConfig;
const applyChemistry = placer.applyChemistry;
const applyChemistryWithConfig = placer.applyChemistryWithConfig;
const protonation = @import("protonation.zig");
const lookup = placer.lookup;
const findAtomPos = lookup.findAtomPos;
const padName = lookup.padName;
test "place hydrogens on ALA" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    const initial_count = mdl.atoms.items.len;
    try testing.expectEqual(@as(usize, 5), initial_count);

    const result = try addHydrogens(&mdl, null, null);

    // ALA should get H atoms added
    try testing.expect(result.n_placed > 0);
    try testing.expect(mdl.atoms.items.len > initial_count);
    try testing.expectEqual(@as(u32, 1), result.n_residues);

    // Find HA atom and check bond length to CA (~1.10 A)
    const ca_pos = mdl.atoms.items[1].pos; // CA is index 1
    var found_ha = false;
    for (mdl.atoms.items) |atom| {
        if (std.mem.eql(u8, atom.nameSlice(), "HA")) {
            found_ha = true;
            const dist = atom.pos.distance(ca_pos);
            try testing.expect(dist > 0.8);
            try testing.expect(dist < 1.4);
            break;
        }
    }
    try testing.expect(found_ha);
}

test "xray mode produces shorter C-H bond than neutron" {
    const source = @embedFile("../test_data/tiny.cif");

    // Place in neutron mode (default)
    var mdl_n = try mmcif.parseModel(testing.allocator, source);
    defer mdl_n.deinit();
    _ = try addHydrogensWithConfig(&mdl_n, null, null, .{ .bond_policy = .{ .mode = .neutron } });

    // Place in xray mode
    var mdl_x = try mmcif.parseModel(testing.allocator, source);
    defer mdl_x.deinit();
    _ = try addHydrogensWithConfig(&mdl_x, null, null, .{ .bond_policy = .{ .mode = .xray } });

    // Find HA in both and compare distance to CA
    const ca_pos_n = mdl_n.atoms.items[1].pos;
    const ca_pos_x = mdl_x.atoms.items[1].pos;
    var dist_n: f32 = 0;
    var dist_x: f32 = 0;
    for (mdl_n.atoms.items) |atom| {
        if (std.mem.eql(u8, atom.nameSlice(), "HA")) {
            dist_n = atom.pos.distance(ca_pos_n);
            break;
        }
    }
    for (mdl_x.atoms.items) |atom| {
        if (std.mem.eql(u8, atom.nameSlice(), "HA")) {
            dist_x = atom.pos.distance(ca_pos_x);
            break;
        }
    }
    // Neutron C-H ~1.10, xray C-H ~0.98
    try testing.expect(dist_n > 1.0);
    try testing.expectApproxEqAbs(@as(f32, 0.98), dist_x, 0.02);
    try testing.expect(dist_x < dist_n);
}

test "protonation override fixes HIS tautomer during placement" {
    const source = @embedFile("../test_data/his.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();
    mdl.residues.items[0].auth_seq_id = 1;

    var overrides = try protonation.parseString(testing.allocator,
        \\A:1 HIS HIE
    );
    defer overrides.deinit();

    applyChemistryWithConfig(&mdl, .{ .protonation = &overrides });
    _ = try addHydrogensWithConfig(&mdl, null, null, .{ .protonation = &overrides });

    try testing.expect(findAddedAtomIdx(&mdl, 0, "HE2") != null);
    try testing.expect(findAddedAtomIdx(&mdl, 0, "HD1") == null);
}

test "protonation override adds ASP sidechain proton" {
    const source =
        \\data_ASP
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\_atom_site.auth_asym_id
        \\_atom_site.auth_seq_id
        \\ATOM 1 N N ASP A 1 0.0 0.0 0.0 1.00 10.0 . A 1
        \\ATOM 2 C CA ASP A 1 1.5 0.0 0.0 1.00 10.0 . A 1
        \\ATOM 3 C C ASP A 1 2.1 1.4 0.0 1.00 10.0 . A 1
        \\ATOM 4 O O ASP A 1 3.3 1.6 0.0 1.00 10.0 . A 1
        \\ATOM 5 C CB ASP A 1 2.0 -0.8 1.2 1.00 10.0 . A 1
        \\ATOM 6 C CG ASP A 1 3.4 -0.4 1.4 1.00 10.0 . A 1
        \\ATOM 7 O OD1 ASP A 1 4.2 0.3 0.7 1.00 10.0 . A 1
        \\ATOM 8 O OD2 ASP A 1 3.8 -0.8 2.6 1.00 10.0 . A 1
        \\#
    ;
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    var overrides = try protonation.parseString(testing.allocator,
        \\A:1 ASP OD2
    );
    defer overrides.deinit();

    applyChemistryWithConfig(&mdl, .{ .protonation = &overrides });
    _ = try addHydrogensWithConfig(&mdl, null, null, .{ .protonation = &overrides });

    try testing.expect(findAddedAtomIdx(&mdl, 0, "HD2") != null);
}

test "protonation override adds GLU sidechain proton" {
    const source =
        \\data_GLU
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\_atom_site.auth_asym_id
        \\_atom_site.auth_seq_id
        \\ATOM 1 N N GLU A 1 0.0 0.0 0.0 1.00 10.0 . A 1
        \\ATOM 2 C CA GLU A 1 1.5 0.0 0.0 1.00 10.0 . A 1
        \\ATOM 3 C C GLU A 1 2.1 1.4 0.0 1.00 10.0 . A 1
        \\ATOM 4 O O GLU A 1 3.3 1.6 0.0 1.00 10.0 . A 1
        \\ATOM 5 C CB GLU A 1 2.0 -0.8 1.2 1.00 10.0 . A 1
        \\ATOM 6 C CG GLU A 1 3.4 -0.4 1.4 1.00 10.0 . A 1
        \\ATOM 7 C CD GLU A 1 4.3 -1.3 0.6 1.00 10.0 . A 1
        \\ATOM 8 O OE1 GLU A 1 5.5 -0.9 0.4 1.00 10.0 . A 1
        \\ATOM 9 O OE2 GLU A 1 3.8 -2.4 0.2 1.00 10.0 . A 1
        \\#
    ;
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    var overrides = try protonation.parseString(testing.allocator,
        \\A:1 GLU OE2
    );
    defer overrides.deinit();

    applyChemistryWithConfig(&mdl, .{ .protonation = &overrides });
    _ = try addHydrogensWithConfig(&mdl, null, null, .{ .protonation = &overrides });

    try testing.expect(findAddedAtomIdx(&mdl, 0, "HE2") != null);
    try testing.expect(findAddedAtomIdx(&mdl, 0, "HE1") == null);
}

test "protonation override LYS neutral skips HZ3" {
    const source =
        \\data_LYS
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\_atom_site.auth_asym_id
        \\_atom_site.auth_seq_id
        \\ATOM  1 N N   LYS A 1  0.000  0.000  0.000 1.00 10.0 . A 1
        \\ATOM  2 C CA  LYS A 1  1.458  0.000  0.000 1.00 10.0 . A 1
        \\ATOM  3 C C   LYS A 1  2.009  1.420  0.000 1.00 10.0 . A 1
        \\ATOM  4 O O   LYS A 1  3.200  1.600  0.000 1.00 10.0 . A 1
        \\ATOM  5 C CB  LYS A 1  1.986 -0.760  1.220 1.00 10.0 . A 1
        \\ATOM  6 C CG  LYS A 1  3.500 -0.800  1.220 1.00 10.0 . A 1
        \\ATOM  7 C CD  LYS A 1  4.028 -1.560  2.440 1.00 10.0 . A 1
        \\ATOM  8 C CE  LYS A 1  5.542 -1.600  2.440 1.00 10.0 . A 1
        \\ATOM  9 N NZ  LYS A 1  6.070 -2.360  3.660 1.00 10.0 . A 1
        \\#
    ;
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    var overrides = try protonation.parseString(testing.allocator,
        \\A:1 LYS NEUTRAL
    );
    defer overrides.deinit();

    applyChemistryWithConfig(&mdl, .{ .protonation = &overrides });
    _ = try addHydrogensWithConfig(&mdl, null, null, .{ .protonation = &overrides });

    try testing.expect(findAddedAtomIdx(&mdl, 0, "HZ1") != null);
    try testing.expect(findAddedAtomIdx(&mdl, 0, "HZ2") != null);
    try testing.expect(findAddedAtomIdx(&mdl, 0, "HZ3") == null);
}

test "protonation override CYS thiolate skips HG" {
    const source =
        \\data_CYS
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\_atom_site.auth_asym_id
        \\_atom_site.auth_seq_id
        \\ATOM 1 N N   CYS A 1  0.000  0.000  0.000 1.00 10.0 . A 1
        \\ATOM 2 C CA  CYS A 1  1.458  0.000  0.000 1.00 10.0 . A 1
        \\ATOM 3 C C   CYS A 1  2.009  1.420  0.000 1.00 10.0 . A 1
        \\ATOM 4 O O   CYS A 1  3.200  1.600  0.000 1.00 10.0 . A 1
        \\ATOM 5 C CB  CYS A 1  1.986 -0.760  1.220 1.00 10.0 . A 1
        \\ATOM 6 S SG  CYS A 1  1.300 -0.200  2.800 1.00 10.0 . A 1
        \\#
    ;
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    var overrides = try protonation.parseString(testing.allocator,
        \\A:1 CYS THIOLATE
    );
    defer overrides.deinit();

    applyChemistryWithConfig(&mdl, .{ .protonation = &overrides });
    _ = try addHydrogensWithConfig(&mdl, null, null, .{ .protonation = &overrides });

    try testing.expect(findAddedAtomIdx(&mdl, 0, "HG") == null);
}

fn findAddedAtomIdx(mdl: *const Model, residue_idx: u32, name: []const u8) ?u32 {
    for (mdl.atoms.items, 0..) |atom, idx| {
        if (atom.residue_idx != residue_idx) continue;
        if (!atom.is_added) continue;
        if (std.mem.eql(u8, atom.nameSlice(), name)) return @intCast(idx);
    }
    return null;
}

test "placed atoms have correct metadata" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    _ = try addHydrogens(&mdl, null, null);

    // Check that newly added atoms have correct flags
    for (mdl.atoms.items[5..]) |atom| {
        try testing.expect(atom.is_hydrogen);
        try testing.expect(atom.is_added);
        try testing.expectEqual(@as(u32, 0), atom.residue_idx);
    }
}

test "PlacementResult tracks counts" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    const result = try addHydrogens(&mdl, null, null);

    // ALA has 5 plans; backbone H skipped on N-term but NH3+ (H1,H2,H3) added = 4+3=7
    try testing.expectEqual(@as(u32, 7), result.n_placed + result.totalSkipped());
}

test "placed H inherits parent atom metadata" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    // tiny.cif atoms have occupancy=1.0, b_factor=10.0, altloc=' '
    _ = try addHydrogens(&mdl, null, null);

    // All placed H atoms should inherit b_factor=10.0 from parent
    for (mdl.atoms.items) |atom| {
        if (atom.is_added and atom.is_hydrogen) {
            try testing.expectEqual(@as(f32, 10.0), atom.b_factor);
            try testing.expectEqual(@as(f32, 1.0), atom.occupancy);
        }
    }
}

test "duplicate H atoms are not placed" {
    const source = @embedFile("../test_data/ala_with_h.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    _ = try addHydrogens(&mdl, null, null);

    // Count HA atoms — should be exactly 1 (the pre-existing one)
    var ha_count: u32 = 0;
    for (mdl.atoms.items) |atom| {
        if (std.mem.eql(u8, atom.nameSlice(), "HA")) ha_count += 1;
    }
    try testing.expectEqual(@as(u32, 1), ha_count);

    // The original HA should NOT be overwritten (b_factor should remain 12.0)
    for (mdl.atoms.items) |atom| {
        if (std.mem.eql(u8, atom.nameSlice(), "HA")) {
            try testing.expectEqual(@as(f32, 12.0), atom.b_factor);
            break;
        }
    }
}

test "PlacementResult counts duplicates as skipped" {
    const source = @embedFile("../test_data/ala_with_h.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    const result = try addHydrogens(&mdl, null, null);

    // HA was pre-existing so should be counted as skipped (existing_h)
    // Total plans attempted should still be the same as clean ALA
    try testing.expect(result.n_skipped_existing >= 1);
    try testing.expectEqual(@as(u32, 1), result.n_residues);
}

test "placement succeeds on stretched geometry with bond topology" {
    const source = @embedFile("../test_data/ala_stretched.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    _ = try addHydrogens(&mdl, null, null);

    // With bond topology, HA should be placed even though CB is >1.9A from CA
    // (HA placement type is hxr3 which needs the 3rd neighbor = CB)
    var found_ha = false;
    for (mdl.atoms.items) |atom| {
        if (std.mem.eql(u8, atom.nameSlice(), "HA")) {
            found_ha = true;
            break;
        }
    }
    try testing.expect(found_ha);
}

test "applyChemistry sets backbone C to C_eq_O" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    // Before: C has generic element type
    // tiny.cif atoms: N(0), CA(1), C(2), O(3), CB(4)
    try testing.expectEqual(element.AtomType.C, mdl.atoms.items[2].element_type);

    applyChemistry(&mdl);

    // After: C has carbonyl type
    try testing.expectEqual(element.AtomType.C_eq_O, mdl.atoms.items[2].element_type);
    try testing.expectApproxEqAbs(@as(f32, 1.65), mdl.atoms.items[2].vdw_radius, 1e-6);
}

test "applyChemistry sets backbone O acceptor flag" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    try testing.expect(!mdl.atoms.items[3].flags.acceptor);
    applyChemistry(&mdl);
    try testing.expect(mdl.atoms.items[3].flags.acceptor);
}

test "applyChemistry sets backbone N donor flag" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    try testing.expect(!mdl.atoms.items[0].flags.donor);
    applyChemistry(&mdl);
    try testing.expect(mdl.atoms.items[0].flags.donor);
}

test "placed H atoms have correct flags from element table" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);
    _ = try addHydrogens(&mdl, null, null);

    for (mdl.atoms.items) |atom| {
        if (atom.is_added and atom.is_hydrogen) {
            const type_flags = atom.element_type.info().flags;
            // Hpol atoms should have donor flag
            if (atom.element_type == .Hpol) {
                try testing.expect(atom.flags.donor);
            }
            // All placed H flags should match their element_type flags
            try testing.expectEqual(type_flags.donor, atom.flags.donor);
            try testing.expectEqual(type_flags.aromatic, atom.flags.aromatic);
        }
    }
}

test "applyChemistry adds positive flag to N-terminal N" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);

    // tiny.cif has single ALA — both N-term and C-term
    // N atom (index 0) should have donor (standard) + positive (terminal)
    try testing.expect(mdl.atoms.items[0].flags.donor);
    try testing.expect(mdl.atoms.items[0].flags.positive);
}

test "applyChemistry adds negative flag to C-terminal O" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);

    // O atom (index 3) should have acceptor (standard) + negative (terminal)
    try testing.expect(mdl.atoms.items[3].flags.acceptor);
    try testing.expect(mdl.atoms.items[3].flags.negative);
}

test "applyChemistry annotates OXT as negative acceptor" {
    const source = @embedFile("../test_data/ala_cterm.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);

    // Find OXT atom
    var oxt_found = false;
    for (mdl.atoms.items) |atom| {
        if (std.mem.eql(u8, atom.nameSlice(), "OXT")) {
            oxt_found = true;
            try testing.expectEqual(element.AtomType.O, atom.element_type);
            try testing.expect(atom.flags.negative);
            try testing.expect(atom.flags.acceptor);
            break;
        }
    }
    try testing.expect(oxt_found);
}

test "multi-chain terminal detection is correct" {
    const source = @embedFile("../test_data/multi_chain.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);

    // Chain A: ALA(N-term, res 0) + GLY(C-term, res 1)
    // Chain B: VAL(both N-term and C-term, res 2)

    // ALA N (atom 0): N-terminal → positive + donor
    try testing.expect(mdl.atoms.items[0].flags.positive);
    try testing.expect(mdl.atoms.items[0].flags.donor);

    // GLY N (atom 4): internal-ish (C-terminal residue, but N is not annotated for C-term)
    try testing.expect(!mdl.atoms.items[4].flags.positive);

    // GLY O (atom 7): C-terminal → negative + acceptor
    try testing.expect(mdl.atoms.items[7].flags.negative);
    try testing.expect(mdl.atoms.items[7].flags.acceptor);

    // VAL N (atom 8): N-terminal of chain B → positive + donor
    try testing.expect(mdl.atoms.items[8].flags.positive);
    try testing.expect(mdl.atoms.items[8].flags.donor);
}

test "OXT does not receive hydrogen atoms" {
    const source = @embedFile("../test_data/ala_cterm.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);
    _ = try addHydrogens(&mdl, null, null);

    // No hydrogen should be bonded to OXT
    for (mdl.atoms.items) |atom| {
        if (atom.is_added and atom.is_hydrogen) {
            const name = atom.nameSlice();
            try testing.expect(!std.mem.eql(u8, name, "HOXT"));
        }
    }
}

test "multi-conformer residue places H per conformer" {
    const source = @embedFile("../test_data/ala_altloc.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);
    const result = try addHydrogens(&mdl, null, null);

    // Should place H for both conformers A and B
    var ha_a_count: u32 = 0;
    var ha_b_count: u32 = 0;
    for (mdl.atoms.items) |atom| {
        if (atom.is_added and std.mem.eql(u8, atom.nameSlice(), "HA")) {
            if (atom.altloc == 'A') ha_a_count += 1;
            if (atom.altloc == 'B') ha_b_count += 1;
        }
    }
    try testing.expectEqual(@as(u32, 1), ha_a_count);
    try testing.expectEqual(@as(u32, 1), ha_b_count);
    try testing.expect(result.n_placed > 0);
}

test "conformer A and B H atoms have different positions" {
    const source = @embedFile("../test_data/ala_altloc.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);
    _ = try addHydrogens(&mdl, null, null);

    var ha_a_pos: ?Vec3f32 = null;
    var ha_b_pos: ?Vec3f32 = null;
    for (mdl.atoms.items) |atom| {
        if (atom.is_added and std.mem.eql(u8, atom.nameSlice(), "HA")) {
            if (atom.altloc == 'A') ha_a_pos = atom.pos;
            if (atom.altloc == 'B') ha_b_pos = atom.pos;
        }
    }
    try testing.expect(ha_a_pos != null);
    try testing.expect(ha_b_pos != null);
    const diff = ha_a_pos.?.distance(ha_b_pos.?);
    try testing.expect(diff > 0.01);
}

test "placed H atoms have correct mover_hint" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);
    _ = try addHydrogens(&mdl, null, null);

    // ALA methyl H (HB1/HB2/HB3) should have rotate_methyl hint
    var methyl_count: u32 = 0;
    for (mdl.atoms.items) |atom| {
        if (atom.is_added and atom.mover_hint == .rotate_methyl) {
            methyl_count += 1;
        }
    }
    // ALA has 3 methyl H on CB
    try testing.expectEqual(@as(u32, 3), methyl_count);
}

test "chain break residue keeps single backbone H" {
    const source = @embedFile("../test_data/gap_chain.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);
    _ = try addHydrogens(&mdl, null, null);

    // Second residue (seq_id 3, index 1) is after a chain break.
    // It should keep the single backbone amide H, not gain NH3+.
    var h_count: u32 = 0;
    var h1_count: u32 = 0;
    var h2_count: u32 = 0;
    var h3_count: u32 = 0;
    for (mdl.atoms.items) |atom| {
        if (!atom.is_added or atom.residue_idx != 1) continue;
        if (std.mem.eql(u8, atom.nameSlice(), "H")) h_count += 1;
        if (std.mem.eql(u8, atom.nameSlice(), "H1")) h1_count += 1;
        if (std.mem.eql(u8, atom.nameSlice(), "H2")) h2_count += 1;
        if (std.mem.eql(u8, atom.nameSlice(), "H3")) h3_count += 1;
    }
    try testing.expectEqual(@as(u32, 1), h_count);
    try testing.expectEqual(@as(u32, 0), h1_count);
    try testing.expectEqual(@as(u32, 0), h2_count);
    try testing.expectEqual(@as(u32, 0), h3_count);
}

test "chain break residue is not annotated as positively charged N-terminus" {
    const source = @embedFile("../test_data/gap_chain.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);

    const res1 = mdl.residues.items[1];
    const atoms1 = mdl.atoms.items[res1.atom_start..res1.atom_end];
    var n_positive = false;
    for (atoms1) |atom| {
        if (std.mem.eql(u8, atom.nameSlice(), "N")) {
            n_positive = atom.flags.positive;
        }
    }
    try testing.expect(!n_positive);
}

test "first observed residue with N-terminal disorder gets NH3+" {
    const source = @embedFile("../test_data/nterm_disorder.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    // First observed residue (seq_id 2, index 0) has is_chain_break_before
    // because seq_id 1 is unobserved. It is still the physical N-terminus
    // and should receive NH3+ (H1, H2, H3), not a single amide H.
    try testing.expect(mdl.residues.items[0].is_chain_break_before);

    applyChemistry(&mdl);
    _ = try addHydrogens(&mdl, null, null);

    var h_count: u32 = 0;
    var h1_count: u32 = 0;
    var h2_count: u32 = 0;
    var h3_count: u32 = 0;
    for (mdl.atoms.items) |atom| {
        if (!atom.is_added or atom.residue_idx != 0) continue;
        if (std.mem.eql(u8, atom.nameSlice(), "H")) h_count += 1;
        if (std.mem.eql(u8, atom.nameSlice(), "H1")) h1_count += 1;
        if (std.mem.eql(u8, atom.nameSlice(), "H2")) h2_count += 1;
        if (std.mem.eql(u8, atom.nameSlice(), "H3")) h3_count += 1;
    }
    // NH3+: no single backbone H, but H1/H2/H3 present
    try testing.expectEqual(@as(u32, 0), h_count);
    try testing.expectEqual(@as(u32, 1), h1_count);
    try testing.expectEqual(@as(u32, 1), h2_count);
    try testing.expectEqual(@as(u32, 1), h3_count);
}

test "first observed residue with N-terminal disorder gets positive charge" {
    const source = @embedFile("../test_data/nterm_disorder.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);

    const res0 = mdl.residues.items[0];
    const atoms0 = mdl.atoms.items[res0.atom_start..res0.atom_end];
    var n_positive = false;
    for (atoms0) |atom| {
        if (std.mem.eql(u8, atom.nameSlice(), "N")) {
            n_positive = atom.flags.positive;
        }
    }
    try testing.expect(n_positive);
}

test "residue before chain break gets C-terminal charge" {
    const source = @embedFile("../test_data/gap_chain.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);

    // First residue (seq_id 1, index 0): before gap -> C-terminal -> O gets negative
    const res0 = mdl.residues.items[0];
    const atoms0 = mdl.atoms.items[res0.atom_start..res0.atom_end];
    var o_negative = false;
    for (atoms0) |atom| {
        if (std.mem.eql(u8, atom.nameSlice(), "O")) {
            o_negative = atom.flags.negative;
        }
    }
    try testing.expect(o_negative);
}

test "addHydrogens skips H on bonded_inter_residue atom" {
    const mmcif_mod = @import("../mmcif.zig");
    const source = @embedFile("../test_data/disulfide.cif");
    var mdl = try mmcif_mod.parseModel(testing.allocator, source);
    defer mdl.deinit();

    // Apply chemistry first (it overwrites flags), then set bonded_inter_residue
    applyChemistry(&mdl);

    // Manually set bonded_inter_residue on SG atoms (index 5 and 11)
    mdl.atoms.items[5].flags.bonded_inter_residue = true;
    mdl.atoms.items[11].flags.bonded_inter_residue = true;

    const result = try addHydrogens(&mdl, null, null);

    // Verify no HG was placed on either CYS SG
    for (mdl.atoms.items) |atom| {
        if (atom.is_added) {
            const name = atom.nameSlice();
            // SG should not have HG placed
            try std.testing.expect(!std.mem.eql(u8, name, "HG"));
        }
    }
    _ = result;
}

test "water placement adds two hydrogens when enabled" {
    const water_cif =
        \\data_WATER
        \\#
        \\loop_
        \\_entity.id
        \\_entity.type
        \\1 polymer
        \\2 water
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_entity_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\ATOM   1  N  N   GLY A 1 1  0.000 0.000 0.000 1.00 10.0 .
        \\ATOM   2  O  O   GLY A 1 1  1.200 1.100 0.000 1.00 10.0 .
        \\ATOM   3  O  O   HOH B 2 1  2.700 0.200 0.000 1.00 12.0 .
        \\#
    ;

    var mdl = try mmcif.parseModel(testing.allocator, water_cif);
    defer mdl.deinit();

    const result = try addHydrogensWithConfig(&mdl, null, null, .{
        .water = .{ .enabled = true },
    });

    try testing.expectEqual(@as(u32, 2), result.n_placed);

    const oxygen_pos = lookup.findAtomPos(&mdl, mdl.residues.items[1], padName("O"), ' ') orelse unreachable;
    var h_positions: [2]Vec3f32 = undefined;
    var h_count: u32 = 0;
    for (mdl.atoms.items) |atom| {
        if (atom.is_added and atom.residue_idx == 1) {
            try testing.expect(atom.is_hydrogen);
            try testing.expectApproxEqAbs(@as(f32, 1.0), atom.occupancy, 1e-6);
            try testing.expectApproxEqAbs(@as(f32, 12.0), atom.b_factor, 1e-6);
            if (h_count < 2) h_positions[h_count] = atom.pos;
            h_count += 1;
        }
    }
    try testing.expectEqual(@as(u32, 2), h_count);

    // Verify O-H bond length (~0.97 A)
    try testing.expectApproxEqAbs(@as(f32, 0.97), oxygen_pos.distance(h_positions[0]), 0.02);
    try testing.expectApproxEqAbs(@as(f32, 0.97), oxygen_pos.distance(h_positions[1]), 0.02);

    // Verify H-O-H angle (~104.5 degrees)
    const hoh_angle = math_mod.angle(f32, h_positions[0], oxygen_pos, h_positions[1]);
    try testing.expectApproxEqAbs(@as(f32, 104.5), hoh_angle, 1.0);
}

test "water placement respects occupancy cutoff" {
    const water_cif =
        \\data_WATER_OCC
        \\#
        \\loop_
        \\_entity.id
        \\_entity.type
        \\1 polymer
        \\2 water
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_entity_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\ATOM   1  N  N   GLY A 1 1  0.000 0.000 0.000 1.00 10.0 .
        \\ATOM   2  O  O   GLY A 1 1  1.200 1.100 0.000 1.00 10.0 .
        \\ATOM   3  O  O   HOH B 2 1  2.700 0.200 0.000 0.50 12.0 .
        \\#
    ;

    var mdl = try mmcif.parseModel(testing.allocator, water_cif);
    defer mdl.deinit();

    const result = try addHydrogensWithConfig(&mdl, null, null, .{
        .water = .{
            .enabled = true,
            .occupancy_cutoff = 0.66,
        },
    });

    try testing.expectEqual(@as(u32, 0), result.n_placed);
    try testing.expectEqual(@as(u32, 1), result.n_skipped_quality_filter);
}

test "water placement skips coordinated water oxygen" {
    const water_cif =
        \\data_WATER_METAL
        \\#
        \\loop_
        \\_entity.id
        \\_entity.type
        \\1 polymer
        \\2 water
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_entity_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\ATOM   1  N  N   GLY A 1 1  0.000 0.000 0.000 1.00 10.0 .
        \\ATOM   2  O  O   GLY A 1 1  1.200 1.100 0.000 1.00 10.0 .
        \\ATOM   3  O  O   HOH B 2 1  2.700 0.200 0.000 1.00 12.0 .
        \\#
    ;

    var mdl = try mmcif.parseModel(testing.allocator, water_cif);
    defer mdl.deinit();
    mdl.atoms.items[2].flags.bonded_inter_residue = true;

    const result = try addHydrogensWithConfig(&mdl, null, null, .{
        .water = .{ .enabled = true },
    });

    try testing.expectEqual(@as(u32, 0), result.n_placed);
    try testing.expect(result.n_skipped_inter_residue > 0);
}

test "water phantom mode places zero-occupancy hydrogens for isolated water" {
    const water_cif =
        \\data_WATER_PHANTOM
        \\#
        \\loop_
        \\_entity.id
        \\_entity.type
        \\1 water
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_entity_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\HETATM 1  O  O   HOH A 1 1  0.000 0.000 0.000 1.00 12.0 .
        \\#
    ;

    var mdl = try mmcif.parseModel(testing.allocator, water_cif);
    defer mdl.deinit();

    const result = try addHydrogensWithConfig(&mdl, null, null, .{
        .water = .{
            .enabled = true,
            .phantom = true,
        },
    });

    try testing.expectEqual(@as(u32, 2), result.n_placed);
    var zero_occ: u32 = 0;
    for (mdl.atoms.items) |atom| {
        if (atom.is_added) {
            try testing.expectApproxEqAbs(@as(f32, 0.0), atom.occupancy, 1e-6);
            zero_occ += 1;
        }
    }
    try testing.expectEqual(@as(u32, 2), zero_occ);
}

test "water placement skips water near metal by distance" {
    const water_cif =
        \\data_WATER_METAL_DIST
        \\#
        \\loop_
        \\_entity.id
        \\_entity.type
        \\1 water
        \\2 non-polymer
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_entity_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\HETATM 1  O  O   HOH A 1 1  0.000 0.000 0.000 1.00 12.0 .
        \\HETATM 2  ZN ZN   ZN B 2 1  2.500 0.000 0.000 1.00 10.0 .
        \\#
    ;

    var mdl = try mmcif.parseModel(testing.allocator, water_cif);
    defer mdl.deinit();

    const result = try addHydrogensWithConfig(&mdl, null, null, .{
        .water = .{ .enabled = true },
    });

    try testing.expectEqual(@as(u32, 0), result.n_placed);
    try testing.expect(result.n_skipped_inter_residue > 0);
}

test "water placement respects B-factor cutoff" {
    const water_cif =
        \\data_WATER_BFAC
        \\#
        \\loop_
        \\_entity.id
        \\_entity.type
        \\1 polymer
        \\2 water
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_entity_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\ATOM   1  N  N   GLY A 1 1  0.000 0.000 0.000 1.00 10.0 .
        \\ATOM   2  O  O   GLY A 1 1  1.200 1.100 0.000 1.00 10.0 .
        \\ATOM   3  O  O   HOH B 2 1  2.700 0.200 0.000 1.00 50.0 .
        \\#
    ;

    var mdl = try mmcif.parseModel(testing.allocator, water_cif);
    defer mdl.deinit();

    const result = try addHydrogensWithConfig(&mdl, null, null, .{
        .water = .{
            .enabled = true,
            .b_factor_cutoff = 40.0,
        },
    });

    try testing.expectEqual(@as(u32, 0), result.n_placed);
    try testing.expectEqual(@as(u32, 1), result.n_skipped_quality_filter);
}

test "water placement skips water with existing H atoms" {
    const water_cif =
        \\data_WATER_EXIST_H
        \\#
        \\loop_
        \\_entity.id
        \\_entity.type
        \\1 polymer
        \\2 water
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_entity_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\ATOM   1  N  N   GLY A 1 1  0.000 0.000 0.000 1.00 10.0 .
        \\ATOM   2  O  O   GLY A 1 1  1.200 1.100 0.000 1.00 10.0 .
        \\ATOM   3  O  O   HOH B 2 1  2.700 0.200 0.000 1.00 12.0 .
        \\ATOM   4  H  H1  HOH B 2 1  3.100 0.800 0.500 1.00 12.0 .
        \\#
    ;

    var mdl = try mmcif.parseModel(testing.allocator, water_cif);
    defer mdl.deinit();

    const result = try addHydrogensWithConfig(&mdl, null, null, .{
        .water = .{ .enabled = true },
    });

    try testing.expectEqual(@as(u32, 0), result.n_placed);
    try testing.expectEqual(@as(u32, 1), result.n_skipped_existing);
}

test "backbone NH placed in peptide plane using C(i-1)" {
    // Two-residue ALA-ALA with realistic geometry.
    // ALA 1: N-term (gets NH3+, no backbone H)
    // ALA 2: should get backbone H in the C(1)-N(2)-CA(2) peptide plane.
    const two_ala_cif =
        \\data_TWO_ALA
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\ATOM 1  N  N   ALA A 1  -1.200  0.000  0.000 1.00 10.0 .
        \\ATOM 2  C  CA  ALA A 1   0.000  0.000  0.000 1.00 10.0 .
        \\ATOM 3  C  C   ALA A 1   0.550  1.420  0.000 1.00 10.0 .
        \\ATOM 4  O  O   ALA A 1   1.720  1.600  0.000 1.00 10.0 .
        \\ATOM 5  C  CB  ALA A 1   0.550 -0.760  1.200 1.00 10.0 .
        \\ATOM 6  N  N   ALA A 2   -0.100  2.500  0.000 1.00 10.0 .
        \\ATOM 7  C  CA  ALA A 2   0.400  3.870  0.000 1.00 10.0 .
        \\ATOM 8  C  C   ALA A 2   1.920  3.900  0.000 1.00 10.0 .
        \\ATOM 9  O  O   ALA A 2   2.500  4.980  0.000 1.00 10.0 .
        \\ATOM 10 C  CB  ALA A 2  -0.100  4.600  1.200 1.00 10.0 .
        \\#
    ;

    const mmcif_mod = @import("../mmcif.zig");
    var mdl = try mmcif_mod.parseModel(testing.allocator, two_ala_cif);
    defer mdl.deinit();

    applyChemistry(&mdl);
    _ = try addHydrogens(&mdl, null, null);

    // Find the backbone H on residue 2 (ALA 2)
    var backbone_h_pos: ?Vec3f32 = null;
    for (mdl.atoms.items) |atom| {
        if (atom.is_added and atom.residue_idx == 1 and
            std.mem.eql(u8, atom.nameSlice(), "H"))
        {
            backbone_h_pos = atom.pos;
        }
    }
    // Backbone H must exist on residue 2
    const h_pos = backbone_h_pos orelse return error.TestUnexpectedResult;

    // Get reference atoms
    const n2 = lookup.findAtomPos(&mdl, mdl.residues.items[1], .{ ' ', 'N', ' ', ' ' }, ' ') orelse unreachable;
    const ca2 = lookup.findAtomPos(&mdl, mdl.residues.items[1], .{ ' ', 'C', 'A', ' ' }, ' ') orelse unreachable;
    const c1 = lookup.findAtomPos(&mdl, mdl.residues.items[0], .{ ' ', 'C', ' ', ' ' }, ' ') orelse unreachable;

    // Check H-N bond length (~1.02 A)
    try testing.expectApproxEqAbs(h_pos.distance(n2), 1.02, 0.05);

    // Check C(i-1)-N-H and CA-N-H angles are approximately equal (bisector)
    const cn_h = math_mod.angle(f32, c1, n2, h_pos);
    const ca_n_h = math_mod.angle(f32, ca2, n2, h_pos);
    try testing.expectApproxEqAbs(cn_h, ca_n_h, 5.0);

    // Both angles should be roughly 119° (peptide plane bisector)
    try testing.expect(cn_h > 110.0 and cn_h < 130.0);

    // H should lie approximately in the C(i-1)-N-CA plane (z ≈ 0 for this fixture)
    try testing.expectApproxEqAbs(h_pos.z, 0.0, 0.1);
}

// ---------------------------------------------------------------------------
// Nucleotide placement tests
// ---------------------------------------------------------------------------

test "place hydrogens on DC (DNA cytidine)" {
    const source = @embedFile("../test_data/dc_residue.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);
    const result = try addHydrogens(&mdl, null, null);

    // DC has 11 H atoms: 7 sugar + 4 base (H5, H6, H41, H42)
    // Some sugar H may be missing_ref if neighbors aren't all bonded.
    // At minimum H5, H6, H41, H42 on the base ring should be placed.
    try testing.expect(result.n_placed >= 4);
    try testing.expect(result.n_residues >= 1);

    // Verify bond lengths are reasonable (0.8 – 1.2 Å)
    for (mdl.atoms.items) |atom| {
        if (!atom.is_added or !atom.is_hydrogen) continue;
        // Find nearest heavy atom in the same residue
        var min_dist: f32 = std.math.inf(f32);
        for (mdl.atoms.items) |other| {
            if (other.is_hydrogen) continue;
            if (other.residue_idx != atom.residue_idx) continue;
            const d = atom.pos.distance(other.pos);
            if (d < min_dist) min_dist = d;
        }
        try testing.expect(min_dist > 0.8 and min_dist < 1.2);
    }
}

test "place hydrogens on RNA adenosine (A)" {
    // Use an inline CIF with idealized coordinates for adenosine (A).
    // We construct just the ribose + base heavy atoms needed to place H.
    const source =
        \\data_RNA_A
        \\#
        \\_entry.id RNA_A
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\_atom_site.auth_asym_id
        \\ATOM 1  P  P    A A 1  -2.397  0.829  1.300 1.00 10.0 . A
        \\ATOM 2  O  OP1  A A 1  -3.498  1.740  1.058 1.00 10.0 . A
        \\ATOM 3  O  OP2  A A 1  -2.703 -0.629  1.400 1.00 10.0 . A
        \\ATOM 4  O  O5'  A A 1  -1.279  1.195  0.255 1.00 10.0 . A
        \\ATOM 5  C  C5'  A A 1   0.017  0.620  0.236 1.00 10.0 . A
        \\ATOM 6  C  C4'  A A 1   0.965  1.534 -0.510 1.00 10.0 . A
        \\ATOM 7  O  O4'  A A 1   0.489  2.877 -0.513 1.00 10.0 . A
        \\ATOM 8  C  C3'  A A 1   1.232  1.052 -1.941 1.00 10.0 . A
        \\ATOM 9  O  O3'  A A 1   2.578  0.613 -2.089 1.00 10.0 . A
        \\ATOM 10 C  C2'  A A 1   0.858  2.282 -2.778 1.00 10.0 . A
        \\ATOM 11 O  O2'  A A 1   1.903  3.222 -2.982 1.00 10.0 . A
        \\ATOM 12 C  C1'  A A 1   0.817  3.372 -1.718 1.00 10.0 . A
        \\ATOM 13 N  N9   A A 1  -0.518  3.933 -1.469 1.00 10.0 . A
        \\ATOM 14 C  C8   A A 1  -0.949  5.219 -1.612 1.00 10.0 . A
        \\ATOM 15 N  N7   A A 1  -2.198  5.393 -1.366 1.00 10.0 . A
        \\ATOM 16 C  C5   A A 1  -2.600  4.163 -1.064 1.00 10.0 . A
        \\ATOM 17 C  C6   A A 1  -3.891  3.671 -0.763 1.00 10.0 . A
        \\ATOM 18 N  N6   A A 1  -4.963  4.417 -0.740 1.00 10.0 . A
        \\ATOM 19 N  N1   A A 1  -4.148  2.382 -0.529 1.00 10.0 . A
        \\ATOM 20 C  C2   A A 1  -3.128  1.467 -0.495 1.00 10.0 . A
        \\ATOM 21 N  N3   A A 1  -1.894  1.755 -0.719 1.00 10.0 . A
        \\ATOM 22 C  C4   A A 1  -1.600  3.070 -0.975 1.00 10.0 . A
        \\#
    ;
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);
    const result = try addHydrogens(&mdl, null, null);

    // RNA A has 11 H atoms: 7 sugar (H1', H2', HO2', H3', H4', H5', H5'') + 4 base (H2, H8, H61, H62)
    // In a single isolated residue some sugar H may be missing_ref (no 3rd neighbor for HXR3).
    // The base H (H2, H8, H61, H62) should all succeed.
    try testing.expect(result.n_placed >= 4);

    // Verify bond lengths for all placed H
    var n_h: u32 = 0;
    for (mdl.atoms.items) |atom| {
        if (!atom.is_added or !atom.is_hydrogen) continue;
        n_h += 1;
        var min_dist: f32 = std.math.inf(f32);
        for (mdl.atoms.items) |other| {
            if (other.is_hydrogen) continue;
            if (other.residue_idx != atom.residue_idx) continue;
            const d = atom.pos.distance(other.pos);
            if (d < min_dist) min_dist = d;
        }
        try testing.expect(min_dist > 0.8 and min_dist < 1.2);
    }
    try testing.expect(n_h >= 4);
}

// ── N-terminal mode tests (issue #251) ───────────────────────────────────────

const BackboneHCounts = struct { h: u32, h1: u32, h2: u32, h3: u32 };

fn countAddedBackboneHOnResidue(mdl: *const Model, residue_idx: u32) BackboneHCounts {
    var counts = BackboneHCounts{ .h = 0, .h1 = 0, .h2 = 0, .h3 = 0 };
    for (mdl.atoms.items) |atom| {
        if (atom.residue_idx != residue_idx) continue;
        if (!atom.is_added) continue;
        const name = atom.nameSlice();
        if (std.mem.eql(u8, name, "H")) counts.h += 1;
        if (std.mem.eql(u8, name, "H1")) counts.h1 += 1;
        if (std.mem.eql(u8, name, "H2")) counts.h2 += 1;
        if (std.mem.eql(u8, name, "H3")) counts.h3 += 1;
    }
    return counts;
}

/// Find the backbone N atom of the given residue and assert its donor and
/// positive charge flags match the expected values. Fails the test on missing N.
fn expectBackboneNFlags(mdl: *const Model, residue_idx: u32, expect_donor: bool, expect_positive: bool) !void {
    for (mdl.atoms.items) |atom| {
        if (atom.residue_idx != residue_idx) continue;
        if (atom.is_hydrogen) continue;
        if (!std.mem.eql(u8, atom.nameSlice(), "N")) continue;
        try testing.expectEqual(expect_donor, atom.flags.donor);
        try testing.expectEqual(expect_positive, atom.flags.positive);
        return;
    }
    return error.BackboneNotFound;
}

test "nterm auto (default) places NH3+ on real N-term only, break-amide on gap" {
    const source = @embedFile("../test_data/gap_chain.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistryWithConfig(&mdl, .{});
    _ = try addHydrogensWithConfig(&mdl, null, null, .{});

    // Residue 0: real N-term → H1/H2/H3 (NH3+), no single H, backbone N positive
    const r0 = countAddedBackboneHOnResidue(&mdl, 0);
    try testing.expectEqual(@as(u32, 0), r0.h);
    try testing.expectEqual(@as(u32, 1), r0.h1);
    try testing.expectEqual(@as(u32, 1), r0.h2);
    try testing.expectEqual(@as(u32, 1), r0.h3);
    try expectBackboneNFlags(&mdl, 0, true, true);

    // Residue 1: post-gap break → single H (break-amide), no NH3+, N donor only
    const r1 = countAddedBackboneHOnResidue(&mdl, 1);
    try testing.expectEqual(@as(u32, 1), r1.h);
    try testing.expectEqual(@as(u32, 0), r1.h1);
    try testing.expectEqual(@as(u32, 0), r1.h2);
    try testing.expectEqual(@as(u32, 0), r1.h3);
    try expectBackboneNFlags(&mdl, 1, true, false);
}

test "nterm aggressive places NH3+ on both real N-term and chain-break residue" {
    const source = @embedFile("../test_data/gap_chain.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistryWithConfig(&mdl, .{ .nterm_mode = .aggressive });
    _ = try addHydrogensWithConfig(&mdl, null, null, .{ .nterm_mode = .aggressive });

    const r0 = countAddedBackboneHOnResidue(&mdl, 0);
    try testing.expectEqual(@as(u32, 0), r0.h);
    try testing.expectEqual(@as(u32, 1), r0.h1);
    try testing.expectEqual(@as(u32, 1), r0.h2);
    try testing.expectEqual(@as(u32, 1), r0.h3);
    try expectBackboneNFlags(&mdl, 0, true, true);

    // In aggressive mode, the post-gap residue also gets NH3+ and the positive flag.
    const r1 = countAddedBackboneHOnResidue(&mdl, 1);
    try testing.expectEqual(@as(u32, 0), r1.h);
    try testing.expectEqual(@as(u32, 1), r1.h1);
    try testing.expectEqual(@as(u32, 1), r1.h2);
    try testing.expectEqual(@as(u32, 1), r1.h3);
    try expectBackboneNFlags(&mdl, 1, true, true);
}

test "nterm neutral: real non-PRO N-term gets NH2, break residue keeps single H" {
    const source = @embedFile("../test_data/gap_chain.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistryWithConfig(&mdl, .{ .nterm_mode = .neutral });
    _ = try addHydrogensWithConfig(&mdl, null, null, .{ .nterm_mode = .neutral });

    // Residue 0 (real N-term ALA): only H2/H3, no single H, no H1, donor-only N.
    const r0 = countAddedBackboneHOnResidue(&mdl, 0);
    try testing.expectEqual(@as(u32, 0), r0.h);
    try testing.expectEqual(@as(u32, 0), r0.h1);
    try testing.expectEqual(@as(u32, 1), r0.h2);
    try testing.expectEqual(@as(u32, 1), r0.h3);
    try expectBackboneNFlags(&mdl, 0, true, false);

    // Residue 1 (post-gap ALA): unchanged from auto — single break-amide H,
    // no NH2, donor-only N.
    const r1 = countAddedBackboneHOnResidue(&mdl, 1);
    try testing.expectEqual(@as(u32, 1), r1.h);
    try testing.expectEqual(@as(u32, 0), r1.h1);
    try testing.expectEqual(@as(u32, 0), r1.h2);
    try testing.expectEqual(@as(u32, 0), r1.h3);
    try expectBackboneNFlags(&mdl, 1, true, false);
}

test "nterm neutral on N-terminal PRO keeps NH2+ geometry AND positive flag" {
    // PRO is the documented exception: neutral mode falls back to the NH2Pro
    // secondary-amine geometry, and the positive charge flag must follow the
    // placement to keep chemistry and atom positions consistent.
    const source =
        \\data_PRO_NTERM
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\_atom_site.auth_asym_id
        \\ATOM 1 N N   PRO A 1 1.000 2.000 3.000 1.00 10.0 . A
        \\ATOM 2 C CA  PRO A 1 2.000 3.000 4.000 1.00 10.0 . A
        \\ATOM 3 C C   PRO A 1 3.000 4.000 5.000 1.00 10.0 . A
        \\ATOM 4 O O   PRO A 1 4.000 5.000 6.000 1.00 10.0 . A
        \\ATOM 5 C CB  PRO A 1 2.500 2.500 3.500 1.00 10.0 . A
        \\ATOM 6 C CG  PRO A 1 1.500 2.000 4.500 1.00 10.0 . A
        \\ATOM 7 C CD  PRO A 1 0.500 3.000 3.500 1.00 10.0 . A
        \\#
    ;
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistryWithConfig(&mdl, .{ .nterm_mode = .neutral });
    _ = try addHydrogensWithConfig(&mdl, null, null, .{ .nterm_mode = .neutral });

    // PRO N-term gets H2/H3 (NH2Pro geometry, 2 H, no H1).
    const r0 = countAddedBackboneHOnResidue(&mdl, 0);
    try testing.expectEqual(@as(u32, 0), r0.h);
    try testing.expectEqual(@as(u32, 0), r0.h1);
    try testing.expectEqual(@as(u32, 1), r0.h2);
    try testing.expectEqual(@as(u32, 1), r0.h3);

    // The documented exception: PRO keeps the positive charge flag even in
    // neutral mode because its secondary amine is protonated at physiological
    // pH. Placement and chemistry must agree.
    try expectBackboneNFlags(&mdl, 0, true, true);
}

test "nterm aggressive on N-terminal PRO places NH2+ with positive flag" {
    // PRO at a chain break under aggressive mode — same NH2Pro geometry
    // and positive flag as at a real N-terminus.
    const source =
        \\data_PRO_AGG
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\_atom_site.auth_asym_id
        \\ATOM 1 N N   PRO A 1 1.000 2.000 3.000 1.00 10.0 . A
        \\ATOM 2 C CA  PRO A 1 2.000 3.000 4.000 1.00 10.0 . A
        \\ATOM 3 C C   PRO A 1 3.000 4.000 5.000 1.00 10.0 . A
        \\ATOM 4 O O   PRO A 1 4.000 5.000 6.000 1.00 10.0 . A
        \\ATOM 5 C CB  PRO A 1 2.500 2.500 3.500 1.00 10.0 . A
        \\ATOM 6 C CG  PRO A 1 1.500 2.000 4.500 1.00 10.0 . A
        \\ATOM 7 C CD  PRO A 1 0.500 3.000 3.500 1.00 10.0 . A
        \\#
    ;
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistryWithConfig(&mdl, .{ .nterm_mode = .aggressive });
    _ = try addHydrogensWithConfig(&mdl, null, null, .{ .nterm_mode = .aggressive });

    const r0 = countAddedBackboneHOnResidue(&mdl, 0);
    try testing.expectEqual(@as(u32, 1), r0.h2);
    try testing.expectEqual(@as(u32, 1), r0.h3);
    try expectBackboneNFlags(&mdl, 0, true, true);
}

// ── NtermResult.missing_ref surfacing (issue #254) ───────────────────────────

test "N-terminal placement tallies missing_ref when backbone C is absent" {
    // First residue has no backbone C atom — placeNtermNH3 cannot compute
    // its h3xr dihedral reference and must report the failure through
    // NtermResult.missing_ref so the caller can surface it instead of
    // silently dropping the terminus.
    const source =
        \\data_MISSING_C
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\_atom_site.auth_asym_id
        \\ATOM 1 N N   ALA A 1 1.458 0.000 0.000 1.00 10.0 . A
        \\ATOM 2 C CA  ALA A 1 0.000 0.000 0.000 1.00 10.0 . A
        \\ATOM 3 O O   ALA A 1 -1.742 1.560 0.000 1.00 10.0 . A
        \\ATOM 4 C CB  ALA A 1 -0.700 -1.100 1.800 1.00 10.0 . A
        \\#
    ;
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    const result = try addHydrogensWithConfig(&mdl, null, null, .{});

    // The missing C backbone reference must surface through the summary
    // counter so a caller like main.zig can warn the user.
    try testing.expect(result.n_skipped_missing_ref > 0);

    // No H1/H2/H3 should have been placed for the broken N-terminal residue.
    const r0 = countAddedBackboneHOnResidue(&mdl, 0);
    try testing.expectEqual(@as(u32, 0), r0.h1);
    try testing.expectEqual(@as(u32, 0), r0.h2);
    try testing.expectEqual(@as(u32, 0), r0.h3);
}

test "N-terminal PRO placement tallies missing_ref when CD is absent" {
    // PRO requires both CA and CD to compute the NH2Pro geometry. If CD is
    // missing, placeNtermNH2Pro must surface the skip through missing_ref.
    const source =
        \\data_PRO_NO_CD
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\_atom_site.auth_asym_id
        \\ATOM 1 N N   PRO A 1 1.458 0.000 0.000 1.00 10.0 . A
        \\ATOM 2 C CA  PRO A 1 0.000 0.000 0.000 1.00 10.0 . A
        \\ATOM 3 C C   PRO A 1 -0.523 1.413 0.000 1.00 10.0 . A
        \\ATOM 4 O O   PRO A 1 -1.742 1.560 0.000 1.00 10.0 . A
        \\ATOM 5 C CB  PRO A 1 -0.700 -1.100 1.800 1.00 10.0 . A
        \\#
    ;
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    const result = try addHydrogensWithConfig(&mdl, null, null, .{});

    try testing.expect(result.n_skipped_missing_ref > 0);

    // No H2/H3 should be placed on backbone N without CD.
    const r0 = countAddedBackboneHOnResidue(&mdl, 0);
    try testing.expectEqual(@as(u32, 0), r0.h1);
    try testing.expectEqual(@as(u32, 0), r0.h2);
    try testing.expectEqual(@as(u32, 0), r0.h3);
}

test "distance-derive fallback places H on unknown ligand" {
    const source = @embedFile("../test_data/unknown_ligand.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    // No CCD, no SDF → distance-derive path for "XYZ" residue.
    const result = try addHydrogens(&mdl, null, null);

    // ALA (standard) + XYZ (distance-derived) = 2 residues processed.
    try testing.expectEqual(@as(u32, 2), result.n_residues);
    try testing.expect(result.n_placed > 0);
    try testing.expect(result.n_distance_derived > 0);
}

test "SDF dict overrides distance-derive fallback" {
    const sdf_source =
        \\XYZ
        \\     test
        \\  test
        \\  3  2  0  0  0  0  0  0  0  0999 V2000
        \\   10.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        \\   11.4300    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
        \\   10.0000    1.5400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        \\  1  2  1  0
        \\  1  3  1  0
        \\M  END
        \\$$$$
    ;
    const sdf_mod = @import("../sdf.zig");
    var sdf_dict = try sdf_mod.parseSdf(testing.allocator, sdf_source);
    defer sdf_dict.deinit();

    const cif_source = @embedFile("../test_data/unknown_ligand.cif");
    var mdl = try mmcif.parseModel(testing.allocator, cif_source);
    defer mdl.deinit();

    const result = try addHydrogensWithConfig(&mdl, null, null, .{
        .sdf_dict = &sdf_dict,
    });

    try testing.expectEqual(@as(u32, 2), result.n_residues);
    try testing.expect(result.n_placed > 0);
    // SDF path was used → distance_derived should be 0.
    try testing.expectEqual(@as(u32, 0), result.n_distance_derived);
}

test "applyChemistry gives ligand atoms element-default flags and radii" {
    const source = @embedFile("../test_data/unknown_ligand.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);

    var found_o1 = false;
    for (mdl.atoms.items) |atom| {
        if (!std.mem.eql(u8, atom.nameSlice(), "O1")) continue;
        found_o1 = true;
        // Ligand O must be an H-bond acceptor with the explicit O radius,
        // otherwise rotors never orient toward ligand acceptors.
        try testing.expect(atom.flags.acceptor);
        try testing.expectApproxEqAbs(element.AtomType.O.info().explicit_radius, atom.vdw_radius, 1e-6);
    }
    try testing.expect(found_o1);
}
