//! End-to-end integration tests for the full zhyd pipeline.
//!
//! These tests verify the complete flow: parse → place H → write → re-parse.

const std = @import("std");
const testing = std.testing;
const mmcif = @import("mmcif.zig");
const place = @import("place.zig");
const writer = @import("writer.zig");
const cif = @import("cif.zig");
const optimize = @import("optimize/optimize.zig");

test "end-to-end: tiny.cif placement" {
    const source = @embedFile("test_data/tiny.cif");
    var model = try mmcif.parseModel(testing.allocator, source);
    defer model.deinit();

    // Place hydrogens
    const result = try place.addHydrogens(&model, null, null);

    // ALA should get H atoms placed
    try testing.expect(result.n_placed > 0);
    try testing.expect(model.atoms.items.len > 5); // 5 heavy + H atoms

    // Write output and verify it's valid CIF
    var buf: std.Io.Writer.Allocating = .init(testing.allocator);
    defer buf.deinit();
    try writer.mmcif_writer.write(&buf.writer, &model, "TEST");

    // Re-parse the output (slices in doc point into buf.writer.buffered(), so buf must stay alive)
    var doc = try cif.readString(testing.allocator, buf.writer.buffered());
    defer doc.deinit();
    try testing.expect(doc.blocks.items.len > 0);

    // Verify _atom_site loop exists with more atoms than original
    const block = &doc.blocks.items[0];
    const loop = block.findLoop("_atom_site.Cartn_x");
    try testing.expect(loop != null);
    try testing.expect(loop.?.length() > 5);
}

test "end-to-end: multi-chain placement" {
    const source = @embedFile("test_data/multi_chain.cif");
    var model = try mmcif.parseModel(testing.allocator, source);
    defer model.deinit();

    const initial_atoms = model.atoms.items.len;
    try testing.expectEqual(@as(usize, 11), initial_atoms);

    const result = try place.addHydrogens(&model, null, null);

    // Should place H on multiple residues
    try testing.expect(result.n_residues >= 2);
    try testing.expect(result.n_placed > 0);
    try testing.expect(model.atoms.items.len > initial_atoms);

    // All added atoms should be hydrogen
    for (model.atoms.items[initial_atoms..]) |atom| {
        try testing.expect(atom.is_hydrogen);
        try testing.expect(atom.is_added);
    }
}

test "end-to-end: H bond lengths are physical" {
    // Use coordinates that form a proper tetrahedral geometry around CA so
    // the placement algorithms do not hit degenerate cross-product cases.
    // ALA heavy atoms at roughly correct bond lengths and angles:
    //   N  at (-1.20,  0.00, 0.00)
    //   CA at ( 0.00,  0.00, 0.00)
    //   C  at ( 0.55,  1.42, 0.00)
    //   O  at ( 1.72,  1.60, 0.00)
    //   CB at ( 0.55, -0.76, 1.20)
    const realistic_cif =
        \\data_REAL
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
        \\ATOM 1 N  N  ALA A 1 -1.200  0.000  0.000 1.00 10.0 .
        \\ATOM 2 C  CA ALA A 1  0.000  0.000  0.000 1.00 10.0 .
        \\ATOM 3 C  C  ALA A 1  0.550  1.420  0.000 1.00 10.0 .
        \\ATOM 4 O  O  ALA A 1  1.720  1.600  0.000 1.00 10.0 .
        \\ATOM 5 C  CB ALA A 1  0.550 -0.760  1.200 1.00 10.0 .
        \\#
    ;

    var model = try mmcif.parseModel(testing.allocator, realistic_cif);
    defer model.deinit();

    const initial_count = model.atoms.items.len;
    _ = try place.addHydrogens(&model, null, null);

    // Should have added at least one H
    try testing.expect(model.atoms.items.len > initial_count);

    // Check that all added H atoms have reasonable bond lengths to their
    // parent heavy atoms (0.8 – 1.5 Å for typical C-H / N-H bonds).
    for (model.atoms.items) |atom| {
        if (!atom.is_added) continue;

        // Find the closest heavy atom in the same residue
        var min_dist: f32 = std.math.inf(f32);
        for (model.atoms.items) |other| {
            if (other.is_hydrogen) continue;
            if (other.residue_idx != atom.residue_idx) continue;
            const dist = atom.pos.distance(other.pos);
            if (dist < min_dist) min_dist = dist;
        }

        // Bond length should be physical for C-H / N-H bonds
        try testing.expect(min_dist > 0.8);
        try testing.expect(min_dist < 1.5);
    }
}

test "end-to-end: JSON log output" {
    const source = @embedFile("test_data/tiny.cif");
    var model = try mmcif.parseModel(testing.allocator, source);
    defer model.deinit();

    const result = try place.addHydrogens(&model, null, null);
    const n_added: u32 = @intCast(model.atoms.items.len - 5);

    var buf: std.Io.Writer.Allocating = .init(testing.allocator);
    defer buf.deinit();

    try writer.json_writer.writeLog(
        &buf.writer,
        "0.1.0",
        "tiny.cif",
        n_added,
        .{},
        &.{}, // no movers
        model.residues.items,
        model.chains.items,
    );

    const output = buf.writer.buffered();
    try testing.expect(std.mem.find(u8, output, "\"version\"") != null);
    try testing.expect(std.mem.find(u8, output, "\"hydrogens_added\"") != null);
    _ = result;
}

test "integration: disulfide SG-HG not placed" {
    const source = @embedFile("test_data/disulfide.cif");
    var doc = cif.readString(testing.allocator, source) catch unreachable;
    defer doc.deinit();
    const block = &doc.blocks.items[0];

    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    // Full pipeline: chemistry first, then bond parsing, then placement
    place.applyChemistry(&mdl);
    var lookup = try mmcif.buildAtomLookup(testing.allocator, block);
    defer lookup.deinit();
    try mmcif.parseStructConn(&mdl, block, &lookup);
    _ = try place.addHydrogens(&mdl, null, null);

    // Check no HG atom was placed on bonded SG atoms
    for (mdl.atoms.items) |atom| {
        if (atom.is_added and std.mem.eql(u8, atom.nameSlice(), "HG")) {
            try testing.expect(false); // unexpected HG placement on disulfide SG
        }
    }
}

test "integration: inline dict used for placement" {
    const source = @embedFile("test_data/inline_comp.cif");
    var doc = cif.readString(testing.allocator, source) catch unreachable;
    defer doc.deinit();
    const block = &doc.blocks.items[0];

    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    var inline_dict = try mmcif.parseInlineComponents(testing.allocator, block);
    defer if (inline_dict) |*d| d.deinit();

    place.applyChemistry(&mdl);
    const result = try place.addHydrogens(&mdl, null, if (inline_dict) |*d| d else null);

    // ALA should have hydrogens placed via inline dict
    try testing.expect(result.n_placed > 0);
}

test "end-to-end: HIS sentinel cleanup matches output and JSON count" {
    const source = @embedFile("test_data/his.cif");
    var model = try mmcif.parseModel(testing.allocator, source);
    defer model.deinit();

    place.applyChemistry(&model);
    _ = try place.addHydrogens(&model, null, null);

    const gen_result = try optimize.generateMovers(testing.allocator, &model, false, null, null, null, .neutron);
    const movers = gen_result.movers;
    defer {
        for (0..movers.len) |i| @constCast(&movers[i]).deinit();
        testing.allocator.free(movers);
    }

    _ = try optimize.optimizer.optimize(testing.allocator, movers, &model, .{});

    var final_h_count: u32 = 0;
    for (model.atoms.items) |*atom| {
        if (optimize.mover.isAbsentH(atom.*)) atom.is_added = false;
        if (atom.is_added and atom.is_hydrogen) final_h_count += 1;
    }

    var cif_buf: std.Io.Writer.Allocating = .init(testing.allocator);
    defer cif_buf.deinit();
    try writer.mmcif_writer.write(&cif_buf.writer, &model, "HIS");

    var h_rows: u32 = 0;
    var lines = std.mem.tokenizeScalar(u8, cif_buf.writer.buffered(), '\n');
    while (lines.next()) |line| {
        if (!std.mem.startsWith(u8, line, "ATOM ")) continue;
        if (std.mem.find(u8, line, " H ")) |_| {
            h_rows += 1;
        }
    }

    try testing.expectEqual(final_h_count, h_rows);

    var json_buf: std.Io.Writer.Allocating = .init(testing.allocator);
    defer json_buf.deinit();
    try writer.json_writer.writeLog(
        &json_buf.writer,
        "0.1.0",
        "his.cif",
        final_h_count,
        .{},
        movers,
        model.residues.items,
        model.chains.items,
    );
    try testing.expect(std.mem.find(u8, json_buf.writer.buffered(), "\"hydrogens_added\": 9") != null);
}
