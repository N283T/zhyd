//! Regression tests for issue #233.
//!
//! Covers three areas:
//!   1. Pinned H-count + coordinate baselines for small fixtures (tiny.cif,
//!      ala_with_h.cif, asn.cif, his.cif, multi_chain.cif).
//!   2. Multithreaded vs single-threaded optimization equivalence.
//!   3. CCD ligand / HETATM placement via inline _chem_comp_* dictionary.
//!
//! Baseline values were captured by running the current build once and
//! hard-coded here.  Any change to the placer geometry or bond policies
//! that shifts a coordinate beyond ±0.001 Å must be explicitly reviewed.

const std = @import("std");
const testing = std.testing;
const model_mod = @import("model.zig");
const Model = model_mod.Model;
const mmcif = @import("mmcif.zig");
const place = @import("place.zig");
const optimize = @import("optimize/optimize.zig");
const cif = @import("cif.zig");

// ============================================================
// Helpers
// ============================================================

const HPos = struct { name: [4]u8, x: f32, y: f32, z: f32 };

/// Collect all placed-H positions for a model into a heap-allocated slice.
/// Caller must free.
fn collectPlacedH(allocator: std.mem.Allocator, mdl: *const Model) ![]HPos {
    var list = std.ArrayListUnmanaged(HPos).empty;
    for (mdl.atoms.items) |atom| {
        if (!atom.is_added or !atom.is_hydrogen) continue;
        try list.append(allocator, .{
            .name = atom.name.buf,
            .x = atom.pos.x,
            .y = atom.pos.y,
            .z = atom.pos.z,
        });
    }
    return list.toOwnedSlice(allocator);
}

/// Run place-only (no optimization) on an in-memory CIF string.
fn runPlaceOnly(allocator: std.mem.Allocator, source: []const u8) !struct { mdl: Model, h_pos: []HPos } {
    var mdl = try mmcif.parseModel(allocator, source);
    errdefer mdl.deinit();
    place.applyChemistry(&mdl);
    _ = try place.addHydrogens(&mdl, null, null);
    const h_pos = try collectPlacedH(allocator, &mdl);
    return .{ .mdl = mdl, .h_pos = h_pos };
}

/// Run place-only with inline component dictionary (for CCD ligand test).
fn runPlaceWithInlineDict(
    allocator: std.mem.Allocator,
    source: []const u8,
) !struct { mdl: Model, h_pos: []HPos } {
    var doc = try cif.readString(allocator, source);
    defer doc.deinit();
    const block = &doc.blocks.items[0];

    var mdl = try mmcif.parseModel(allocator, source);
    errdefer mdl.deinit();

    var inline_dict = try mmcif.parseInlineComponents(allocator, block);
    defer if (inline_dict) |*d| d.deinit();

    place.applyChemistry(&mdl);

    var lookup = try mmcif.buildAtomLookup(allocator, block);
    defer lookup.deinit();
    try mmcif.parseStructConn(&mdl, block, &lookup);
    try mmcif.parseBranchLinks(allocator, &mdl, block, &lookup);
    mmcif.flagLeavingAtoms(&mdl, if (inline_dict) |*d| d else null, null);

    _ = try place.addHydrogens(&mdl, null, if (inline_dict) |*d| d else null);

    const h_pos = try collectPlacedH(allocator, &mdl);
    return .{ .mdl = mdl, .h_pos = h_pos };
}

/// Run full pipeline (place + optimize) with given thread count.
fn runPlaceAndOptimize(
    allocator: std.mem.Allocator,
    source: []const u8,
    n_threads: u32,
) !struct { mdl: Model, h_pos: []HPos } {
    var doc = try cif.readString(allocator, source);
    defer doc.deinit();
    const block = &doc.blocks.items[0];

    var mdl = try mmcif.parseModel(allocator, source);
    errdefer mdl.deinit();

    var inline_dict = try mmcif.parseInlineComponents(allocator, block);
    defer if (inline_dict) |*d| d.deinit();

    place.applyChemistry(&mdl);

    var lookup = try mmcif.buildAtomLookup(allocator, block);
    defer lookup.deinit();
    try mmcif.parseStructConn(&mdl, block, &lookup);
    try mmcif.parseBranchLinks(allocator, &mdl, block, &lookup);
    mmcif.flagLeavingAtoms(&mdl, if (inline_dict) |*d| d else null, null);

    _ = try place.addHydrogens(&mdl, null, if (inline_dict) |*d| d else null);

    const gen_result = try optimize.generateMovers(
        allocator,
        &mdl,
        false,
        null,
        if (inline_dict) |*d| d else null,
        null,
        .neutron,
    );
    var movers = gen_result.movers;
    defer {
        for (0..movers.len) |i| movers[i].deinit();
        allocator.free(movers);
    }

    if (movers.len > 0) {
        _ = try optimize.optimizer.optimize(allocator, movers, &mdl, .{ .n_threads = n_threads });
    }

    for (mdl.atoms.items) |*atom| {
        if (optimize.mover.isAbsentH(atom.*)) atom.is_added = false;
    }

    const h_pos = try collectPlacedH(allocator, &mdl);
    return .{ .mdl = mdl, .h_pos = h_pos };
}

/// Assert a single H position matches a pinned baseline within tolerance.
fn assertHPos(h_pos: []const HPos, name: []const u8, ex: f32, ey: f32, ez: f32, tol: f32) !void {
    for (h_pos) |h| {
        if (!std.mem.eql(u8, std.mem.trimEnd(u8, &h.name, " "), name)) continue;
        const dx = @abs(h.x - ex);
        const dy = @abs(h.y - ey);
        const dz = @abs(h.z - ez);
        if (dx > tol or dy > tol or dz > tol) {
            std.debug.print(
                "H {s}: expected ({d:.5},{d:.5},{d:.5}) got ({d:.5},{d:.5},{d:.5})\n",
                .{ name, ex, ey, ez, h.x, h.y, h.z },
            );
            return error.TestUnexpectedResult;
        }
        return;
    }
    std.debug.print("H '{s}' not found in placed H list\n", .{name});
    return error.TestUnexpectedResult;
}

// ============================================================
// 1. Pinned regression: H counts
//
// Baselines were captured by running the placer on the current build.
// Any count change must be explicitly reviewed.
// ============================================================

test "pinned regression: tiny.cif H count is stable" {
    // ALA (treated as N-terminal single residue):
    //   HA + HB1 + HB2 + HB3 + H1 + H2 + H3 = 7 H atoms
    const source = @embedFile("test_data/tiny.cif");
    var r = try runPlaceOnly(testing.allocator, source);
    defer {
        r.mdl.deinit();
        testing.allocator.free(r.h_pos);
    }
    try testing.expectEqual(@as(usize, 7), r.h_pos.len);
}

test "pinned regression: ala_with_h.cif H count is stable" {
    // ala_with_h already has HA → skipped; placed = HB1 + HB2 + HB3 + H1 + H2 + H3 = 6
    const source = @embedFile("test_data/ala_with_h.cif");
    var r = try runPlaceOnly(testing.allocator, source);
    defer {
        r.mdl.deinit();
        testing.allocator.free(r.h_pos);
    }
    try testing.expectEqual(@as(usize, 6), r.h_pos.len);
}

test "pinned regression: asn.cif H count is stable" {
    // ASN N-terminal: HA + HB2 + HB3 + HD21 + HD22 + H1 + H2 + H3 = 8 H atoms
    const source = @embedFile("test_data/asn.cif");
    var r = try runPlaceOnly(testing.allocator, source);
    defer {
        r.mdl.deinit();
        testing.allocator.free(r.h_pos);
    }
    try testing.expectEqual(@as(usize, 8), r.h_pos.len);
}

test "pinned regression: his.cif H count is stable" {
    // HIS N-terminal delta-tautomer:
    //   HA + HB2 + HB3 + HD2 + HE1 + HD1 + HE2 + H1 + H2 + H3 = 10 H atoms
    const source = @embedFile("test_data/his.cif");
    var r = try runPlaceOnly(testing.allocator, source);
    defer {
        r.mdl.deinit();
        testing.allocator.free(r.h_pos);
    }
    try testing.expectEqual(@as(usize, 10), r.h_pos.len);
}

test "pinned regression: multi_chain.cif H count is stable" {
    // ALA(A) + GLY(A) + VAL(B) each as N-terminal single residues:
    //   ALA: H1+H2+H3 = 3; GLY: H+HA2+HA3 = 3; VAL: H1+H2+H3 = 3 → total = 9
    const source = @embedFile("test_data/multi_chain.cif");
    var r = try runPlaceOnly(testing.allocator, source);
    defer {
        r.mdl.deinit();
        testing.allocator.free(r.h_pos);
    }
    try testing.expectEqual(@as(usize, 9), r.h_pos.len);
}

// ============================================================
// 2. Pinned regression: specific H coordinates
//
// The following tests pin representative coordinates for each fixture to
// catch any silent shift in placement geometry.  Tolerance is ±0.001 Å.
// ============================================================

test "pinned regression: tiny.cif HA coordinate within 0.001 Å of baseline" {
    // HA baseline: x=1.37069, y=3.62931, z=4.62931
    const source = @embedFile("test_data/tiny.cif");
    var r = try runPlaceOnly(testing.allocator, source);
    defer {
        r.mdl.deinit();
        testing.allocator.free(r.h_pos);
    }
    try assertHPos(r.h_pos, "HA", 1.37069, 3.62931, 4.62931, 0.001);
}

test "pinned regression: asn.cif HA coordinate within 0.001 Å of baseline" {
    // HA baseline: x=-0.35068, y=-0.49946, z=-0.90314
    const source = @embedFile("test_data/asn.cif");
    var r = try runPlaceOnly(testing.allocator, source);
    defer {
        r.mdl.deinit();
        testing.allocator.free(r.h_pos);
    }
    try assertHPos(r.h_pos, "HA", -0.35068, -0.49946, -0.90314, 0.001);
}

test "pinned regression: asn.cif HD21 coordinate within 0.001 Å of baseline" {
    // HD21 baseline: x=-0.92609, y=-3.87171, z=0.34020
    const source = @embedFile("test_data/asn.cif");
    var r = try runPlaceOnly(testing.allocator, source);
    defer {
        r.mdl.deinit();
        testing.allocator.free(r.h_pos);
    }
    try assertHPos(r.h_pos, "HD21", -0.92609, -3.87171, 0.34020, 0.001);
}

test "pinned regression: his.cif HA coordinate within 0.001 Å of baseline" {
    // HA baseline: x=-0.35068, y=-0.49946, z=-0.90314
    const source = @embedFile("test_data/his.cif");
    var r = try runPlaceOnly(testing.allocator, source);
    defer {
        r.mdl.deinit();
        testing.allocator.free(r.h_pos);
    }
    try assertHPos(r.h_pos, "HA", -0.35068, -0.49946, -0.90314, 0.001);
}

test "pinned regression: his.cif HD1 coordinate within 0.001 Å of baseline" {
    // HD1 baseline: x=1.26163, y=-2.50608, z=2.55485
    const source = @embedFile("test_data/his.cif");
    var r = try runPlaceOnly(testing.allocator, source);
    defer {
        r.mdl.deinit();
        testing.allocator.free(r.h_pos);
    }
    try assertHPos(r.h_pos, "HD1", 1.26163, -2.50608, 2.55485, 0.001);
}

test "pinned regression: his.cif HE1 coordinate within 0.001 Å of baseline" {
    // HE1 baseline: x=1.26984, y=-5.00327, z=1.84322
    const source = @embedFile("test_data/his.cif");
    var r = try runPlaceOnly(testing.allocator, source);
    defer {
        r.mdl.deinit();
        testing.allocator.free(r.h_pos);
    }
    try assertHPos(r.h_pos, "HE1", 1.26984, -5.00327, 1.84322, 0.001);
}

// ============================================================
// 3. Multithreaded optimization equivalence
// ============================================================

test "multithreaded equivalence: his.cif n_threads=1 vs n_threads=4 produce same H count" {
    // HIS has a His-ring flipper; running with 4 threads must yield the same
    // final H count as 1 thread.
    const source = @embedFile("test_data/his.cif");

    var r1 = try runPlaceAndOptimize(testing.allocator, source, 1);
    defer {
        r1.mdl.deinit();
        testing.allocator.free(r1.h_pos);
    }

    var r4 = try runPlaceAndOptimize(testing.allocator, source, 4);
    defer {
        r4.mdl.deinit();
        testing.allocator.free(r4.h_pos);
    }

    try testing.expectEqual(r1.h_pos.len, r4.h_pos.len);
}

test "multithreaded equivalence: asn.cif n_threads=1 vs n_threads=4 H count and max position delta" {
    // ASN has an amide flipper; verify both thread counts yield identical H count
    // and near-identical positions (optimizer is deterministic per-mover).
    const source = @embedFile("test_data/asn.cif");

    var r1 = try runPlaceAndOptimize(testing.allocator, source, 1);
    defer {
        r1.mdl.deinit();
        testing.allocator.free(r1.h_pos);
    }

    var r4 = try runPlaceAndOptimize(testing.allocator, source, 4);
    defer {
        r4.mdl.deinit();
        testing.allocator.free(r4.h_pos);
    }

    // H count must be identical.
    try testing.expectEqual(r1.h_pos.len, r4.h_pos.len);

    // Compute max coordinate delta across all placed H atoms.
    // The optimizer processes each mover independently so results should be
    // bit-identical across thread counts; we allow 0.01 Å for any float-order
    // non-associativity that might arise from different scheduling.
    var max_delta: f32 = 0.0;
    for (r1.h_pos, r4.h_pos) |h1, h4| {
        const dx = @abs(h1.x - h4.x);
        const dy = @abs(h1.y - h4.y);
        const dz = @abs(h1.z - h4.z);
        const d = @max(dx, @max(dy, dz));
        if (d > max_delta) max_delta = d;
    }

    std.log.info("asn.cif 1-vs-4-thread max position delta: {d:.6} Å", .{max_delta});

    try testing.expect(max_delta < 0.01);
}

// ============================================================
// 4. CCD ligand placement via inline _chem_comp_* dictionary
// ============================================================

test "CCD ligand placement: inline dict is parsed for LIG residue" {
    // The ccd_ligand.cif fixture embeds _chem_comp_atom / _chem_comp_bond loops.
    // Verify they are parsed and the LIG component is available.
    const source = @embedFile("test_data/ccd_ligand.cif");

    var doc = try cif.readString(testing.allocator, source);
    defer doc.deinit();
    const block = &doc.blocks.items[0];

    var inline_dict = try mmcif.parseInlineComponents(testing.allocator, block);
    defer if (inline_dict) |*d| d.deinit();

    try testing.expect(inline_dict != null);
    if (inline_dict) |*d| {
        try testing.expect(d.get("LIG") != null);
    }
}

test "CCD ligand placement: H atoms placed on HETATM LIG via inline dict" {
    // Verify that CCD-derived placement actually adds hydrogens to LIG and that
    // all bond lengths are physical (0.7–1.6 Å).
    const source = @embedFile("test_data/ccd_ligand.cif");

    var doc = try cif.readString(testing.allocator, source);
    defer doc.deinit();
    const block = &doc.blocks.items[0];

    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    var inline_dict = try mmcif.parseInlineComponents(testing.allocator, block);
    defer if (inline_dict) |*d| d.deinit();

    try testing.expect(inline_dict != null);

    place.applyChemistry(&mdl);

    var lookup = try mmcif.buildAtomLookup(testing.allocator, block);
    defer lookup.deinit();
    try mmcif.parseStructConn(&mdl, block, &lookup);
    try mmcif.parseBranchLinks(testing.allocator, &mdl, block, &lookup);
    mmcif.flagLeavingAtoms(&mdl, if (inline_dict) |*d| d else null, null);

    const result = try place.addHydrogens(&mdl, null, if (inline_dict) |*d| d else null);

    try testing.expect(result.n_placed > 0);

    // Verify every placed H has a physical bond length to its nearest heavy atom.
    var n_placed_h: u32 = 0;
    for (mdl.atoms.items) |atom| {
        if (!atom.is_added or !atom.is_hydrogen) continue;
        n_placed_h += 1;

        var min_dist: f32 = std.math.inf(f32);
        for (mdl.atoms.items) |other| {
            if (other.is_hydrogen) continue;
            if (other.residue_idx != atom.residue_idx) continue;
            const dist = atom.pos.distance(other.pos);
            if (dist < min_dist) min_dist = dist;
        }

        if (min_dist < 0.7 or min_dist > 1.6) {
            std.debug.print(
                "Bad H bond length for {s}: {d:.4} Å\n",
                .{ atom.nameSlice(), min_dist },
            );
            return error.TestUnexpectedResult;
        }
    }

    try testing.expect(n_placed_h > 0);

    // H on the same heavy atom must not coincide: the fixture has no ideal
    // coordinates, so their dihedrals come from the fallback spread.
    for (mdl.atoms.items, 0..) |a, i| {
        if (!a.is_added or !a.is_hydrogen) continue;
        for (mdl.atoms.items[i + 1 ..]) |b| {
            if (!b.is_added or !b.is_hydrogen) continue;
            if (a.pos.distance(b.pos) < 1.0) {
                std.debug.print("H {s} and {s} are {d:.3} A apart\n", .{ a.nameSlice(), b.nameSlice(), a.pos.distance(b.pos) });
                return error.TestUnexpectedResult;
            }
        }
    }
}

test "pinned regression: ccd_ligand.cif H count is stable" {
    // LIG: C1(2H) + C2(2H) + O1(1H) + N1(2H) = 7 H atoms
    // (C1 has 2 heavy bonds — C2 and N1 — so sp3 C gets 2H)
    const source = @embedFile("test_data/ccd_ligand.cif");
    var r = try runPlaceWithInlineDict(testing.allocator, source);
    defer {
        r.mdl.deinit();
        testing.allocator.free(r.h_pos);
    }
    try testing.expectEqual(@as(usize, 7), r.h_pos.len);
}

test "pinned regression: ccd_ligand.cif HO1 coordinate within 0.001 Å of baseline" {
    // HO1 baseline: x=3.06608, y=1.13283, z=0.00000
    const source = @embedFile("test_data/ccd_ligand.cif");
    var r = try runPlaceWithInlineDict(testing.allocator, source);
    defer {
        r.mdl.deinit();
        testing.allocator.free(r.h_pos);
    }
    try assertHPos(r.h_pos, "HO1", 3.06608, 1.13283, 0.00000, 0.001);
}
