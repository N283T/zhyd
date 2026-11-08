//! End-to-end tests with real PDB structures.
//!
//! Tests chain gap handling, insertion code handling, and full pipeline
//! correctness on deposited PDB structures (examples/data/*.cif.gz).
//!
//! Structures used:
//!   1rqf — 4 chains, mid-chain gaps in B/C (CK2 beta subunit)
//!   3rk2 — 4 chains, internal gaps in B/F (SNARE complex)
//!   6fys — 1 chain, 8 insertion code positions (nanobody)
//!   2cf8 — 3 chains, chymotrypsinogen insertion codes (thrombin)
//!   2hnt — 4 chains, both gaps AND insertion codes (gamma-thrombin)

const std = @import("std");
const testing = std.testing;
const zhyd = @import("root.zig");
const mmcif = zhyd.mmcif;
const place = zhyd.place;
const optimize = zhyd.optimize;
const cif = zhyd.cif;

/// Helper: run the full pipeline on a .cif.gz file, return model and placement stats.
/// Caller owns model and must deinit.
const PipelineResult = struct {
    model: zhyd.model.Model,
    n_placed: u32,
    n_residues: u32,
    n_skipped_existing: u32,
    n_skipped_inter_residue: u32,
    n_skipped_missing_ref: u32,
    n_movers: usize,

    pub fn deinit(self: *PipelineResult) void {
        self.model.deinit();
    }
};

/// Run the full pipeline without an external CCD dictionary.
/// Structures that need CCD-derived placement (non-standard ligands)
/// will have those residues skipped, which is acceptable for these tests.
fn runPipeline(allocator: std.mem.Allocator, path: []const u8) !PipelineResult {
    // 1. Read and decompress
    const source = try zhyd.run.readFile(allocator, path);
    defer allocator.free(source);

    // 2. Parse CIF document
    var doc = try cif.readString(allocator, source);
    defer doc.deinit();
    const block = &doc.blocks.items[0];

    // 3. Parse model
    var mdl = try mmcif.parseModel(allocator, source);
    errdefer mdl.deinit();

    // 4. Inline components
    var inline_dict = try mmcif.parseInlineComponents(allocator, block);
    defer if (inline_dict) |*d| d.deinit();

    // 5. Chemistry + bonds
    place.applyChemistry(&mdl);
    var atom_lookup = try mmcif.buildAtomLookup(allocator, block);
    defer atom_lookup.deinit();
    try mmcif.parseStructConn(&mdl, block, &atom_lookup);
    try mmcif.parseBranchLinks(allocator, &mdl, block, &atom_lookup);
    mmcif.flagLeavingAtoms(&mdl, if (inline_dict) |*d| d else null, null);

    // 6. Place hydrogens
    const place_result = try place.addHydrogens(&mdl, null, if (inline_dict) |*d| d else null);

    // 7. Optimize
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
        _ = try optimize.optimizer.optimize(allocator, movers, &mdl, .{ .n_threads = 1 });
    }

    // 8. Mark absent H
    for (mdl.atoms.items) |*atom| {
        if (optimize.mover.isAbsentH(atom.*)) atom.is_added = false;
    }

    // 9. Validate (should not crash on any real structure)
    var validation = try zhyd.validate.validateModel(allocator, &mdl);
    validation.deinit();

    return .{
        .model = mdl,
        .n_placed = place_result.n_placed,
        .n_residues = place_result.n_residues,
        .n_skipped_existing = place_result.n_skipped_existing,
        .n_skipped_inter_residue = place_result.n_skipped_inter_residue,
        .n_skipped_missing_ref = place_result.n_skipped_missing_ref,
        .n_movers = movers.len,
    };
}

/// Count added hydrogen atoms in the final model.
fn countFinalH(mdl: *const zhyd.model.Model) u32 {
    var count: u32 = 0;
    for (mdl.atoms.items) |atom| {
        if (atom.is_added and atom.is_hydrogen) count += 1;
    }
    return count;
}

/// Count residues with is_chain_break_before = true.
fn countChainBreaks(mdl: *const zhyd.model.Model) u32 {
    var count: u32 = 0;
    for (mdl.residues.items) |res| {
        if (res.is_chain_break_before) count += 1;
    }
    return count;
}

/// Check that the output is valid mmCIF (write and re-parse).
fn verifyOutputRoundTrip(allocator: std.mem.Allocator, mdl: *const zhyd.model.Model) !void {
    var buf: std.Io.Writer.Allocating = .init(allocator);
    defer buf.deinit();
    try zhyd.writer.mmcif_writer.write(&buf.writer, mdl, "TEST");

    // Re-parse to verify valid CIF
    var doc2 = try cif.readString(allocator, buf.writer.buffered());
    defer doc2.deinit();
    try testing.expect(doc2.blocks.items.len > 0);
    const loop = doc2.blocks.items[0].findLoop("_atom_site.Cartn_x");
    try testing.expect(loop != null);
}

/// Verify added H bond lengths. Allows up to 1% outliers for real structures
/// (strained geometry, low-occupancy altlocs, etc.).
fn verifyBondLengths(mdl: *const zhyd.model.Model) !void {
    var n_checked: u32 = 0;
    var n_bad: u32 = 0;

    for (mdl.atoms.items) |atom| {
        if (!atom.is_added or !atom.is_hydrogen) continue;

        // Find nearest heavy atom in same residue as proxy for bond length.
        // Not exact (nearest != bonded parent) but sufficient for smoke testing.
        var min_dist: f32 = std.math.inf(f32);
        for (mdl.atoms.items) |other| {
            if (other.is_hydrogen) continue;
            if (other.residue_idx != atom.residue_idx) continue;
            const dist = atom.pos.distance(other.pos);
            if (dist < min_dist) min_dist = dist;
        }

        n_checked += 1;
        if (min_dist < 0.5 or min_dist > 1.6) {
            n_bad += 1;
        }
    }

    // Allow up to 1% outliers (strained real-world geometries)
    const max_bad = @max(n_checked / 100, 3);
    if (n_bad > max_bad) {
        std.debug.print(
            "\nToo many bad H bond lengths: {d}/{d} (max allowed: {d})\n",
            .{ n_bad, n_checked, max_bad },
        );
        return error.TestUnexpectedResult;
    }
}

// ============================================================
// Chain gap tests
// ============================================================

test "real-file: 1rqf — chain gaps, 4-chain CK2 beta" {
    // Paths are relative to project root (zig build test CWD).
    var result = try runPipeline(testing.allocator, "examples/data/1rqf.cif.gz");
    defer result.deinit();

    // 1rqf: ~1418 residues, ~10467 H placed (measured); bounds ±10%
    try testing.expect(result.n_placed > 9420);
    try testing.expect(result.n_placed < 11514);
    try testing.expect(result.n_residues > 1276);
    try testing.expect(result.n_residues < 1560);
    try testing.expect(result.n_movers > 0);

    // Mid-chain gaps in B/C
    const n_breaks = countChainBreaks(&result.model);
    try testing.expect(n_breaks > 0);

    const final_h = countFinalH(&result.model);
    try testing.expect(final_h > 9420);
    try testing.expect(final_h < 11514);

    try verifyBondLengths(&result.model);
    try verifyOutputRoundTrip(testing.allocator, &result.model);
}

test "real-file: 3rk2 — chain gaps, SNARE complex" {
    var result = try runPipeline(testing.allocator, "examples/data/3rk2.cif.gz");
    defer result.deinit();

    // 3rk2: ~505 residues, ~3734 H placed (measured); bounds ±10%
    try testing.expect(result.n_placed > 3361);
    try testing.expect(result.n_placed < 4107);
    try testing.expect(result.n_residues > 455);
    try testing.expect(result.n_residues < 556);
    try testing.expect(result.n_movers > 0);

    const n_breaks = countChainBreaks(&result.model);
    try testing.expect(n_breaks > 0);

    const final_h = countFinalH(&result.model);
    try testing.expect(final_h > 3361);
    try testing.expect(final_h < 4107);

    try verifyBondLengths(&result.model);
    try verifyOutputRoundTrip(testing.allocator, &result.model);
}

// ============================================================
// Insertion code tests
// ============================================================

test "real-file: 6fys — insertion codes, nanobody" {
    var result = try runPipeline(testing.allocator, "examples/data/6fys.cif.gz");
    defer result.deinit();

    // 6fys: ~1218 residues (with HOH), ~7320 H placed (measured); bounds ±10%
    try testing.expect(result.n_placed > 6588);
    try testing.expect(result.n_placed < 8052);
    try testing.expect(result.n_residues > 1096);
    try testing.expect(result.n_residues < 1340);
    try testing.expect(result.n_movers > 0);

    const final_h = countFinalH(&result.model);
    try testing.expect(final_h > 6588);
    try testing.expect(final_h < 8052);

    try verifyBondLengths(&result.model);
    try verifyOutputRoundTrip(testing.allocator, &result.model);
}

test "real-file: 2cf8 — insertion codes, thrombin 3-chain" {
    var result = try runPipeline(testing.allocator, "examples/data/2cf8.cif.gz");
    defer result.deinit();

    // 2cf8: ~684 residues, ~151 H placed (many HOH skipped without CCD); bounds ±30%
    try testing.expect(result.n_placed > 106);
    try testing.expect(result.n_placed < 196);
    try testing.expect(result.n_residues > 477);
    try testing.expect(result.n_residues < 889);

    // Multi-chain: should have at least 3 chains
    try testing.expect(result.model.chains.items.len >= 3);

    const final_h = countFinalH(&result.model);
    try testing.expect(final_h > 106);
    try testing.expect(final_h < 196);

    try verifyBondLengths(&result.model);
    try verifyOutputRoundTrip(testing.allocator, &result.model);
}

// ============================================================
// Both gaps and insertion codes
// ============================================================

test "real-file: 2hnt — gaps + insertion codes, gamma-thrombin" {
    var result = try runPipeline(testing.allocator, "examples/data/2hnt.cif.gz");
    defer result.deinit();

    // 2hnt: ~448 residues, ~2100 H placed (measured); bounds ±10%
    try testing.expect(result.n_placed > 1890);
    try testing.expect(result.n_placed < 2310);
    try testing.expect(result.n_residues > 403);
    try testing.expect(result.n_residues < 493);
    try testing.expect(result.n_movers > 0);

    // Should have chain breaks
    const n_breaks = countChainBreaks(&result.model);
    try testing.expect(n_breaks > 0);

    // Multi-chain
    try testing.expect(result.model.chains.items.len >= 2);

    const final_h = countFinalH(&result.model);
    try testing.expect(final_h > 1890);
    try testing.expect(final_h < 2310);

    try verifyBondLengths(&result.model);
    try verifyOutputRoundTrip(testing.allocator, &result.model);
}

// ============================================================
// Experimental pKa model on real structures
// ============================================================

fn loadForPka(allocator: std.mem.Allocator, path: []const u8) !zhyd.model.Model {
    const source = try zhyd.run.readFile(allocator, path);
    defer allocator.free(source);
    var mdl = try mmcif.parseModel(allocator, source);
    errdefer mdl.deinit();
    place.applyChemistry(&mdl);
    return mdl;
}

test "real-file pka: 2cf8 disulfide cysteines excluded, values bounded, deterministic" {
    var mdl = try loadForPka(testing.allocator, "examples/data/2cf8.cif.gz");
    defer mdl.deinit();
    var first = try zhyd.pka.predict(testing.allocator, &mdl, .{});
    defer first.deinit();
    var second = try zhyd.pka.predict(testing.allocator, &mdl, .{});
    defer second.deinit();

    // Thrombin heavy+light chains: every cysteine is in a disulfide.
    var n_disulfide: usize = 0;
    for (first.excluded) |site| {
        if (site.reason == .disulfide) n_disulfide += 1;
    }
    try testing.expect(n_disulfide >= 8);
    for (first.sites) |site| try testing.expect(site.kind != .cys);

    try testing.expect(first.sites.len > 50);
    try testing.expectEqual(first.sites.len, second.sites.len);
    for (first.sites, second.sites) |a, b| {
        try testing.expect(a.predicted_pka >= -2.0 and a.predicted_pka <= 16.0);
        try testing.expectEqual(a.predicted_pka, b.predicted_pka);
    }
}

test "real-file pka: crystal waters do not change structure terms" {
    var mdl = try loadForPka(testing.allocator, "examples/data/2cf8.cif.gz");
    defer mdl.deinit();
    var with_water = try zhyd.pka.predict(testing.allocator, &mdl, .{ .include_coupling = false });
    defer with_water.deinit();

    // Move every water far away; burial and H-bond terms must be identical.
    var n_water: usize = 0;
    for (mdl.atoms.items) |*atom| {
        const comp = mdl.residues.items[atom.residue_idx].compIdSlice();
        if (!std.mem.eql(u8, comp, "HOH")) continue;
        atom.pos = .{ .x = 1.0e4, .y = 1.0e4, .z = 1.0e4 };
        n_water += 1;
    }
    try testing.expect(n_water > 0);
    var without_water = try zhyd.pka.predict(testing.allocator, &mdl, .{ .include_coupling = false });
    defer without_water.deinit();

    for (with_water.sites, without_water.sites) |a, b| {
        try testing.expectEqual(a.terms.burial, b.terms.burial);
        try testing.expectEqual(a.terms.hydrogen_bond, b.terms.hydrogen_bond);
        try testing.expectEqual(a.burial_fraction, b.burial_fraction);
    }
}

test "real-file pka: 2cf8 PROPKA engine matches PROPKA 3.5.1" {
    if (!zhyd.pka_propka.available) return error.SkipZigTest;
    var mdl = try loadForPka(testing.allocator, "examples/data/2cf8.cif.gz");
    defer mdl.deinit();
    var result = try zhyd.pka_propka.predict(testing.allocator, &mdl, .{});
    defer result.deinit();

    // Reference: propka3 3.5.1 on the same model (hydrogens and waters
    // removed, PDB written by gemmi). Chain H starts at Ile16, so its N+
    // group also checks TER handling between chains.
    const Expected = struct { chain: []const u8, num: i32, kind: zhyd.pka.SiteKind, pka: f32 };
    const expected = [_]Expected{
        .{ .chain = "H", .num = 194, .kind = .asp, .pka = 4.41 },
        .{ .chain = "H", .num = 192, .kind = .glu, .pka = 4.71 },
        .{ .chain = "H", .num = 57, .kind = .his, .pka = 4.86 },
        .{ .chain = "H", .num = 16, .kind = .n_term, .pka = 7.67 },
        .{ .chain = "I", .num = 2, .kind = .n_term, .pka = 7.49 },
    };
    for (expected) |e| {
        var found = false;
        for (result.sites) |site| {
            const res = mdl.residues.items[site.residue_idx];
            const chain = mdl.chains.items[res.chain_idx];
            if (site.kind != e.kind or res.auth_seq_id != e.num) continue;
            if (!std.mem.eql(u8, chain.authSlice(), e.chain)) continue;
            try testing.expectApproxEqAbs(e.pka, site.predicted_pka, 0.006);
            found = true;
        }
        try testing.expect(found);
    }
    // PROPKA marks all eight cysteines 99.99; zhyd reports them as excluded.
    var n_disulfide: usize = 0;
    for (result.excluded) |site| {
        if (site.reason == .disulfide) n_disulfide += 1;
    }
    try testing.expectEqual(@as(usize, 8), n_disulfide);
    for (result.sites) |site| try testing.expect(site.kind != .cys);
}
