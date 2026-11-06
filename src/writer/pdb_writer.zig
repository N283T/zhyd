//! PDB format writer: outputs a Model with added hydrogen atoms in PDB fixed-column format.
//!
//! The writer iterates a PdbRecord slice:
//!   - raw_line records are written as-is (passthrough for HEADER, REMARK, etc.)
//!   - atom_site markers trigger output of all atoms for the current residue group
//!
//! TER records are inserted between chains, and serial numbers are renumbered
//! sequentially from 1. Absent flipper sentinel H atoms are skipped.

const std = @import("std");
const model_mod = @import("../model.zig");
const Model = model_mod.Model;
const Atom = model_mod.Atom;
const Residue = model_mod.Residue;
const Chain = model_mod.Chain;
const pdb_mod = @import("../pdb.zig");
const PdbRecord = pdb_mod.PdbRecord;
const mover_mod = @import("../optimize/mover.zig");
const place = @import("../place.zig");
const format = @import("format.zig");

const elementSymbol = format.elementSymbol;
const writeFixedFloat3 = format.writeFixedFloat3;
const writeFixedFloat2 = format.writeFixedFloat2;

// ── Atom name PDB formatting ─────────────────────────────────────────────────

/// Format an atom name into a 4-byte PDB-convention name field.
/// PDB convention: 1-char elements start at column 14 (index 1), 2-char elements
/// start at column 13 (index 0).
/// Raw atom name string (trimmed) → 4-char space-padded array.
pub fn formatAtomNamePdb(trimmed: []const u8, elem_sym: []const u8) [4]u8 {
    var out = [4]u8{ ' ', ' ', ' ', ' ' };
    if (trimmed.len == 0) return out;

    // Special: hydrogen names that start with a digit (e.g. "1HB", "2HG1")
    // are written starting at column 13 (index 0) in PDB format.
    const starts_with_digit = std.ascii.isDigit(trimmed[0]);
    if (starts_with_digit) {
        // Place from index 0, up to 4 chars
        const n = @min(trimmed.len, 4);
        for (0..n) |i| out[i] = trimmed[i];
        return out;
    }

    // 2-char element symbol → name starts at column 13 (index 0)
    // 1-char element symbol → name starts at column 14 (index 1)
    const start: usize = if (elem_sym.len >= 2) 0 else 1;

    // If the trimmed name is already 4 chars, write from index 0 regardless
    // (e.g. "HG11") — standard PDB practice for 4-char hydrogen names.
    if (trimmed.len == 4) {
        for (0..4) |i| out[i] = trimmed[i];
        return out;
    }

    const n = @min(trimmed.len, 4 - start);
    for (0..n) |i| out[start + i] = trimmed[i];
    return out;
}

// ── Float-to-fixed-width helpers ─────────────────────────────────────────────

/// Format a coordinate (8.3f) right-justified into an 8-char field.
fn writeCoord(writer: anytype, val: f32) !void {
    var buf: [32]u8 = undefined;
    var w: std.Io.Writer = .fixed(&buf);
    try writeFixedFloat3(&w, val);
    const s = w.buffered();
    if (s.len <= 8) {
        // Right-justify in 8 chars
        for (0..8 - s.len) |_| try writer.writeByte(' ');
        try writer.writeAll(s);
    } else {
        // Overflow: write as many chars as fit (truncate from left padding)
        try writer.writeAll(s[0..@min(s.len, 8)]);
    }
}

/// Format occupancy/b-factor (6.2f) right-justified into a 6-char field.
fn writeOccupancy(writer: anytype, val: f32) !void {
    var buf: [32]u8 = undefined;
    var w: std.Io.Writer = .fixed(&buf);
    try writeFixedFloat2(&w, val);
    const s = w.buffered();
    if (s.len <= 6) {
        // Right-justify in 6 chars
        for (0..6 - s.len) |_| try writer.writeByte(' ');
        try writer.writeAll(s);
    } else {
        // Overflow: write as many chars as fit (truncate from left padding)
        try writer.writeAll(s[0..@min(s.len, 6)]);
    }
}

// ── Single ATOM/HETATM line writer ───────────────────────────────────────────

/// Write a single ATOM or HETATM line.
/// `is_hetatm` controls the record type.
/// `serial` is the sequential serial number (1-based).
fn writeAtomLine(
    writer: anytype,
    atom: Atom,
    res: Residue,
    chain_id: u8,
    serial: u32,
    is_hetatm: bool,
    output_isotope: place.OutputIsotope,
) !void {
    // cols 1-6: record type
    if (is_hetatm) {
        try writer.writeAll("HETATM");
    } else {
        try writer.writeAll("ATOM  ");
    }

    // cols 7-11: serial (right-justified, 5 chars)
    try writer.print("{d:>5}", .{serial});

    // col 12: space
    try writer.writeByte(' ');

    // cols 13-16: atom name (4 chars, PDB convention)
    const raw_name = atom.nameSlice();
    const sym = blk: {
        // For added H in deuterium mode, output "D"
        if (atom.is_added and atom.is_hydrogen and output_isotope == .deuterium) {
            break :blk "D";
        }
        break :blk elementSymbol(atom.element_type);
    };
    const name4 = formatAtomNamePdb(raw_name, sym);
    try writer.writeAll(&name4);

    // col 17: altloc
    try writer.writeByte(if (atom.altloc == ' ') ' ' else atom.altloc);

    // cols 18-20: residue name (right-justified in 3 chars)
    const comp = res.compIdSlice();
    switch (comp.len) {
        0 => try writer.writeAll("   "),
        1 => { try writer.writeAll("  "); try writer.writeByte(comp[0]); },
        2 => { try writer.writeByte(' ');  try writer.writeAll(comp); },
        else => try writer.writeAll(comp[0..3]),
    }

    // col 21: space
    try writer.writeByte(' ');

    // col 22: chain ID
    try writer.writeByte(chain_id);

    // cols 23-26: seq number (right-justified, 4 chars)
    try writer.print("{d:>4}", .{res.auth_seq_id});

    // col 27: insertion code
    try writer.writeByte(if (res.ins_code == 0 or res.ins_code == ' ') ' ' else res.ins_code);

    // cols 28-30: spaces
    try writer.writeAll("   ");

    // cols 31-38: X (8.3f)
    try writeCoord(writer, atom.pos.x);

    // cols 39-46: Y (8.3f)
    try writeCoord(writer, atom.pos.y);

    // cols 47-54: Z (8.3f)
    try writeCoord(writer, atom.pos.z);

    // cols 55-60: occupancy (6.2f)
    try writeOccupancy(writer, atom.occupancy);

    // cols 61-66: B-factor (6.2f)
    try writeOccupancy(writer, atom.b_factor);

    // cols 67-76: spaces (10 chars)
    try writer.writeAll("          ");

    // cols 77-78: element symbol (right-justified in 2 chars)
    if (sym.len == 1) {
        try writer.writeByte(' ');
        try writer.writeAll(sym);
    } else {
        // 2-char symbols: write up to 2 chars
        try writer.writeAll(sym[0..@min(sym.len, 2)]);
    }

    try writer.writeByte('\n');
}

/// Write a TER record.
fn writeTerLine(
    writer: anytype,
    serial: u32,
    res: Residue,
    chain_id: u8,
) !void {
    // cols 1-6: "TER   "
    try writer.writeAll("TER   ");

    // cols 7-11: serial (right-justified, 5 chars)
    try writer.print("{d:>5}", .{serial});

    // cols 12-17: spaces (6 chars)
    try writer.writeAll("      ");

    // cols 18-20: residue name
    const comp = res.compIdSlice();
    switch (comp.len) {
        0 => try writer.writeAll("   "),
        1 => { try writer.writeAll("  "); try writer.writeByte(comp[0]); },
        2 => { try writer.writeByte(' ');  try writer.writeAll(comp); },
        else => try writer.writeAll(comp[0..3]),
    }

    // col 21: space
    try writer.writeByte(' ');

    // col 22: chain ID
    try writer.writeByte(chain_id);

    // cols 23-26: seq number
    try writer.print("{d:>4}", .{res.auth_seq_id});

    // col 27: insertion code
    try writer.writeByte(if (res.ins_code == 0 or res.ins_code == ' ') ' ' else res.ins_code);

    try writer.writeByte('\n');
}

// ── Chain ID extraction ───────────────────────────────────────────────────────

fn chainId(chain: Chain) u8 {
    const s = chain.authSlice();
    if (s.len == 0) return ' ';
    return s[0];
}

// ── HETATM predicate for original atoms ──────────────────────────────────────

fn isHetatm(res: Residue) bool {
    return switch (res.entity_type) {
        .non_polymer, .water, .unknown => true,
        .polymer, .branched => false,
    };
}

// ── Pre-built added-H index ───────────────────────────────────────────────────

const AddedHIndex = struct {
    offsets: []u32, // length n_residues + 1
    indices: []u32, // atom indices of added H atoms
    allocator: std.mem.Allocator,

    fn build(model: *const Model) !AddedHIndex {
        const allocator = model.allocator;
        const n_res = model.residues.items.len;

        const counts = try allocator.alloc(u32, n_res + 1);
        errdefer allocator.free(counts);
        @memset(counts, 0);

        for (model.atoms.items) |a| {
            if (a.is_added) counts[a.residue_idx] += 1;
        }

        const offsets = try allocator.alloc(u32, n_res + 1);
        errdefer allocator.free(offsets);
        offsets[0] = 0;
        for (0..n_res) |r| offsets[r + 1] = offsets[r] + counts[r];

        const total = offsets[n_res];
        const indices = try allocator.alloc(u32, total);
        errdefer allocator.free(indices);

        @memset(counts, 0); // reuse as cursor
        for (model.atoms.items, 0..) |a, ai| {
            if (!a.is_added) continue;
            const r = a.residue_idx;
            indices[offsets[r] + counts[r]] = @intCast(ai);
            counts[r] += 1;
        }

        allocator.free(counts);
        return .{ .offsets = offsets, .indices = indices, .allocator = allocator };
    }

    fn deinit(self: *AddedHIndex) void {
        self.allocator.free(self.offsets);
        self.allocator.free(self.indices);
    }

    fn slice(self: *const AddedHIndex, res_idx: usize) []const u32 {
        return self.indices[self.offsets[res_idx]..self.offsets[res_idx + 1]];
    }
};

// ── Public API ────────────────────────────────────────────────────────────────

/// Write a Model to PDB format, interleaving passthrough records from the
/// original parse and inserting added hydrogen atoms residue-by-residue.
///
/// Write multiple models with MODEL/ENDMDL records.
/// Header records are written first, then each model is wrapped with MODEL/ENDMDL.
pub fn writeMultiModel(
    writer: anytype,
    entries: []const pdb_mod.PdbModelEntry,
    header_records: []const PdbRecord,
    output_isotope: place.OutputIsotope,
) !void {
    // Write header records (HEADER, CRYST1, etc.)
    for (header_records) |rec| {
        switch (rec) {
            .raw_line => |line| {
                // Skip END record in header — we'll write our own at the end
                if (line.len >= 3 and std.mem.eql(u8, std.mem.trimEnd(u8, line, " "), "END")) continue;
                try writer.print("{s}\n", .{line});
            },
            .atom_site => {},
        }
    }

    const multi = entries.len > 1;
    for (entries) |entry| {
        if (multi) {
            try writer.print("MODEL     {d:>4}\n", .{entry.model_num});
        }
        // Use writeModel for each entry's per-model records (atoms + TER).
        try writeModel(writer, &entry.model, entry.records.items, output_isotope);
        if (multi) {
            try writer.writeAll("ENDMDL\n");
        }
    }

    try writer.writeAll("END\n");
}

/// `records` must match the original parse order: `atom_site` tags mark where
/// coordinate atoms were parsed; `raw_line` entries are written verbatim.
/// TER records in `records` are replaced with fresh TER lines at chain boundaries.
/// An END record is appended at the end if one was present in the passthrough.
pub fn writeModel(
    writer: anytype,
    model: *const Model,
    records: []const PdbRecord,
    output_isotope: place.OutputIsotope,
) !void {
    var serial: u32 = 1;

    // Build added-H index for O(1) per-residue lookup.
    var h_idx = try AddedHIndex.build(model);
    defer h_idx.deinit();

    // We walk through the records list. For each run of consecutive atom_site
    // markers we need to emit one complete residue group (original heavy atoms
    // + added H). We track the original atom index (for original atoms) and
    // the residue index implied by each atom.
    var orig_atom_cursor: u32 = 0; // next original atom to output
    var last_chain_idx: u32 = std.math.maxInt(u32); // sentinel: no chain yet

    // Track the last residue we visited, for TER line data.
    var last_res_idx: u32 = 0;

    // Scan for a contiguous block of atom_site records representing one residue.
    // We iterate records; on each atom_site we look at model.atoms[orig_atom_cursor]
    // to know which residue we're in, then output all atoms for that residue.

    // We need to detect residue boundaries within atom_site runs. Because the
    // parser guarantees atoms are in order and residue.atom_start/atom_end are set,
    // we can track when orig_atom_cursor crosses into a new residue.

    var rec_idx: usize = 0;
    while (rec_idx < records.len) : (rec_idx += 1) {
        const rec = records[rec_idx];
        switch (rec) {
            .raw_line => |line| {
                // Skip original TER lines — we regenerate them ourselves.
                if (line.len >= 3 and std.mem.eql(u8, line[0..3], "TER")) continue;
                // Write other raw lines verbatim.
                try writer.writeAll(line);
                try writer.writeByte('\n');
            },
            .atom_site => {
                // This atom_site corresponds to model.atoms[orig_atom_cursor].
                // If we're at the start of a new residue, output the full residue
                // group (original atoms + added H). If we're mid-residue we already
                // output everything; skip.

                if (orig_atom_cursor >= model.original_atom_count) {
                    // No more original atoms to emit (shouldn't happen in a well-formed
                    // parse result, but guard anyway).
                    orig_atom_cursor += 1;
                    continue;
                }

                const atom = model.atoms.items[orig_atom_cursor];
                const res_idx = atom.residue_idx;
                const res = model.residues.items[res_idx];
                const chain = model.chains.items[res.chain_idx];
                const cid = chainId(chain);

                // Check if we just entered a new chain → emit TER for the previous chain.
                if (res.chain_idx != last_chain_idx and last_chain_idx != std.math.maxInt(u32)) {
                    // Emit TER for the last residue of the previous chain.
                    const prev_chain = model.chains.items[last_chain_idx];
                    if (prev_chain.residue_end > prev_chain.residue_start) {
                        const last_res = model.residues.items[prev_chain.residue_end - 1];
                        const prev_cid = chainId(prev_chain);
                        try writeTerLine(writer, serial, last_res, prev_cid);
                        serial += 1;
                    }
                }
                last_chain_idx = res.chain_idx;
                last_res_idx = res_idx;

                // Only emit this residue group when we hit its first atom.
                // Subsequent atom_site records for the same residue are skipped.
                if (orig_atom_cursor == res.atom_start) {
                    const hetatm = isHetatm(res);

                    // Original heavy atoms
                    for (res.atom_start..res.atom_end) |ai| {
                        const a = model.atoms.items[ai];
                        try writeAtomLine(writer, a, res, cid, serial, hetatm, output_isotope);
                        serial += 1;
                    }

                    // Added H atoms for this residue
                    for (h_idx.slice(res_idx)) |ai| {
                        const a = model.atoms.items[ai];
                        if (mover_mod.isAbsentH(a)) continue;
                        // Added H atoms always use ATOM record type (not HETATM)
                        try writeAtomLine(writer, a, res, cid, serial, false, output_isotope);
                        serial += 1;
                    }
                }

                orig_atom_cursor += 1;
            },
        }
    }

    // After all records, emit TER for the final chain if there were any atoms.
    if (last_chain_idx != std.math.maxInt(u32)) {
        const last_chain = model.chains.items[last_chain_idx];
        if (last_chain.residue_end > last_chain.residue_start) {
            const last_res = model.residues.items[last_chain.residue_end - 1];
            const cid = chainId(last_chain);
            try writeTerLine(writer, serial, last_res, cid);
            serial += 1;
        }
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

const testing = std.testing;

test "formatAtomNamePdb" {
    // 1-char element, short name -> starts at col 14
    try testing.expectEqualStrings(" CA ", &formatAtomNamePdb("CA", "C"));
    try testing.expectEqualStrings(" N  ", &formatAtomNamePdb("N", "N"));
    // 2-char element -> starts at col 13
    try testing.expectEqualStrings("FE  ", &formatAtomNamePdb("FE", "Fe"));
    // 4-char name -> fills all columns
    try testing.expectEqualStrings("HG11", &formatAtomNamePdb("HG11", "H"));
    // Digit-prefixed hydrogen
    try testing.expectEqualStrings("1HB ", &formatAtomNamePdb("1HB", "H"));
}

test "write tiny PDB round-trip" {
    const pdb_parse_mod = @import("../pdb.zig");
    const source = @embedFile("../test_data/tiny.pdb");
    var result = try pdb_parse_mod.parse(testing.allocator, source);
    defer result.deinit(testing.allocator);

    var buf: std.Io.Writer.Allocating = .init(testing.allocator);
    defer buf.deinit();
    const w = &buf.writer;
    try writeModel(w, &result.model, result.records.items, .hydrogen);
    const output = buf.writer.buffered();

    // Should contain HEADER passthrough
    try testing.expect(std.mem.find(u8, output, "HEADER") != null);

    // Should contain 5 ATOM lines
    var atom_count: usize = 0;
    var line_it = std.mem.splitScalar(u8, output, '\n');
    while (line_it.next()) |line| {
        if (line.len >= 4 and std.mem.eql(u8, line[0..4], "ATOM")) atom_count += 1;
    }
    try testing.expectEqual(@as(usize, 5), atom_count);

    // Should end with END
    try testing.expect(std.mem.find(u8, output, "END") != null);
}

test "PDB writer preserves coordinates" {
    const pdb_parse_mod = @import("../pdb.zig");
    const source = @embedFile("../test_data/tiny.pdb");
    var result = try pdb_parse_mod.parse(testing.allocator, source);
    defer result.deinit(testing.allocator);

    var buf: std.Io.Writer.Allocating = .init(testing.allocator);
    defer buf.deinit();
    const w = &buf.writer;
    try writeModel(w, &result.model, result.records.items, .hydrogen);

    // Re-parse the output
    var mdl2 = try pdb_parse_mod.parseModel(testing.allocator, buf.writer.buffered());
    defer mdl2.deinit();

    try testing.expectEqual(result.model.atoms.items.len, mdl2.atoms.items.len);
    for (result.model.atoms.items, mdl2.atoms.items) |a1, a2| {
        try testing.expectApproxEqAbs(a1.pos.x, a2.pos.x, 1e-3);
        try testing.expectApproxEqAbs(a1.pos.y, a2.pos.y, 1e-3);
        try testing.expectApproxEqAbs(a1.pos.z, a2.pos.z, 1e-3);
    }
}

test "PDB writer multi-chain with TER records" {
    const pdb_parse_mod = @import("../pdb.zig");
    const source = @embedFile("../test_data/multi_chain.pdb");
    var result = try pdb_parse_mod.parse(testing.allocator, source);
    defer result.deinit(testing.allocator);

    var buf: std.Io.Writer.Allocating = .init(testing.allocator);
    defer buf.deinit();
    const w = &buf.writer;
    try writeModel(w, &result.model, result.records.items, .hydrogen);
    const output = buf.writer.buffered();

    // Count ATOM lines: 8 chain A + 4 chain B = 12
    var atom_count: usize = 0;
    var ter_count: usize = 0;
    var line_it = std.mem.splitScalar(u8, output, '\n');
    while (line_it.next()) |line| {
        if (line.len >= 4 and std.mem.eql(u8, line[0..4], "ATOM")) atom_count += 1;
        if (line.len >= 3 and std.mem.eql(u8, line[0..3], "TER")) ter_count += 1;
    }
    try testing.expectEqual(@as(usize, 12), atom_count);
    try testing.expectEqual(@as(usize, 2), ter_count);
}

test "PDB writer HETATM records" {
    const pdb_parse_mod = @import("../pdb.zig");
    const source = @embedFile("../test_data/hetatm.pdb");
    var result = try pdb_parse_mod.parse(testing.allocator, source);
    defer result.deinit(testing.allocator);

    var buf: std.Io.Writer.Allocating = .init(testing.allocator);
    defer buf.deinit();
    const w = &buf.writer;
    try writeModel(w, &result.model, result.records.items, .hydrogen);
    const output = buf.writer.buffered();

    var atom_count: usize = 0;
    var hetatm_count: usize = 0;
    var line_it = std.mem.splitScalar(u8, output, '\n');
    while (line_it.next()) |line| {
        if (line.len >= 4 and std.mem.eql(u8, line[0..4], "ATOM")) atom_count += 1;
        if (line.len >= 6 and std.mem.eql(u8, line[0..6], "HETATM")) hetatm_count += 1;
    }
    // 5 polymer ATOM + 4 non_polymer/water HETATM
    try testing.expectEqual(@as(usize, 5), atom_count);
    try testing.expectEqual(@as(usize, 4), hetatm_count);
}
