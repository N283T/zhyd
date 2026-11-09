//! mmCIF writer: outputs model atoms plus a custom _zhyd_log category.

const std = @import("std");
const model_mod = @import("../model.zig");
const Model = model_mod.Model;
const Atom = model_mod.Atom;
const Residue = model_mod.Residue;
const Chain = model_mod.Chain;
const mover_mod = @import("../optimize/mover.zig");
const element = @import("../element.zig");
const AtomType = element.AtomType;
const cif = @import("../cif.zig");
const mmcif = @import("../mmcif.zig");
const place = @import("../place.zig");
const format = @import("format.zig");

const elementSymbol = format.elementSymbol;
const writeFixedFloat3 = format.writeFixedFloat3;
const writeFixedFloat2 = format.writeFixedFloat2;
const writeAtomName = format.writeAtomName;
const writeAltId = format.writeAltId;
const writeCifValue = format.writeCifValue;
const writeCifValueInLoop = format.writeCifValueInLoop;
const writePairCifValue = format.writePairCifValue;
const writePairValue = format.writePairValue;

const cif_types = cif.types;
const Document = cif_types.Document;
const Block = cif_types.Block;
const Loop = cif_types.Loop;
const Pair = cif_types.Pair;
const Item = cif_types.Item;

// ── Shared types ──────────────────────────────────────────────────────────────

/// Maps _atom_site column tags to their indices in a Loop.
const ColMap = struct {
    group_pdb: ?usize = null,
    id: ?usize = null,
    type_symbol: ?usize = null,
    label_atom_id: ?usize = null,
    label_alt_id: ?usize = null,
    label_comp_id: ?usize = null,
    label_asym_id: ?usize = null,
    label_entity_id: ?usize = null,
    label_seq_id: ?usize = null,
    cartn_x: ?usize = null,
    cartn_y: ?usize = null,
    cartn_z: ?usize = null,
    occupancy: ?usize = null,
    b_factor: ?usize = null,
    auth_seq_id: ?usize = null,
    auth_comp_id: ?usize = null,
    auth_asym_id: ?usize = null,
    auth_atom_id: ?usize = null,
    pdb_model_num: ?usize = null,
};

/// Build a ColMap from the tags of an _atom_site loop.
fn buildColMap(orig_loop: *const Loop) ColMap {
    var cm = ColMap{};
    for (orig_loop.tags.items, 0..) |tag, i| {
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.group_PDB")) cm.group_pdb = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.id")) cm.id = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.type_symbol")) cm.type_symbol = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.label_atom_id")) cm.label_atom_id = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.label_alt_id")) cm.label_alt_id = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.label_comp_id")) cm.label_comp_id = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.label_asym_id")) cm.label_asym_id = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.label_entity_id")) cm.label_entity_id = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.label_seq_id")) cm.label_seq_id = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.Cartn_x")) cm.cartn_x = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.Cartn_y")) cm.cartn_y = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.Cartn_z")) cm.cartn_z = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.occupancy")) cm.occupancy = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.B_iso_or_equiv")) cm.b_factor = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.auth_seq_id")) cm.auth_seq_id = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.auth_comp_id")) cm.auth_comp_id = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.auth_asym_id")) cm.auth_asym_id = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.auth_atom_id")) cm.auth_atom_id = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.pdbx_PDB_model_num")) cm.pdb_model_num = i;
    }
    return cm;
}

/// Result of buildAddedHIndex — caller must call deinit to free slices.
const AddedHIndex = struct {
    offsets: []u32,
    indices: []u32,

    pub fn deinit(self: AddedHIndex, allocator: std.mem.Allocator) void {
        allocator.free(self.offsets);
        allocator.free(self.indices);
    }
};

/// Build a counting-sort index of added-H atoms by residue.
/// `offsets[r]..offsets[r+1]` gives the slice of `indices` for residue r.
fn buildAddedHIndex(model: *const Model) !AddedHIndex {
    const allocator = model.allocator;
    const n_residues = model.residues.items.len;

    // Phase 1: count added atoms per residue.
    const counts = try allocator.alloc(u32, n_residues + 1);
    defer allocator.free(counts);
    @memset(counts, 0);
    for (model.atoms.items) |atom| {
        if (!atom.is_added) continue;
        counts[atom.residue_idx] += 1;
    }

    // Phase 2: build prefix-sum offsets (length n_residues + 1).
    const offsets = try allocator.alloc(u32, n_residues + 1);
    errdefer allocator.free(offsets);
    offsets[0] = 0;
    for (0..n_residues) |r| {
        offsets[r + 1] = offsets[r] + counts[r];
    }
    const total_added = offsets[n_residues];

    // Phase 3: fill sorted atom indices (reuse counts as cursor).
    const indices = try allocator.alloc(u32, total_added);
    @memset(counts, 0);
    for (model.atoms.items, 0..) |atom, atom_idx| {
        if (!atom.is_added) continue;
        const r = atom.residue_idx;
        indices[offsets[r] + counts[r]] = @intCast(atom_idx);
        counts[r] += 1;
    }

    return .{ .offsets = offsets, .indices = indices };
}

/// Write one added-H atom row using the preserved column layout.
/// `first_orig_row` is the loop row index of the residue's first heavy atom,
/// used to inherit entity_id / auth_seq_id / unknown columns from the original.
fn writeAddedHRow(
    writer: anytype,
    model: *const Model,
    atom: Atom,
    res: Residue,
    serial: u32,
    cm: ColMap,
    col_widths: []const usize,
    orig_loop: *const Loop,
    first_orig_row: usize,
    output_isotope: place.OutputIsotope,
) !void {
    const w = col_widths.len;
    const chain = model.chains.items[res.chain_idx];
    for (0..w) |col| {
        if (col > 0) try writer.writeByte(' ');
        const cw = col_widths[col];
        if (cm.group_pdb != null and col == cm.group_pdb.?) {
            // Same record type as the residue's own atoms (HETATM for ligands).
            if (first_orig_row < orig_loop.length()) {
                try writePaddedCifValue(writer, orig_loop.val(first_orig_row, col) orelse "ATOM", cw);
            } else {
                try writePaddedCell(writer, "ATOM", cw);
            }
        } else if (cm.id != null and col == cm.id.?) {
            try writePaddedInt(writer, serial, cw);
        } else if (cm.type_symbol != null and col == cm.type_symbol.?) {
            try writePaddedCell(writer, atomTypeSymbol(atom, output_isotope), cw);
        } else if (cm.label_atom_id != null and col == cm.label_atom_id.?) {
            try writePaddedCell(writer, atom.nameSlice(), cw);
        } else if (cm.label_alt_id != null and col == cm.label_alt_id.?) {
            if (atom.altloc == ' ') {
                try writePaddedCell(writer, ".", cw);
            } else {
                try writePaddedCell(writer, &[_]u8{atom.altloc}, cw);
            }
        } else if (cm.label_comp_id != null and col == cm.label_comp_id.?) {
            try writePaddedCell(writer, res.compIdSlice(), cw);
        } else if (cm.label_asym_id != null and col == cm.label_asym_id.?) {
            try writePaddedCell(writer, chain.labelSlice(), cw);
        } else if (cm.label_entity_id != null and col == cm.label_entity_id.?) {
            if (first_orig_row < orig_loop.length()) {
                try writePaddedCifValue(writer, orig_loop.val(first_orig_row, col) orelse ".", cw);
            } else {
                try writePaddedCell(writer, ".", cw);
            }
        } else if (cm.label_seq_id != null and col == cm.label_seq_id.?) {
            // Copied from the residue: non-polymer residues have "." here.
            if (first_orig_row < orig_loop.length()) {
                try writePaddedCifValue(writer, orig_loop.val(first_orig_row, col) orelse ".", cw);
            } else {
                try writePaddedInt(writer, res.seq_id, cw);
            }
        } else if (cm.cartn_x != null and col == cm.cartn_x.?) {
            try writePaddedFloat3(writer, atom.pos.x, cw);
        } else if (cm.cartn_y != null and col == cm.cartn_y.?) {
            try writePaddedFloat3(writer, atom.pos.y, cw);
        } else if (cm.cartn_z != null and col == cm.cartn_z.?) {
            try writePaddedFloat3(writer, atom.pos.z, cw);
        } else if (cm.occupancy != null and col == cm.occupancy.?) {
            try writePaddedFloat2(writer, atom.occupancy, cw);
        } else if (cm.b_factor != null and col == cm.b_factor.?) {
            try writePaddedFloat2(writer, atom.b_factor, cw);
        } else if (cm.auth_seq_id != null and col == cm.auth_seq_id.?) {
            if (first_orig_row < orig_loop.length()) {
                try writePaddedCifValue(writer, orig_loop.val(first_orig_row, col) orelse ".", cw);
            } else {
                try writePaddedInt(writer, res.seq_id, cw);
            }
        } else if (cm.auth_comp_id != null and col == cm.auth_comp_id.?) {
            try writePaddedCell(writer, res.compIdSlice(), cw);
        } else if (cm.auth_asym_id != null and col == cm.auth_asym_id.?) {
            try writePaddedCell(writer, chain.authSlice(), cw);
        } else if (cm.auth_atom_id != null and col == cm.auth_atom_id.?) {
            try writePaddedCell(writer, atom.nameSlice(), cw);
        } else if (cm.pdb_model_num != null and col == cm.pdb_model_num.?) {
            try writePaddedInt(writer, model.model_num, cw);
        } else {
            if (first_orig_row < orig_loop.length()) {
                try writePaddedCifValue(writer, orig_loop.val(first_orig_row, col) orelse ".", cw);
            } else {
                try writePaddedCell(writer, ".", cw);
            }
        }
    }
    try writer.writeByte('\n');
}

// ── Public API ────────────────────────────────────────────────────────────────

/// Write a model as mmCIF, preserving all original categories from the source document.
/// The _atom_site loop is replaced with the model's atoms (including added H).
/// If doc is null, writes atom_site only (atom-site-only mode).
pub fn writeWithDocument(writer: anytype, model: *const Model, doc: ?*const Document) !void {
    return writeWithDocumentWithPolicy(writer, model, doc, .{});
}

/// Write a model as mmCIF with explicit bond/output policy.
pub fn writeWithDocumentWithPolicy(writer: anytype, model: *const Model, doc: ?*const Document, bond_policy: place.BondPolicy) !void {
    if (doc) |d| {
        if (d.blocks.items.len > 0) {
            const block = &d.blocks.items[0];
            try writer.print("data_{s}\n", .{block.name});

            // Write all items, replacing _atom_site loop.
            // Each item is preceded by a single '#' separator line.
            var atom_site_written = false;
            for (block.items.items) |item| {
                try writer.writeAll("#\n");
                switch (item) {
                    .pair => |p| {
                        try writePairValue(writer, p.tag, p.value);
                    },
                    .loop => |l| {
                        if (isAtomSiteLoop(&l)) {
                            try writeAtomSitePreserving(writer, model, &l, bond_policy.output_isotope);
                            atom_site_written = true;
                        } else {
                            try writeLoopBody(writer, &l);
                        }
                    },
                }
            }
            try writer.writeAll("#\n");

            // If no atom_site loop was found in original, append it
            if (!atom_site_written) {
                try writeAtomSite(writer, model, bond_policy.output_isotope);
                try writer.writeAll("#\n");
            }
            return;
        }
    }

    // Fallback: atom-site-only mode
    try writer.writeAll("data_ZHYD\n#\n");
    try writeAtomSite(writer, model, bond_policy.output_isotope);
    try writer.writeAll("#\n");
}

/// Write atom-site-only output (no original document preservation).
pub fn write(writer: anytype, model: *const Model, block_name: []const u8) !void {
    return writeWithPolicy(writer, model, block_name, .{});
}

/// Write atom-site-only output with explicit bond/output policy.
pub fn writeWithPolicy(writer: anytype, model: *const Model, block_name: []const u8, bond_policy: place.BondPolicy) !void {
    try writer.print("data_{s}\n#\n", .{block_name});
    try writeAtomSite(writer, model, bond_policy.output_isotope);
    try writer.writeAll("#\n");
}

fn isAtomSiteLoop(loop: *const Loop) bool {
    for (loop.tags.items) |tag| {
        if (std.ascii.startsWithIgnoreCase(tag, "_atom_site.")) return true;
    }
    return false;
}

/// Write a multi-model mmCIF, preserving all original categories from the source document.
/// All models are written into a single _atom_site loop with correct model numbers.
pub fn writeMultiModelWithDocumentWithPolicy(
    writer: anytype,
    entries: []const mmcif.ModelEntry,
    doc: ?*const Document,
    bond_policy: place.BondPolicy,
) !void {
    if (doc) |d| {
        if (d.blocks.items.len > 0) {
            const block = &d.blocks.items[0];
            try writer.print("data_{s}\n", .{block.name});

            var atom_site_written = false;
            for (block.items.items) |item| {
                try writer.writeAll("#\n");
                switch (item) {
                    .pair => |p| {
                        try writePairValue(writer, p.tag, p.value);
                    },
                    .loop => |l| {
                        if (isAtomSiteLoop(&l)) {
                            try writeMultiModelAtomSitePreserving(writer, entries, &l, bond_policy.output_isotope);
                            atom_site_written = true;
                        } else {
                            try writeLoopBody(writer, &l);
                        }
                    },
                }
            }
            try writer.writeAll("#\n");

            if (!atom_site_written) {
                // Fallback: write without preservation
                for (entries) |entry| {
                    try writeAtomSite(writer, &entry.model, bond_policy.output_isotope);
                }
                try writer.writeAll("#\n");
            }
            return;
        }
    }

    // Fallback: atom-site-only mode
    try writer.writeAll("data_ZHYD\n#\n");
    for (entries) |entry| {
        try writeAtomSite(writer, &entry.model, bond_policy.output_isotope);
    }
    try writer.writeAll("#\n");
}

// ── Internal writers ──────────────────────────────────────────────────────────

/// Write _atom_site loop preserving the original column structure.
/// Heavy atom rows are written from the original loop data.
/// Added H atom rows fill known columns and use '.' for the rest.
fn writeAtomSitePreserving(writer: anytype, model: *const Model, orig_loop: *const Loop, output_isotope: place.OutputIsotope) !void {
    const w = orig_loop.width();
    if (w == 0) return;

    // Write loop header with original tags
    try writer.writeAll("loop_\n");
    for (orig_loop.tags.items) |tag| {
        try writer.print("{s}\n", .{tag});
    }

    const cm = buildColMap(orig_loop);

    // Compute column widths from original data for PDBx-styled alignment.
    const col_widths = try model.allocator.alloc(usize, w);
    defer model.allocator.free(col_widths);
    @memset(col_widths, 0);
    for (0..orig_loop.length()) |row| {
        for (0..w) |col| {
            const val = orig_loop.val(row, col) orelse ".";
            if (val.len > col_widths[col]) col_widths[col] = val.len;
        }
    }
    // Ensure minimum widths for columns we generate (serial, coords, etc.)
    // Serial number width: estimate from total atom count.
    const total_atoms = model.atoms.items.len;
    const serial_width = blk: {
        var digits: usize = 1;
        var n = total_atoms;
        while (n >= 10) : (n /= 10) digits += 1;
        break :blk digits;
    };
    if (cm.id) |c| col_widths[c] = @max(col_widths[c], serial_width);
    // Coordinate columns: up to 8 chars (e.g. "-123.456")
    if (cm.cartn_x) |c| col_widths[c] = @max(col_widths[c], 8);
    if (cm.cartn_y) |c| col_widths[c] = @max(col_widths[c], 8);
    if (cm.cartn_z) |c| col_widths[c] = @max(col_widths[c], 8);
    // Occupancy/B-factor: up to 6 chars (e.g. "100.00")
    if (cm.occupancy) |c| col_widths[c] = @max(col_widths[c], 6);
    if (cm.b_factor) |c| col_widths[c] = @max(col_widths[c], 6);
    // Account for added H atom names which may be wider than original names
    // (e.g. "HG11" = 4 chars vs original max of 3 chars like "CG1").
    if (cm.label_atom_id != null or cm.auth_atom_id != null) {
        for (model.atoms.items) |atom| {
            if (!atom.is_added) continue;
            const name_len = atom.nameSlice().len;
            if (cm.label_atom_id) |c| col_widths[c] = @max(col_widths[c], name_len);
            if (cm.auth_atom_id) |c| col_widths[c] = @max(col_widths[c], name_len);
        }
    }

    // Pre-index added H atoms by residue_idx to avoid O(N*R) inner scan.
    const h_index = try buildAddedHIndex(model);
    defer h_index.deinit(model.allocator);

    // Track serial for renumbering
    var serial: u32 = 1;

    // Write rows grouped by residue: original heavy atoms then added H
    for (model.residues.items, 0..) |res, res_idx| {
        // Original heavy atoms: write rows from original loop data,
        // but use model atom coordinates (which may have been updated by
        // optimizers, e.g. amide flip swapping O/N positions).
        for (res.atom_start..res.atom_end) |orig_row_idx| {
            if (orig_row_idx < orig_loop.length()) {
                const atom = model.atoms.items[orig_row_idx];
                for (0..w) |col| {
                    if (col > 0) try writer.writeByte(' ');
                    if (cm.id != null and col == cm.id.?) {
                        try writePaddedInt(writer, serial, col_widths[col]);
                    } else if (cm.cartn_x != null and col == cm.cartn_x.?) {
                        try writePaddedFloat3(writer, atom.pos.x, col_widths[col]);
                    } else if (cm.cartn_y != null and col == cm.cartn_y.?) {
                        try writePaddedFloat3(writer, atom.pos.y, col_widths[col]);
                    } else if (cm.cartn_z != null and col == cm.cartn_z.?) {
                        try writePaddedFloat3(writer, atom.pos.z, col_widths[col]);
                    } else {
                        const val = orig_loop.val(orig_row_idx, col) orelse ".";
                        try writePaddedCifValue(writer, val, col_widths[col]);
                    }
                }
                try writer.writeByte('\n');
                serial += 1;
            }
        }

        // Added H atoms for this residue — look up via pre-built index.
        const h_start = h_index.offsets[res_idx];
        const h_end = h_index.offsets[res_idx + 1];
        for (h_index.indices[h_start..h_end]) |atom_idx| {
            const atom = model.atoms.items[atom_idx];
            // Skip absent H atoms (flipper sentinel)
            if (mover_mod.isAbsentH(atom)) continue;

            const first_orig_row = res.atom_start;
            try writeAddedHRow(writer, model, atom, res, serial, cm, col_widths, orig_loop, first_orig_row, output_isotope);
            serial += 1;
        }
    }
}

/// Write _atom_site loop preserving the original column structure for multiple models.
/// All models are written into a single loop with correct model numbers.
fn writeMultiModelAtomSitePreserving(
    writer: anytype,
    entries: []const mmcif.ModelEntry,
    orig_loop: *const Loop,
    output_isotope: place.OutputIsotope,
) !void {
    if (entries.len == 0) return;

    const w = orig_loop.width();
    if (w == 0) return;

    // Use the first entry's model allocator for temporary allocations
    const allocator = entries[0].model.allocator;

    // Write loop header
    try writer.writeAll("loop_\n");
    for (orig_loop.tags.items) |tag| {
        try writer.print("{s}\n", .{tag});
    }

    const cm = buildColMap(orig_loop);

    // Compute column widths from original data
    const col_widths = try allocator.alloc(usize, w);
    defer allocator.free(col_widths);
    @memset(col_widths, 0);
    for (0..orig_loop.length()) |row| {
        for (0..w) |col| {
            const val = orig_loop.val(row, col) orelse ".";
            if (val.len > col_widths[col]) col_widths[col] = val.len;
        }
    }

    // Estimate total atoms across all models for serial width
    var total_atoms: usize = 0;
    for (entries) |entry| total_atoms += entry.model.atoms.items.len;
    const serial_width = blk: {
        var digits: usize = 1;
        var n = total_atoms;
        while (n >= 10) : (n /= 10) digits += 1;
        break :blk digits;
    };
    if (cm.id) |c| col_widths[c] = @max(col_widths[c], serial_width);
    if (cm.cartn_x) |c| col_widths[c] = @max(col_widths[c], 8);
    if (cm.cartn_y) |c| col_widths[c] = @max(col_widths[c], 8);
    if (cm.cartn_z) |c| col_widths[c] = @max(col_widths[c], 8);
    if (cm.occupancy) |c| col_widths[c] = @max(col_widths[c], 6);
    if (cm.b_factor) |c| col_widths[c] = @max(col_widths[c], 6);
    if (cm.label_atom_id != null or cm.auth_atom_id != null) {
        for (entries) |entry| {
            for (entry.model.atoms.items) |atom| {
                if (!atom.is_added) continue;
                const name_len = atom.nameSlice().len;
                if (cm.label_atom_id) |c| col_widths[c] = @max(col_widths[c], name_len);
                if (cm.auth_atom_id) |c| col_widths[c] = @max(col_widths[c], name_len);
            }
        }
    }

    var serial: u32 = 1;

    // Write each model's atoms
    for (entries) |entry| {
        const model = &entry.model;
        const cif_row_offset = entry.cif_row_start;

        // Pre-index added H atoms for this model
        const h_index = try buildAddedHIndex(model);
        defer h_index.deinit(model.allocator);

        // Write rows grouped by residue
        for (model.residues.items, 0..) |res, res_idx| {
            // Original heavy atoms — use cif_row_offset for orig_loop access
            for (res.atom_start..res.atom_end) |local_idx| {
                const orig_row_idx = cif_row_offset + @as(u32, @intCast(local_idx));
                if (orig_row_idx < orig_loop.length()) {
                    const atom = model.atoms.items[local_idx];
                    for (0..w) |col| {
                        if (col > 0) try writer.writeByte(' ');
                        if (cm.id != null and col == cm.id.?) {
                            try writePaddedInt(writer, serial, col_widths[col]);
                        } else if (cm.cartn_x != null and col == cm.cartn_x.?) {
                            try writePaddedFloat3(writer, atom.pos.x, col_widths[col]);
                        } else if (cm.cartn_y != null and col == cm.cartn_y.?) {
                            try writePaddedFloat3(writer, atom.pos.y, col_widths[col]);
                        } else if (cm.cartn_z != null and col == cm.cartn_z.?) {
                            try writePaddedFloat3(writer, atom.pos.z, col_widths[col]);
                        } else {
                            const val = orig_loop.val(orig_row_idx, col) orelse ".";
                            try writePaddedCifValue(writer, val, col_widths[col]);
                        }
                    }
                    try writer.writeByte('\n');
                    serial += 1;
                }
            }

            // Added H atoms for this residue
            const h_start = h_index.offsets[res_idx];
            const h_end = h_index.offsets[res_idx + 1];
            for (h_index.indices[h_start..h_end]) |atom_idx| {
                const atom = model.atoms.items[atom_idx];
                if (mover_mod.isAbsentH(atom)) continue;

                const first_orig_row = cif_row_offset + res.atom_start;
                try writeAddedHRow(writer, model, atom, res, serial, cm, col_widths, orig_loop, first_orig_row, output_isotope);
                serial += 1;
            }
        }
    }
}

/// Write a CIF loop body (tags + data rows) without surrounding '#' separators.
/// Preserves PDBx-styled column alignment by computing and applying column widths.
fn writeLoopBody(writer: anytype, loop: *const Loop) !void {
    try writer.writeAll("loop_\n");
    for (loop.tags.items) |tag| {
        try writer.print("{s}\n", .{tag});
    }
    const w = loop.width();
    if (w == 0) return;
    const n_rows = loop.length();
    if (n_rows == 0) return;

    // Compute column widths for alignment.
    var col_widths_buf: [64]usize = undefined;
    const col_widths = col_widths_buf[0..@min(w, 64)];
    @memset(col_widths, 0);
    for (0..n_rows) |row| {
        for (0..col_widths.len) |col| {
            const val = loop.val(row, col) orelse ".";
            if (val.len > col_widths[col]) col_widths[col] = val.len;
        }
    }

    for (0..n_rows) |row| {
        for (0..w) |col| {
            if (col > 0) try writer.writeByte(' ');
            const val = loop.val(row, col) orelse ".";
            if (col < col_widths.len) {
                try writePaddedCifValue(writer, val, col_widths[col]);
            } else {
                try writeCifValueInLoop(writer, val);
            }
        }
        try writer.writeByte('\n');
    }
}

/// Write the _atom_site loop with all atoms.
fn writeAtomSite(writer: anytype, model: *const Model, output_isotope: place.OutputIsotope) !void {
    try writer.writeAll(
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
        \\_atom_site.pdbx_PDB_model_num
        \\
    );

    // Write atoms grouped by residue: heavy atoms first, then added H atoms.
    // The placer appends H atoms at the end of model.atoms, so we need to
    // reorder output so each residue's H atoms follow its heavy atoms.
    const h_index = try buildAddedHIndex(model);
    defer h_index.deinit(model.allocator);

    var serial: u32 = 1;
    for (model.residues.items, 0..) |res, res_idx| {
        // Original heavy atoms in this residue's range
        for (model.atoms.items[res.atom_start..res.atom_end]) |atom| {
            try writeAtomRow(writer, model, atom, res, serial, output_isotope);
            serial += 1;
        }
        // Added H atoms for this residue — look up via pre-built index.
        const h_start = h_index.offsets[res_idx];
        const h_end = h_index.offsets[res_idx + 1];
        for (h_index.indices[h_start..h_end]) |atom_idx| {
            const atom = model.atoms.items[atom_idx];
            // Skip absent H atoms (flipper sentinel position)
            if (mover_mod.isAbsentH(atom)) continue;
            try writeAtomRow(writer, model, atom, res, serial, output_isotope);
            serial += 1;
        }
    }
}

fn writeAtomRow(writer: anytype, model: *const Model, atom: Atom, res: Residue, serial: u32, output_isotope: place.OutputIsotope) !void {
    const chain = model.chains.items[res.chain_idx];
    const group = switch (res.entity_type) {
        .water, .non_polymer, .branched => "HETATM",
        .polymer, .unknown => "ATOM",
    };
    const elem_str = atomTypeSymbol(atom, output_isotope);
    try writer.print(
        "{s} {d} {s} {s} {s} {s} {d} ",
        .{
            group,
            serial,
            elem_str,
            atom.nameSlice(),
            res.compIdSlice(),
            chain.labelSlice(),
            res.seq_id,
        },
    );
    try writeFixedFloat3(writer, atom.pos.x);
    try writer.writeByte(' ');
    try writeFixedFloat3(writer, atom.pos.y);
    try writer.writeByte(' ');
    try writeFixedFloat3(writer, atom.pos.z);
    try writer.writeByte(' ');
    try writeFixedFloat2(writer, atom.occupancy);
    try writer.writeByte(' ');
    try writeFixedFloat2(writer, atom.b_factor);
    try writer.writeByte(' ');
    try writeAltId(writer, atom.altloc);
    try writer.print(" {d}", .{model.model_num});
    try writer.writeByte('\n');
}

fn atomTypeSymbol(atom: Atom, output_isotope: place.OutputIsotope) []const u8 {
    if (output_isotope == .deuterium and atom.is_added and atom.is_hydrogen) return "D";
    return elementSymbol(atom.element_type);
}

// ── Padding helpers ───────────────────────────────────────────────────────────

/// Write a cell value padded to min_width with trailing spaces.
fn writePaddedCell(writer: anytype, val: []const u8, min_width: usize) !void {
    try writer.writeAll(val);
    var i: usize = val.len;
    while (i < min_width) : (i += 1) {
        try writer.writeByte(' ');
    }
}

/// Write a CIF value with quoting, then pad to min_width.
fn writePaddedCifValue(writer: anytype, val: []const u8, min_width: usize) !void {
    var buf: [256]u8 = undefined;
    var w: std.Io.Writer = .fixed(&buf);
    writeCifValueInLoop(&w, val) catch {
        // Value too large for buffer — write directly without padding
        try writeCifValueInLoop(writer, val);
        return;
    };
    try writePaddedCell(writer, w.buffered(), min_width);
}

/// Write an integer value padded to min_width.
fn writePaddedInt(writer: anytype, val: anytype, min_width: usize) !void {
    var buf: [20]u8 = undefined;
    var w: std.Io.Writer = .fixed(&buf);
    w.print("{d}", .{val}) catch unreachable;
    try writePaddedCell(writer, w.buffered(), min_width);
}

/// Write a float3 value padded to min_width.
fn writePaddedFloat3(writer: anytype, val: f32, min_width: usize) !void {
    var buf: [64]u8 = undefined;
    var w: std.Io.Writer = .fixed(&buf);
    format.writeFixedFloat3(&w, val) catch unreachable;
    try writePaddedCell(writer, w.buffered(), min_width);
}

/// Write a float2 value padded to min_width.
fn writePaddedFloat2(writer: anytype, val: f32, min_width: usize) !void {
    var buf: [64]u8 = undefined;
    var w: std.Io.Writer = .fixed(&buf);
    format.writeFixedFloat2(&w, val) catch unreachable;
    try writePaddedCell(writer, w.buffered(), min_width);
}

// ── Optimization log ──────────────────────────────────────────────────────────

/// Write optimization log as a custom mmCIF category.
pub fn writeZhydLog(
    writer: anytype,
    movers: []const mover_mod.Mover,
    residues: []const Residue,
    chains: []const Chain,
) !void {
    try writer.writeAll(
        \\loop_
        \\_zhyd_log.residue_id
        \\_zhyd_log.action
        \\_zhyd_log.orientation
        \\_zhyd_log.score
        \\
    );
    for (movers) |m| {
        const res = residues[m.residue_idx];
        const chain = chains[res.chain_idx];
        const action = switch (m.kind) {
            .single_h_rotator, .nh3_rotator, .nh2_rotator => "rotate",
            .amide_flip => "flip_amide",
            .his_flip => "flip_his",
        };
        try writer.print("{s}.{s}.{d} {s} {d} .\n", .{
            chain.labelSlice(),
            res.compIdSlice(),
            res.seq_id,
            action,
            m.best_orientation,
        });
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

const testing = std.testing;

test "write mmCIF round-trip" {
    var mdl = Model.init(testing.allocator);
    defer mdl.deinit();

    // Add one chain
    var chain = Chain{};
    chain.setLabelAsymId("A");
    try mdl.chains.append(mdl.allocator, chain);

    // Add one residue
    var res = Residue{};
    res.setCompId("ALA");
    res.seq_id = 1;
    res.chain_idx = 0;
    res.atom_start = 0;
    res.atom_end = 2;
    try mdl.residues.append(mdl.allocator, res);

    // Add two atoms: N and CA
    var n_atom = Atom{ .pos = .{ .x = 1.0, .y = 2.0, .z = 3.0 } };
    n_atom.setName("N");
    n_atom.element_type = .N;
    n_atom.residue_idx = 0;
    try mdl.atoms.append(mdl.allocator, n_atom);

    var ca_atom = Atom{ .pos = .{ .x = 2.5, .y = 3.5, .z = 4.5 } };
    ca_atom.setName("CA");
    ca_atom.element_type = .C;
    ca_atom.residue_idx = 0;
    try mdl.atoms.append(mdl.allocator, ca_atom);

    var buf: std.Io.Writer.Allocating = .init(testing.allocator);
    defer buf.deinit();

    try write(&buf.writer, &mdl, "TEST");

    const output = buf.writer.buffered();
    try testing.expect(std.mem.find(u8, output, "data_TEST") != null);
    try testing.expect(std.mem.find(u8, output, "_atom_site.Cartn_x") != null);
    try testing.expect(std.mem.find(u8, output, "ALA") != null);
    try testing.expect(std.mem.find(u8, output, "ATOM") != null);
}

test "writeWithDocument preserves multiline pair and loop values" {
    const source =
        \\data_SAMPLE
        \\_struct.entry_id SAMPLE
        \\_struct.title
        \\;first line
        \\second line
        \\;
        \\#
        \\loop_
        \\_pdbx_data_usage.details
        \\_pdbx_data_usage.id
        \\;loop line 1
        \\loop line 2
        \\;
        \\1
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
        \\ATOM 1 N N ALA A 1 0.0 0.0 0.0 1.0 10.0 .
        \\ATOM 2 C CA ALA A 1 1.5 0.0 0.0 1.0 10.0 .
        \\ATOM 3 C C ALA A 1 2.0 1.4 0.0 1.0 10.0 .
        \\ATOM 4 O O ALA A 1 3.2 1.5 0.0 1.0 10.0 .
        \\ATOM 5 C CB ALA A 1 1.9 -0.8 1.2 1.0 10.0 .
        \\#
    ;

    var doc = try cif.readString(testing.allocator, source);
    defer doc.deinit();

    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();
    _ = try place.addHydrogens(&mdl, null, null);

    var buf: std.Io.Writer.Allocating = .init(testing.allocator);
    defer buf.deinit();
    try writeWithDocument(&buf.writer, &mdl, &doc);

    var reparsed = try cif.readString(testing.allocator, buf.writer.buffered());
    defer reparsed.deinit();

    const block = &reparsed.blocks.items[0];
    try testing.expectEqualStrings("SAMPLE", block.name);
    try testing.expectEqualStrings("first line\nsecond line\n", block.findValue("_struct.title").?);

    const data_usage = block.findLoop("_pdbx_data_usage.details").?;
    try testing.expectEqual(@as(usize, 1), data_usage.length());
    try testing.expectEqualStrings("loop line 1 loop line 2 ", data_usage.val(0, 0).?);

    const atom_site = block.findLoop("_atom_site.Cartn_x").?;
    try testing.expect(atom_site.length() > 5);
}

test "writeWithDocument preserves added hydrogen altloc" {
    const source = @embedFile("../test_data/ala_altloc.cif");

    var doc = try cif.readString(testing.allocator, source);
    defer doc.deinit();

    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();
    place.applyChemistry(&mdl);
    _ = try place.addHydrogens(&mdl, null, null);

    var buf: std.Io.Writer.Allocating = .init(testing.allocator);
    defer buf.deinit();
    try writeWithDocument(&buf.writer, &mdl, &doc);

    var reparsed = try cif.readString(testing.allocator, buf.writer.buffered());
    defer reparsed.deinit();

    const block = &reparsed.blocks.items[0];
    const atom_site = block.findLoop("_atom_site.label_atom_id").?;
    const atom_name_col = atom_site.findTag("_atom_site.label_atom_id").?;
    const alt_id_col = atom_site.findTag("_atom_site.label_alt_id").?;

    var ha_a_count: u32 = 0;
    var ha_b_count: u32 = 0;
    for (0..atom_site.length()) |row| {
        const atom_name = atom_site.val(row, atom_name_col) orelse continue;
        if (!std.mem.eql(u8, atom_name, "HA")) continue;
        const alt_id = atom_site.val(row, alt_id_col) orelse continue;
        if (std.mem.eql(u8, alt_id, "A")) ha_a_count += 1;
        if (std.mem.eql(u8, alt_id, "B")) ha_b_count += 1;
    }

    try testing.expectEqual(@as(u32, 1), ha_a_count);
    try testing.expectEqual(@as(u32, 1), ha_b_count);
}

test "writeWithPolicy outputs D for added hydrogens when isotope is deuterium" {
    const source = @embedFile("../test_data/tiny.cif");

    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();
    _ = try place.addHydrogens(&mdl, null, null);

    var out: std.Io.Writer.Allocating = .init(testing.allocator);
    defer out.deinit();
    try writeWithPolicy(&out.writer, &mdl, "TEST", .{
        .output_isotope = .deuterium,
    });

    var doc = try cif.readString(testing.allocator, out.writer.buffered());
    defer doc.deinit();

    const block = &doc.blocks.items[0];
    const atom_site = block.findLoop("_atom_site.label_atom_id").?;
    const atom_name_col = atom_site.findTag("_atom_site.label_atom_id").?;
    const type_symbol_col = atom_site.findTag("_atom_site.type_symbol").?;

    var found_added_h: bool = false;
    for (0..atom_site.length()) |row| {
        const atom_name = atom_site.val(row, atom_name_col) orelse continue;
        if (!std.mem.startsWith(u8, atom_name, "H")) continue;
        found_added_h = true;
        const type_symbol = atom_site.val(row, type_symbol_col) orelse continue;
        try testing.expectEqualStrings("D", type_symbol);
    }
    try testing.expect(found_added_h);
}

test "writeWithDocumentWithPolicy outputs D in preserving mode" {
    const source = @embedFile("../test_data/tiny.cif");

    var doc = try cif.readString(testing.allocator, source);
    defer doc.deinit();

    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();
    _ = try place.addHydrogens(&mdl, null, null);

    var out: std.Io.Writer.Allocating = .init(testing.allocator);
    defer out.deinit();
    try writeWithDocumentWithPolicy(&out.writer, &mdl, &doc, .{
        .output_isotope = .deuterium,
    });

    var parsed = try cif.readString(testing.allocator, out.writer.buffered());
    defer parsed.deinit();

    const block = &parsed.blocks.items[0];
    const atom_site = block.findLoop("_atom_site.label_atom_id").?;
    const atom_name_col = atom_site.findTag("_atom_site.label_atom_id").?;
    const type_symbol_col = atom_site.findTag("_atom_site.type_symbol").?;

    var found_d: bool = false;
    for (0..atom_site.length()) |row| {
        const atom_name = atom_site.val(row, atom_name_col) orelse continue;
        const type_symbol = atom_site.val(row, type_symbol_col) orelse continue;
        if (std.mem.startsWith(u8, atom_name, "H")) {
            // Added H atoms should have type_symbol "D"
            try testing.expectEqualStrings("D", type_symbol);
            found_d = true;
        }
    }
    try testing.expect(found_d);
}

test "deuterium mode preserves existing H atoms as H" {
    const source = @embedFile("../test_data/ala_with_h.cif");

    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();
    _ = try place.addHydrogens(&mdl, null, null);

    var out: std.Io.Writer.Allocating = .init(testing.allocator);
    defer out.deinit();
    try writeWithPolicy(&out.writer, &mdl, "TEST", .{
        .output_isotope = .deuterium,
    });

    var doc = try cif.readString(testing.allocator, out.writer.buffered());
    defer doc.deinit();

    const block = &doc.blocks.items[0];
    const atom_site = block.findLoop("_atom_site.label_atom_id").?;
    const type_symbol_col = atom_site.findTag("_atom_site.type_symbol").?;
    const atom_name_col = atom_site.findTag("_atom_site.label_atom_id").?;

    var found_original_h = false;
    var found_added_d = false;
    for (0..atom_site.length()) |row| {
        const atom_name = atom_site.val(row, atom_name_col) orelse continue;
        const type_symbol = atom_site.val(row, type_symbol_col) orelse continue;
        if (!std.mem.startsWith(u8, atom_name, "H")) continue;

        // Check if this is an original atom (in the first residue's range) or added
        // ala_with_h.cif has pre-existing H atoms; added H will have type_symbol D
        if (std.mem.eql(u8, type_symbol, "H")) {
            found_original_h = true;
        } else if (std.mem.eql(u8, type_symbol, "D")) {
            found_added_d = true;
        }
    }
    // ala_with_h.cif has pre-existing H, so we should see both H and D
    try testing.expect(found_original_h);
}

test "writeAtomSitePreserving serial numbers are contiguous even with out-of-range rows" {
    // Construct a document with an atom_site loop that has fewer rows than
    // the model expects (simulating orig_row_idx >= orig_loop.length()).
    // The serial counter must only increment for actually-written rows.
    const source =
        \\data_TEST
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
        \\ATOM 1 N N ALA A 1 0.0 0.0 0.0 1.0 10.0 .
        \\ATOM 2 C CA ALA A 1 1.5 0.0 0.0 1.0 10.0 .
        \\#
    ;

    var doc = try cif.readString(testing.allocator, source);
    defer doc.deinit();

    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    // Extend the residue's atom range beyond loop length to simulate gap
    mdl.residues.items[0].atom_end = 4; // only 2 rows exist in loop

    var buf: std.Io.Writer.Allocating = .init(testing.allocator);
    defer buf.deinit();
    try writeWithDocument(&buf.writer, &mdl, &doc);

    // Parse output and verify serial numbers are 1 and 2 (no gap)
    var reparsed = try cif.readString(testing.allocator, buf.writer.buffered());
    defer reparsed.deinit();
    const block = &reparsed.blocks.items[0];
    const atom_site = block.findLoop("_atom_site.id").?;
    const id_col = atom_site.findTag("_atom_site.id").?;
    try testing.expectEqual(@as(usize, 2), atom_site.length());
    try testing.expectEqualStrings("1", atom_site.val(0, id_col).?);
    try testing.expectEqualStrings("2", atom_site.val(1, id_col).?);
}

test "writeWithDocument preserves bare null alt ids on original rows" {
    const source = @embedFile("../test_data/tiny.cif");

    var doc = try cif.readString(testing.allocator, source);
    defer doc.deinit();

    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    var buf: std.Io.Writer.Allocating = .init(testing.allocator);
    defer buf.deinit();
    try writeWithDocument(&buf.writer, &mdl, &doc);

    try testing.expect(std.mem.find(u8, buf.writer.buffered(), "'.'") == null);
}

test "writeZhydLog produces expected output" {
    const allocator = testing.allocator;

    var chains = [_]Chain{blk: {
        var c = Chain{};
        c.setLabelAsymId("A");
        break :blk c;
    }};
    var residues = [_]Residue{blk: {
        var r = Residue{};
        r.setCompId("ASN");
        r.seq_id = 5;
        r.chain_idx = 0;
        break :blk r;
    }};

    const positions = try allocator.alloc(@import("../math.zig").Vec3(f32), 1);
    defer allocator.free(positions);
    positions[0] = .{ .x = 0, .y = 0, .z = 0 };

    const orientations = try allocator.alloc(mover_mod.Orientation, 1);
    defer allocator.free(orientations);
    orientations[0] = .{ .positions = positions };

    const atom_indices = try allocator.alloc(u32, 1);
    defer allocator.free(atom_indices);
    atom_indices[0] = 0;

    const movers = [_]mover_mod.Mover{.{
        .kind = .amide_flip,
        .residue_idx = 0,
        .atom_indices = atom_indices,
        .orientations = orientations,
        .best_orientation = 1,
        .current_orientation = 0,
        .allocator = allocator,
    }};

    var buf: std.Io.Writer.Allocating = .init(allocator);
    defer buf.deinit();

    try writeZhydLog(&buf.writer, &movers, &residues, &chains);

    const output = buf.writer.buffered();
    try testing.expect(std.mem.find(u8, output, "_zhyd_log.residue_id") != null);
    try testing.expect(std.mem.find(u8, output, "flip_amide") != null);
    try testing.expect(std.mem.find(u8, output, "A.ASN.5") != null);
}

test "writeMultiModelWithDocumentWithPolicy round-trip preserves model numbers" {
    const allocator = std.testing.allocator;
    const source = @embedFile("../test_data/multi_model.cif");

    // Parse document and models
    var doc = try cif.readString(allocator, source);
    defer doc.deinit();

    var entries = try mmcif.parseModelsFromBlock(allocator, &doc.blocks.items[0], .all);
    defer {
        for (entries.items) |*e| e.model.deinit();
        entries.deinit(allocator);
    }

    // Write multi-model output
    var buf: std.Io.Writer.Allocating = .init(allocator);
    defer buf.deinit();

    try writeMultiModelWithDocumentWithPolicy(&buf.writer, entries.items, &doc, .{});

    const output = buf.writer.buffered();

    // Verify both model numbers appear in output
    const model1_count = blk: {
        var count: usize = 0;
        var pos: usize = 0;
        while (std.mem.findPos(u8, output, pos, " 1\n")) |idx| {
            count += 1;
            pos = idx + 1;
        }
        break :blk count;
    };
    // Model 1 should have at least 4 atom rows ending with " 1\n"
    try testing.expect(model1_count >= 4);

    // Re-parse and verify atom count
    var doc2 = try cif.readString(allocator, output);
    defer doc2.deinit();

    var entries2 = try mmcif.parseModelsFromBlock(allocator, &doc2.blocks.items[0], .all);
    defer {
        for (entries2.items) |*e| e.model.deinit();
        entries2.deinit(allocator);
    }

    // Should have 2 models with 4 atoms each
    try testing.expectEqual(@as(usize, 2), entries2.items.len);
    try testing.expectEqual(@as(usize, 4), entries2.items[0].model.atoms.items.len);
    try testing.expectEqual(@as(usize, 4), entries2.items[1].model.atoms.items.len);

    // Verify model numbers
    try testing.expectEqual(@as(u32, 1), entries2.items[0].model_num);
    try testing.expectEqual(@as(u32, 2), entries2.items[1].model_num);

    // Verify coordinates preserved
    try testing.expectApproxEqAbs(@as(f32, 1.0), entries2.items[0].model.atoms.items[0].pos.x, 0.01);
    try testing.expectApproxEqAbs(@as(f32, 11.0), entries2.items[1].model.atoms.items[0].pos.x, 0.01);
}
