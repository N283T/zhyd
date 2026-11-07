//! mmCIF extraction: parses _atom_site loop data into a Model.

const std = @import("std");
const Allocator = std.mem.Allocator;

const cif = @import("cif.zig");
const model_mod = @import("model.zig");
const element = @import("element.zig");
const bond_mod = @import("model/bond.zig");
const ccd_mod = @import("ccd.zig");

const Model = model_mod.Model;
const Atom = model_mod.Atom;
const Residue = model_mod.Residue;
const Chain = model_mod.Chain;

pub const ModelFilter = union(enum) {
    all,
    first,
    specific: u32,
};

pub const ModelEntry = struct {
    model_num: u32,
    model: Model,
    /// Row range in the _atom_site loop for this model's original atoms.
    /// Used by the preserving writer and buildAtomLookupForRange.
    cif_row_start: u32,
    cif_row_end: u32,

    pub fn deinit(self: *ModelEntry) void {
        self.model.deinit();
    }
};

pub const MmcifError = error{
    NoAtomSiteLoop,
    MissingCoordinateField,
    MissingRequiredField,
    InvalidCoordinateValue,
    CifParseError,
    OutOfMemory,
};

/// Column indices for _atom_site loop fields.
const AtomSiteColumns = struct {
    cartn_x: ?usize = null,
    cartn_y: ?usize = null,
    cartn_z: ?usize = null,
    type_symbol: ?usize = null,
    label_atom_id: ?usize = null,
    label_comp_id: ?usize = null,
    label_asym_id: ?usize = null,
    label_seq_id: ?usize = null,
    auth_asym_id: ?usize = null,
    label_alt_id: ?usize = null,
    occupancy: ?usize = null,
    b_iso_or_equiv: ?usize = null,
    group_pdb: ?usize = null,
    id: ?usize = null,
    pdbx_PDB_ins_code: ?usize = null,
    label_entity_id: ?usize = null,
    auth_seq_id: ?usize = null,
    pdbx_PDB_model_num: ?usize = null,
};

/// Compare a Chain's label_asym_id with a string from the CIF data.
fn chainAsymMatches(chain: Chain, asym: []const u8) bool {
    return std.mem.eql(u8, chain.labelSlice(), asym);
}

/// Scan _pdbx_poly_seq_scheme to detect sequence gaps and mark residues
/// that follow an unobserved residue (auth_seq_num = '?' or '.').
fn detectChainBreaks(
    mdl: *Model,
    pss: *const cif.Loop,
    col_asym: usize,
    col_seq: usize,
    col_auth_seq: usize,
) void {
    const nrows = pss.length();
    var prev_asym: []const u8 = "";
    var gap_pending = false;
    // Cached chain index for the current asym_id; null = not yet found.
    var cached_chain_idx: ?usize = null;

    for (0..nrows) |row| {
        const asym = cif.asString(pss.val(row, col_asym) orelse continue);
        const seq_str = pss.val(row, col_seq) orelse continue;
        const auth_seq = cif.asString(pss.val(row, col_auth_seq) orelse "?");

        // Reset on chain change
        if (!std.mem.eql(u8, asym, prev_asym)) {
            gap_pending = false;
            prev_asym = asym;
            // Find and cache the chain index for this asym_id.
            cached_chain_idx = null;
            for (mdl.chains.items, 0..) |chain, ci| {
                if (chainAsymMatches(chain, asym)) {
                    cached_chain_idx = ci;
                    break;
                }
            }
        }

        // Unobserved residue (asString converts '?' and '.' to empty)
        if (auth_seq.len == 0) {
            gap_pending = true;
            continue;
        }

        // Observed residue after gap -> find and mark within the matching chain only.
        if (gap_pending) {
            const seq_id = cif.value.asIntOr(i32, seq_str, 0);
            if (cached_chain_idx) |ci| {
                const chain = mdl.chains.items[ci];
                for (mdl.residues.items[chain.residue_start..chain.residue_end]) |*res| {
                    if (res.seq_id == seq_id) {
                        res.is_chain_break_before = true;
                        break;
                    }
                }
            }
            gap_pending = false;
        }
    }
}

const EntityType = model_mod.residue.EntityType;

fn entityTypeFromString(s: []const u8) EntityType {
    if (std.ascii.eqlIgnoreCase(s, "polymer")) return .polymer;
    if (std.ascii.eqlIgnoreCase(s, "non-polymer")) return .non_polymer;
    if (std.ascii.eqlIgnoreCase(s, "branched")) return .branched;
    if (std.ascii.eqlIgnoreCase(s, "water")) return .water;
    return .unknown;
}

fn applyEntityTypes(mdl: *Model, block: *const cif.Block) void {
    const ent = block.findLoop("_entity.id") orelse {
        // Single-entity files usually write _entity as key-value pairs.
        const eid = cif.asString(block.findValue("_entity.id") orelse return);
        const etype = entityTypeFromString(cif.asString(block.findValue("_entity.type") orelse return));
        if (etype == .unknown) return;
        for (mdl.chains.items) |chain| {
            if (!std.mem.eql(u8, chain.entityIdSlice(), eid)) continue;
            for (mdl.residues.items[chain.residue_start..chain.residue_end]) |*res| res.entity_type = etype;
        }
        return;
    };
    const col_id = ent.findTag("_entity.id") orelse return;
    const col_type = ent.findTag("_entity.type") orelse return;

    // Iterate chains (C << R) — all residues in a chain share the same entity_id.
    for (mdl.chains.items) |chain| {
        const entity_id = chain.entityIdSlice();
        if (entity_id.len == 0) continue;

        // Find entity type for this chain's entity_id (E is typically 3-5 rows).
        var etype: EntityType = .unknown;
        for (0..ent.length()) |row| {
            const eid = cif.asString(ent.val(row, col_id) orelse continue);
            if (std.mem.eql(u8, eid, entity_id)) {
                const etype_str = cif.asString(ent.val(row, col_type) orelse continue);
                etype = entityTypeFromString(etype_str);
                break;
            }
        }
        if (etype == .unknown) continue;

        // Apply to all residues in this chain.
        const res_slice = mdl.residues.items[chain.residue_start..chain.residue_end];
        for (res_slice) |*res| {
            res.entity_type = etype;
        }
    }
}

/// Parse an mmCIF source string and extract the first model's _atom_site records.
/// Convenience wrapper around parseModels with filter=.first.
pub fn parseModel(allocator: Allocator, source: []const u8) MmcifError!Model {
    var entries = try parseModels(allocator, source, .first);
    defer {
        // Free all entries except index 0 (which we return)
        if (entries.items.len > 1) {
            for (entries.items[1..]) |*e| e.model.deinit();
        }
        entries.deinit(allocator);
    }
    if (entries.items.len == 0) return MmcifError.NoAtomSiteLoop;
    return entries.items[0].model;
}

/// Parse an mmCIF source string and extract _atom_site records into Models.
/// Convenience wrapper that parses the CIF document internally.
pub fn parseModels(allocator: Allocator, source: []const u8, filter: ModelFilter) MmcifError!std.ArrayListUnmanaged(ModelEntry) {
    var doc = cif.readString(allocator, source) catch |err| switch (err) {
        error.OutOfMemory => return MmcifError.OutOfMemory,
        else => return MmcifError.CifParseError,
    };
    defer doc.deinit();

    if (doc.blocks.items.len == 0) return MmcifError.NoAtomSiteLoop;
    return parseModelsFromBlock(allocator, &doc.blocks.items[0], filter);
}

/// Extract _atom_site records from a pre-parsed CIF block into Models.
/// Returns one ModelEntry per model found (filtered by ModelFilter).
/// Row indices in returned entries correspond to the block's _atom_site loop.
pub fn parseModelsFromBlock(allocator: Allocator, block: *const cif.Block, filter: ModelFilter) MmcifError!std.ArrayListUnmanaged(ModelEntry) {
    const loop = block.findLoop("_atom_site.Cartn_x") orelse return MmcifError.NoAtomSiteLoop;

    // Map column indices
    var cols = AtomSiteColumns{};
    cols.cartn_x = loop.findTag("_atom_site.Cartn_x");
    cols.cartn_y = loop.findTag("_atom_site.Cartn_y");
    cols.cartn_z = loop.findTag("_atom_site.Cartn_z");
    cols.type_symbol = loop.findTag("_atom_site.type_symbol");
    cols.label_atom_id = loop.findTag("_atom_site.label_atom_id");
    cols.label_comp_id = loop.findTag("_atom_site.label_comp_id");
    cols.label_asym_id = loop.findTag("_atom_site.label_asym_id");
    cols.label_seq_id = loop.findTag("_atom_site.label_seq_id");
    cols.auth_asym_id = loop.findTag("_atom_site.auth_asym_id");
    cols.label_alt_id = loop.findTag("_atom_site.label_alt_id");
    cols.occupancy = loop.findTag("_atom_site.occupancy");
    cols.b_iso_or_equiv = loop.findTag("_atom_site.B_iso_or_equiv");
    cols.group_pdb = loop.findTag("_atom_site.group_PDB");
    cols.id = loop.findTag("_atom_site.id");
    cols.pdbx_PDB_ins_code = loop.findTag("_atom_site.pdbx_PDB_ins_code");
    cols.label_entity_id = loop.findTag("_atom_site.label_entity_id");
    cols.auth_seq_id = loop.findTag("_atom_site.auth_seq_id");
    cols.pdbx_PDB_model_num = loop.findTag("_atom_site.pdbx_PDB_model_num");

    // Require x, y, z
    if (cols.cartn_x == null or cols.cartn_y == null or cols.cartn_z == null) {
        return MmcifError.MissingCoordinateField;
    }

    var entries = std.ArrayListUnmanaged(ModelEntry).empty;
    errdefer {
        for (entries.items) |*e| e.model.deinit();
        entries.deinit(allocator);
    }

    var mdl = Model.init(allocator);
    errdefer mdl.deinit();

    const nrows = loop.length();

    // State tracking for chain/residue boundaries (reset per model)
    var cur_label_asym_id: []const u8 = "";
    var cur_seq_id: []const u8 = "";
    var cur_comp_id: []const u8 = "";
    var cur_ins_code: []const u8 = "";
    var cur_auth_seq: []const u8 = "";
    var res_atom_start: u32 = 0;
    var in_chain = false;
    var in_residue = false;

    // Multi-model tracking
    var cur_model_str: []const u8 = "";
    var cur_model_row_start: u32 = 0;
    var model_started = false;

    for (0..nrows) |row| {
        // Determine model number for this row
        var row_model_str: []const u8 = "";
        if (cols.pdbx_PDB_model_num) |model_col| {
            const raw = loop.val(row, model_col);
            row_model_str = if (raw) |r| cif.asString(r) else "";
            if (row_model_str.len == 0) continue; // skip CIF-null model numbers
        }

        // Detect model boundary
        const model_changed = model_started and cols.pdbx_PDB_model_num != null and
            !std.mem.eql(u8, row_model_str, cur_model_str);

        if (model_changed) {
            // Finalize and push current model only if it has atoms
            if (mdl.atoms.items.len > 0) {
                finalizeModel(&mdl, in_residue, in_chain);

                const model_num = std.fmt.parseInt(u32, cur_model_str, 10) catch 1;
                mdl.model_num = model_num;

                try entries.append(allocator, .{
                    .model_num = model_num,
                    .model = mdl,
                    .cif_row_start = cur_model_row_start,
                    .cif_row_end = @intCast(row),
                });

                mdl = Model.init(allocator);
            } else {
                // Discard empty model (skipped by filter)
                mdl.deinit();
                mdl = Model.init(allocator);
            }
            in_residue = false;
            in_chain = false;

            // For .first filter, stop after the first model
            if (filter == .first) break;
            // For .specific filter, stop if we already found the target
            switch (filter) {
                .specific => {
                    if (entries.items.len > 0) break;
                },
                else => {},
            }

            // Reset state for new model
            cur_label_asym_id = "";
            cur_seq_id = "";
            cur_comp_id = "";
            cur_ins_code = "";
            cur_auth_seq = "";
            cur_model_row_start = @intCast(row);
        }

        if (!model_started) {
            model_started = true;
            cur_model_row_start = @intCast(row);
        }
        cur_model_str = row_model_str;

        // Apply filter: skip rows for non-matching specific model
        switch (filter) {
            .specific => |target| {
                const row_num = std.fmt.parseInt(u32, row_model_str, 10) catch 0;
                if (cols.pdbx_PDB_model_num != null and row_num != target) {
                    // If we already collected atoms for the target, we can stop
                    // (models are contiguous per mmCIF convention)
                    if (mdl.atoms.items.len > 0) break;
                    continue;
                }
            },
            else => {},
        }

        const x_str = loop.val(row, cols.cartn_x.?) orelse return MmcifError.InvalidCoordinateValue;
        const y_str = loop.val(row, cols.cartn_y.?) orelse return MmcifError.InvalidCoordinateValue;
        const z_str = loop.val(row, cols.cartn_z.?) orelse return MmcifError.InvalidCoordinateValue;

        const x = cif.asFloat(x_str) orelse return MmcifError.InvalidCoordinateValue;
        const y = cif.asFloat(y_str) orelse return MmcifError.InvalidCoordinateValue;
        const z = cif.asFloat(z_str) orelse return MmcifError.InvalidCoordinateValue;

        // Get string fields
        const type_sym = if (cols.type_symbol) |c| cif.asString(loop.val(row, c) orelse ".") else "";
        const atom_name = if (cols.label_atom_id) |c| cif.asString(loop.val(row, c) orelse ".") else "";
        const comp_id = if (cols.label_comp_id) |c| cif.asString(loop.val(row, c) orelse ".") else "";
        const label_asym = if (cols.label_asym_id) |c| cif.asString(loop.val(row, c) orelse ".") else "";
        const seq_id = if (cols.label_seq_id) |c| cif.asString(loop.val(row, c) orelse ".") else "";
        const auth_seq = if (cols.auth_seq_id) |c| cif.asString(loop.val(row, c) orelse ".") else "";
        const auth_asym = if (cols.auth_asym_id) |c| cif.asString(loop.val(row, c) orelse ".") else "";
        const alt_loc_str = if (cols.label_alt_id) |c| loop.val(row, c) orelse "." else ".";
        const occ_str = if (cols.occupancy) |c| loop.val(row, c) orelse "1.0" else "1.0";
        const bfac_str = if (cols.b_iso_or_equiv) |c| loop.val(row, c) orelse "0.0" else "0.0";
        const serial_str = if (cols.id) |c| loop.val(row, c) orelse "0" else "0";
        const ins_code_str = if (cols.pdbx_PDB_ins_code) |c| cif.asString(loop.val(row, c) orelse ".") else "";

        // Chain boundary detection
        if (!in_chain or !std.mem.eql(u8, label_asym, cur_label_asym_id)) {
            // Close previous residue
            if (in_residue) {
                const atom_end: u32 = @intCast(mdl.atoms.items.len);
                const res_idx = mdl.residues.items.len - 1;
                mdl.residues.items[res_idx].atom_end = atom_end;
                in_residue = false;
            }
            // Close previous chain
            if (in_chain) {
                const res_end: u32 = @intCast(mdl.residues.items.len);
                const chain_idx = mdl.chains.items.len - 1;
                mdl.chains.items[chain_idx].residue_end = res_end;
            }
            // Start new chain
            var new_chain = Chain{};
            new_chain.setLabelAsymId(label_asym);
            new_chain.setAuthAsymId(auth_asym);
            const entity_id_str = if (cols.label_entity_id) |c| cif.asString(loop.val(row, c) orelse ".") else "";
            new_chain.setEntityId(entity_id_str);
            new_chain.residue_start = @intCast(mdl.residues.items.len);
            try mdl.chains.append(mdl.allocator, new_chain);
            cur_label_asym_id = label_asym;
            in_chain = true;
            // Force new residue
            cur_seq_id = "";
            cur_comp_id = "";
            cur_ins_code = "";
            cur_auth_seq = "";
        }

        // Residue boundary detection (includes insertion code and auth_seq for non-polymer)
        const new_residue = !in_residue or
            !std.mem.eql(u8, seq_id, cur_seq_id) or
            !std.mem.eql(u8, comp_id, cur_comp_id) or
            !std.mem.eql(u8, ins_code_str, cur_ins_code) or
            !std.mem.eql(u8, auth_seq, cur_auth_seq);

        if (new_residue) {
            // Close previous residue
            if (in_residue) {
                const atom_end: u32 = @intCast(mdl.atoms.items.len);
                const res_idx = mdl.residues.items.len - 1;
                mdl.residues.items[res_idx].atom_end = atom_end;
            }
            // Start new residue
            res_atom_start = @intCast(mdl.atoms.items.len);
            var new_res = Residue{};
            new_res.setCompId(comp_id);
            const label_seq_int = cif.value.asIntOr(i32, seq_id, 0);
            const auth_seq_int = cif.value.asIntOr(i32, auth_seq, 0);
            new_res.seq_id = if (seq_id.len == 0) auth_seq_int else label_seq_int;
            new_res.auth_seq_id = auth_seq_int;
            new_res.ins_code = if (ins_code_str.len > 0) ins_code_str[0] else ' ';
            new_res.chain_idx = @intCast(mdl.chains.items.len - 1);
            new_res.atom_start = res_atom_start;
            new_res.atom_end = res_atom_start; // will be updated
            try mdl.residues.append(mdl.allocator, new_res);
            cur_seq_id = seq_id;
            cur_comp_id = comp_id;
            cur_ins_code = ins_code_str;
            cur_auth_seq = auth_seq;
            in_residue = true;
        }

        // Build atom
        var atom = Atom{
            .pos = .{ .x = x, .y = y, .z = z },
        };
        atom.setName(atom_name);
        atom.element_type = element.elementFromSymbol(type_sym);
        atom.is_hydrogen = switch (atom.element_type) {
            .H, .Har, .Hpol, .Ha_p, .HOd => true,
            else => false,
        };
        atom.occupancy = cif.asFloatOr(occ_str, 1.0);
        atom.b_factor = cif.asFloatOr(bfac_str, 0.0);
        atom.serial = @intCast(cif.value.asIntOr(u32, serial_str, 0));
        atom.residue_idx = @intCast(mdl.residues.items.len - 1);

        // Alt location: first char if not '.'
        const alt_raw = cif.asString(alt_loc_str);
        atom.altloc = if (alt_raw.len > 0) alt_raw[0] else ' ';

        try mdl.atoms.append(mdl.allocator, atom);
    }

    // Finalize the last model (if it has atoms)
    if (mdl.atoms.items.len > 0) {
        finalizeModel(&mdl, in_residue, in_chain);

        const model_num = if (cur_model_str.len > 0)
            std.fmt.parseInt(u32, cur_model_str, 10) catch 1
        else
            1;
        mdl.model_num = model_num;

        try entries.append(allocator, .{
            .model_num = model_num,
            .model = mdl,
            .cif_row_start = cur_model_row_start,
            .cif_row_end = @intCast(nrows),
        });
        // mdl ownership transferred to entries; prevent errdefer double-free
        mdl = Model.init(allocator);
    }
    mdl.deinit(); // free the empty spare model

    // Apply per-model post-processing
    for (entries.items) |*entry| {
        // Parse _pdbx_poly_seq_scheme for chain-break detection (optional)
        if (block.findLoop("_pdbx_poly_seq_scheme.seq_id")) |pss| {
            const col_asym = pss.findTag("_pdbx_poly_seq_scheme.asym_id");
            const col_seq = pss.findTag("_pdbx_poly_seq_scheme.seq_id");
            const col_auth_seq = pss.findTag("_pdbx_poly_seq_scheme.auth_seq_num");

            if (col_asym != null and col_seq != null and col_auth_seq != null) {
                detectChainBreaks(&entry.model, pss, col_asym.?, col_seq.?, col_auth_seq.?);
            }
        }

        // Parse _pdbx_unobs_or_zero_occ_atoms count (optional, diagnostic)
        if (block.findLoop("_pdbx_unobs_or_zero_occ_atoms.label_atom_id")) |unobs| {
            entry.model.n_unobs_atoms = @intCast(unobs.length());
        }

        // Parse _entity loop for entity types (optional)
        applyEntityTypes(&entry.model, block);
    }

    return entries;
}

/// Finalize a model: close the last residue and chain, set original_atom_count.
fn finalizeModel(mdl: *Model, in_residue: bool, in_chain: bool) void {
    if (in_residue) {
        const atom_end: u32 = @intCast(mdl.atoms.items.len);
        const res_idx = mdl.residues.items.len - 1;
        mdl.residues.items[res_idx].atom_end = atom_end;
    }
    mdl.original_atom_count = @intCast(mdl.atoms.items.len);
    if (in_chain) {
        const res_end: u32 = @intCast(mdl.residues.items.len);
        const chain_idx = mdl.chains.items.len - 1;
        mdl.chains.items[chain_idx].residue_end = res_end;
    }
}

// ── Sub-module re-exports ──────────────────────────────────────────────────────

const conn_mod = @import("mmcif/conn.zig");
const inline_comp_mod = @import("mmcif/inline_comp.zig");

pub const AtomLookupKey = conn_mod.AtomLookupKey;
pub const AtomLookupContext = conn_mod.AtomLookupContext;
pub const AtomLookup = conn_mod.AtomLookup;
pub const buildAtomLookup = conn_mod.buildAtomLookup;
pub const buildAtomLookupForRange = conn_mod.buildAtomLookupForRange;
pub const isCovalentConnType = conn_mod.isCovalentConnType;
pub const parseStructConn = conn_mod.parseStructConn;
pub const parseBranchLinks = conn_mod.parseBranchLinks;

pub const parseInlineComponents = inline_comp_mod.parseInlineComponents;
pub const flagLeavingAtoms = inline_comp_mod.flagLeavingAtoms;

// ── Tests ─────────────────────────────────────────────────────────────────────

const testing = std.testing;

test "parse tiny mmCIF" {
    const source = @embedFile("test_data/tiny.cif");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    try testing.expectEqual(@as(usize, 5), mdl.atoms.items.len);
    try testing.expectEqual(@as(usize, 1), mdl.residues.items.len);
    try testing.expectEqual(@as(usize, 1), mdl.chains.items.len);
    try testing.expectApproxEqAbs(mdl.atoms.items[0].pos.x, 1.0, 1e-3);
    try testing.expectApproxEqAbs(mdl.atoms.items[0].pos.y, 2.0, 1e-3);
    try testing.expectEqualStrings("ALA", mdl.residues.items[0].compIdSlice());
}

test "parse multi-chain multi-residue mmCIF" {
    const source = @embedFile("test_data/multi_chain.cif");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    // 11 atoms total
    try testing.expectEqual(@as(usize, 11), mdl.atoms.items.len);
    // 3 residues: ALA(A/1), GLY(A/2), VAL(B/1)
    try testing.expectEqual(@as(usize, 3), mdl.residues.items.len);
    // 2 chains: A and B
    try testing.expectEqual(@as(usize, 2), mdl.chains.items.len);

    // Chain A: residues 0-1, atoms 0-7
    try testing.expectEqualStrings("A", mdl.chains.items[0].labelSlice());
    try testing.expectEqual(@as(u32, 0), mdl.chains.items[0].residue_start);
    try testing.expectEqual(@as(u32, 2), mdl.chains.items[0].residue_end);

    // Chain B: residue 2, atoms 8-10
    try testing.expectEqualStrings("B", mdl.chains.items[1].labelSlice());
    try testing.expectEqual(@as(u32, 2), mdl.chains.items[1].residue_start);
    try testing.expectEqual(@as(u32, 3), mdl.chains.items[1].residue_end);

    // Residue 0: ALA, atoms 0-3
    try testing.expectEqualStrings("ALA", mdl.residues.items[0].compIdSlice());
    try testing.expectEqual(@as(u32, 0), mdl.residues.items[0].atom_start);
    try testing.expectEqual(@as(u32, 4), mdl.residues.items[0].atom_end);

    // Residue 1: GLY, atoms 4-7
    try testing.expectEqualStrings("GLY", mdl.residues.items[1].compIdSlice());
    try testing.expectEqual(@as(u32, 4), mdl.residues.items[1].atom_start);
    try testing.expectEqual(@as(u32, 8), mdl.residues.items[1].atom_end);

    // Residue 2: VAL, atoms 8-10
    try testing.expectEqualStrings("VAL", mdl.residues.items[2].compIdSlice());
    try testing.expectEqual(@as(u32, 8), mdl.residues.items[2].atom_start);
    try testing.expectEqual(@as(u32, 11), mdl.residues.items[2].atom_end);

    // Last atom coordinate check
    try testing.expectApproxEqAbs(mdl.atoms.items[10].pos.x, 13.0, 1e-3);
}

test "error: no atom_site loop" {
    const source = "data_EMPTY\n_entry.id EMPTY\n";
    const result = parseModel(testing.allocator, source);
    try testing.expectError(MmcifError.NoAtomSiteLoop, result);
}

test "error: missing coordinate field" {
    const source =
        \\data_BAD
        \\loop_
        \\_atom_site.label_atom_id
        \\CA
    ;
    const result = parseModel(testing.allocator, source);
    // No Cartn_x -> NoAtomSiteLoop (findLoop won't find it)
    try testing.expectError(MmcifError.NoAtomSiteLoop, result);
}

test "error: invalid coordinate value" {
    const source =
        \\data_BAD
        \\loop_
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\1.0 2.0 abc
    ;
    const result = parseModel(testing.allocator, source);
    try testing.expectError(MmcifError.InvalidCoordinateValue, result);
}

test "parse chain break from pdbx_poly_seq_scheme" {
    const source = @embedFile("test_data/gap_chain.cif");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    try testing.expectEqual(@as(usize, 2), mdl.residues.items.len);
    try testing.expectEqual(@as(usize, 1), mdl.chains.items.len);

    // First residue: no chain break before
    try testing.expect(!mdl.residues.items[0].is_chain_break_before);
    // Second residue: chain break before (seq_id 2 is unobserved)
    try testing.expect(mdl.residues.items[1].is_chain_break_before);
}

test "parse without pdbx_poly_seq_scheme is backward compatible" {
    const source = @embedFile("test_data/tiny.cif");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    try testing.expect(!mdl.residues.items[0].is_chain_break_before);
}

test "model reports zero unobs atoms for tiny.cif" {
    const source = @embedFile("test_data/tiny.cif");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    try testing.expectEqual(@as(u32, 0), mdl.n_unobs_atoms);
}

test "insertion code splits residues with same seq_id" {
    const source = @embedFile("test_data/ins_code.cif");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    // Same seq_id=1 but different ins_codes (blank vs 'A') -> 2 residues
    try testing.expectEqual(@as(usize, 2), mdl.residues.items.len);
    try testing.expectEqual(@as(usize, 1), mdl.chains.items.len);

    // First residue: ALA with ins_code=' '
    try testing.expectEqualStrings("ALA", mdl.residues.items[0].compIdSlice());
    try testing.expectEqual(@as(u8, ' '), mdl.residues.items[0].ins_code);

    // Second residue: GLY with ins_code='A'
    try testing.expectEqualStrings("GLY", mdl.residues.items[1].compIdSlice());
    try testing.expectEqual(@as(u8, 'A'), mdl.residues.items[1].ins_code);
}

test "parseModel stores auth_seq_id and falls back" {
    const source = @embedFile("test_data/entity_type.cif");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    // Polymer ALA: label_seq_id=1, auth_seq_id=1
    try testing.expectEqual(@as(i32, 1), mdl.residues.items[0].seq_id);
    try testing.expectEqual(@as(i32, 1), mdl.residues.items[0].auth_seq_id);

    // Non-polymer EDO: label_seq_id="." -> fallback to auth_seq_id=1
    try testing.expectEqual(@as(i32, 1), mdl.residues.items[1].seq_id);
    try testing.expectEqual(@as(i32, 1), mdl.residues.items[1].auth_seq_id);

    // Water HOH: label_seq_id="." -> fallback to auth_seq_id=1
    try testing.expectEqual(@as(i32, 1), mdl.residues.items[2].seq_id);

    // Branched NAG residue 1: label_seq_id="." -> fallback to auth_seq_id=1
    try testing.expectEqual(@as(i32, 1), mdl.residues.items[3].seq_id);
    try testing.expectEqual(@as(i32, 1), mdl.residues.items[3].auth_seq_id);

    // Branched NAG residue 2: label_seq_id="." -> fallback to auth_seq_id=2
    try testing.expectEqual(@as(i32, 2), mdl.residues.items[4].seq_id);
    try testing.expectEqual(@as(i32, 2), mdl.residues.items[4].auth_seq_id);
}

test "parseModel sets entity_type from _entity key-value pairs" {
    const source =
        \\data_PAIR
        \\_entity.id   1
        \\_entity.type polymer
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
        \\_atom_site.auth_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\_atom_site.pdbx_PDB_ins_code
        \\ATOM 1 N N  ALA A 1 1 1 1.0 0.0 0.0 1.0 0.0 . .
        \\ATOM 2 C CA ALA A 1 1 1 2.0 0.0 0.0 1.0 0.0 . .
        \\#
        \\
    ;
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();
    try testing.expectEqual(@as(usize, 1), mdl.residues.items.len);
    try testing.expectEqual(EntityType.polymer, mdl.residues.items[0].entity_type);
}

test "parseModel sets entity_type from _entity loop" {
    const source = @embedFile("test_data/entity_type.cif");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    // 5 residues: ALA, EDO, HOH, NAG(1), NAG(2)
    try testing.expectEqual(@as(usize, 5), mdl.residues.items.len);

    // ALA (entity 1) -> polymer
    try testing.expectEqual(EntityType.polymer, mdl.residues.items[0].entity_type);

    // EDO (entity 2) -> non_polymer
    try testing.expectEqual(EntityType.non_polymer, mdl.residues.items[1].entity_type);

    // HOH (entity 3) -> water
    try testing.expectEqual(EntityType.water, mdl.residues.items[2].entity_type);

    // NAG (entity 4) -> branched
    try testing.expectEqual(EntityType.branched, mdl.residues.items[3].entity_type);
    try testing.expectEqual(EntityType.branched, mdl.residues.items[4].entity_type);
}

test "parseModel without _entity loop keeps entity_type unknown" {
    const source = @embedFile("test_data/tiny.cif");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    // tiny.cif has no _entity loop -- all residues should stay .unknown
    for (mdl.residues.items) |res| {
        try testing.expectEqual(EntityType.unknown, res.entity_type);
    }
}

test "parseModel extracts only the first model from multi-model mmCIF" {
    const source = @embedFile("test_data/multi_model.cif");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    // Only 4 atoms from model 1, not 8
    try testing.expectEqual(@as(usize, 4), mdl.atoms.items.len);

    // Verify coordinates are from model 1
    try testing.expectApproxEqAbs(@as(f32, 1.0), mdl.atoms.items[0].pos.x, 0.01);
    try testing.expectApproxEqAbs(@as(f32, 2.0), mdl.atoms.items[1].pos.x, 0.01);

    // model_num should be parsed from the column
    try testing.expectEqual(@as(u32, 1), mdl.model_num);
}

test "parseModel sets model_num to 1 when column is absent" {
    const source = @embedFile("test_data/tiny.cif");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();
    try testing.expectEqual(@as(u32, 1), mdl.model_num);
}

test "parseModel skips CIF-null model numbers and filters by first real model" {
    const source = @embedFile("test_data/multi_model_null.cif");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    // Row 1: model '?' -> skipped (CIF null)
    // Row 2: model '1' -> first_model_num = "1", included
    // Row 3: model '1' -> matches, included
    // Row 4: model '2' -> break
    try testing.expectEqual(@as(usize, 2), mdl.atoms.items.len);
    try testing.expectApproxEqAbs(@as(f32, 2.0), mdl.atoms.items[0].pos.x, 0.01);
    try testing.expectApproxEqAbs(@as(f32, 3.0), mdl.atoms.items[1].pos.x, 0.01);
}

test "parseModels with filter .all returns both models" {
    const source = @embedFile("test_data/multi_model.cif");
    var entries = try parseModels(testing.allocator, source, .all);
    defer {
        for (entries.items) |*e| e.model.deinit();
        entries.deinit(testing.allocator);
    }

    try testing.expectEqual(@as(usize, 2), entries.items.len);

    // Model 1
    try testing.expectEqual(@as(u32, 1), entries.items[0].model_num);
    try testing.expectEqual(@as(usize, 4), entries.items[0].model.atoms.items.len);
    try testing.expectApproxEqAbs(@as(f32, 1.0), entries.items[0].model.atoms.items[0].pos.x, 0.01);
    try testing.expectEqual(@as(u32, 0), entries.items[0].cif_row_start);
    try testing.expectEqual(@as(u32, 4), entries.items[0].cif_row_end);

    // Model 2
    try testing.expectEqual(@as(u32, 2), entries.items[1].model_num);
    try testing.expectEqual(@as(usize, 4), entries.items[1].model.atoms.items.len);
    try testing.expectApproxEqAbs(@as(f32, 11.0), entries.items[1].model.atoms.items[0].pos.x, 0.01);
    try testing.expectEqual(@as(u32, 4), entries.items[1].cif_row_start);
    try testing.expectEqual(@as(u32, 8), entries.items[1].cif_row_end);
}

test "parseModels with filter .specific selects one model" {
    const source = @embedFile("test_data/multi_model.cif");
    var entries = try parseModels(testing.allocator, source, .{ .specific = 2 });
    defer {
        for (entries.items) |*e| e.model.deinit();
        entries.deinit(testing.allocator);
    }

    try testing.expectEqual(@as(usize, 1), entries.items.len);
    try testing.expectEqual(@as(u32, 2), entries.items[0].model_num);
    try testing.expectApproxEqAbs(@as(f32, 11.0), entries.items[0].model.atoms.items[0].pos.x, 0.01);
}

test "parseModels with non-existent model returns empty" {
    const source = @embedFile("test_data/multi_model.cif");
    var entries = try parseModels(testing.allocator, source, .{ .specific = 99 });
    defer {
        for (entries.items) |*e| e.model.deinit();
        entries.deinit(testing.allocator);
    }

    try testing.expectEqual(@as(usize, 0), entries.items.len);
}

test "parseModels without model column returns single entry" {
    const source = @embedFile("test_data/tiny.cif");
    var entries = try parseModels(testing.allocator, source, .all);
    defer {
        for (entries.items) |*e| e.model.deinit();
        entries.deinit(testing.allocator);
    }

    try testing.expectEqual(@as(usize, 1), entries.items.len);
    try testing.expectEqual(@as(u32, 1), entries.items[0].model_num);
}

test "parseModelsFromBlock works with pre-parsed document" {
    const source = @embedFile("test_data/multi_model.cif");
    var doc = try cif.readString(testing.allocator, source);
    defer doc.deinit();

    var entries = try parseModelsFromBlock(testing.allocator, &doc.blocks.items[0], .all);
    defer {
        for (entries.items) |*e| e.model.deinit();
        entries.deinit(testing.allocator);
    }

    try testing.expectEqual(@as(usize, 2), entries.items.len);
    try testing.expectEqual(@as(u32, 1), entries.items[0].model_num);
    try testing.expectEqual(@as(u32, 2), entries.items[1].model_num);
}

test "parseModelsFromBlock with strip_h produces consistent row indices" {
    const source = @embedFile("test_data/multi_model_with_h.cif");
    var doc = try cif.readString(testing.allocator, source);
    defer doc.deinit();

    // Strip H from doc first (like processFileMmcif does)
    const block = &doc.blocks.items[0];
    const loop = block.findLoopMut("_atom_site.type_symbol").?;
    const type_col = loop.findTag("_atom_site.type_symbol").?;
    const w = loop.width();
    const n_rows = loop.length();
    var write: usize = 0;
    var read: usize = 0;
    while (read < n_rows) : (read += 1) {
        const type_sym = loop.val(read, type_col) orelse {
            if (write != read) {
                for (0..w) |col| {
                    loop.values.items[write * w + col] = loop.values.items[read * w + col];
                }
            }
            write += 1;
            continue;
        };
        if (std.ascii.eqlIgnoreCase(type_sym, "H") or std.ascii.eqlIgnoreCase(type_sym, "D")) continue;
        if (write != read) {
            for (0..w) |col| {
                loop.values.items[write * w + col] = loop.values.items[read * w + col];
            }
        }
        write += 1;
    }
    loop.values.items.len = write * w;

    // Now parse from the stripped block
    var entries = try parseModelsFromBlock(testing.allocator, block, .all);
    defer {
        for (entries.items) |*e| e.model.deinit();
        entries.deinit(testing.allocator);
    }

    // 6 rows total, 2 H stripped → 4 heavy atoms: 2 per model
    try testing.expectEqual(@as(usize, 2), entries.items.len);
    try testing.expectEqual(@as(usize, 2), entries.items[0].model.atoms.items.len);
    try testing.expectEqual(@as(usize, 2), entries.items[1].model.atoms.items.len);

    // cif_row indices should be consistent with the stripped loop (4 rows total)
    try testing.expectEqual(@as(u32, 0), entries.items[0].cif_row_start);
    try testing.expectEqual(@as(u32, 2), entries.items[0].cif_row_end);
    try testing.expectEqual(@as(u32, 2), entries.items[1].cif_row_start);
    try testing.expectEqual(@as(u32, 4), entries.items[1].cif_row_end);

    // Coordinates should be from the correct models
    try testing.expectApproxEqAbs(@as(f32, 1.0), entries.items[0].model.atoms.items[0].pos.x, 0.01);
    try testing.expectApproxEqAbs(@as(f32, 11.0), entries.items[1].model.atoms.items[0].pos.x, 0.01);
}

test "entityTypeFromString maps correctly" {
    try testing.expectEqual(EntityType.polymer, entityTypeFromString("polymer"));
    try testing.expectEqual(EntityType.polymer, entityTypeFromString("POLYMER"));
    try testing.expectEqual(EntityType.non_polymer, entityTypeFromString("non-polymer"));
    try testing.expectEqual(EntityType.branched, entityTypeFromString("branched"));
    try testing.expectEqual(EntityType.water, entityTypeFromString("water"));
    try testing.expectEqual(EntityType.unknown, entityTypeFromString("macrolide"));
    try testing.expectEqual(EntityType.unknown, entityTypeFromString(""));
}

test {
    // Pull in sub-module tests
    _ = @import("mmcif/conn.zig");
    _ = @import("mmcif/inline_comp.zig");
}
