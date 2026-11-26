//! PDB format parser: parses ATOM/HETATM records into a Model.
//!
//! Supports multi-model PDB files via parseAll(). The single-model
//! parse() function extracts only the first model.
//! CONECT and MASTER records are dropped entirely.
//! All other non-coordinate records are kept as raw_line passthrough.

const std = @import("std");
const Allocator = std.mem.Allocator;

const model_mod = @import("model.zig");
const element = @import("element.zig");
const math = @import("math.zig");
const mmcif_mod = @import("mmcif.zig");
pub const ModelFilter = mmcif_mod.ModelFilter;

const Model = model_mod.Model;
const Atom = model_mod.Atom;
const Residue = model_mod.Residue;
const Chain = model_mod.Chain;

pub const PdbError = error{
    OutOfMemory,
    InvalidCoordinateValue,
};

/// A parsed PDB record: either a coordinate atom or a raw passthrough line.
pub const PdbRecord = union(enum) {
    atom_site,
    raw_line: []const u8,
};

/// Full parse result including the model and passthrough records.
pub const PdbParseResult = struct {
    model: Model,
    records: std.ArrayListUnmanaged(PdbRecord),
    source: []const u8,

    pub fn deinit(self: *PdbParseResult, allocator: Allocator) void {
        self.model.deinit();
        self.records.deinit(allocator);
    }
};

/// A single model entry from a multi-model PDB file.
pub const PdbModelEntry = struct {
    model_num: u32,
    model: Model,
    /// Per-model records (atom_site markers). Header records are stored
    /// separately in PdbMultiModelResult.header_records.
    records: std.ArrayListUnmanaged(PdbRecord),

    pub fn deinit(self: *PdbModelEntry, allocator: Allocator) void {
        self.model.deinit();
        self.records.deinit(allocator);
    }
};

/// Result from parsing a multi-model PDB file.
pub const PdbMultiModelResult = struct {
    entries: std.ArrayListUnmanaged(PdbModelEntry),
    /// Records that appear before the first MODEL (HEADER, CRYST1, etc.)
    /// and after the last ENDMDL (END, etc.) — shared across models.
    header_records: std.ArrayListUnmanaged(PdbRecord),
    source: []const u8,

    pub fn deinit(self: *PdbMultiModelResult, allocator: Allocator) void {
        for (self.entries.items) |*e| e.deinit(allocator);
        self.entries.deinit(allocator);
        self.header_records.deinit(allocator);
    }
};

// ---------------------------------------------------------------------------
// Standard polymer residue names (20 AA + special + nucleotides)
// ---------------------------------------------------------------------------

const polymer_comp_ids = std.StaticStringMap(void).initComptime(.{
    // 20 standard amino acids
    .{ "ALA", {} }, .{ "ARG", {} }, .{ "ASN", {} }, .{ "ASP", {} },
    .{ "CYS", {} }, .{ "GLN", {} }, .{ "GLU", {} }, .{ "GLY", {} },
    .{ "HIS", {} }, .{ "ILE", {} }, .{ "LEU", {} }, .{ "LYS", {} },
    .{ "MET", {} }, .{ "PHE", {} }, .{ "PRO", {} }, .{ "SER", {} },
    .{ "THR", {} }, .{ "TRP", {} }, .{ "TYR", {} }, .{ "VAL", {} },
    // Special amino acids
    .{ "MSE", {} }, .{ "SEC", {} }, .{ "PYL", {} },
    // Common modified amino acids
    .{ "TPO", {} },
    .{ "SEP", {} }, .{ "PTR", {} }, .{ "MLY", {} }, .{ "CSO", {} },
    .{ "HYP", {} }, .{ "CSS", {} }, .{ "CME", {} }, .{ "CSD", {} },
    .{ "OCS", {} }, .{ "KCX", {} }, .{ "LLP", {} }, .{ "M3L", {} },
    .{ "ALY", {} }, .{ "CGU", {} },
    // Standard nucleotides (RNA)
    .{ "A", {} },   .{ "C", {} },
    .{ "G", {} },   .{ "U", {} },
    // Standard nucleotides (DNA)
      .{ "DA", {} },  .{ "DC", {} },
    .{ "DG", {} },  .{ "DT", {} },  .{ "DI", {} },
});

fn entityTypeFromCompId(comp_id: []const u8) model_mod.residue.EntityType {
    if (std.mem.eql(u8, comp_id, "HOH") or
        std.mem.eql(u8, comp_id, "WAT") or
        std.mem.eql(u8, comp_id, "H2O"))
    {
        return .water;
    }
    if (polymer_comp_ids.has(comp_id)) return .polymer;
    return .non_polymer;
}

// ---------------------------------------------------------------------------
// Line parsing helpers
// ---------------------------------------------------------------------------

/// Return line[start..end] if the line is long enough, else return "".
fn safeSlice(line: []const u8, start: usize, end: usize) []const u8 {
    if (start >= line.len) return "";
    const e = @min(end, line.len);
    return line[start..e];
}

/// Parse a float strictly — returns error for empty or malformed values.
fn parseFloatStrict(s: []const u8) PdbError!f32 {
    const trimmed = std.mem.trim(u8, s, " ");
    if (trimmed.len == 0) return PdbError.InvalidCoordinateValue;
    return std.fmt.parseFloat(f32, trimmed) catch return PdbError.InvalidCoordinateValue;
}

/// Parse a float with a default for empty fields. Non-empty malformed values return error.
fn parseFloatOpt(s: []const u8, default: f32) PdbError!f32 {
    const trimmed = std.mem.trim(u8, s, " ");
    if (trimmed.len == 0) return default;
    return std.fmt.parseFloat(f32, trimmed) catch return PdbError.InvalidCoordinateValue;
}

/// Parse an integer from a slice, trimming whitespace. Returns 0 on failure.
fn parseInt(comptime T: type, s: []const u8) T {
    const trimmed = std.mem.trim(u8, s, " ");
    if (trimmed.len == 0) return 0;
    return std.fmt.parseInt(T, trimmed, 10) catch 0;
}

/// Infer element from PDB atom name when columns 76-77 are absent.
/// PDB atom names are 4 chars (cols 12-15). PDB convention:
///   - Names starting with a space (col 13): 1-char element at col 14 (index 1).
///   - Names starting with a digit (col 13): hydrogen, element is H.
///   - Names starting with a letter (col 13): 2-char element at cols 13-14 (e.g. "FE  " -> Fe).
fn inferElementFromName(raw_name: []const u8) element.AtomType {
    // raw_name is the 4-char field from cols 12-15 (may have leading/trailing spaces)
    if (raw_name.len == 0) return .unknown;

    // Digit-prefixed names (e.g. "1HB ", "2HG1"): these are always hydrogen.
    // Use the untrimmed first char to detect the leading-digit convention.
    // Also trim to handle names shorter than 4 chars.
    const first = raw_name[0];
    if (std.ascii.isDigit(first)) {
        // Hydrogen with numeric prefix
        const trimmed = std.mem.trim(u8, raw_name, " ");
        // Strip the leading digit(s) and take 1 char for element
        var i: usize = 0;
        while (i < trimmed.len and std.ascii.isDigit(trimmed[i])) : (i += 1) {}
        if (i >= trimmed.len) return .unknown;
        return element.elementFromSymbol(trimmed[i .. i + 1]);
    }

    if (first == ' ') {
        // Space in column 13: 1-char element at column 14 (index 1)
        if (raw_name.len < 2) return .unknown;
        const ch = raw_name[1];
        if (ch == ' ') return .unknown;
        return element.elementFromSymbol(raw_name[1..2]);
    }

    // Non-space, non-digit first char: 2-char element starting at col 13
    // (e.g. "FE  " -> Fe, "CL  " -> Cl)
    const end = if (raw_name.len >= 2 and raw_name[1] != ' ') @as(usize, 2) else @as(usize, 1);
    return element.elementFromSymbol(raw_name[0..end]);
}

// ---------------------------------------------------------------------------
// Core parser state
// ---------------------------------------------------------------------------

const ParseState = struct {
    mdl: Model,
    allocator: Allocator,

    // Current chain tracking
    cur_chain_id: u8 = 0,
    cur_chain_idx: u32 = 0,
    in_chain: bool = false,
    next_entity_id: u32 = 1,

    // Current residue tracking
    cur_seq_id: i32 = 0,
    cur_ins_code: u8 = 0,
    cur_comp_id: [5]u8 = undefined,
    cur_comp_id_len: usize = 0,
    in_residue: bool = false,
    res_atom_start: u32 = 0,

    fn init(allocator: Allocator) ParseState {
        return .{
            .mdl = Model.init(allocator),
            .allocator = allocator,
        };
    }

    fn deinit(self: *ParseState) void {
        self.mdl.deinit();
    }

    fn curCompIdSlice(self: *const ParseState) []const u8 {
        return self.cur_comp_id[0..self.cur_comp_id_len];
    }

    fn setCompId(self: *ParseState, s: []const u8) void {
        const n = @min(s.len, 5);
        self.cur_comp_id_len = n;
        for (0..n) |i| self.cur_comp_id[i] = s[i];
    }

    /// Flush the in-progress residue (set atom_end) and append it.
    fn flushResidue(self: *ParseState) !void {
        if (!self.in_residue) return;
        const atom_end: u32 = @intCast(self.mdl.atoms.items.len);
        var res = Residue{
            .chain_idx = self.cur_chain_idx,
            .seq_id = self.cur_seq_id,
            .auth_seq_id = self.cur_seq_id,
            .ins_code = self.cur_ins_code,
            .atom_start = self.res_atom_start,
            .atom_end = atom_end,
            .entity_type = entityTypeFromCompId(self.curCompIdSlice()),
        };
        res.setCompId(self.curCompIdSlice());
        try self.mdl.residues.append(self.allocator, res);
        self.in_residue = false;
    }

    /// Flush the in-progress chain (set residue_end) and append it.
    fn flushChain(self: *ParseState) !void {
        if (!self.in_chain) return;
        try self.flushResidue();
        const res_end: u32 = @intCast(self.mdl.residues.items.len);
        self.mdl.chains.items[self.cur_chain_idx].residue_end = res_end;
        self.in_chain = false;
    }

    /// Open a new chain for `chain_id`. Does not flush the old chain; caller
    /// must call flushChain() first when transitioning between chains.
    fn openChain(self: *ParseState, chain_id: u8) !void {
        const chain_str = [_]u8{chain_id};
        var entity_buf: [8]u8 = undefined;
        const entity_str = std.fmt.bufPrint(&entity_buf, "{d}", .{self.next_entity_id}) catch "1";
        self.next_entity_id += 1;

        var ch = Chain{
            .residue_start = @intCast(self.mdl.residues.items.len),
            .residue_end = @intCast(self.mdl.residues.items.len),
        };
        ch.setLabelAsymId(&chain_str);
        ch.setAuthAsymId(&chain_str);
        ch.setEntityId(entity_str);

        self.cur_chain_idx = @intCast(self.mdl.chains.items.len);
        self.cur_chain_id = chain_id;
        self.in_chain = true;
        try self.mdl.chains.append(self.allocator, ch);
    }

    /// Process a single ATOM/HETATM line.
    fn processAtomLine(self: *ParseState, line: []const u8) !void {
        // Minimum length to reach end of Z coordinate (column 54).
        if (line.len < 54) return PdbError.InvalidCoordinateValue;

        // --- Extract fixed-width fields ---
        const raw_name = safeSlice(line, 12, 16); // 4 chars, space-padded
        const altloc_ch = if (line.len > 16) line[16] else ' ';
        const comp_id_raw = std.mem.trim(u8, safeSlice(line, 17, 20), " ");
        const chain_id_ch = if (line.len > 21) line[21] else ' ';
        const seq_id_str = safeSlice(line, 22, 26);
        const ins_code_ch: u8 = if (line.len > 26 and line[26] != ' ') line[26] else ' ';
        const x_str = safeSlice(line, 30, 38);
        const y_str = safeSlice(line, 38, 46);
        const z_str = safeSlice(line, 46, 54);
        const occ_str = safeSlice(line, 54, 60);
        const bfac_str = safeSlice(line, 60, 66);
        const elem_raw = safeSlice(line, 76, 78);

        const seq_id = parseInt(i32, seq_id_str);
        const x = try parseFloatStrict(x_str);
        const y = try parseFloatStrict(y_str);
        const z = try parseFloatStrict(z_str);

        // --- Chain transition ---
        if (!self.in_chain or chain_id_ch != self.cur_chain_id) {
            try self.flushChain();
            try self.openChain(chain_id_ch);
            self.in_residue = false;
        }

        // --- Residue transition ---
        const same_residue = self.in_residue and
            self.cur_seq_id == seq_id and
            self.cur_ins_code == ins_code_ch and
            std.mem.eql(u8, self.curCompIdSlice(), comp_id_raw);

        if (!same_residue) {
            try self.flushResidue();
            self.cur_seq_id = seq_id;
            self.cur_ins_code = ins_code_ch;
            self.setCompId(comp_id_raw);
            self.res_atom_start = @intCast(self.mdl.atoms.items.len);
            self.in_residue = true;
        }

        // --- Determine element ---
        const elem_trimmed = std.mem.trim(u8, elem_raw, " ");
        const atom_type = if (elem_trimmed.len > 0)
            element.elementFromSymbol(elem_trimmed)
        else
            inferElementFromName(raw_name);

        const is_h = switch (atom_type) {
            .H, .Har, .Hpol, .Ha_p, .HOd => true,
            else => false,
        };
        const vdw = atom_type.info().explicit_radius;

        // --- Trim atom name ---
        const name_trimmed = std.mem.trim(u8, raw_name, " ");

        var atm = Atom{
            .pos = .{ .x = x, .y = y, .z = z },
            .element_type = atom_type,
            .residue_idx = @intCast(self.mdl.residues.items.len),
            .altloc = altloc_ch,
            .occupancy = try parseFloatOpt(occ_str, 1.0),
            .b_factor = try parseFloatOpt(bfac_str, 0.0),
            .is_hydrogen = is_h,
            .vdw_radius = vdw,
        };
        atm.setName(name_trimmed);

        try self.mdl.atoms.append(self.allocator, atm);
    }

    /// Finalise after all lines are processed.
    fn finalize(self: *ParseState) !void {
        try self.flushChain();
        self.mdl.original_atom_count = @intCast(self.mdl.atoms.items.len);
        // Apply chain-break detection based on seq_id gaps
        try detectChainBreaks(&self.mdl);
    }
};

/// Detect chain breaks: for consecutive polymer residues in the same chain,
/// mark is_chain_break_before when seq_id gap > 1.
fn detectChainBreaks(mdl: *Model) !void {
    if (mdl.residues.items.len < 2) return;
    for (mdl.chains.items) |chain| {
        const res_slice = mdl.residues.items[chain.residue_start..chain.residue_end];
        if (res_slice.len < 2) continue;
        var prev_seq: i32 = res_slice[0].seq_id;
        var prev_is_polymer = res_slice[0].entity_type == .polymer;
        for (res_slice[1..]) |*res| {
            const is_poly = res.entity_type == .polymer;
            if (is_poly and prev_is_polymer and (res.seq_id - prev_seq) > 1) {
                res.is_chain_break_before = true;
            }
            prev_seq = res.seq_id;
            prev_is_polymer = is_poly;
        }
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Strip DOS line ending from a line slice (remove trailing \r).
fn stripCR(line: []const u8) []const u8 {
    if (line.len > 0 and line[line.len - 1] == '\r') return line[0 .. line.len - 1];
    return line;
}

/// Parse PDB source, keeping passthrough records. Returns a PdbParseResult.
pub fn parse(allocator: Allocator, source: []const u8) PdbError!PdbParseResult {
    var state = ParseState.init(allocator);
    errdefer state.deinit();

    var records = std.ArrayListUnmanaged(PdbRecord).empty;
    errdefer records.deinit(allocator);

    var in_model_block = false;
    var model_done = false;

    var iter = std.mem.splitScalar(u8, source, '\n');
    while (iter.next()) |raw_line| {
        const line = stripCR(raw_line);
        if (line.len == 0) continue;

        // Record type is cols 0-5 (up to 6 chars), left-justified
        const rec_type = safeSlice(line, 0, 6);

        // Handle MODEL/ENDMDL — only parse MODEL 1
        if (std.mem.startsWith(u8, rec_type, "MODEL ") or std.mem.eql(u8, std.mem.trim(u8, rec_type, " "), "MODEL")) {
            if (!model_done and !in_model_block) {
                in_model_block = true;
                // Parse model serial number from columns 10-14
                if (line.len >= 14) {
                    const num_str = std.mem.trim(u8, safeSlice(line, 10, 14), " ");
                    state.mdl.model_num = std.fmt.parseInt(u32, num_str, 10) catch 1;
                }
            } else {
                // Skip subsequent models; still passthrough the line
                try records.append(allocator, .{ .raw_line = line });
            }
            continue;
        }
        if (std.mem.startsWith(u8, rec_type, "ENDMDL")) {
            if (in_model_block) {
                in_model_block = false;
                model_done = true;
            }
            try records.append(allocator, .{ .raw_line = line });
            continue;
        }

        // If we've finished MODEL 1 and there are more atoms, skip them
        if (model_done and
            (std.mem.startsWith(u8, rec_type, "ATOM  ") or
                std.mem.startsWith(u8, rec_type, "HETATM")))
        {
            continue;
        }

        // Drop CONECT and MASTER records entirely
        if (std.mem.startsWith(u8, rec_type, "CONECT") or
            std.mem.startsWith(u8, rec_type, "MASTER"))
        {
            continue;
        }

        // TER record: flush current chain
        if (std.mem.startsWith(u8, rec_type, "TER")) {
            try state.flushChain();
            try records.append(allocator, .{ .raw_line = line });
            continue;
        }

        // Coordinate records
        if (std.mem.startsWith(u8, rec_type, "ATOM  ") or
            std.mem.startsWith(u8, rec_type, "HETATM"))
        {
            try state.processAtomLine(line);
            try records.append(allocator, .atom_site);
            continue;
        }

        // Everything else: passthrough
        try records.append(allocator, .{ .raw_line = line });
    }

    try state.finalize();

    return PdbParseResult{
        .model = state.mdl,
        .records = records,
        .source = source,
    };
}

/// Parse PDB source into a Model, discarding passthrough records.
pub fn parseModel(allocator: Allocator, source: []const u8) PdbError!Model {
    var result = try parse(allocator, source);
    result.records.deinit(allocator);
    return result.model;
}

/// Parse a multi-model PDB file, returning one PdbModelEntry per model.
/// For single-model files (no MODEL records), returns one entry.
pub fn parseAll(allocator: Allocator, source: []const u8, filter: ModelFilter) PdbError!PdbMultiModelResult {
    var entries = std.ArrayListUnmanaged(PdbModelEntry).empty;
    errdefer {
        for (entries.items) |*e| e.deinit(allocator);
        entries.deinit(allocator);
    }

    var header_records = std.ArrayListUnmanaged(PdbRecord).empty;
    errdefer header_records.deinit(allocator);

    var state = ParseState.init(allocator);
    errdefer state.deinit();

    var cur_records = std.ArrayListUnmanaged(PdbRecord).empty;
    errdefer cur_records.deinit(allocator);

    var in_model_block = false;
    var has_model_records = false; // true if any MODEL record was seen
    var cur_model_num: u32 = 1;

    var iter = std.mem.splitScalar(u8, source, '\n');
    while (iter.next()) |raw_line| {
        const line = stripCR(raw_line);
        if (line.len == 0) continue;

        const rec_type = safeSlice(line, 0, 6);

        // MODEL record: start a new model block
        if (std.mem.startsWith(u8, rec_type, "MODEL ") or std.mem.eql(u8, std.mem.trim(u8, rec_type, " "), "MODEL")) {
            has_model_records = true;

            // Parse model serial number from columns 10-14
            var model_num: u32 = 1;
            if (line.len >= 14) {
                const num_str = std.mem.trim(u8, safeSlice(line, 10, 14), " ");
                model_num = std.fmt.parseInt(u32, num_str, 10) catch 1;
            }

            // Check filter: should we process this model?
            const want = switch (filter) {
                .all => true,
                .first => entries.items.len == 0,
                .specific => |target| model_num == target,
            };

            if (want) {
                // Free previous state/records before overwriting
                // (handles consecutive MODEL without ENDMDL)
                state.deinit();
                cur_records.deinit(allocator);
                in_model_block = true;
                cur_model_num = model_num;
                state = ParseState.init(allocator);
                cur_records = .empty;
            }
            continue;
        }

        // ENDMDL: finalize current model
        if (std.mem.startsWith(u8, rec_type, "ENDMDL")) {
            if (in_model_block) {
                try state.finalize();
                state.mdl.model_num = cur_model_num;

                try entries.append(allocator, .{
                    .model_num = cur_model_num,
                    .model = state.mdl,
                    .records = cur_records,
                });

                // Reset for next model (don't deinit — ownership transferred)
                state = ParseState.init(allocator);
                cur_records = .empty;
                in_model_block = false;

                // For .first, stop after one model
                if (filter == .first) break;
                // For .specific, stop after finding the target
                switch (filter) {
                    .specific => {
                        if (entries.items.len > 0) break;
                    },
                    else => {},
                }
            }
            continue;
        }

        // Drop CONECT and MASTER records
        if (std.mem.startsWith(u8, rec_type, "CONECT") or
            std.mem.startsWith(u8, rec_type, "MASTER"))
        {
            continue;
        }

        // Coordinate records
        if (std.mem.startsWith(u8, rec_type, "ATOM  ") or
            std.mem.startsWith(u8, rec_type, "HETATM"))
        {
            if (in_model_block or !has_model_records) {
                try state.processAtomLine(line);
                try cur_records.append(allocator, .atom_site);
            }
            continue;
        }

        // TER record
        if (std.mem.startsWith(u8, rec_type, "TER")) {
            if (in_model_block or !has_model_records) {
                try state.flushChain();
                try cur_records.append(allocator, .{ .raw_line = line });
            }
            continue;
        }

        // Non-coordinate records before first MODEL → header
        if (!has_model_records and !in_model_block) {
            try header_records.append(allocator, .{ .raw_line = line });
        }
    }

    // For single-model files (no MODEL records), finalize and add the one model
    if (!has_model_records and state.mdl.atoms.items.len > 0) {
        try state.finalize();
        try entries.append(allocator, .{
            .model_num = 1,
            .model = state.mdl,
            .records = cur_records,
        });
        state = ParseState.init(allocator);
        cur_records = .empty;
    }

    // Clean up spare state
    state.deinit();
    cur_records.deinit(allocator);

    return PdbMultiModelResult{
        .entries = entries,
        .header_records = header_records,
        .source = source,
    };
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

const testing = std.testing;

test "parse tiny PDB" {
    const source = @embedFile("test_data/tiny.pdb");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();
    try testing.expectEqual(@as(usize, 5), mdl.atoms.items.len);
    try testing.expectEqual(@as(usize, 1), mdl.residues.items.len);
    try testing.expectEqual(@as(usize, 1), mdl.chains.items.len);
    try testing.expectApproxEqAbs(mdl.atoms.items[0].pos.x, 1.0, 1e-3);
    try testing.expectApproxEqAbs(mdl.atoms.items[0].pos.y, 2.0, 1e-3);
    try testing.expectEqualStrings("ALA", mdl.residues.items[0].compIdSlice());
    try testing.expectEqualStrings("N", mdl.atoms.items[0].nameSlice());
}

test "parse multi-chain PDB" {
    const source = @embedFile("test_data/multi_chain.pdb");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();
    try testing.expectEqual(@as(usize, 12), mdl.atoms.items.len);
    try testing.expectEqual(@as(usize, 3), mdl.residues.items.len);
    try testing.expectEqual(@as(usize, 2), mdl.chains.items.len);
    try testing.expectEqualStrings("A", mdl.chains.items[0].labelSlice());
    try testing.expectEqual(@as(u32, 0), mdl.chains.items[0].residue_start);
    try testing.expectEqual(@as(u32, 2), mdl.chains.items[0].residue_end);
    try testing.expectEqualStrings("B", mdl.chains.items[1].labelSlice());
    try testing.expectEqual(@as(u32, 2), mdl.chains.items[1].residue_start);
    try testing.expectEqual(@as(u32, 3), mdl.chains.items[1].residue_end);
    try testing.expectApproxEqAbs(mdl.atoms.items[11].pos.x, 13.0, 1e-3);
}

test "parse HETATM and entity types" {
    const source = @embedFile("test_data/hetatm.pdb");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();
    try testing.expectEqual(@as(usize, 9), mdl.atoms.items.len);
    try testing.expectEqual(@as(usize, 3), mdl.residues.items.len);
    // ALA -> polymer
    try testing.expectEqual(model_mod.residue.EntityType.polymer, mdl.residues.items[0].entity_type);
    // EDO -> non_polymer
    try testing.expectEqual(model_mod.residue.EntityType.non_polymer, mdl.residues.items[1].entity_type);
    // HOH -> water
    try testing.expectEqual(model_mod.residue.EntityType.water, mdl.residues.items[2].entity_type);
}

test "parse MODEL 1 only from multi-model PDB" {
    const source =
        \\MODEL        1
        \\ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 10.00           N
        \\ATOM      2  CA  ALA A   1       2.000   3.000   4.000  1.00 10.00           C
        \\ENDMDL
        \\MODEL        2
        \\ATOM      1  N   ALA A   1       9.000   9.000   9.000  1.00 10.00           N
        \\ATOM      2  CA  ALA A   1       8.000   8.000   8.000  1.00 10.00           C
        \\ENDMDL
        \\END
    ;
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();
    try testing.expectEqual(@as(usize, 2), mdl.atoms.items.len);
    // Coordinates should be from MODEL 1, not MODEL 2
    try testing.expectApproxEqAbs(mdl.atoms.items[0].pos.x, 1.0, 1e-3);
}

test "parseAll returns all models" {
    const source =
        \\MODEL        1
        \\ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 10.00           N
        \\ATOM      2  CA  ALA A   1       2.000   3.000   4.000  1.00 10.00           C
        \\ENDMDL
        \\MODEL        2
        \\ATOM      1  N   ALA A   1       9.000   9.000   9.000  1.00 10.00           N
        \\ATOM      2  CA  ALA A   1       8.000   8.000   8.000  1.00 10.00           C
        \\ENDMDL
        \\END
    ;
    var result = try parseAll(testing.allocator, source, .all);
    defer result.deinit(testing.allocator);

    try testing.expectEqual(@as(usize, 2), result.entries.items.len);
    try testing.expectEqual(@as(u32, 1), result.entries.items[0].model_num);
    try testing.expectEqual(@as(u32, 2), result.entries.items[1].model_num);
    try testing.expectEqual(@as(usize, 2), result.entries.items[0].model.atoms.items.len);
    try testing.expectEqual(@as(usize, 2), result.entries.items[1].model.atoms.items.len);
    try testing.expectApproxEqAbs(@as(f32, 1.0), result.entries.items[0].model.atoms.items[0].pos.x, 0.01);
    try testing.expectApproxEqAbs(@as(f32, 9.0), result.entries.items[1].model.atoms.items[0].pos.x, 0.01);
}

test "parseAll with specific filter" {
    const source =
        \\MODEL        1
        \\ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 10.00           N
        \\ENDMDL
        \\MODEL        2
        \\ATOM      1  N   ALA A   1       9.000   9.000   9.000  1.00 10.00           N
        \\ENDMDL
    ;
    var result = try parseAll(testing.allocator, source, .{ .specific = 2 });
    defer result.deinit(testing.allocator);

    try testing.expectEqual(@as(usize, 1), result.entries.items.len);
    try testing.expectEqual(@as(u32, 2), result.entries.items[0].model_num);
    try testing.expectApproxEqAbs(@as(f32, 9.0), result.entries.items[0].model.atoms.items[0].pos.x, 0.01);
}

test "parseAll single-model PDB (no MODEL records)" {
    const source =
        \\ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 10.00           N
        \\END
    ;
    var result = try parseAll(testing.allocator, source, .all);
    defer result.deinit(testing.allocator);

    try testing.expectEqual(@as(usize, 1), result.entries.items.len);
    try testing.expectEqual(@as(u32, 1), result.entries.items[0].model_num);
}

test "parse altloc" {
    const source =
        \\ATOM      1  CA AALA A   1       1.000   2.000   3.000  1.00 10.00           C
        \\ATOM      2  CA BALA A   1       1.500   2.500   3.500  0.50 10.00           C
        \\END
    ;
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();
    try testing.expectEqual(@as(usize, 2), mdl.atoms.items.len);
    try testing.expectEqual(@as(u8, 'A'), mdl.atoms.items[0].altloc);
    try testing.expectEqual(@as(u8, 'B'), mdl.atoms.items[1].altloc);
}

test "parse insertion code splits residues" {
    const source =
        \\ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 10.00           N
        \\ATOM      2  N   GLY A   1A      2.000   3.000   4.000  1.00 10.00           N
        \\END
    ;
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();
    try testing.expectEqual(@as(usize, 2), mdl.residues.items.len);
    try testing.expectEqual(@as(u8, ' '), mdl.residues.items[0].ins_code);
    try testing.expectEqual(@as(u8, 'A'), mdl.residues.items[1].ins_code);
}

test "chain break detection from seq_id gap" {
    const source =
        \\ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 10.00           N
        \\ATOM      2  N   GLY A   5       2.000   3.000   4.000  1.00 10.00           N
        \\END
    ;
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();
    try testing.expectEqual(@as(usize, 2), mdl.residues.items.len);
    try testing.expect(!mdl.residues.items[0].is_chain_break_before);
    try testing.expect(mdl.residues.items[1].is_chain_break_before);
}

test "inferElementFromName" {
    // Standard: " CA " -> C
    try testing.expectEqual(element.AtomType.C, inferElementFromName(" CA "));
    // Hydrogen with digit prefix: "1HB " -> H
    try testing.expectEqual(element.AtomType.H, inferElementFromName("1HB "));
    // Iron: "FE  " -> Fe
    try testing.expectEqual(element.AtomType.Fe, inferElementFromName("FE  "));
    // Single letter: " N  " -> N
    try testing.expectEqual(element.AtomType.N, inferElementFromName(" N  "));
}

test "invalid coordinate returns error" {
    const source =
        \\ATOM      1  N   ALA A   1       abc     2.000   3.000  1.00 10.00           N
    ;
    const result = parseModel(testing.allocator, source);
    try testing.expectError(PdbError.InvalidCoordinateValue, result);
}

test "short ATOM line returns error" {
    const source = "ATOM      1  N   ALA A   1\n";
    const result = parseModel(testing.allocator, source);
    try testing.expectError(PdbError.InvalidCoordinateValue, result);
}
