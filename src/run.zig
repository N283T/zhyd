//! Single-file processing pipeline: parse, place, optimize, validate, write.

const std = @import("std");
const Allocator = std.mem.Allocator;

fn defaultIo() std.Io {
    return std.Io.Threaded.global_single_threaded.io();
}
const zreduce = @import("root.zig");

pub const ProcessConfig = struct {
    input_path: []const u8,
    output_path: ?[]const u8 = null, // must be non-null when called from batch workers (stdout is not thread-safe)
    dict: ?*const zreduce.ccd.ComponentDict = null, // shared read-only; caller owns lifetime
    sdf_dict: ?*const zreduce.ccd.ComponentDict = null, // SDF-derived topology; caller owns lifetime
    json_path: ?[]const u8 = null,
    state_json_path: ?[]const u8 = null,
    json_version: []const u8 = "", // version string for JSON log (passed from main)
    no_opt: bool = false,
    no_flip: bool = false,
    validate_flag: bool = false,
    opt_threads: u32 = 0, // 0 = auto; batch sets to 1
    quiet: bool = false, // suppress diagnostic prints (batch mode)
    water: zreduce.place.WaterConfig = .{},
    bond_policy: zreduce.place.BondPolicy = .{},
    nterm_mode: zreduce.place.NtermMode = .auto,
    protonation_path: ?[]const u8 = null,
    fix_path: ?[]const u8 = null,
    dump_movers_path: ?[]const u8 = null,
    strip_h: bool = false,
    format: InputFormat = .mmcif,
    model_filter: ModelFilter = .all,
};

pub const ProcessResult = struct {
    n_placed: u32,
    n_residues: u32,
    n_skipped_existing: u32,
    n_skipped_inter_residue: u32,
    n_skipped_missing_ref: u32,
    n_skipped_quality_filter: u32 = 0,
    n_distance_derived: u32 = 0,
    n_movers: u32 = 0,
    n_singletons: u32 = 0,
    n_brute_force: u32 = 0,
    n_vertex_cut: u32 = 0,

    pub fn totalSkipped(self: ProcessResult) u32 {
        return self.n_skipped_existing + self.n_skipped_inter_residue + self.n_skipped_missing_ref + self.n_skipped_quality_filter;
    }
};

pub const InputFormat = enum {
    mmcif,
    pdb,
};

pub const ModelFilter = zreduce.mmcif.ModelFilter;

/// Detect input format from file extension.
pub fn detectFormat(path: []const u8) InputFormat {
    if (std.mem.endsWith(u8, path, ".pdb") or
        std.mem.endsWith(u8, path, ".pdb.gz") or
        std.mem.endsWith(u8, path, ".ent") or
        std.mem.endsWith(u8, path, ".ent.gz"))
    {
        return .pdb;
    }
    return .mmcif;
}

pub fn readFile(allocator: Allocator, path: []const u8) ![]const u8 {
    const io = defaultIo();
    if (std.mem.endsWith(u8, path, ".gz")) {
        return zreduce.gzip.readGzip(allocator, path);
    }
    return try std.Io.Dir.cwd().readFileAlloc(io, path, allocator, .limited(1024 * 1024 * 1024));
}

/// Remove hydrogen rows from the CIF document's _atom_site loop.
/// This keeps the document in sync with the model after stripHydrogens().
pub fn stripDocumentHydrogens(block: *zreduce.cif.Block) void {
    const loop = block.findLoopMut("_atom_site.type_symbol") orelse return;
    const type_col = loop.findTag("_atom_site.type_symbol") orelse return;
    const w = loop.width();
    const n_rows = loop.length();

    var write: usize = 0;
    var read: usize = 0;
    while (read < n_rows) : (read += 1) {
        const type_sym = loop.val(read, type_col) orelse {
            // Missing type_symbol — not hydrogen, keep the row
            if (write != read) {
                for (0..w) |col| {
                    loop.values.items[write * w + col] = loop.values.items[read * w + col];
                }
            }
            write += 1;
            continue;
        };
        const is_h = std.ascii.eqlIgnoreCase(type_sym, "H") or
            std.ascii.eqlIgnoreCase(type_sym, "D");
        if (is_h) continue;

        if (write != read) {
            for (0..w) |col| {
                loop.values.items[write * w + col] = loop.values.items[read * w + col];
            }
        }
        write += 1;
    }
    loop.values.items.len = write * w;
}

/// Shared per-model processing: chemistry, place, movers/optimize, markAbsent, validate, snapshot.
/// Called by both mmCIF and PDB pipelines. Format-specific steps (bond parsing, strip_h) are
/// performed by the caller before invoking this function.
///
/// Parameters:
///   mdl                   - the model to process (atoms already stripped if needed)
///   model_num             - 1-based model number (for diagnostics and snapshots)
///   n_models              - total model count in file (to decide whether to print progress)
///   inline_dict           - CIF inline component dictionary (null for PDB)
///   protonation_overrides - optional protonation overrides (caller owns, may be null)
///   fix_overrides         - optional fix overrides (caller owns, may be null)
///   mover_snapshots       - accumulator list for JSON log snapshots (appended to)
///   result                - accumulator for ProcessResult fields (incremented in place)
///   config                - full ProcessConfig for opt/validate/dump flags
fn processModelShared(
    allocator: Allocator,
    mdl: *zreduce.model.Model,
    model_num: u32,
    n_models: usize,
    inline_dict: ?*const zreduce.ccd.ComponentDict,
    protonation_overrides: ?*zreduce.place.ProtonationOverrides,
    fix_overrides: ?*zreduce.optimize.fix.FixOverrides,
    mover_snapshots: *std.ArrayListUnmanaged(zreduce.writer.json_writer.MoverSnapshot),
    result: *ProcessResult,
    config: ProcessConfig,
) !void {
    if (!config.quiet and n_models > 1) {
        std.debug.print("  Processing model {d} ({d} atoms)\n", .{ model_num, mdl.atoms.items.len });
    }

    // Apply chemistry annotations
    zreduce.place.applyChemistryWithConfig(mdl, .{
        .protonation = protonation_overrides,
        .nterm_mode = config.nterm_mode,
    });

    // Place hydrogens
    const place_result = try zreduce.place.addHydrogensWithConfig(
        mdl,
        config.dict,
        inline_dict,
        .{
            .water = config.water,
            .bond_policy = config.bond_policy,
            .protonation = protonation_overrides,
            .nterm_mode = config.nterm_mode,
            .sdf_dict = config.sdf_dict,
        },
    );

    result.n_placed += place_result.n_placed;
    result.n_residues += place_result.n_residues;
    result.n_skipped_existing += place_result.n_skipped_existing;
    result.n_skipped_inter_residue += place_result.n_skipped_inter_residue;
    result.n_skipped_missing_ref += place_result.n_skipped_missing_ref;
    result.n_skipped_quality_filter += place_result.n_skipped_quality_filter;
    result.n_distance_derived += place_result.n_distance_derived;

    // Generate movers and optimize
    const needs_movers = !config.no_opt or config.fix_path != null or config.dump_movers_path != null;
    if (needs_movers) {
        const gen_result = try zreduce.optimize.generateMovers(
            allocator,
            mdl,
            config.no_flip,
            config.dict,
            inline_dict,
            protonation_overrides,
            config.bond_policy.mode,
        );
        var movers = gen_result.movers;
        defer {
            for (0..movers.len) |i| movers[i].deinit();
            allocator.free(movers);
        }
        result.n_movers += @intCast(movers.len);

        if (fix_overrides) |ov| {
            try zreduce.optimize.fix.applyFixes(ov, mdl, movers);
            for (movers) |*m| {
                if (m.is_fixed) m.applyOrientation(mdl.atoms.items, m.best_orientation);
            }
            if (!config.quiet) ov.warnUnmatched(mdl, movers);
        }

        if (config.dump_movers_path) |dump_path| {
            const io = defaultIo();
            var dump_buf: [4096]u8 = undefined;
            const dump_file = try std.Io.Dir.cwd().createFile(io, dump_path, .{});
            defer dump_file.close(io);
            var dump_fw = dump_file.writer(io, &dump_buf);
            try zreduce.optimize.fix.dumpMovers(&dump_fw.interface, mdl, movers);
            try dump_fw.interface.flush();
        }

        if (!config.no_opt and movers.len > 0) {
            const opt_result = try zreduce.optimize.optimizer.optimize(
                allocator,
                movers,
                mdl,
                .{ .n_threads = config.opt_threads },
            );
            result.n_singletons += opt_result.n_singletons;
            result.n_brute_force += opt_result.n_brute_force;
            result.n_vertex_cut += opt_result.n_vertex_cut;
        }

        // Capture mover snapshots for JSON log before movers are freed
        if (config.json_path != null) {
            try mover_snapshots.ensureUnusedCapacity(allocator, movers.len);
            for (movers) |m| {
                mover_snapshots.appendAssumeCapacity(
                    zreduce.writer.json_writer.MoverSnapshot.capture(m, mdl.residues.items, mdl.chains.items, model_num),
                );
            }
        }
    }

    // Mark absent H atoms
    markAbsentHydrogens(mdl);

    // Validate
    {
        var validation = try zreduce.validate.validateModel(allocator, mdl);
        defer validation.deinit();

        if (!validation.ok()) {
            if (!config.quiet) std.debug.print("  Model {d}: {d} validation issue(s)\n", .{ model_num, validation.issues.len });
            if (config.validate_flag) {
                zreduce.validate.reportIssues(validation.issues, mdl);
            }
        }
    }
}

fn markAbsentHydrogens(mdl: *zreduce.model.Model) void {
    for (mdl.atoms.items) |*atom| {
        if (zreduce.optimize.mover.isAbsentH(atom.*)) {
            atom.is_added = false;
        }
    }
}

fn countAddedHydrogens(mdl: *const zreduce.model.Model) u32 {
    var count: u32 = 0;
    for (mdl.atoms.items) |atom| {
        if (atom.is_added and atom.is_hydrogen) count += 1;
    }
    return count;
}

/// Process a single structure file (mmCIF or PDB) through the full pipeline.
/// Errors are returned (not exit'd) -- the caller decides how to handle them.
pub fn processFile(allocator: Allocator, config: ProcessConfig) !ProcessResult {
    // 1. Read input
    const source = try readFile(allocator, config.input_path);
    defer allocator.free(source);

    return switch (config.format) {
        .mmcif => processFileMmcif(allocator, config, source),
        .pdb => processFilePdb(allocator, config, source),
    };
}

/// mmCIF pipeline: supports multi-model processing.
fn processFileMmcif(allocator: Allocator, config: ProcessConfig, source: []const u8) !ProcessResult {
    // 2. Parse CIF document (for preserving non-atom_site categories in output)
    var doc = try zreduce.cif.readString(allocator, source);
    defer doc.deinit();

    // 2a. Strip existing hydrogens from the document BEFORE parsing models.
    // This ensures cif_row_start/end in ModelEntry are consistent with the
    // (possibly stripped) document loop used by the writer and atom lookup.
    if (doc.blocks.items.len == 0) {
        if (!config.quiet) std.debug.print("  No CIF data blocks found\n", .{});
        return ProcessResult{
            .n_placed = 0,
            .n_residues = 0,
            .n_skipped_existing = 0,
            .n_skipped_inter_residue = 0,
            .n_skipped_missing_ref = 0,
        };
    }

    if (config.strip_h) {
        stripDocumentHydrogens(&doc.blocks.items[0]);
    }

    // 3. Extract models from CIF block (shares the same doc — no double parse).
    // When strip_h is active, models are parsed from the already-stripped loop,
    // so they contain no H atoms and row indices match the stripped document.
    var entries = try zreduce.mmcif.parseModelsFromBlock(allocator, &doc.blocks.items[0], config.model_filter);
    defer {
        for (entries.items) |*e| e.model.deinit();
        entries.deinit(allocator);
    }

    if (entries.items.len == 0) {
        if (!config.quiet) std.debug.print("  No models found matching filter\n", .{});
        return ProcessResult{
            .n_placed = 0,
            .n_residues = 0,
            .n_skipped_existing = 0,
            .n_skipped_inter_residue = 0,
            .n_skipped_missing_ref = 0,
        };
    }

    if (config.strip_h and !config.quiet) {
        // Count how many H atoms were in the original source (before strip)
        // by comparing original_atom_count with a hypothetical full parse.
        // Since we stripped the doc before parsing, the models already lack H.
        std.debug.print("  Stripped existing H atoms from document\n", .{});
    }

    // 3b. Parse inline components (once — model-independent)
    var inline_dict: ?zreduce.ccd.ComponentDict = null;
    defer if (inline_dict) |*d| d.deinit();
    const block = &doc.blocks.items[0];
    inline_dict = try zreduce.mmcif.parseInlineComponents(allocator, block);

    // Parse overrides (once — shared across models)
    var protonation_overrides: ?zreduce.place.ProtonationOverrides = null;
    defer if (protonation_overrides) |*ov| ov.deinit();
    if (config.protonation_path) |path| {
        protonation_overrides = zreduce.place.protonation.parseFile(allocator, path) catch |err| {
            std.debug.print("Error: failed to load protonation override file '{s}': {s}\n", .{ path, @errorName(err) });
            return err;
        };
    }

    var fix_overrides: ?zreduce.optimize.fix.FixOverrides = null;
    defer if (fix_overrides) |*ov| ov.deinit();
    if (config.fix_path) |path| {
        fix_overrides = zreduce.optimize.fix.parseFile(allocator, path) catch |err| {
            std.debug.print("Error: failed to load fix override file '{s}': {s}\n", .{ path, @errorName(err) });
            return err;
        };
    }

    // Process each model independently
    var result = ProcessResult{
        .n_placed = 0,
        .n_residues = 0,
        .n_skipped_existing = 0,
        .n_skipped_inter_residue = 0,
        .n_skipped_missing_ref = 0,
    };

    // Accumulate mover snapshots for JSON log (movers are freed per model)
    var mover_snapshots = std.ArrayListUnmanaged(zreduce.writer.json_writer.MoverSnapshot).empty;
    defer mover_snapshots.deinit(allocator);

    for (entries.items) |*entry| {
        const mdl = &entry.model;

        // Build per-model atom lookup for bond parsing (mmCIF-specific)
        var atom_lookup = try zreduce.mmcif.buildAtomLookupForRange(allocator, block, entry.cif_row_start, entry.cif_row_end);
        defer atom_lookup.deinit();

        // mmCIF-specific: parse bonds from struct_conn/branch_link and flag leaving atoms
        try zreduce.mmcif.parseStructConn(mdl, block, &atom_lookup);
        try zreduce.mmcif.parseBranchLinks(allocator, mdl, block, &atom_lookup);
        zreduce.mmcif.flagLeavingAtoms(mdl, if (inline_dict) |*d| d else null, config.dict);

        // Shared per-model pipeline: chemistry → place → movers/opt → markAbsent → validate → snapshot
        try processModelShared(
            allocator,
            mdl,
            entry.model_num,
            entries.items.len,
            if (inline_dict) |*d| d else null,
            if (protonation_overrides) |*ov| ov else null,
            if (fix_overrides) |*ov| ov else null,
            &mover_snapshots,
            &result,
            config,
        );
    }

    // Warn about unmatched overrides (after all models processed)
    if (protonation_overrides) |*ov| {
        if (!config.quiet) {
            // Use first model for name resolution (overrides are file-level, not per-model)
            ov.warnUnmatched(&entries.items[0].model);
        }
    }

    // 9. Write output (all models into one file)
    const io = defaultIo();
    var out_buf: [4096]u8 = undefined;
    if (config.output_path) |out_path| {
        if (std.mem.endsWith(u8, out_path, ".gz")) {
            var gw = try zreduce.gzip.GzipWriter.init(allocator, out_path);
            errdefer gw.close() catch {};
            try zreduce.writer.mmcif_writer.writeMultiModelWithDocumentWithPolicy(gw.writer(), entries.items, &doc, config.bond_policy);
            try gw.close();
        } else {
            const file = try std.Io.Dir.cwd().createFile(io, out_path, .{});
            defer file.close(io);
            var fw = file.writer(io, &out_buf);
            try zreduce.writer.mmcif_writer.writeMultiModelWithDocumentWithPolicy(&fw.interface, entries.items, &doc, config.bond_policy);
            try fw.interface.flush();
        }
    } else {
        const stdout = std.Io.File.stdout();
        var sw = stdout.writer(io, &out_buf);
        try zreduce.writer.mmcif_writer.writeMultiModelWithDocumentWithPolicy(&sw.interface, entries.items, &doc, config.bond_policy);
        try sw.interface.flush();
    }

    // 10. Write JSON log (optional, all models)
    if (config.json_path) |json_path| {
        var json_buf: [4096]u8 = undefined;
        const file = try std.Io.Dir.cwd().createFile(io, json_path, .{});
        defer file.close(io);
        var jw = file.writer(io, &json_buf);
        var total_added: u32 = 0;
        for (entries.items) |entry| total_added += countAddedHydrogens(&entry.model);
        try zreduce.writer.json_writer.writeMultiModelLog(
            &jw.interface,
            config.json_version,
            config.input_path,
            total_added,
            config.bond_policy,
            mover_snapshots.items,
        );
        try jw.interface.flush();
    }

    if (config.state_json_path) |state_path| {
        const reports = try allocator.alloc(zreduce.writer.protonation_json.StateReportModel, entries.items.len);
        defer allocator.free(reports);
        for (entries.items, 0..) |*entry, i| reports[i] = .{ .model_num = entry.model_num, .model = &entry.model };
        try writeStateReportFile(
            state_path,
            config.input_path,
            config.output_path,
            reports,
            if (protonation_overrides) |*ov| ov else null,
            config.nterm_mode,
        );
    }

    return result;
}

/// PDB pipeline: supports multi-model processing.
fn processFilePdb(allocator: Allocator, config: ProcessConfig, source: []const u8) !ProcessResult {
    var pdb_result = try zreduce.pdb.parseAll(allocator, source, config.model_filter);
    defer pdb_result.deinit(allocator);

    if (pdb_result.entries.items.len == 0) {
        if (!config.quiet) std.debug.print("  No models found matching filter\n", .{});
        return ProcessResult{
            .n_placed = 0,
            .n_residues = 0,
            .n_skipped_existing = 0,
            .n_skipped_inter_residue = 0,
            .n_skipped_missing_ref = 0,
        };
    }

    // Parse overrides (once — shared across models)
    var protonation_overrides: ?zreduce.place.ProtonationOverrides = null;
    defer if (protonation_overrides) |*ov| ov.deinit();
    if (config.protonation_path) |path| {
        protonation_overrides = zreduce.place.protonation.parseFile(allocator, path) catch |err| {
            std.debug.print("Error: failed to load protonation override file '{s}': {s}\n", .{ path, @errorName(err) });
            return err;
        };
    }

    var fix_overrides: ?zreduce.optimize.fix.FixOverrides = null;
    defer if (fix_overrides) |*ov| ov.deinit();
    if (config.fix_path) |path| {
        fix_overrides = zreduce.optimize.fix.parseFile(allocator, path) catch |err| {
            std.debug.print("Error: failed to load fix override file '{s}': {s}\n", .{ path, @errorName(err) });
            return err;
        };
    }

    var result = ProcessResult{
        .n_placed = 0,
        .n_residues = 0,
        .n_skipped_existing = 0,
        .n_skipped_inter_residue = 0,
        .n_skipped_missing_ref = 0,
    };

    var mover_snapshots = std.ArrayListUnmanaged(zreduce.writer.json_writer.MoverSnapshot).empty;
    defer mover_snapshots.deinit(allocator);

    for (pdb_result.entries.items) |*entry| {
        const mdl = &entry.model;

        // PDB-specific: strip existing H per model (mmCIF strips at document level)
        if (config.strip_h) {
            const n_stripped = mdl.stripHydrogens();
            if (!config.quiet and n_stripped > 0) {
                std.debug.print("  Model {d}: stripped {d} existing H atoms\n", .{ entry.model_num, n_stripped });
            }
        }

        // Shared per-model pipeline: chemistry → place → movers/opt → markAbsent → validate → snapshot
        try processModelShared(
            allocator,
            mdl,
            entry.model_num,
            pdb_result.entries.items.len,
            null, // no inline_dict for PDB
            if (protonation_overrides) |*ov| ov else null,
            if (fix_overrides) |*ov| ov else null,
            &mover_snapshots,
            &result,
            config,
        );
    }

    if (protonation_overrides) |*ov| {
        if (!config.quiet) ov.warnUnmatched(&pdb_result.entries.items[0].model);
    }

    // Write output
    const io = defaultIo();
    var out_buf: [4096]u8 = undefined;
    if (config.output_path) |out_path| {
        if (std.mem.endsWith(u8, out_path, ".gz")) {
            var gw = try zreduce.gzip.GzipWriter.init(allocator, out_path);
            errdefer gw.close() catch {};
            try zreduce.writer.pdb_writer.writeMultiModel(gw.writer(), pdb_result.entries.items, pdb_result.header_records.items, config.bond_policy.output_isotope);
            try gw.close();
        } else {
            const file = try std.Io.Dir.cwd().createFile(io, out_path, .{});
            defer file.close(io);
            var fw = file.writer(io, &out_buf);
            try zreduce.writer.pdb_writer.writeMultiModel(&fw.interface, pdb_result.entries.items, pdb_result.header_records.items, config.bond_policy.output_isotope);
            try fw.interface.flush();
        }
    } else {
        const stdout = std.Io.File.stdout();
        var sw = stdout.writer(io, &out_buf);
        try zreduce.writer.pdb_writer.writeMultiModel(&sw.interface, pdb_result.entries.items, pdb_result.header_records.items, config.bond_policy.output_isotope);
        try sw.interface.flush();
    }

    if (config.json_path) |json_path| {
        var json_buf: [4096]u8 = undefined;
        const file = try std.Io.Dir.cwd().createFile(io, json_path, .{});
        defer file.close(io);
        var jw = file.writer(io, &json_buf);
        var total_added: u32 = 0;
        for (pdb_result.entries.items) |entry| total_added += countAddedHydrogens(&entry.model);
        try zreduce.writer.json_writer.writeMultiModelLog(
            &jw.interface,
            config.json_version,
            config.input_path,
            total_added,
            config.bond_policy,
            mover_snapshots.items,
        );
        try jw.interface.flush();
    }

    if (config.state_json_path) |state_path| {
        const reports = try allocator.alloc(zreduce.writer.protonation_json.StateReportModel, pdb_result.entries.items.len);
        defer allocator.free(reports);
        for (pdb_result.entries.items, 0..) |*entry, i| reports[i] = .{ .model_num = entry.model_num, .model = &entry.model };
        try writeStateReportFile(
            state_path,
            config.input_path,
            config.output_path,
            reports,
            if (protonation_overrides) |*ov| ov else null,
            config.nterm_mode,
        );
    }

    return result;
}

fn writeStateReportFile(
    path: []const u8,
    input_path: []const u8,
    coordinate_path: ?[]const u8,
    reports: []const zreduce.writer.protonation_json.StateReportModel,
    overrides: ?*const zreduce.place.ProtonationOverrides,
    nterm_mode: zreduce.place.NtermMode,
) !void {
    const io = defaultIo();
    var buffer: [4096]u8 = undefined;
    const file = try std.Io.Dir.cwd().createFile(io, path, .{});
    defer file.close(io);
    var fw = file.writer(io, &buffer);
    try zreduce.writer.protonation_json.writeStateReport(&fw.interface, input_path, coordinate_path, reports, overrides, nterm_mode);
    try fw.interface.flush();
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "detectFormat identifies PDB extensions" {
    try std.testing.expectEqual(InputFormat.pdb, detectFormat("input.pdb"));
    try std.testing.expectEqual(InputFormat.pdb, detectFormat("input.pdb.gz"));
    try std.testing.expectEqual(InputFormat.pdb, detectFormat("input.ent"));
    try std.testing.expectEqual(InputFormat.pdb, detectFormat("input.ent.gz"));
    try std.testing.expectEqual(InputFormat.mmcif, detectFormat("input.cif"));
    try std.testing.expectEqual(InputFormat.mmcif, detectFormat("input.cif.gz"));
    try std.testing.expectEqual(InputFormat.mmcif, detectFormat("input.mmcif"));
}

test "detectFormat defaults to mmcif for unknown/ambiguous extensions" {
    // Ambiguous extension: no path match → mmcif default
    try std.testing.expectEqual(InputFormat.mmcif, detectFormat("structure.txt"));
    // No extension at all
    try std.testing.expectEqual(InputFormat.mmcif, detectFormat("noextension"));
    // .mmcif extension
    try std.testing.expectEqual(InputFormat.mmcif, detectFormat("data.mmcif"));
    // Uppercase: detectFormat is case-sensitive per std.mem.endsWith, so ".PDB" is mmcif
    try std.testing.expectEqual(InputFormat.mmcif, detectFormat("input.PDB"));
    // Partial match: "something.pdb.bak" is not a PDB extension
    try std.testing.expectEqual(InputFormat.mmcif, detectFormat("input.pdb.bak"));
}

test "ProcessResult.totalSkipped aggregates all skip categories" {
    const r = ProcessResult{
        .n_placed = 100,
        .n_residues = 50,
        .n_skipped_existing = 3,
        .n_skipped_inter_residue = 7,
        .n_skipped_missing_ref = 5,
        .n_skipped_quality_filter = 2,
    };
    try std.testing.expectEqual(@as(u32, 17), r.totalSkipped());
}

test "ProcessResult.totalSkipped is zero when nothing skipped" {
    const r = ProcessResult{
        .n_placed = 10,
        .n_residues = 2,
        .n_skipped_existing = 0,
        .n_skipped_inter_residue = 0,
        .n_skipped_missing_ref = 0,
    };
    try std.testing.expectEqual(@as(u32, 0), r.totalSkipped());
}

test "stripDocumentHydrogens removes H and D rows from _atom_site loop" {
    // Build a minimal CIF block with a 3-column _atom_site loop:
    //   type_symbol, label_atom_id, id
    // Rows: N(1), H(2), CA(3), D(4), C(5)
    // After strip: N(1), CA(3), C(5) — 3 rows remain.
    const allocator = std.testing.allocator;

    var block = zreduce.cif.Block{ .name = "TEST", .items = .empty };
    defer block.deinit(allocator);

    // Append a loop item
    var loop = zreduce.cif.Loop{};

    try loop.tags.append(allocator, "_atom_site.type_symbol");
    try loop.tags.append(allocator, "_atom_site.label_atom_id");
    try loop.tags.append(allocator, "_atom_site.id");

    // Row 1: N
    try loop.values.append(allocator, "N");
    try loop.values.append(allocator, "N");
    try loop.values.append(allocator, "1");
    // Row 2: H
    try loop.values.append(allocator, "H");
    try loop.values.append(allocator, "H");
    try loop.values.append(allocator, "2");
    // Row 3: C (heavy — keep)
    try loop.values.append(allocator, "C");
    try loop.values.append(allocator, "CA");
    try loop.values.append(allocator, "3");
    // Row 4: D (deuterium — remove)
    try loop.values.append(allocator, "D");
    try loop.values.append(allocator, "D");
    try loop.values.append(allocator, "4");
    // Row 5: C (heavy — keep)
    try loop.values.append(allocator, "C");
    try loop.values.append(allocator, "C");
    try loop.values.append(allocator, "5");

    try block.items.append(allocator, .{ .loop = loop });

    stripDocumentHydrogens(&block);

    const result_loop = block.findLoop("_atom_site.type_symbol").?;
    try std.testing.expectEqual(@as(usize, 3), result_loop.length());

    // Verify the surviving rows are the heavy atoms in original order
    const type_col = result_loop.findTag("_atom_site.type_symbol").?;
    const id_col = result_loop.findTag("_atom_site.id").?;
    try std.testing.expectEqualStrings("N", result_loop.val(0, type_col).?);
    try std.testing.expectEqualStrings("1", result_loop.val(0, id_col).?);
    try std.testing.expectEqualStrings("C", result_loop.val(1, type_col).?);
    try std.testing.expectEqualStrings("3", result_loop.val(1, id_col).?);
    try std.testing.expectEqualStrings("C", result_loop.val(2, type_col).?);
    try std.testing.expectEqualStrings("5", result_loop.val(2, id_col).?);
}

test "stripDocumentHydrogens is no-op when loop has no H atoms" {
    const allocator = std.testing.allocator;

    var block = zreduce.cif.Block{ .name = "TEST", .items = .empty };
    defer block.deinit(allocator);

    var loop = zreduce.cif.Loop{};
    try loop.tags.append(allocator, "_atom_site.type_symbol");
    try loop.tags.append(allocator, "_atom_site.id");

    try loop.values.append(allocator, "N");
    try loop.values.append(allocator, "1");
    try loop.values.append(allocator, "C");
    try loop.values.append(allocator, "2");

    try block.items.append(allocator, .{ .loop = loop });

    stripDocumentHydrogens(&block);

    const result_loop = block.findLoop("_atom_site.type_symbol").?;
    try std.testing.expectEqual(@as(usize, 2), result_loop.length());
}

test "processFile processes tiny.cif via mmcif path" {
    const allocator = std.testing.allocator;

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    // Write the tiny.cif fixture to a temp file
    const cif_content = @embedFile("test_data/tiny.cif");
    const io = std.testing.io;
    try tmp_dir.dir.writeFile(io, .{ .sub_path = "tiny.cif", .data = cif_content });

    const input_path = try tmp_dir.dir.realPathFileAlloc(io, "tiny.cif", allocator);
    defer allocator.free(input_path);
    const output_path = try tmp_dir.dir.realPathFileAlloc(io, ".", allocator);
    defer allocator.free(output_path);
    const out_cif = try std.fmt.allocPrint(allocator, "{s}/out.cif", .{output_path});
    defer allocator.free(out_cif);

    const result = try processFile(allocator, .{
        .input_path = input_path,
        .output_path = out_cif,
        .format = .mmcif,
        .no_opt = true,
        .quiet = true,
    });

    // ALA has heavy atoms → some H should be placed
    try std.testing.expect(result.n_placed > 0);
    try std.testing.expectEqual(@as(u32, 1), result.n_residues);

    // Output file must exist
    const out_file = try tmp_dir.dir.openFile(io, "out.cif", .{});
    out_file.close(io);
}

test "processFile processes tiny.pdb via pdb path" {
    const allocator = std.testing.allocator;

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    const pdb_content = @embedFile("test_data/tiny.pdb");
    const io = std.testing.io;
    try tmp_dir.dir.writeFile(io, .{ .sub_path = "tiny.pdb", .data = pdb_content });

    const input_path = try tmp_dir.dir.realPathFileAlloc(io, "tiny.pdb", allocator);
    defer allocator.free(input_path);
    const output_path = try tmp_dir.dir.realPathFileAlloc(io, ".", allocator);
    defer allocator.free(output_path);
    const out_pdb = try std.fmt.allocPrint(allocator, "{s}/out.pdb", .{output_path});
    defer allocator.free(out_pdb);

    const result = try processFile(allocator, .{
        .input_path = input_path,
        .output_path = out_pdb,
        .format = .pdb,
        .no_opt = true,
        .quiet = true,
    });

    try std.testing.expect(result.n_placed > 0);
    try std.testing.expectEqual(@as(u32, 1), result.n_residues);

    const out_file = try tmp_dir.dir.openFile(io, "out.pdb", .{});
    out_file.close(io);
}

test "processFile strip_h removes existing hydrogens before placement" {
    const allocator = std.testing.allocator;

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    // ala_with_h.cif already contains a pre-placed HA — strip_h should remove it
    const cif_content = @embedFile("test_data/ala_with_h.cif");
    const io = std.testing.io;
    try tmp_dir.dir.writeFile(io, .{ .sub_path = "ala_with_h.cif", .data = cif_content });

    const input_path = try tmp_dir.dir.realPathFileAlloc(io, "ala_with_h.cif", allocator);
    defer allocator.free(input_path);
    const out_dir = try tmp_dir.dir.realPathFileAlloc(io, ".", allocator);
    defer allocator.free(out_dir);
    const out_path = try std.fmt.allocPrint(allocator, "{s}/out.cif", .{out_dir});
    defer allocator.free(out_path);

    // Without strip_h: the existing H should be counted as skipped_existing
    const result_keep = try processFile(allocator, .{
        .input_path = input_path,
        .output_path = out_path,
        .format = .mmcif,
        .no_opt = true,
        .quiet = true,
    });
    try std.testing.expect(result_keep.n_skipped_existing > 0);

    // With strip_h: the existing H is removed first, so n_skipped_existing drops to 0
    const result_strip = try processFile(allocator, .{
        .input_path = input_path,
        .output_path = out_path,
        .format = .mmcif,
        .no_opt = true,
        .quiet = true,
        .strip_h = true,
    });
    try std.testing.expectEqual(@as(u32, 0), result_strip.n_skipped_existing);
    // And fresh H atoms should be placed
    try std.testing.expect(result_strip.n_placed > 0);
}
