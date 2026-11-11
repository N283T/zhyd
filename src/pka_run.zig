//! Per-file pKa prediction pipeline shared by `zhyd pka <file>` and
//! `zhyd pka <dir>`: parse, annotate chemistry, predict, write the JSON
//! report. Errors are returned with the stage that failed; the caller decides
//! whether to exit (single file) or record the failure and continue (batch).

const std = @import("std");
const Allocator = std.mem.Allocator;
const zhyd = @import("root.zig");
const timing = zhyd.timing;

fn defaultIo() std.Io {
    return std.Io.Threaded.global_single_threaded.io();
}

pub const Engine = enum {
    builtin,
    propka,
    empirical,

    /// PROPKA when compiled in (zig build -Dpropka=true), else the MIT
    /// empirical engine.
    pub fn default() Engine {
        return if (zhyd.pka_propka.available) .propka else .empirical;
    }
};

/// Fully resolved prediction options. Borrowed data (parameter overrides)
/// must outlive every call that uses these options; it is only read, so one
/// `Options` may be shared by all batch workers.
pub const Options = struct {
    engine: Engine,
    model_filter: zhyd.run.ModelFilter = .all,
    prediction_config: zhyd.pka.Config = .{},
    propka_options: zhyd.pka_propka.Options = .{},
    empirical_options: zhyd.pka_empirical.Options = .{},
};

/// Stages reported by `zhyd pka --timing`, in pipeline order.
pub const stages = [_]timing.Stage{ .load, .read, .parse, .chemistry, .pka, .write };

/// Where a file failed.
pub const FailStage = enum {
    read,
    parse,
    no_models,
    predict,
    write,

    pub fn name(self: FailStage) []const u8 {
        return switch (self) {
            .read => "read",
            .parse => "parse",
            .no_models => "parse",
            .predict => "predict",
            .write => "write",
        };
    }
};

pub const Failure = struct {
    stage: FailStage = .read,
    err: anyerror = error.Unknown,
};

pub const FileSummary = struct {
    /// Predicted sites summed over the reported models.
    n_sites: u32 = 0,
};

/// Read the input structure (gzip-aware). The caller frees the result.
pub fn readInput(allocator: Allocator, path: []const u8) ![]const u8 {
    return zhyd.run.readFile(allocator, path);
}

/// Predict and write the report for an already-read input. `input_path` is
/// echoed verbatim into the report's "input" field. `output_path` null
/// writes to stdout (single-file mode only; not thread-safe).
pub fn processSource(
    allocator: Allocator,
    source: []const u8,
    input_path: []const u8,
    output_path: ?[]const u8,
    options: *const Options,
    failure: *Failure,
    times: *timing.StageTimes,
) anyerror!FileSummary {
    var sw = timing.Stopwatch.start();
    if (zhyd.run.detectFormat(input_path) == .mmcif) {
        var entries = zhyd.mmcif.parseModels(allocator, source, options.model_filter) catch |err| {
            failure.* = .{ .stage = .parse, .err = err };
            return err;
        };
        defer {
            for (entries.items) |*entry| entry.model.deinit();
            entries.deinit(allocator);
        }
        sw.lapInto(times, .parse);
        return predictEntries(allocator, entries.items, input_path, output_path, options, failure, times, &sw);
    } else {
        var parsed = zhyd.pdb.parseAll(allocator, source, options.model_filter) catch |err| {
            failure.* = .{ .stage = .parse, .err = err };
            return err;
        };
        defer parsed.deinit(allocator);
        sw.lapInto(times, .parse);
        return predictEntries(allocator, parsed.entries.items, input_path, output_path, options, failure, times, &sw);
    }
}

/// Read, predict and write one file (batch worker entry point).
pub fn processFile(
    allocator: Allocator,
    input_path: []const u8,
    output_path: []const u8,
    options: *const Options,
    failure: *Failure,
    times: *timing.StageTimes,
) anyerror!FileSummary {
    var sw = timing.Stopwatch.start();
    const source = readInput(allocator, input_path) catch |err| {
        failure.* = .{ .stage = .read, .err = err };
        return err;
    };
    defer allocator.free(source);
    sw.lapInto(times, .read);
    return processSource(allocator, source, input_path, output_path, options, failure, times);
}

fn predictEntries(
    allocator: Allocator,
    entries: anytype,
    input_path: []const u8,
    output_path: ?[]const u8,
    options: *const Options,
    failure: *Failure,
    times: *timing.StageTimes,
    sw: *timing.Stopwatch,
) anyerror!FileSummary {
    if (entries.len == 0) {
        failure.* = .{ .stage = .no_models, .err = error.NoMatchingModels };
        return error.NoMatchingModels;
    }
    const predictions = allocator.alloc(zhyd.pka.Prediction, entries.len) catch |err| {
        failure.* = .{ .stage = .predict, .err = err };
        return err;
    };
    defer allocator.free(predictions);
    var predicted: usize = 0;
    defer for (predictions[0..predicted]) |*prediction| prediction.deinit();
    const reports = allocator.alloc(zhyd.writer.protonation_json.ReportModel, entries.len) catch |err| {
        failure.* = .{ .stage = .predict, .err = err };
        return err;
    };
    defer allocator.free(reports);
    var summary = FileSummary{};
    for (entries, 0..) |*entry, i| {
        zhyd.place.applyChemistry(&entry.model);
        sw.lapInto(times, .chemistry);
        predictions[i] = runEngine(allocator, options, &entry.model) catch |err| {
            failure.* = .{ .stage = .predict, .err = err };
            return err;
        };
        sw.lapInto(times, .pka);
        predicted += 1;
        summary.n_sites += @intCast(predictions[i].sites.len);
        reports[i] = .{ .model_num = entry.model_num, .model = &entry.model, .prediction = &predictions[i] };
    }
    writeReport(output_path, input_path, reports) catch |err| {
        failure.* = .{ .stage = .write, .err = err };
        return err;
    };
    sw.lapInto(times, .write);
    return summary;
}

fn runEngine(allocator: Allocator, options: *const Options, mdl: *const zhyd.model.Model) !zhyd.pka.Prediction {
    return switch (options.engine) {
        .builtin => zhyd.pka.predict(allocator, mdl, options.prediction_config),
        .propka => zhyd.pka_propka.predict(allocator, mdl, options.propka_options),
        .empirical => zhyd.pka_empirical.predict(allocator, mdl, options.empirical_options),
    };
}

fn writeReport(
    output_path: ?[]const u8,
    input_path: []const u8,
    reports: []const zhyd.writer.protonation_json.ReportModel,
) !void {
    const io = defaultIo();
    var buffer: [4096]u8 = undefined;
    if (output_path) |path| {
        const file = try std.Io.Dir.cwd().createFile(io, path, .{});
        defer file.close(io);
        var fw = file.writer(io, &buffer);
        try zhyd.writer.protonation_json.writePredictionReport(&fw.interface, input_path, reports);
        try fw.interface.flush();
    } else {
        const stdout = std.Io.File.stdout();
        var sw = stdout.writer(io, &buffer);
        try zhyd.writer.protonation_json.writePredictionReport(&sw.interface, input_path, reports);
        try sw.interface.flush();
    }
}

const testing = std.testing;

test "processFile reports parse failures with their stage" {
    var failure = Failure{};
    var times = timing.StageTimes{};
    const options = Options{ .engine = .builtin };
    try testing.expectError(error.FileNotFound, processFile(testing.allocator, "/nonexistent/x.cif", "/dev/null", &options, &failure, &times));
    try testing.expectEqual(FailStage.read, failure.stage);
}

test "processFile predicts sites for his.cif" {
    var failure = Failure{};
    var times = timing.StageTimes{};
    const options = Options{ .engine = .builtin };
    const summary = try processFile(testing.allocator, "src/test_data/his.cif", "/dev/null", &options, &failure, &times);
    try testing.expectEqual(@as(u32, 3), summary.n_sites); // HIS plus both termini
}
