//! `zhyd pka <dir>`: pKa reports for every structure file in a directory.
//!
//! Each `<dir>/<name>` gets `<outdir>/<stem>.json`, byte-identical to
//! `zhyd pka <dir>/<name> -o <outdir>/<stem>.json` with the same options
//! (the report echoes the input path as `<dir>/<name>`). Files are processed
//! in parallel; one failing file does not stop the batch.

const std = @import("std");
const Allocator = std.mem.Allocator;
const zhyd = @import("root.zig");
const batch = zhyd.batch;
const pka_run = zhyd.pka_run;
const thread_pool = @import("thread_pool.zig");
const timing = zhyd.timing;

fn defaultIo() std.Io {
    return std.Io.Threaded.global_single_threaded.io();
}

pub const Config = struct {
    input_dir: []const u8,
    /// Default: `<input_dir>_pka/`.
    output_dir: ?[]const u8 = null,
    jsonl_path: ?[]const u8 = null,
    /// 0 = auto-detect CPU count.
    n_threads: u32 = 0,
    quiet: bool = false,
    show_timing: bool = false,
    options: pka_run.Options,
    /// Time spent loading parameter files before `run` (reported as "load").
    load_ns: u64 = 0,
};

pub const FileResult = struct {
    input_path: []const u8,
    output_path: []const u8,
    status: Status = .err,
    n_sites: u32 = 0,
    time_ns: u64 = 0,
    /// Failing stage and error; meaningful when status == .err.
    failure: pka_run.Failure = .{},
    /// Set when the output name collides with an earlier file's.
    collides_with: ?[]const u8 = null,
    timing: timing.StageTimes = .{},

    pub const Status = enum { ok, err };
};

/// `<stem>.json` for a structure file name: drop `.gz`, then the
/// `.cif`/`.pdb`/`.ent` extension.
pub fn outputName(allocator: Allocator, filename: []const u8) ![]u8 {
    var base = std.fs.path.basename(filename);
    if (std.mem.endsWith(u8, base, ".gz")) base = base[0 .. base.len - 3];
    for ([_][]const u8{ ".cif", ".pdb", ".ent" }) |ext| {
        if (std.mem.endsWith(u8, base, ext)) {
            base = base[0 .. base.len - ext.len];
            break;
        }
    }
    return std.fmt.allocPrint(allocator, "{s}.json", .{base});
}

fn writeJsonlLine(writer: anytype, r: FileResult) !void {
    try writer.writeAll("{\"input\":");
    try batch.writeJsonString(writer, r.input_path);
    try writer.writeAll(",\"output\":");
    try batch.writeJsonString(writer, r.output_path);
    const seconds = @as(f64, @floatFromInt(r.time_ns)) / 1e9;
    switch (r.status) {
        .ok => try writer.print(",\"status\":\"ok\",\"n_sites\":{d},\"seconds\":{d:.4}}}\n", .{ r.n_sites, seconds }),
        .err => {
            try writer.print(",\"status\":\"error\",\"seconds\":{d:.4},\"error\":", .{seconds});
            var msg_buf: [256]u8 = undefined;
            const msg = errorMessage(&msg_buf, r);
            try batch.writeJsonString(writer, msg);
            try writer.writeAll("}\n");
        },
    }
}

fn errorMessage(buf: []u8, r: FileResult) []const u8 {
    if (r.collides_with) |other| {
        return std.fmt.bufPrint(buf, "output name collides with '{s}'", .{other}) catch "output name collision";
    }
    return std.fmt.bufPrint(buf, "{s}: {s}", .{ r.failure.stage.name(), @errorName(r.failure.err) }) catch @errorName(r.failure.err);
}

const Context = struct {
    results: []FileResult,
    options: *const pka_run.Options,
    processed_count: std.atomic.Value(usize) = .init(0),
    jsonl_stream: ?*batch.JsonlStreamWriter,
};

fn processOne(ctx: *Context, idx: usize) void {
    const r = &ctx.results[idx];
    if (r.collides_with == null) {
        // Per-file arena: engines allocate freely; nothing outlives the file.
        var arena = std.heap.ArenaAllocator.init(std.heap.smp_allocator);
        defer arena.deinit();
        var sw = timing.Stopwatch.start();
        if (pka_run.processFile(arena.allocator(), r.input_path, r.output_path, ctx.options, &r.failure, &r.timing)) |summary| {
            r.status = .ok;
            r.n_sites = summary.n_sites;
        } else |_| {
            r.status = .err;
        }
        r.time_ns = sw.read();
    }
    if (ctx.jsonl_stream) |stream| {
        var line_buf: [8192]u8 = undefined;
        var w: std.Io.Writer = .fixed(&line_buf);
        if (writeJsonlLine(&w, r.*)) {
            stream.writeLine(w.buffered());
        } else |_| {
            _ = stream.write_errors.fetchAdd(1, .monotonic);
        }
    }
    _ = ctx.processed_count.fetchAdd(1, .monotonic);
}

/// Thread-pool work function: one file per chunk.
fn worker(ctx: *Context, start: usize, end: usize) void {
    for (start..end) |idx| processOne(ctx, idx);
}

const ProgressArgs = struct {
    ctx: *Context,
    done: *std.atomic.Value(bool),
};

fn progressLoop(pa: *ProgressArgs) void {
    const total = pa.ctx.results.len;
    while (!pa.done.load(.acquire)) {
        const processed = pa.ctx.processed_count.load(.acquire);
        std.debug.print("\rProcessing: {d}/{d}", .{ processed, total });
        std.Io.sleep(defaultIo(), .fromNanoseconds(100 * std.time.ns_per_ms), .awake) catch {};
    }
    std.debug.print("\rProcessing: {d}/{d}\n", .{ total, total });
}

pub fn run(allocator: Allocator, config: Config) !void {
    const io = defaultIo();
    var wall = timing.Stopwatch.start();

    const files = try batch.scanDirectory(allocator, config.input_dir);
    defer {
        for (files) |f| allocator.free(f);
        allocator.free(files);
    }
    if (files.len == 0) {
        std.debug.print("Error: no structure files (.cif, .pdb, .ent, or .gz variants) found in '{s}'\n", .{config.input_dir});
        return error.NoFilesFound;
    }
    std.debug.print("Found {d} structure files in '{s}'\n", .{ files.len, config.input_dir });

    const output_dir: []const u8 = config.output_dir orelse try std.fmt.allocPrint(allocator, "{s}_pka", .{config.input_dir});
    defer if (config.output_dir == null) allocator.free(output_dir);
    std.Io.Dir.cwd().createDirPath(io, output_dir) catch |err| {
        std.debug.print("Error: cannot create output directory '{s}': {s}\n", .{ output_dir, @errorName(err) });
        return err;
    };

    // Paths are built up front so workers never allocate shared memory.
    var path_arena = std.heap.ArenaAllocator.init(allocator);
    defer path_arena.deinit();
    const pa = path_arena.allocator();
    const results = try allocator.alloc(FileResult, files.len);
    defer allocator.free(results);
    var seen: std.StringHashMapUnmanaged([]const u8) = .empty;
    defer seen.deinit(allocator);
    for (files, results) |name, *r| {
        const out_name = try outputName(pa, name);
        r.* = .{
            .input_path = try std.fs.path.join(pa, &.{ config.input_dir, name }),
            .output_path = try std.fs.path.join(pa, &.{ output_dir, out_name }),
        };
        // Files are sorted, so the first of `x.cif` / `x.cif.gz` wins.
        const gop = try seen.getOrPut(allocator, out_name);
        if (gop.found_existing) r.collides_with = gop.value_ptr.* else gop.value_ptr.* = name;
    }

    var jsonl_file: ?std.Io.File = null;
    if (config.jsonl_path) |path| jsonl_file = try std.Io.Dir.cwd().createFile(io, path, .{});
    defer if (jsonl_file) |f| f.close(io);
    var jsonl_storage: batch.JsonlStreamWriter = undefined;
    if (jsonl_file) |f| jsonl_storage = .{ .file = f };

    var ctx = Context{
        .results = results,
        .options = &config.options,
        .jsonl_stream = if (jsonl_file != null) &jsonl_storage else null,
    };

    var progress_done = std.atomic.Value(bool).init(false);
    var progress_args = ProgressArgs{ .ctx = &ctx, .done = &progress_done };
    const progress_thread: ?std.Thread = if (!config.quiet)
        std.Thread.spawn(.{}, progressLoop, .{&progress_args}) catch null
    else
        null;

    const n_threads = batch.resolveThreadCount(config.n_threads, files.len);
    if (n_threads <= 1) {
        worker(&ctx, 0, files.len);
    } else {
        var pool = try thread_pool.ThreadPool(*Context, void).init(allocator, n_threads, worker, &ctx, files.len, 1);
        defer pool.deinit();
        pool.run();
    }

    progress_done.store(true, .release);
    if (progress_thread) |t| t.join();

    var successful: u32 = 0;
    var failed: u32 = 0;
    var n_sites: u64 = 0;
    var file_time_ns: u64 = 0;
    var stages = timing.StageTimes{};
    for (results) |r| {
        file_time_ns += r.time_ns;
        switch (r.status) {
            .ok => {
                successful += 1;
                n_sites += r.n_sites;
                stages.merge(r.timing);
            },
            .err => failed += 1,
        }
    }

    if (config.jsonl_path) |path| {
        std.debug.print("JSONL log written to '{s}'", .{path});
        const errs = jsonl_storage.write_errors.load(.acquire);
        if (errs > 0) std.debug.print(" ({d} write error(s))", .{errs});
        std.debug.print("\n", .{});
    }

    std.debug.print("\n--- pKa Batch Summary ---\n", .{});
    std.debug.print("  Total files : {d}\n", .{files.len});
    std.debug.print("  Successful  : {d}\n", .{successful});
    std.debug.print("  Failed      : {d}\n", .{failed});
    std.debug.print("  Sites       : {d}\n", .{n_sites});
    std.debug.print("  Total time  : {d:.1} ms\n", .{@as(f64, @floatFromInt(file_time_ns)) / 1e6});
    std.debug.print("  Output dir  : {s}\n", .{output_dir});
    if (failed > 0) {
        std.debug.print("\nFailed files:\n", .{});
        const max_listed = 20;
        var listed: usize = 0;
        for (results) |r| {
            if (r.status != .err) continue;
            if (listed == max_listed) {
                std.debug.print("  ... and {d} more\n", .{failed - max_listed});
                break;
            }
            var buf: [256]u8 = undefined;
            std.debug.print("  {s}: {s}\n", .{ r.input_path, errorMessage(&buf, r) });
            listed += 1;
        }
    }

    if (config.show_timing) {
        stages.addNs(.load, config.load_ns);
        timing.printStageProfile("PKA_BATCH", &pka_run.stages, stages, wall.read() + config.load_ns, .{
            .files = files.len,
            .successful = successful,
            .threads = n_threads,
        });
    }

    if (failed > 0) return error.SomeFilesFailed;
}

const testing = std.testing;

test "outputName strips gzip and structure extensions" {
    const cases = [_][2][]const u8{
        .{ "1abc.cif", "1abc.json" },
        .{ "1abc.cif.gz", "1abc.json" },
        .{ "x.pdb.gz", "x.json" },
        .{ "pdb1abc.ent", "pdb1abc.json" },
        .{ "dir/AF-P1-F1-model_v6.cif", "AF-P1-F1-model_v6.json" },
    };
    for (cases) |c| {
        const got = try outputName(testing.allocator, c[0]);
        defer testing.allocator.free(got);
        try testing.expectEqualStrings(c[1], got);
    }
}

test "writeJsonlLine reports ok and error results" {
    var buf: std.Io.Writer.Allocating = .init(testing.allocator);
    defer buf.deinit();
    try writeJsonlLine(&buf.writer, .{ .input_path = "d/a.cif", .output_path = "o/a.json", .status = .ok, .n_sites = 7, .time_ns = 1_500_000 });
    try writeJsonlLine(&buf.writer, .{ .input_path = "d/b.cif", .output_path = "o/b.json", .failure = .{ .stage = .parse, .err = error.CifParseError } });
    try writeJsonlLine(&buf.writer, .{ .input_path = "d/a.cif.gz", .output_path = "o/a.json", .collides_with = "a.cif" });
    try testing.expectEqualStrings(
        \\{"input":"d/a.cif","output":"o/a.json","status":"ok","n_sites":7,"seconds":0.0015}
        \\{"input":"d/b.cif","output":"o/b.json","status":"error","seconds":0.0000,"error":"parse: CifParseError"}
        \\{"input":"d/a.cif.gz","output":"o/a.json","status":"error","seconds":0.0000,"error":"output name collides with 'a.cif'"}
        \\
    , buf.written());
}
