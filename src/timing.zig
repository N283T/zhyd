//! Per-stage wall-clock accounting for `--timing`.
//!
//! Stage times are always collected (a clock read per stage is negligible);
//! they are only printed when the user asks for `--timing`.
//! Report layout adapted from zsasa (MIT, same author).

const std = @import("std");

fn defaultIo() std.Io {
    return std.Io.Threaded.global_single_threaded.io();
}

pub const Stage = enum {
    /// Loading shared inputs once: CCD/SDF dictionaries, parameter files.
    load,
    /// Reading the input file, including gzip decompression.
    read,
    /// CIF/PDB parsing, model extraction and bond records.
    parse,
    /// Chemistry annotations (donor/acceptor/charge flags).
    chemistry,
    /// Hydrogen placement.
    place,
    /// Mover generation and optimization.
    optimize,
    /// Post-placement model validation.
    validate,
    /// pKa prediction.
    pka,
    /// Writing outputs (coordinates, JSON reports and logs).
    write,

    pub fn label(self: Stage) []const u8 {
        return switch (self) {
            .load => "Load dictionaries/parameters",
            .read => "Read + decompress",
            .parse => "Parse",
            .chemistry => "Chemistry",
            .place => "Place H",
            .optimize => "Optimize",
            .validate => "Validate",
            .pka => "pKa predict",
            .write => "Write output",
        };
    }

    /// Key for the machine-readable `<PREFIX>_<KEY>_TIME_MS:` lines.
    pub fn key(self: Stage) []const u8 {
        return switch (self) {
            .load => "LOAD",
            .read => "READ",
            .parse => "PARSE",
            .chemistry => "CHEMISTRY",
            .place => "PLACE",
            .optimize => "OPTIMIZE",
            .validate => "VALIDATE",
            .pka => "PKA",
            .write => "WRITE",
        };
    }
};

const n_stages = std.meta.fields(Stage).len;

pub const StageTimes = struct {
    ns: [n_stages]u64 = @splat(0),

    pub fn addNs(self: *StageTimes, stage: Stage, ns: u64) void {
        self.ns[@intFromEnum(stage)] += ns;
    }

    pub fn get(self: StageTimes, stage: Stage) u64 {
        return self.ns[@intFromEnum(stage)];
    }

    pub fn merge(self: *StageTimes, other: StageTimes) void {
        for (&self.ns, other.ns) |*a, b| a.* += b;
    }

    pub fn sum(self: StageTimes) u64 {
        var total: u64 = 0;
        for (self.ns) |v| total += v;
        return total;
    }
};

/// Monotonic stopwatch. `lap` returns the nanoseconds since the previous lap
/// (or since `start`) and restarts the measurement.
pub const Stopwatch = struct {
    last: std.Io.Timestamp,

    pub fn start() Stopwatch {
        return .{ .last = std.Io.Clock.awake.now(defaultIo()) };
    }

    pub fn read(self: Stopwatch) u64 {
        const ns = self.last.untilNow(defaultIo(), .awake).toNanoseconds();
        return if (ns > 0) @intCast(ns) else 0;
    }

    pub fn lap(self: *Stopwatch) u64 {
        const now = std.Io.Clock.awake.now(defaultIo());
        const ns = self.last.durationTo(now).toNanoseconds();
        self.last = now;
        return if (ns > 0) @intCast(ns) else 0;
    }

    /// Add the lap to `stage` of `times`.
    pub fn lapInto(self: *Stopwatch, times: *StageTimes, stage: Stage) void {
        times.addNs(stage, self.lap());
    }
};

pub const BatchInfo = struct {
    /// Files processed (successful and failed).
    files: usize,
    successful: usize,
    threads: usize,
};

const ns_per_ms = 1_000_000.0;

fn ms(ns: u64) f64 {
    return @as(f64, @floatFromInt(ns)) / ns_per_ms;
}

/// Print the stage breakdown to stderr: a human-readable table, then
/// `<prefix>_<STAGE>_TIME_MS:<value>` lines for scripts. `stages` lists the
/// stages that exist for the command, in pipeline order. With `batch`, stage
/// times are sums over files (CPU-time-like across threads) and the wall time
/// is followed by the file count and throughput.
pub fn printStageProfile(
    prefix: []const u8,
    stages: []const Stage,
    times: StageTimes,
    wall_ns: u64,
    batch: ?BatchInfo,
) void {
    if (batch) |b| {
        std.debug.print("\nTiming breakdown (stage sums over {d} file(s), {d} thread(s)):\n", .{ b.files, b.threads });
    } else {
        std.debug.print("\nTiming breakdown:\n", .{});
    }
    var stage_sum: u64 = 0;
    for (stages) |stage| {
        const v = times.get(stage);
        stage_sum += v;
        std.debug.print("  {s:<30} {d:>10.2} ms\n", .{ stage.label(), ms(v) });
    }
    std.debug.print("  {s:<30} {d:>10.2} ms\n", .{ "Sum of stages", ms(stage_sum) });
    std.debug.print("  {s:<30} {d:>10.2} ms\n", .{ "Wall time", ms(wall_ns) });
    if (batch) |b| {
        const secs = @as(f64, @floatFromInt(wall_ns)) / 1e9;
        const throughput = if (secs > 0) @as(f64, @floatFromInt(b.files)) / secs else 0.0;
        std.debug.print("  {s:<30} {d:>10}\n", .{ "Files", b.files });
        std.debug.print("  {s:<30} {d:>10.1} files/sec\n", .{ "Throughput", throughput });
    }

    for (stages) |stage| {
        std.debug.print("{s}_{s}_TIME_MS:{d:.2}\n", .{ prefix, stage.key(), ms(times.get(stage)) });
    }
    std.debug.print("{s}_WALL_TIME_MS:{d:.2}\n", .{ prefix, ms(wall_ns) });
    if (batch) |b| {
        std.debug.print("{s}_FILES:{d}\n", .{ prefix, b.files });
        std.debug.print("{s}_SUCCESS:{d}\n", .{ prefix, b.successful });
        std.debug.print("{s}_THREADS:{d}\n", .{ prefix, b.threads });
    }
}

test "StageTimes accumulates and merges" {
    var a = StageTimes{};
    a.addNs(.read, 10);
    a.addNs(.read, 5);
    a.addNs(.write, 1);
    var b = StageTimes{};
    b.addNs(.read, 100);
    b.addNs(.pka, 7);
    a.merge(b);
    try std.testing.expectEqual(@as(u64, 115), a.get(.read));
    try std.testing.expectEqual(@as(u64, 7), a.get(.pka));
    try std.testing.expectEqual(@as(u64, 123), a.sum());
}

test "Stopwatch laps are monotonic" {
    var sw = Stopwatch.start();
    var times = StageTimes{};
    sw.lapInto(&times, .parse);
    sw.lapInto(&times, .parse);
    try std.testing.expect(times.get(.parse) < std.time.ns_per_s * 60);
}
