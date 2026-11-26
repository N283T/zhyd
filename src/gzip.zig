//! Gzip I/O via Zig std.compress.flate.

const std = @import("std");
const Allocator = std.mem.Allocator;

fn defaultIo() std.Io {
    return std.Io.Threaded.global_single_threaded.io();
}

/// Default max decompressed size: 1 GB (matches plain file read limit in run.zig)
pub const DEFAULT_MAX_SIZE: usize = 1024 * 1024 * 1024;

/// Decompress a gzip file. Caller owns the returned slice.
pub fn readGzip(allocator: Allocator, path: []const u8) ![]u8 {
    return readGzipLimit(allocator, path, DEFAULT_MAX_SIZE);
}

/// Decompress a gzip file with a custom size limit. Caller owns the returned slice.
/// Non-gzip files are returned as raw bytes for compatibility with historical transparent reads.
pub fn readGzipLimit(allocator: Allocator, path: []const u8, max_size: usize) ![]u8 {
    const io = defaultIo();
    const data = try std.Io.Dir.cwd().readFileAlloc(io, path, allocator, .limited(max_size));
    errdefer allocator.free(data);

    if (data.len < 2 or data[0] != 0x1f or data[1] != 0x8b) {
        return data;
    }

    var source_reader: std.Io.Reader = .fixed(data);
    var window: [std.compress.flate.max_window_len]u8 = undefined;
    var decompress: std.compress.flate.Decompress = .init(&source_reader, .gzip, &window);
    var out_writer: std.Io.Writer.Allocating = .init(allocator);
    errdefer out_writer.deinit();
    _ = try decompress.reader.streamRemaining(&out_writer.writer);
    if (out_writer.writer.buffered().len > max_size) return error.StreamTooLong;
    const out = try out_writer.toOwnedSlice();
    allocator.free(data);
    return out;
}

/// A writer that buffers uncompressed output and writes gzip-compressed data on close.
pub const GzipWriter = struct {
    allocator: Allocator,
    path: []const u8,
    buffer: std.Io.Writer.Allocating,
    closed: bool = false,

    pub fn init(allocator: Allocator, path: []const u8) !GzipWriter {
        return .{
            .allocator = allocator,
            .path = path,
            .buffer = .init(allocator),
        };
    }

    pub fn writer(self: *GzipWriter) *std.Io.Writer {
        return &self.buffer.writer;
    }

    pub fn close(self: *GzipWriter) !void {
        if (self.closed) return;
        self.closed = true;

        const data = try self.buffer.toOwnedSlice();
        defer self.allocator.free(data);

        const io = defaultIo();
        var file = try std.Io.Dir.cwd().createFile(io, self.path, .{});
        defer file.close(io);

        var file_buf: [4096]u8 = undefined;
        var file_writer = file.writer(io, &file_buf);
        const flate_buf = try self.allocator.alloc(u8, std.compress.flate.max_window_len);
        defer self.allocator.free(flate_buf);

        var compressor = try std.compress.flate.Compress.init(
            &file_writer.interface,
            flate_buf,
            .gzip,
            .default,
        );
        try compressor.writer.writeAll(data);
        try compressor.finish();
        try file_writer.interface.flush();
    }

    pub fn deinit(self: *GzipWriter) void {
        if (!self.closed) self.buffer.deinit();
        self.* = undefined;
    }
};

test "readGzip decompresses gzip store block" {
    const allocator = std.testing.allocator;
    const io = std.testing.io;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    const gz_data = [_]u8{
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03,
        0x01, 0x0c, 0x00, 0xf3, 0xff, 'H',  'e',  'l',  'l',  'o',
        ' ',  'w',  'o',  'r',  'l',  'd',  '\n', 0xd5, 0xe0, 0x39,
        0xb7, 0x0c, 0x00, 0x00, 0x00,
    };
    try tmp_dir.dir.writeFile(io, .{ .sub_path = "test.gz", .data = &gz_data });

    const tmp_path = try tmp_dir.dir.realPathFileAlloc(io, "test.gz", allocator);
    defer allocator.free(tmp_path);

    const out = try readGzip(allocator, tmp_path);
    defer allocator.free(out);
    try std.testing.expectEqualStrings("Hello world\n", out);
}

test "readGzipLimit rejects decompressed data above limit" {
    const allocator = std.testing.allocator;
    const io = std.testing.io;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    const gz_data = [_]u8{
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03,
        0x01, 0x0c, 0x00, 0xf3, 0xff, 'H',  'e',  'l',  'l',  'o',
        ' ',  'w',  'o',  'r',  'l',  'd',  '\n', 0xd5, 0xe0, 0x39,
        0xb7, 0x0c, 0x00, 0x00, 0x00,
    };
    try tmp_dir.dir.writeFile(io, .{ .sub_path = "test.gz", .data = &gz_data });

    const tmp_path = try tmp_dir.dir.realPathFileAlloc(io, "test.gz", allocator);
    defer allocator.free(tmp_path);

    const result = readGzipLimit(allocator, tmp_path, 5);
    try std.testing.expectError(error.StreamTooLong, result);
}

test "readGzip reads non-gzip file transparently" {
    const allocator = std.testing.allocator;
    const io = std.testing.io;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    try tmp_dir.dir.writeFile(io, .{ .sub_path = "plain.gz", .data = "raw content" });
    const tmp_path = try tmp_dir.dir.realPathFileAlloc(io, "plain.gz", allocator);
    defer allocator.free(tmp_path);

    const out = try readGzip(allocator, tmp_path);
    defer allocator.free(out);
    try std.testing.expectEqualStrings("raw content", out);
}

test "readGzip empty file returns empty" {
    const allocator = std.testing.allocator;
    const io = std.testing.io;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    try tmp_dir.dir.writeFile(io, .{ .sub_path = "empty.gz", .data = "" });
    const tmp_path = try tmp_dir.dir.realPathFileAlloc(io, "empty.gz", allocator);
    defer allocator.free(tmp_path);

    const out = try readGzip(allocator, tmp_path);
    defer allocator.free(out);
    try std.testing.expectEqualStrings("", out);
}

test "GzipWriter writes gzip output" {
    const allocator = std.testing.allocator;
    const io = std.testing.io;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    const tmp_path = try tmp_dir.dir.realPathFileAlloc(io, ".", allocator);
    defer allocator.free(tmp_path);
    const out_path = try allocator.print("{s}/out.gz", .{tmp_path});
    defer allocator.free(out_path);

    var gw = try GzipWriter.init(allocator, out_path);
    try gw.writer().writeAll("Hello, gzip writer!\nLine two.\n");
    try gw.close();
    gw.deinit();

    const out = try readGzip(allocator, out_path);
    defer allocator.free(out);
    try std.testing.expectEqualStrings("Hello, gzip writer!\nLine two.\n", out);
}
