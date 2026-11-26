//! Character classification lookup table for CIF format parsing.

pub const CharType = enum(u2) {
    whitespace = 0,
    ordinary = 1,
    special = 2,
};

pub const char_table: [128]CharType = blk: {
    var table: [128]CharType = @splat(.ordinary);
    // Control characters (0-32) → whitespace
    for (0..33) |i| table[i] = .whitespace;
    table[127] = .whitespace; // DEL
    // Explicit whitespace (redundant but clear)
    table[' '] = .whitespace;
    table['\t'] = .whitespace;
    table['\n'] = .whitespace;
    table['\r'] = .whitespace;
    // Special characters
    table['#'] = .special;
    table['$'] = .special;
    table['\''] = .special;
    table['"'] = .special;
    table['_'] = .special;
    table[';'] = .special;
    break :blk table;
};

pub fn isWhitespace(c: u8) bool {
    if (c > 127) return false;
    return char_table[c] == .whitespace;
}

pub fn isOrdinary(c: u8) bool {
    if (c > 127) return true; // Non-ASCII treated as ordinary
    return char_table[c] == .ordinary;
}

// --- Tests ---

const std = @import("std");

test "whitespace detection" {
    try std.testing.expect(isWhitespace(' '));
    try std.testing.expect(isWhitespace('\t'));
    try std.testing.expect(isWhitespace('\n'));
    try std.testing.expect(isWhitespace('\r'));
    try std.testing.expect(isWhitespace(0));
    try std.testing.expect(isWhitespace(127));
    try std.testing.expect(!isWhitespace('a'));
    try std.testing.expect(!isWhitespace('_'));
}

test "ordinary detection" {
    try std.testing.expect(isOrdinary('a'));
    try std.testing.expect(isOrdinary('Z'));
    try std.testing.expect(isOrdinary('5'));
    try std.testing.expect(!isOrdinary('_'));
    try std.testing.expect(!isOrdinary('#'));
    try std.testing.expect(!isOrdinary(' '));
    // Non-ASCII treated as ordinary
    try std.testing.expect(isOrdinary(200));
}
