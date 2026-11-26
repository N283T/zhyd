//! Generic fixed-capacity string type for space-padded PDB-style identifiers.
//!
//! FixedString(N) bundles a [N]u8 buffer with a compact length field and
//! provides set/slice helpers.  It replaces the repeated
//!   buf: [N]u8, buf_len: uX
//! pattern found across Atom, Residue, and Chain.

const std = @import("std");

/// A fixed-capacity string that stores up to N bytes plus a compact length.
/// The buffer is zero-initialised to spaces (' ') so it is compatible with
/// PDB space-padded atom/chain names.
pub fn FixedString(comptime N: comptime_int) type {
    comptime std.debug.assert(N > 0 and N <= 128);
    // Minimum unsigned int wide enough to hold values 0..N.
    const LenInt = @Int(.unsigned, std.math.log2_int_ceil(u8, N + 1));

    return struct {
        buf: [N]u8 = @splat(' '),
        len: LenInt = 0,

        const Self = @This();

        /// Return the active slice (first `len` bytes of the buffer).
        pub fn slice(self: *const Self) []const u8 {
            return self.buf[0..self.len];
        }

        /// Copy `value` into the buffer, truncating to N characters.
        pub fn set(self: *Self, value: []const u8) void {
            const n: LenInt = @intCast(@min(value.len, N));
            self.buf = @splat(' ');
            for (0..n) |i| self.buf[i] = value[i];
            self.len = n;
        }
    };
}

// ── Tests ─────────────────────────────────────────────────────────────────────

test "FixedString(4) set and slice" {
    var s: FixedString(4) = .{};
    s.set("CA");
    try std.testing.expectEqualStrings("CA", s.slice());
    try std.testing.expectEqual(@as(u3, 2), s.len);
}

test "FixedString(4) truncates at capacity" {
    var s: FixedString(4) = .{};
    s.set("ABCDE");
    try std.testing.expectEqualStrings("ABCD", s.slice());
    try std.testing.expectEqual(@as(u3, 4), s.len);
}

test "FixedString(5) set and slice" {
    var s: FixedString(5) = .{};
    s.set("BGLAN");
    try std.testing.expectEqualStrings("BGLAN", s.slice());
    try std.testing.expectEqual(@as(u3, 5), s.len);
}

test "FixedString default value is spaces" {
    const s: FixedString(4) = .{};
    try std.testing.expectEqualStrings("", s.slice());
    try std.testing.expectEqualSlices(u8, "    ", &s.buf);
}
