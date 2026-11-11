const std = @import("std");

fn defaultIo() std.Io {
    return std.Io.Threaded.global_single_threaded.io();
}
const model_mod = @import("../model.zig");
const mover_mod = @import("mover.zig");
const fixed_string = @import("../model/fixed_string.zig");

const FixedString = fixed_string.FixedString;

const Model = model_mod.Model;
const Mover = mover_mod.Mover;
const MoverKind = mover_mod.MoverKind;

pub const Selector = struct {
    chain_id: []const u8,
    auth_seq_id: i32,
    ins_code: u8 = ' ',

    pub fn matches(self: Selector, mdl: *const Model, mover: *const Mover) bool {
        const res = mdl.residues.items[mover.residue_idx];
        const chain = mdl.chains.items[res.chain_idx];
        if (res.auth_seq_id != self.auth_seq_id) return false;
        if (res.ins_code != self.ins_code) return false;
        return std.mem.eql(u8, chain.authSlice(), self.chain_id) or
            std.mem.eql(u8, chain.labelSlice(), self.chain_id);
    }
};

pub const Entry = struct {
    selector: Selector,
    comp_id: FixedString(5) = .{},
    target: []const u8,
    value: Value,

    pub fn compIdSlice(self: *const Entry) []const u8 {
        return self.comp_id.slice();
    }
};

pub const Value = union(enum) {
    orientation: u16,
    amide: enum { original, flip },
    his: enum { hie, hid, hie_flip, hid_flip },
};

pub const FixOverrides = struct {
    allocator: std.mem.Allocator,
    entries: []Entry,

    pub fn deinit(self: *FixOverrides) void {
        for (self.entries) |entry| {
            self.allocator.free(entry.selector.chain_id);
            self.allocator.free(entry.target);
        }
        self.allocator.free(self.entries);
        self.entries = &.{};
    }

    pub fn warnUnmatched(self: *const FixOverrides, mdl: *const Model, movers: []const Mover) void {
        for (self.entries) |entry| {
            var matched = false;
            for (movers) |*m| {
                if (!entry.selector.matches(mdl, m)) continue;
                const res = mdl.residues.items[m.residue_idx];
                if (!std.mem.eql(u8, entry.compIdSlice(), res.compIdSlice())) continue;
                const target = moverTargetName(mdl, m);
                if (std.mem.eql(u8, target, entry.target)) {
                    matched = true;
                    break;
                }
            }
            if (!matched) {
                std.log.warn("fix override for {s}:{d} {s} {s} did not match any mover", .{
                    entry.selector.chain_id,
                    entry.selector.auth_seq_id,
                    entry.compIdSlice(),
                    entry.target,
                });
            }
        }
    }
};

pub fn parseFile(allocator: std.mem.Allocator, path: []const u8) !FixOverrides {
    const io = defaultIo();
    const file = try std.Io.Dir.cwd().openFile(io, path, .{});
    defer file.close(io);
    var read_buf: [4096]u8 = undefined;
    var reader = file.reader(io, &read_buf);
    const source = try reader.interface.readAlloc(allocator, 1024 * 1024);
    defer allocator.free(source);
    return parseString(allocator, source);
}

pub fn parseString(allocator: std.mem.Allocator, source: []const u8) !FixOverrides {
    var entries = std.ArrayListUnmanaged(Entry).empty;
    errdefer {
        for (entries.items) |entry| {
            allocator.free(entry.selector.chain_id);
            allocator.free(entry.target);
        }
        entries.deinit(allocator);
    }

    var lines = std.mem.splitScalar(u8, source, '\n');
    var line_num: u32 = 0;
    while (lines.next()) |line_raw| {
        line_num += 1;
        const line = std.mem.trim(u8, line_raw, " \t\r");
        if (line.len == 0 or line[0] == '#') continue;

        var toks = std.mem.tokenizeAny(u8, line, " \t");
        const selector_tok = toks.next() orelse {
            std.log.warn("fix override line {d}: expected 'chain:seq comp target value'", .{line_num});
            return error.InvalidFixOverride;
        };
        const comp_tok = toks.next() orelse {
            std.log.warn("fix override line {d}: missing comp_id", .{line_num});
            return error.InvalidFixOverride;
        };
        const target_tok = toks.next() orelse {
            std.log.warn("fix override line {d}: missing target", .{line_num});
            return error.InvalidFixOverride;
        };
        const value_tok = toks.next() orelse {
            std.log.warn("fix override line {d}: missing value", .{line_num});
            return error.InvalidFixOverride;
        };
        if (toks.next() != null) {
            std.log.warn("fix override line {d}: unexpected extra tokens", .{line_num});
            return error.InvalidFixOverride;
        }
        if (comp_tok.len == 0 or comp_tok.len > 5) {
            std.log.warn("fix override line {d}: invalid comp_id '{s}'", .{ line_num, comp_tok });
            return error.InvalidFixOverride;
        }

        // Validate value before allocating (parseValue is pure)
        const value = parseValue(target_tok, value_tok) catch {
            std.log.warn("fix override line {d}: invalid value '{s}' for target '{s}'", .{ line_num, value_tok, target_tok });
            return error.InvalidFixOverride;
        };
        const selector = parseSelector(allocator, selector_tok) catch {
            std.log.warn("fix override line {d}: invalid selector '{s}' (expected chain:seq[:ins])", .{ line_num, selector_tok });
            return error.InvalidFixOverride;
        };
        errdefer allocator.free(selector.chain_id);
        const target = try allocator.dupe(u8, target_tok);

        var comp_id: FixedString(5) = .{};
        const comp_id_len: u3 = @intCast(comp_tok.len);
        comp_id.len = comp_id_len;
        for (0..comp_id_len) |i| comp_id.buf[i] = std.ascii.toUpper(comp_tok[i]);

        try entries.append(allocator, .{
            .selector = selector,
            .comp_id = comp_id,
            .target = target,
            .value = value,
        });
    }

    return .{ .allocator = allocator, .entries = try entries.toOwnedSlice(allocator) };
}

/// Lock each mover named by a fix entry. Entries that match no mover (a typo,
/// a group that has no mover such as a methyl, or a residue absent from this
/// model) are reported with a warning and otherwise ignored.
pub fn applyFixes(overrides: *const FixOverrides, mdl: *const Model, movers: []Mover) !void {
    const matched = try overrides.allocator.alloc(bool, overrides.entries.len);
    defer overrides.allocator.free(matched);
    @memset(matched, false);

    for (movers) |*m| {
        const res = mdl.residues.items[m.residue_idx];
        var i = overrides.entries.len;
        while (i > 0) {
            i -= 1;
            const entry = &overrides.entries[i];
            if (!entry.selector.matches(mdl, m)) continue;
            if (!std.mem.eql(u8, entry.compIdSlice(), res.compIdSlice())) continue;
            if (!std.mem.eql(u8, moverTargetName(mdl, m), entry.target)) continue;
            matched[i] = true;
            if (entryOrientation(m, entry.value)) |idx| {
                m.lockToOrientation(idx);
            } else {
                std.log.warn("fix override for {s}:{d} {s} {s}: value is not valid for this mover (has {d} orientations)", .{
                    entry.selector.chain_id,
                    entry.selector.auth_seq_id,
                    entry.compIdSlice(),
                    entry.target,
                    m.nOrientations(),
                });
                return error.InvalidFixOverride;
            }
            break;
        }
    }

    for (overrides.entries, matched) |*entry, hit| {
        if (hit) continue;
        std.log.warn("fix override for {s}:{d} {s} {s} matches no mover; ignored", .{
            entry.selector.chain_id,
            entry.selector.auth_seq_id,
            entry.compIdSlice(),
            entry.target,
        });
    }
}

pub fn dumpMovers(writer: anytype, mdl: *const Model, movers: []const Mover) !void {
    try writer.writeAll(
        "# selector comp_id target current states\n",
    );
    for (movers) |*m| {
        const res = mdl.residues.items[m.residue_idx];
        const chain = mdl.chains.items[res.chain_idx];
        const target = moverTargetName(mdl, m);
        try writer.print("{s}:{d}", .{ chain.authSlice(), res.auth_seq_id });
        if (res.ins_code != ' ') try writer.print(":{c}", .{res.ins_code});
        try writer.print(" {s} {s} {d} ", .{ res.compIdSlice(), target, m.best_orientation });
        try dumpStates(writer, m);
        try writer.writeByte('\n');
    }
}

fn dumpStates(writer: anytype, m: *const Mover) !void {
    switch (m.kind) {
        .amide_flip => try writer.writeAll("ORIGINAL|FLIP"),
        // HIP (2 orientations) falls through to indices: 0 original, 1 flipped.
        .his_flip => if (m.nOrientations() == 4) try writer.writeAll("HIE|HID|HIE_FLIP|HID_FLIP") else try writer.writeAll("0|1"),
        else => try writer.print("0..{d}", .{m.nOrientations() - 1}),
    }
}

pub fn moverTargetName(mdl: *const Model, m: *const Mover) []const u8 {
    return switch (m.kind) {
        .amide_flip => "amide",
        .his_flip => "his",
        .carboxyl_flip => "carboxyl",
        .single_h_rotator, .nh3_rotator, .nh2_rotator => blk: {
            const center_idx = m.center_idx orelse break :blk "unknown";
            break :blk mdl.atoms.items[center_idx].nameSlice();
        },
    };
}

fn entryOrientation(m: *const Mover, value: Value) ?u16 {
    return switch (value) {
        .orientation => |idx| if (idx < m.nOrientations()) idx else null,
        .amide => |amide| switch (m.kind) {
            .amide_flip => switch (amide) {
                .original => 0,
                .flip => 1,
            },
            else => null,
        },
        .his => |his| switch (m.kind) {
            // Tautomer names apply to the neutral 4-orientation mover only.
            .his_flip => if (m.nOrientations() != 4) null else switch (his) {
                .hie => 0,
                .hid => 1,
                .hie_flip => 2,
                .hid_flip => 3,
            },
            else => null,
        },
    };
}

fn parseSelector(allocator: std.mem.Allocator, token: []const u8) !Selector {
    var parts = std.mem.splitScalar(u8, token, ':');
    const chain_tok = parts.next() orelse return error.InvalidFixOverride;
    const seq_tok = parts.next() orelse return error.InvalidFixOverride;
    const ins_tok = parts.next();
    if (parts.next() != null) return error.InvalidFixOverride;
    const auth_seq_id = std.fmt.parseInt(i32, seq_tok, 10) catch return error.InvalidFixOverride;
    var ins_code: u8 = ' ';
    if (ins_tok) |tok| {
        if (tok.len == 1) ins_code = tok[0] else if (!std.mem.eql(u8, tok, ".")) return error.InvalidFixOverride;
    }
    return .{
        .chain_id = try allocator.dupe(u8, chain_tok),
        .auth_seq_id = auth_seq_id,
        .ins_code = ins_code,
    };
}

fn parseValue(target_tok: []const u8, value_tok: []const u8) !Value {
    if (std.ascii.eqlIgnoreCase(target_tok, "amide")) {
        if (std.ascii.eqlIgnoreCase(value_tok, "ORIGINAL")) return .{ .amide = .original };
        if (std.ascii.eqlIgnoreCase(value_tok, "FLIP")) return .{ .amide = .flip };
        return error.InvalidFixOverride;
    }
    if (std.ascii.eqlIgnoreCase(target_tok, "his")) {
        if (std.ascii.eqlIgnoreCase(value_tok, "HIE")) return .{ .his = .hie };
        if (std.ascii.eqlIgnoreCase(value_tok, "HID")) return .{ .his = .hid };
        if (std.ascii.eqlIgnoreCase(value_tok, "HIE_FLIP")) return .{ .his = .hie_flip };
        if (std.ascii.eqlIgnoreCase(value_tok, "HID_FLIP")) return .{ .his = .hid_flip };
        return error.InvalidFixOverride;
    }
    const idx = std.fmt.parseInt(u16, value_tok, 10) catch return error.InvalidFixOverride;
    return .{ .orientation = idx };
}

const testing = std.testing;
const mmcif = @import("../mmcif.zig");
const place = @import("../place.zig");
const optimize = @import("optimize.zig");

test "parse fix override file" {
    var overrides = try parseString(testing.allocator,
        \\A:1 ASN amide FLIP
        \\A:2 HIS his HID_FLIP
        \\B:9 SER OG 6
    );
    defer overrides.deinit();

    try testing.expectEqual(@as(usize, 3), overrides.entries.len);
    try testing.expectEqualStrings("amide", overrides.entries[0].target);
    try testing.expect(overrides.entries[0].value == .amide);
    try testing.expectEqual(@as(u16, 6), overrides.entries[2].value.orientation);
}

test "parser rejects invalid input" {
    // Missing tokens
    try testing.expectError(error.InvalidFixOverride, parseString(testing.allocator, "A:1 ASN amide"));
    // Extra tokens
    try testing.expectError(error.InvalidFixOverride, parseString(testing.allocator, "A:1 ASN amide FLIP extra"));
    // Invalid comp_id (too long - more than 5 chars)
    try testing.expectError(error.InvalidFixOverride, parseString(testing.allocator, "A:1 ASNXYZ amide FLIP"));
    // Invalid selector
    try testing.expectError(error.InvalidFixOverride, parseString(testing.allocator, "A ASN amide FLIP"));
    // Invalid value for target
    try testing.expectError(error.InvalidFixOverride, parseString(testing.allocator, "A:1 ASN amide BADVALUE"));
    // Invalid numeric value for rotator
    try testing.expectError(error.InvalidFixOverride, parseString(testing.allocator, "A:1 SER OG abc"));
}

test "apply fix locks his mover state" {
    const source = @embedFile("../test_data/his.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();
    mdl.residues.items[0].auth_seq_id = 1;

    place.applyChemistry(&mdl);
    _ = try place.addHydrogens(&mdl, null, null);
    const gen = try optimize.generateMovers(testing.allocator, &mdl, false, null, null, null, .neutron);
    const movers = gen.movers;
    defer {
        for (movers) |*m| m.deinit();
        testing.allocator.free(movers);
    }

    var overrides = try parseString(testing.allocator,
        \\A:1 HIS his HID_FLIP
    );
    defer overrides.deinit();
    try applyFixes(&overrides, &mdl, movers);

    // HIS N-terminal: his_flip (movers[0]) + nh3_rotator (movers[1])
    try testing.expectEqual(@as(usize, 2), movers.len);
    try testing.expect(movers[0].is_fixed);
    try testing.expectEqual(@as(u16, 3), movers[0].best_orientation);
}

test "apply fix locks amide and rotator movers" {
    {
        const source = @embedFile("../test_data/asn.cif");
        var mdl = try mmcif.parseModel(testing.allocator, source);
        defer mdl.deinit();
        mdl.residues.items[0].auth_seq_id = 1;

        place.applyChemistry(&mdl);
        _ = try place.addHydrogens(&mdl, null, null);
        const gen = try optimize.generateMovers(testing.allocator, &mdl, false, null, null, null, .neutron);
        const movers = gen.movers;
        defer {
            for (movers) |*m| m.deinit();
            testing.allocator.free(movers);
        }

        var overrides = try parseString(testing.allocator,
            \\A:1 ASN amide FLIP
        );
        defer overrides.deinit();
        try applyFixes(&overrides, &mdl, movers);

        try testing.expect(movers[0].is_fixed);
        try testing.expectEqual(@as(u16, 1), movers[0].best_orientation);
    }

    {
        const source = @embedFile("../test_data/tiny.cif");
        var mdl = try mmcif.parseModel(testing.allocator, source);
        defer mdl.deinit();
        mdl.residues.items[0].auth_seq_id = 1;

        place.applyChemistry(&mdl);
        _ = try place.addHydrogens(&mdl, null, null);
        const gen = try optimize.generateMovers(testing.allocator, &mdl, false, null, null, null, .neutron);
        const movers = gen.movers;
        defer {
            for (movers) |*m| m.deinit();
            testing.allocator.free(movers);
        }

        var overrides = try parseString(testing.allocator,
            \\A:1 ALA N 2
        );
        defer overrides.deinit();
        try applyFixes(&overrides, &mdl, movers);

        try testing.expect(movers[0].is_fixed);
        try testing.expectEqual(@as(u16, 2), movers[0].best_orientation);
    }
}

test "dump movers lists symbolic states" {
    const source = @embedFile("../test_data/his.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();
    mdl.residues.items[0].auth_seq_id = 1;

    place.applyChemistry(&mdl);
    _ = try place.addHydrogens(&mdl, null, null);
    const gen = try optimize.generateMovers(testing.allocator, &mdl, false, null, null, null, .neutron);
    const movers = gen.movers;
    defer {
        for (movers) |*m| m.deinit();
        testing.allocator.free(movers);
    }

    var buf: std.Io.Writer.Allocating = .init(testing.allocator);
    defer buf.deinit();
    try dumpMovers(&buf.writer, &mdl, movers);
    try testing.expect(std.mem.find(u8, buf.writer.buffered(), "A:1 HIS his 0 HIE|HID|HIE_FLIP|HID_FLIP") != null);
}

test "applyFixes rejects out-of-range orientation" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();
    mdl.residues.items[0].auth_seq_id = 1;

    place.applyChemistry(&mdl);
    _ = try place.addHydrogens(&mdl, null, null);
    const gen = try optimize.generateMovers(testing.allocator, &mdl, false, null, null, null, .neutron);
    const movers = gen.movers;
    defer {
        for (movers) |*m| m.deinit();
        testing.allocator.free(movers);
    }

    var overrides = try parseString(testing.allocator,
        \\A:1 ALA N 999
    );
    defer overrides.deinit();
    try testing.expectError(error.InvalidFixOverride, applyFixes(&overrides, &mdl, movers));
    try testing.expect(!movers[0].is_fixed);
}

test "fixed mover preserved through optimize()" {
    const source = @embedFile("../test_data/his.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();
    mdl.residues.items[0].auth_seq_id = 1;

    place.applyChemistry(&mdl);
    _ = try place.addHydrogens(&mdl, null, null);
    const gen = try optimize.generateMovers(testing.allocator, &mdl, false, null, null, null, .neutron);
    const movers = gen.movers;
    defer {
        for (movers) |*m| m.deinit();
        testing.allocator.free(movers);
    }

    // Lock the His flipper to HID_FLIP (orientation 3)
    movers[0].lockToOrientation(3);
    movers[0].applyOrientation(mdl.atoms.items, 3);

    const opt_result = try @import("optimizer.zig").optimize(testing.allocator, movers, &mdl, .{});
    _ = opt_result;

    // Orientation must be preserved after optimization
    try testing.expect(movers[0].is_fixed);
    try testing.expectEqual(@as(u16, 3), movers[0].best_orientation);
}

test "applyFixes ignores an entry that matches no mover" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();
    mdl.residues.items[0].auth_seq_id = 1;

    place.applyChemistry(&mdl);
    _ = try place.addHydrogens(&mdl, null, null);
    const gen = try optimize.generateMovers(testing.allocator, &mdl, false, null, null, null, .neutron);
    const movers = gen.movers;
    defer {
        for (movers) |*m| m.deinit();
        testing.allocator.free(movers);
    }

    // The ALA CB methyl has no mover; the N-terminal NH3+ does.
    var overrides = try parseString(testing.allocator,
        \\A:1 ALA CB 1
        \\A:1 ALA N 2
    );
    defer overrides.deinit();
    try applyFixes(&overrides, &mdl, movers);

    try testing.expect(movers[0].is_fixed);
    try testing.expectEqual(@as(u16, 2), movers[0].best_orientation);
}
