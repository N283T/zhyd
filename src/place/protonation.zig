const std = @import("std");

fn defaultIo() std.Io {
    return std.Io.Threaded.global_single_threaded.io();
}
const model_mod = @import("../model.zig");
const fixed_string = @import("../model/fixed_string.zig");

const FixedString = fixed_string.FixedString;

const Model = model_mod.Model;

pub const HisState = enum { auto, hid, hie, hip };
pub const AcidState = enum { deprotonated, atom1, atom2 };
pub const LysState = enum { charged, neutral };
pub const CysState = enum { thiol, thiolate };

pub const ResidueState = union(enum) {
    his: HisState,
    asp: AcidState,
    glu: AcidState,
    lys: LysState,
    cys: CysState,
};

pub const ResidueSelector = struct {
    chain_id: []const u8,
    auth_seq_id: i32,
    ins_code: u8 = ' ',

    pub fn matches(self: ResidueSelector, mdl: *const Model, res_idx: usize) bool {
        const res = mdl.residues.items[res_idx];
        const chain = mdl.chains.items[res.chain_idx];
        if (res.auth_seq_id != self.auth_seq_id) return false;
        if (res.ins_code != self.ins_code) return false;
        return std.mem.eql(u8, chain.authSlice(), self.chain_id) or
            std.mem.eql(u8, chain.labelSlice(), self.chain_id);
    }
};

pub const Entry = struct {
    selector: ResidueSelector,
    comp_id: FixedString(5) = .{},
    state: ResidueState,

    pub fn compIdSlice(self: *const Entry) []const u8 {
        return self.comp_id.slice();
    }
};

pub const ProtonationOverrides = struct {
    allocator: std.mem.Allocator,
    entries: []Entry,

    pub fn deinit(self: *ProtonationOverrides) void {
        for (self.entries) |entry| self.allocator.free(entry.selector.chain_id);
        self.allocator.free(self.entries);
        self.entries = &.{};
    }

    /// Find the protonation state override for a residue.
    /// When multiple entries match, the last entry in the file wins.
    pub fn find(self: *const ProtonationOverrides, mdl: *const Model, res_idx: usize) ?ResidueState {
        const res = mdl.residues.items[res_idx];
        var i = self.entries.len;
        while (i > 0) {
            i -= 1;
            const entry = &self.entries[i];
            if (!std.mem.eql(u8, entry.compIdSlice(), res.compIdSlice())) continue;
            if (entry.selector.matches(mdl, res_idx)) return entry.state;
        }
        return null;
    }

    /// Warn about override entries that did not match any residue.
    pub fn warnUnmatched(self: *const ProtonationOverrides, mdl: *const Model) void {
        for (self.entries) |entry| {
            var matched = false;
            for (0..mdl.residues.items.len) |res_idx| {
                if (!entry.selector.matches(mdl, res_idx)) continue;
                const res = mdl.residues.items[res_idx];
                if (std.mem.eql(u8, entry.compIdSlice(), res.compIdSlice())) {
                    matched = true;
                    break;
                } else {
                    std.log.warn("protonation override for {s}:{d} expects {s} but residue is {s}", .{
                        entry.selector.chain_id,
                        entry.selector.auth_seq_id,
                        entry.compIdSlice(),
                        res.compIdSlice(),
                    });
                    matched = true; // selector matched but comp_id mismatch
                    break;
                }
            }
            if (!matched) {
                std.log.warn("protonation override for {s}:{d} {s} did not match any residue", .{
                    entry.selector.chain_id,
                    entry.selector.auth_seq_id,
                    entry.compIdSlice(),
                });
            }
        }
    }
};

pub fn parseFile(allocator: std.mem.Allocator, path: []const u8) !ProtonationOverrides {
    const io = defaultIo();
    // readAlloc on a streaming reader reports EndOfStream for ordinary short
    // files under Zig 0.16. A bounded whole-file read has the intended
    // semantics for this deliberately small control format.
    const source = try std.Io.Dir.cwd().readFileAlloc(io, path, allocator, .limited(1024 * 1024));
    defer allocator.free(source);
    return parseString(allocator, source);
}

pub fn parseString(allocator: std.mem.Allocator, source: []const u8) !ProtonationOverrides {
    var entries = std.ArrayListUnmanaged(Entry).empty;
    errdefer {
        for (entries.items) |entry| allocator.free(entry.selector.chain_id);
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
            std.log.warn("protonation override line {d}: expected 'chain:seq comp state'", .{line_num});
            return error.InvalidProtonationOverride;
        };
        const comp_tok = toks.next() orelse {
            std.log.warn("protonation override line {d}: missing comp_id and state", .{line_num});
            return error.InvalidProtonationOverride;
        };
        const state_tok = toks.next() orelse {
            std.log.warn("protonation override line {d}: missing state", .{line_num});
            return error.InvalidProtonationOverride;
        };
        if (toks.next() != null) {
            std.log.warn("protonation override line {d}: unexpected extra tokens", .{line_num});
            return error.InvalidProtonationOverride;
        }
        if (comp_tok.len > 5) {
            std.log.warn("protonation override line {d}: comp_id '{s}' exceeds 5 characters", .{ line_num, comp_tok });
            return error.InvalidProtonationOverride;
        }

        const selector = parseSelector(allocator, selector_tok) catch {
            std.log.warn("protonation override line {d}: invalid selector '{s}' (expected chain:seq[:ins])", .{ line_num, selector_tok });
            return error.InvalidProtonationOverride;
        };
        errdefer allocator.free(selector.chain_id);
        const state = parseState(comp_tok, state_tok) catch {
            std.log.warn("protonation override line {d}: invalid state '{s}' for {s}", .{ line_num, state_tok, comp_tok });
            return error.InvalidProtonationOverride;
        };

        var comp_id: FixedString(5) = .{};
        const comp_id_len: u3 = @intCast(comp_tok.len);
        comp_id.len = comp_id_len;
        for (0..comp_id_len) |i| comp_id.buf[i] = std.ascii.toUpper(comp_tok[i]);

        entries.append(allocator, .{
            .selector = selector,
            .comp_id = comp_id,
            .state = state,
        }) catch |err| {
            allocator.free(selector.chain_id);
            return err;
        };
    }

    return .{
        .allocator = allocator,
        .entries = try entries.toOwnedSlice(allocator),
    };
}

fn parseSelector(allocator: std.mem.Allocator, token: []const u8) !ResidueSelector {
    var parts = std.mem.splitScalar(u8, token, ':');
    const chain_tok = parts.next() orelse return error.InvalidProtonationOverride;
    const seq_tok = parts.next() orelse return error.InvalidProtonationOverride;
    const ins_tok = parts.next();
    if (parts.next() != null) return error.InvalidProtonationOverride;
    if (chain_tok.len == 0) return error.InvalidProtonationOverride;

    const auth_seq_id = std.fmt.parseInt(i32, seq_tok, 10) catch return error.InvalidProtonationOverride;
    var ins_code: u8 = ' ';
    if (ins_tok) |tok| {
        if (tok.len == 0 or std.mem.eql(u8, tok, ".")) {
            ins_code = ' ';
        } else {
            if (tok.len != 1) return error.InvalidProtonationOverride;
            ins_code = tok[0];
        }
    }

    return .{
        .chain_id = try allocator.dupe(u8, chain_tok),
        .auth_seq_id = auth_seq_id,
        .ins_code = ins_code,
    };
}

fn parseState(comp_tok: []const u8, state_tok: []const u8) !ResidueState {
    var comp_buf: [3]u8 = .{ ' ', ' ', ' ' };
    for (comp_tok[0..@min(comp_tok.len, 3)], 0..) |c, i| comp_buf[i] = std.ascii.toUpper(c);
    const comp_id = std.mem.trimEnd(u8, &comp_buf, " ");

    if (std.mem.eql(u8, comp_id, "HIS")) {
        if (std.ascii.eqlIgnoreCase(state_tok, "AUTO")) return .{ .his = .auto };
        if (std.ascii.eqlIgnoreCase(state_tok, "HID")) return .{ .his = .hid };
        if (std.ascii.eqlIgnoreCase(state_tok, "HIE")) return .{ .his = .hie };
        if (std.ascii.eqlIgnoreCase(state_tok, "HIP")) return .{ .his = .hip };
        return error.InvalidProtonationOverride;
    }
    if (std.mem.eql(u8, comp_id, "ASP")) {
        if (std.ascii.eqlIgnoreCase(state_tok, "DEPROTONATED")) return .{ .asp = .deprotonated };
        if (std.ascii.eqlIgnoreCase(state_tok, "OD1")) return .{ .asp = .atom1 };
        if (std.ascii.eqlIgnoreCase(state_tok, "OD2")) return .{ .asp = .atom2 };
        return error.InvalidProtonationOverride;
    }
    if (std.mem.eql(u8, comp_id, "GLU")) {
        if (std.ascii.eqlIgnoreCase(state_tok, "DEPROTONATED")) return .{ .glu = .deprotonated };
        if (std.ascii.eqlIgnoreCase(state_tok, "OE1")) return .{ .glu = .atom1 };
        if (std.ascii.eqlIgnoreCase(state_tok, "OE2")) return .{ .glu = .atom2 };
        return error.InvalidProtonationOverride;
    }
    if (std.mem.eql(u8, comp_id, "LYS")) {
        if (std.ascii.eqlIgnoreCase(state_tok, "CHARGED")) return .{ .lys = .charged };
        if (std.ascii.eqlIgnoreCase(state_tok, "NEUTRAL")) return .{ .lys = .neutral };
        return error.InvalidProtonationOverride;
    }
    if (std.mem.eql(u8, comp_id, "CYS")) {
        if (std.ascii.eqlIgnoreCase(state_tok, "THIOL")) return .{ .cys = .thiol };
        if (std.ascii.eqlIgnoreCase(state_tok, "THIOLATE")) return .{ .cys = .thiolate };
        return error.InvalidProtonationOverride;
    }
    return error.InvalidProtonationOverride;
}

const testing = std.testing;

test "parse protonation override file" {
    var overrides = try parseString(testing.allocator,
        \\# chain:auth_seq[:ins] comp state
        \\A:42 HIS HIE
        \\B:7:A ASP OD2
        \\C:9 LYS neutral
    );
    defer overrides.deinit();

    try testing.expectEqual(@as(usize, 3), overrides.entries.len);
    try testing.expectEqualStrings("A", overrides.entries[0].selector.chain_id);
    try testing.expectEqual(@as(i32, 42), overrides.entries[0].selector.auth_seq_id);
    try testing.expect(overrides.entries[0].state == .his);
    try testing.expectEqual(HisState.hie, overrides.entries[0].state.his);
    try testing.expectEqual(@as(u8, 'A'), overrides.entries[1].selector.ins_code);
    try testing.expectEqual(LysState.neutral, overrides.entries[2].state.lys);
}

test "selector matches residue by chain and auth seq" {
    var mdl = model_mod.Model.init(testing.allocator);
    defer mdl.deinit();

    var chain = model_mod.Chain{};
    chain.setLabelAsymId("A");
    chain.setAuthAsymId("A");
    chain.residue_start = 0;
    chain.residue_end = 1;
    try mdl.chains.append(testing.allocator, chain);

    var residue = model_mod.Residue{};
    residue.setCompId("HIS");
    residue.chain_idx = 0;
    residue.auth_seq_id = 1;
    residue.atom_start = 0;
    residue.atom_end = 0;
    try mdl.residues.append(testing.allocator, residue);

    var overrides = try parseString(testing.allocator,
        \\A:1 HIS HID
    );
    defer overrides.deinit();

    try testing.expect(overrides.entries[0].selector.matches(&mdl, 0));
    try testing.expectEqual(HisState.hid, overrides.entries[0].state.his);
}

test "parser rejects invalid input" {
    // Extra token
    try testing.expectError(error.InvalidProtonationOverride, parseString(testing.allocator, "A:1 HIS HIE extra"));
    // Missing state
    try testing.expectError(error.InvalidProtonationOverride, parseString(testing.allocator, "A:1 HIS"));
    // Unknown state for known residue
    try testing.expectError(error.InvalidProtonationOverride, parseString(testing.allocator, "A:1 HIS BADSTATE"));
    // Unknown residue type
    try testing.expectError(error.InvalidProtonationOverride, parseString(testing.allocator, "A:1 ALA SOMETHING"));
    // No colon in selector
    try testing.expectError(error.InvalidProtonationOverride, parseString(testing.allocator, "A HIS HIE"));
    // Non-numeric seq_id
    try testing.expectError(error.InvalidProtonationOverride, parseString(testing.allocator, "A:XYZ HIS HIE"));
    // comp_id > 5 characters
    try testing.expectError(error.InvalidProtonationOverride, parseString(testing.allocator, "A:1 HISTXY HIE"));
}
