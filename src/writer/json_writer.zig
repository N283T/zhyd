//! JSON log writer: outputs optimization results as JSON.

const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;
const model_mod = @import("../model.zig");
const mover_mod = @import("../optimize/mover.zig");
const bond_policy_mod = @import("../place/bond_policy.zig");

/// Snapshot of a mover's result, captured before the mover is freed.
/// Used to accumulate mover info across per-model loops for JSON log.
pub const MoverSnapshot = struct {
    chain_label: [4]u8 = .{ ' ', ' ', ' ', ' ' },
    chain_label_len: u8 = 0,
    comp_id: [4]u8 = .{ ' ', ' ', ' ', ' ' },
    comp_id_len: u8 = 0,
    seq_id: i32 = 0,
    model_num: u32 = 1,
    kind: mover_mod.MoverKind = .single_h_rotator,
    best_orientation: u16 = 0,

    pub fn capture(
        m: mover_mod.Mover,
        residues: []const model_mod.Residue,
        chains: []const model_mod.Chain,
        model_num: u32,
    ) MoverSnapshot {
        const res = residues[m.residue_idx];
        const chain = chains[res.chain_idx];
        var snap = MoverSnapshot{
            .seq_id = res.seq_id,
            .model_num = model_num,
            .kind = m.kind,
            .best_orientation = m.best_orientation,
        };
        const label = chain.labelSlice();
        const llen: u8 = @intCast(@min(label.len, 4));
        @memcpy(snap.chain_label[0..llen], label[0..llen]);
        snap.chain_label_len = llen;
        const comp = res.compIdSlice();
        const clen: u8 = @intCast(@min(comp.len, 4));
        @memcpy(snap.comp_id[0..clen], comp[0..clen]);
        snap.comp_id_len = clen;
        return snap;
    }

    // Both accessors take a pointer: a by-value `self` is a temporary copy, and the slice
    // would point into it after the call returns.
    pub fn chainLabel(self: *const MoverSnapshot) []const u8 {
        return self.chain_label[0..self.chain_label_len];
    }

    pub fn compId(self: *const MoverSnapshot) []const u8 {
        return self.comp_id[0..self.comp_id_len];
    }
};

/// Write JSON log of optimization results.
pub fn writeLog(
    out_writer: anytype,
    version: []const u8,
    input_file: []const u8,
    n_hydrogens: u32,
    bond_policy: bond_policy_mod.BondPolicy,
    movers: []const mover_mod.Mover,
    residues: []const model_mod.Residue,
    chains: []const model_mod.Chain,
) !void {
    try out_writer.writeAll("{\n");
    try out_writer.writeAll("  \"version\": \"");
    try writeJsonString(out_writer, version);
    try out_writer.writeAll("\",\n");
    try out_writer.writeAll("  \"input\": \"");
    try writeJsonString(out_writer, input_file);
    try out_writer.writeAll("\",\n");
    try out_writer.print("  \"hydrogens_added\": {d},\n", .{n_hydrogens});
    try out_writer.print("  \"bond_mode\": \"{s}\",\n", .{bondModeStr(bond_policy.mode)});
    try out_writer.print("  \"output_isotope\": \"{s}\",\n", .{outputIsotopeStr(bond_policy.output_isotope)});
    try out_writer.writeAll("  \"movers\": [\n");

    for (movers, 0..) |m, i| {
        const res = residues[m.residue_idx];
        const chain = chains[res.chain_idx];
        const kind_str = moverKindStr(m.kind);

        try out_writer.print("    {{\"residue\": \"{s}.{s}.{d}\", \"type\": \"{s}\", \"orientation\": {d}}}", .{
            chain.labelSlice(),
            res.compIdSlice(),
            res.seq_id,
            kind_str,
            m.best_orientation,
        });

        if (i < movers.len - 1) try out_writer.writeAll(",");
        try out_writer.writeAll("\n");
    }

    try out_writer.writeAll("  ]\n}\n");
}

/// Write JSON log for multi-model processing using captured mover snapshots.
pub fn writeMultiModelLog(
    out_writer: anytype,
    version: []const u8,
    input_file: []const u8,
    n_hydrogens: u32,
    bond_policy: bond_policy_mod.BondPolicy,
    snapshots: []const MoverSnapshot,
) !void {
    try out_writer.writeAll("{\n");
    try out_writer.writeAll("  \"version\": \"");
    try writeJsonString(out_writer, version);
    try out_writer.writeAll("\",\n");
    try out_writer.writeAll("  \"input\": \"");
    try writeJsonString(out_writer, input_file);
    try out_writer.writeAll("\",\n");
    try out_writer.print("  \"hydrogens_added\": {d},\n", .{n_hydrogens});
    try out_writer.print("  \"bond_mode\": \"{s}\",\n", .{bondModeStr(bond_policy.mode)});
    try out_writer.print("  \"output_isotope\": \"{s}\",\n", .{outputIsotopeStr(bond_policy.output_isotope)});
    try out_writer.writeAll("  \"movers\": [\n");

    for (snapshots, 0..) |*snap, i| {
        try out_writer.print("    {{\"model\": {d}, \"residue\": \"{s}.{s}.{d}\", \"type\": \"{s}\", \"orientation\": {d}}}", .{
            snap.model_num,
            snap.chainLabel(),
            snap.compId(),
            snap.seq_id,
            moverKindStr(snap.kind),
            snap.best_orientation,
        });

        if (i < snapshots.len - 1) try out_writer.writeAll(",");
        try out_writer.writeAll("\n");
    }

    try out_writer.writeAll("  ]\n}\n");
}

/// Write a JSON-escaped string (without surrounding quotes).
fn writeJsonString(writer: anytype, s: []const u8) !void {
    for (s) |c| {
        switch (c) {
            '"' => try writer.writeAll("\\\""),
            '\\' => try writer.writeAll("\\\\"),
            '\n' => try writer.writeAll("\\n"),
            '\r' => try writer.writeAll("\\r"),
            '\t' => try writer.writeAll("\\t"),
            else => {
                if (c < 0x20) {
                    try writer.print("\\u{x:0>4}", .{c});
                } else {
                    try writer.writeByte(c);
                }
            },
        }
    }
}

fn moverKindStr(kind: mover_mod.MoverKind) []const u8 {
    return switch (kind) {
        .single_h_rotator => "single_h_rotator",
        .nh3_rotator => "nh3_rotator",
        .nh2_rotator => "nh2_rotator",
        .amide_flip => "amide_flip",
        .his_flip => "his_flip",
        .carboxyl_flip => "carboxyl_flip",
    };
}

fn bondModeStr(mode: bond_policy_mod.BondLengthMode) []const u8 {
    return switch (mode) {
        .xray => "xray",
        .neutron => "neutron",
    };
}

fn outputIsotopeStr(isotope: bond_policy_mod.OutputIsotope) []const u8 {
    return switch (isotope) {
        .hydrogen => "hydrogen",
        .deuterium => "deuterium",
    };
}

test "write JSON log" {
    var buf: std.Io.Writer.Allocating = .init(testing.allocator);
    defer buf.deinit();

    try writeLog(&buf.writer, "0.1.0", "test.cif", 42, .{}, &.{}, &.{}, &.{});

    const output = buf.writer.buffered();
    try testing.expect(std.mem.find(u8, output, "\"version\": \"0.1.0\"") != null);
    try testing.expect(std.mem.find(u8, output, "\"hydrogens_added\": 42") != null);
    try testing.expect(std.mem.find(u8, output, "\"bond_mode\": \"neutron\"") != null);
    try testing.expect(std.mem.find(u8, output, "\"output_isotope\": \"hydrogen\"") != null);
}

test "multi-model JSON log keeps chain and residue names" {
    var buf: std.Io.Writer.Allocating = .init(testing.allocator);
    defer buf.deinit();

    // Many snapshots, so the writer runs long enough to reuse the stack slot
    // a by-value accessor would have returned a slice into.
    var snaps: [200]MoverSnapshot = undefined;
    for (&snaps, 0..) |*s, i| {
        s.* = .{ .chain_label_len = 1, .comp_id_len = 3, .seq_id = @intCast(i) };
        s.chain_label[0] = 'B';
        @memcpy(s.comp_id[0..3], "LEU");
    }
    try writeMultiModelLog(&buf.writer, "0.1.0", "test.cif", 0, .{}, &snaps);

    const output = buf.writer.buffered();
    try testing.expect(std.mem.indexOfScalar(u8, output, 0) == null);
    try testing.expectEqual(@as(usize, snaps.len), std.mem.count(u8, output, "\"residue\": \"B.LEU."));
}
