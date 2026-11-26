//! Post-placement model validation.
//!
//! Checks the model for common issues after hydrogen placement and optimization:
//! - Sentinel positions (absent H atoms that should have been removed)
//! - Atoms with NaN/Inf coordinates

const std = @import("std");
const model_mod = @import("model.zig");
const Model = model_mod.Model;
const Atom = model_mod.Atom;
const mover_mod = @import("optimize/mover.zig");

pub const Issue = struct {
    kind: Kind,
    atom_idx: u32,
    detail: [64]u8 = @splat(0),
    detail_len: u8 = 0,

    pub const Kind = enum {
        sentinel_position,
        nan_coordinate,
        inf_coordinate,
    };

    pub fn detailSlice(self: *const Issue) []const u8 {
        return self.detail[0..self.detail_len];
    }
};

pub const ValidationResult = struct {
    issues: []Issue,
    allocator: std.mem.Allocator,

    pub fn deinit(self: *ValidationResult) void {
        self.allocator.free(self.issues);
    }

    pub fn ok(self: *const ValidationResult) bool {
        return self.issues.len == 0;
    }
};

/// Validate a model after placement and optimization.
/// Returns a list of issues found. Caller must call deinit on the result.
pub fn validateModel(allocator: std.mem.Allocator, mdl: *const Model) !ValidationResult {
    var issues = std.ArrayListUnmanaged(Issue).empty;
    errdefer issues.deinit(allocator);

    for (mdl.atoms.items, 0..) |atom, idx| {
        // Check for sentinel positions (absent H from flipper)
        if (mover_mod.isAbsentH(atom)) {
            var issue = Issue{
                .kind = .sentinel_position,
                .atom_idx = @intCast(idx),
            };
            const name = atom.nameSlice();
            const len: u8 = @intCast(@min(name.len, 64));
            @memcpy(issue.detail[0..len], name[0..len]);
            issue.detail_len = len;
            try issues.append(allocator, issue);
        }

        // Check for NaN coordinates
        if (std.math.isNan(atom.pos.x) or std.math.isNan(atom.pos.y) or std.math.isNan(atom.pos.z)) {
            try issues.append(allocator, .{
                .kind = .nan_coordinate,
                .atom_idx = @intCast(idx),
            });
        }

        // Check for Inf coordinates
        if (std.math.isInf(atom.pos.x) or std.math.isInf(atom.pos.y) or std.math.isInf(atom.pos.z)) {
            try issues.append(allocator, .{
                .kind = .inf_coordinate,
                .atom_idx = @intCast(idx),
            });
        }
    }

    const result_issues = try allocator.dupe(Issue, issues.items);
    issues.deinit(allocator);

    return .{
        .issues = result_issues,
        .allocator = allocator,
    };
}

/// Print validation issues to stderr.
pub fn reportIssues(issues: []const Issue, mdl: *const Model) void {
    for (issues) |issue| {
        const atom = mdl.atoms.items[issue.atom_idx];
        const res = mdl.residues.items[atom.residue_idx];
        const comp_id = res.compIdSlice();
        const name = atom.nameSlice();
        switch (issue.kind) {
            .sentinel_position => {
                std.debug.print("WARN: sentinel position on {s} {s} (res {d}) — absent H not removed\n", .{ comp_id, name, atom.residue_idx });
            },
            .nan_coordinate => {
                std.debug.print("ERROR: NaN coordinate on {s} {s} (res {d})\n", .{ comp_id, name, atom.residue_idx });
            },
            .inf_coordinate => {
                std.debug.print("ERROR: Inf coordinate on {s} {s} (res {d})\n", .{ comp_id, name, atom.residue_idx });
            },
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

const testing = std.testing;

test "clean model passes validation" {
    const mmcif = @import("mmcif.zig");
    const source = @embedFile("test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    var result = try validateModel(testing.allocator, &mdl);
    defer result.deinit();

    try testing.expect(result.ok());
}

test "sentinel position detected" {
    var mdl = Model.init(testing.allocator);
    defer mdl.deinit();

    // Add a normal atom
    try mdl.atoms.append(testing.allocator, .{
        .pos = .{ .x = 1.0, .y = 2.0, .z = 3.0 },
    });

    // Add a sentinel atom
    var sentinel = Atom{ .pos = .{ .x = 1000.0, .y = 1000.0, .z = 1000.0 } };
    sentinel.is_added = true;
    sentinel.is_hydrogen = true;
    sentinel.setName("HD1");
    try mdl.atoms.append(testing.allocator, sentinel);

    // Need at least one residue for reporting
    try mdl.residues.append(testing.allocator, .{});

    var result = try validateModel(testing.allocator, &mdl);
    defer result.deinit();

    try testing.expect(!result.ok());
    try testing.expectEqual(@as(usize, 1), result.issues.len);
    try testing.expectEqual(Issue.Kind.sentinel_position, result.issues[0].kind);
}

test "NaN coordinate detected" {
    var mdl = Model.init(testing.allocator);
    defer mdl.deinit();
    try mdl.atoms.append(testing.allocator, .{
        .pos = .{ .x = std.math.nan(f32), .y = 2.0, .z = 3.0 },
    });
    try mdl.residues.append(testing.allocator, .{});
    var result = try validateModel(testing.allocator, &mdl);
    defer result.deinit();
    try testing.expect(!result.ok());
    try testing.expectEqual(Issue.Kind.nan_coordinate, result.issues[0].kind);
}

test "Inf coordinate detected" {
    var mdl = Model.init(testing.allocator);
    defer mdl.deinit();
    try mdl.atoms.append(testing.allocator, .{
        .pos = .{ .x = std.math.inf(f32), .y = 2.0, .z = 3.0 },
    });
    try mdl.residues.append(testing.allocator, .{});
    var result = try validateModel(testing.allocator, &mdl);
    defer result.deinit();
    try testing.expect(!result.ok());
    try testing.expectEqual(Issue.Kind.inf_coordinate, result.issues[0].kind);
}
