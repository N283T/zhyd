//! Mover interface and types for hydrogen optimization.
//! A Mover controls one or more atoms and provides a set of discrete orientations.

const std = @import("std");
const math_mod = @import("../math.zig");
const model_mod = @import("../model.zig");
const Vec3 = math_mod.Vec3;
const Atom = model_mod.Atom;
const element = @import("../element.zig");

pub const MoverKind = enum(u4) {
    single_h_rotator, // OH, SH, SeH — 12 orientations at 30° intervals
    nh3_rotator, // NH3+ — 3 orientations (120° spacing)
    methyl_rotator, // CH3 — 3 orientations (60° stagger from initial)
    aromatic_methyl, // Aromatic CH3 — same as methyl
    nh2_rotator, // NH2 (non-PRO neutral N-term, PRO NH2+) — 12 orientations at 30°
    amide_flip, // Asn/Gln — 2 orientations (original + O↔N swap)
    his_flip, // His — 6 orientations (2 flip × 3 protonation)
};

pub const Orientation = struct {
    positions: []Vec3(f32), // atom positions for this orientation
    flags: ?[]element.AtomFlags = null, // per-atom flags (for flip state chemistry updates)
    penalty: f32 = 0.0,
};

pub const Mover = struct {
    kind: MoverKind,
    residue_idx: u32,
    atom_indices: []u32, // indices into model.atoms that this mover controls
    orientations: []Orientation,
    best_orientation: u16 = 0,
    current_orientation: u16 = 0,
    is_fixed: bool = false,
    allocator: std.mem.Allocator,

    // Rotation axis geometry (for fine search, rotators only)
    center_idx: ?u32 = null,
    axis_idx: ?u32 = null,

    /// Optional backing allocation for all orientation positions in one contiguous block.
    /// When set, `deinit` frees this single block instead of each orientation's positions
    /// slice individually. Each orientation's `positions` field points into this block.
    positions_backing: ?[]Vec3(f32) = null,

    pub fn deinit(self: *Mover) void {
        if (self.positions_backing) |backing| {
            // Positions are sub-slices of the backing block; free the block once.
            self.allocator.free(backing);
        } else {
            for (self.orientations) |o| {
                self.allocator.free(o.positions);
            }
        }
        for (self.orientations) |o| {
            if (o.flags) |f| self.allocator.free(f);
        }
        self.allocator.free(self.orientations);
        self.allocator.free(self.atom_indices);
    }

    /// Sentinel position used by flipper for absent H atoms.
    pub const ABSENT_H_POS = Vec3(f32){ .x = 1000.0, .y = 1000.0, .z = 1000.0 };

    /// Apply the given orientation: update atom positions and optionally flags.
    pub fn applyOrientation(self: *const Mover, atoms: []Atom, idx: u16) void {
        const orient = self.orientations[idx];
        for (self.atom_indices, 0..) |ai, i| {
            atoms[ai].pos = orient.positions[i];
        }
        if (orient.flags) |flags| {
            for (self.atom_indices, 0..) |ai, i| {
                atoms[ai].flags = flags[i];
            }
        }
    }

    /// Get the penalty for a given orientation.
    pub fn orientationPenalty(self: *const Mover, idx: u16) f32 {
        return self.orientations[idx].penalty;
    }

    pub fn nOrientations(self: *const Mover) u16 {
        return @intCast(self.orientations.len);
    }

    pub fn effectiveOrientations(self: *const Mover) u16 {
        return if (self.is_fixed) 1 else self.nOrientations();
    }

    pub fn lockToOrientation(self: *Mover, idx: u16) void {
        self.best_orientation = idx;
        self.current_orientation = idx;
        self.is_fixed = true;
    }
};

/// Check if an atom is at the sentinel (absent) position.
pub fn isAbsentH(atom: Atom) bool {
    return atom.is_added and isAbsentPos(atom.pos);
}

/// Whether a position is the `ABSENT_H_POS` sentinel. Absent H atoms of
/// different movers share this position, so scoring and the interaction graph
/// must skip them or they bump into each other.
pub fn isAbsentPos(pos: Vec3(f32)) bool {
    return pos.x > 999.0 and pos.y > 999.0 and pos.z > 999.0;
}

// ── Tests ─────────────────────────────────────────────────────────────────────

const testing = std.testing;

test "Mover applyOrientation updates atom positions" {
    const allocator = testing.allocator;

    var atoms = [_]Atom{
        .{ .pos = .{ .x = 2.0, .y = 0.0, .z = 0.0 } },
    };

    const positions = try allocator.alloc(Vec3(f32), 1);
    positions[0] = .{ .x = 0.0, .y = 2.0, .z = 0.0 };

    const orientations = try allocator.alloc(Orientation, 1);
    orientations[0] = .{ .positions = positions, .penalty = 0.0 };

    const atom_indices = try allocator.alloc(u32, 1);
    atom_indices[0] = 0;

    var mover = Mover{
        .kind = .single_h_rotator,
        .residue_idx = 0,
        .atom_indices = atom_indices,
        .orientations = orientations,
        .allocator = allocator,
    };
    defer mover.deinit();

    mover.applyOrientation(&atoms, 0);
    try testing.expectApproxEqAbs(@as(f32, 0.0), atoms[0].pos.x, 0.001);
    try testing.expectApproxEqAbs(@as(f32, 2.0), atoms[0].pos.y, 0.001);
}

test "Mover nOrientations" {
    const allocator = testing.allocator;

    const pos0 = try allocator.alloc(Vec3(f32), 1);
    pos0[0] = .{ .x = 1.0, .y = 0.0, .z = 0.0 };
    const pos1 = try allocator.alloc(Vec3(f32), 1);
    pos1[0] = .{ .x = 0.0, .y = 1.0, .z = 0.0 };

    const orientations = try allocator.alloc(Orientation, 2);
    orientations[0] = .{ .positions = pos0, .penalty = 0.0 };
    orientations[1] = .{ .positions = pos1, .penalty = 0.5 };

    const atom_indices = try allocator.alloc(u32, 1);
    atom_indices[0] = 0;

    var mover = Mover{
        .kind = .single_h_rotator,
        .residue_idx = 0,
        .atom_indices = atom_indices,
        .orientations = orientations,
        .allocator = allocator,
    };
    defer mover.deinit();

    try testing.expectEqual(@as(u16, 2), mover.nOrientations());
    try testing.expectApproxEqAbs(@as(f32, 0.5), mover.orientationPenalty(1), 0.001);
}
