//! Desolvation and burial weights ("Desolvation").

const std = @import("std");
const params_mod = @import("params.zig");
const sites = @import("sites.zig");

const Parameters = params_mod.Parameters;

pub const Burial = struct {
    /// Atom count N within the count radius.
    count: f64,
    /// Volume sum V.
    volume: f64,
};

/// Count and volume sum around `center`, skipping atoms whose (residue
/// number, chain) equal those of `group_atom` (insertion code ignored). With
/// Q1 only atoms of the group atom's own residue are skipped. Virtual atoms
/// (Q3 C-terminal oxygen) never count.
pub fn burial(conf: *const sites.Conformation, group_atom: u32, center: sites.Vec, params: Parameters) Burial {
    const ga = conf.atoms.items[group_atom];
    const min4 = std.math.pow(f64, params.desolv_min_distance, 4);
    var count: f64 = 0;
    var volume: f64 = 0;
    for (conf.atoms.items) |a| {
        if (!a.isHeavy() or a.virtual) continue;
        if (params.fixes.q1) {
            if (a.residue_idx == ga.residue_idx) continue;
        } else if (a.res_num == ga.res_num and a.chain_key == ga.chain_key) continue;
        const d = a.pos.distance(center);
        if (d < params.desolv_volume_radius) {
            const d4 = d * d * d * d;
            volume += params_mod.relativeVolume(sites.volumeClass(a)) / @max(min4, d4);
        }
        if (d < params.desolv_count_radius) count += 1;
    }
    return .{ .count = count, .volume = volume };
}

/// Burial weight w = clamp((N − 280) / (560 − 280), 0, 1).
pub fn weight(count: f64, params: Parameters) f64 {
    return std.math.clamp((count - params.buried_count_min) / (params.buried_count_max - params.buried_count_min), 0.0, 1.0);
}

/// Pair weight w₁₂ = clamp((N₁ + N₂ − 560) / (1120 − 560), 0, 1).
pub fn pairWeight(n1: f64, n2: f64, params: Parameters) f64 {
    const lo = 2.0 * params.buried_count_min;
    const hi = 2.0 * params.buried_count_max;
    return std.math.clamp((n1 + n2 - lo) / (hi - lo), 0.0, 1.0);
}

/// ΔpKa_desolv = q · (−13.0) · max(0, V − 0.0) · s, s = 1 − (1 − 0.25)(1 − w).
pub fn shift(q: f64, volume: f64, w: f64, params: Parameters) f64 {
    const s = 1.0 - (1.0 - params.surface_scale) * (1.0 - w);
    return q * params.desolv_prefactor * @max(0.0, volume - params.desolv_volume_offset) * s;
}

/// Fill count, weight and desolvation of every titratable group and ion. The
/// desolvation shift is multiplied by the site kind's scale factor (1.0
/// unless overridden).
pub fn apply(conf: *sites.Conformation, params: Parameters) void {
    for (conf.groups.items) |*g| {
        if (!(g.isTitratable() or g.gtype == .ion)) continue;
        const b = burial(conf, g.group_atom, g.center, params);
        g.count = b.count;
        g.weight = weight(b.count, params);
        if (g.isTitratable()) g.desolvation = shift(g.q, b.volume, g.weight, params) * params.desolvation_scale.get(g.kind.?);
    }
}

const testing = std.testing;

test "burial weight and pair weight clamp linearly" {
    const p = Parameters{};
    try testing.expectEqual(@as(f64, 0), weight(100, p));
    try testing.expectEqual(@as(f64, 0.5), weight(420, p));
    try testing.expectEqual(@as(f64, 1), weight(900, p));
    // (300 + 400 − 560) / 560 = 0.25
    try testing.expectEqual(@as(f64, 0.25), pairWeight(300, 400, p));
}

test "desolvation shift by hand" {
    const p = Parameters{};
    // Acid, V = 0.1, w = 0.5: s = 1 − 0.75·0.5 = 0.625; shift = −1·−13·0.1·0.625
    try testing.expectApproxEqAbs(@as(f64, 0.8125), shift(-1, 0.1, 0.5, p), 1e-12);
    // Base on the surface: s = 0.25; shift = 1·−13·0.2·0.25 = −0.65
    try testing.expectApproxEqAbs(@as(f64, -0.65), shift(1, 0.2, 0, p), 1e-12);
}

test "Q1: same-residue exclusion compares the full residue identity" {
    var conf = sites.Conformation{ .allocator = testing.allocator };
    defer conf.deinit();
    // Group atom in residue 0 (number 1); an oxygen 3 Å away in residue 1
    // with the same number and chain (an insertion-code variant).
    try conf.atoms.append(testing.allocator, .{ .pos = .{ .x = 0, .y = 0, .z = 0 }, .element = .o, .residue_idx = 0, .chain_key = 0, .res_num = 1, .order = 0 });
    try conf.atoms.append(testing.allocator, .{ .pos = .{ .x = 3, .y = 0, .z = 0 }, .element = .o, .residue_idx = 1, .chain_key = 0, .res_num = 1, .order = 1 });
    const center = sites.Vec{ .x = 0, .y = 0, .z = 0 };
    // Parity: the variant shares (number, chain) and is skipped.
    try testing.expectEqual(@as(f64, 0), burial(&conf, 0, center, .{}).count);
    // Q1: it is another residue and counts, with volume 1.00 / 3⁴.
    const b = burial(&conf, 0, center, .{ .fixes = .{ .q1 = true } });
    try testing.expectEqual(@as(f64, 1), b.count);
    try testing.expectApproxEqAbs(1.0 / 81.0, b.volume, 1e-15);
}

test "volume sum uses relative volumes and the 2.75 Å floor" {
    var conf = sites.Conformation{ .allocator = testing.allocator };
    defer conf.deinit();
    const p = Parameters{};
    // Group atom (residue 1), a CA 2.0 Å away (floored to 2.75), a side-chain
    // carbon 4 Å away, an oxygen in the group's own residue number (skipped)
    // and a nitrogen 16 Å away (volume only).
    const atoms = [_]sites.Atom{
        .{ .pos = .{ .x = 0, .y = 0, .z = 0 }, .element = .o, .residue_idx = 0, .chain_key = 0, .res_num = 1, .order = 0 },
        .{ .pos = .{ .x = 2, .y = 0, .z = 0 }, .element = .c, .name = .{ 'C', 'A', ' ', ' ' }, .name_len = 2, .residue_idx = 1, .chain_key = 0, .res_num = 2, .order = 1 },
        .{ .pos = .{ .x = 0, .y = 4, .z = 0 }, .element = .c, .name = .{ 'C', 'B', ' ', ' ' }, .name_len = 2, .residue_idx = 1, .chain_key = 0, .res_num = 2, .order = 2 },
        .{ .pos = .{ .x = 0, .y = 0, .z = 1.5 }, .element = .o, .residue_idx = 0, .chain_key = 0, .res_num = 1, .order = 3 },
        .{ .pos = .{ .x = 16, .y = 0, .z = 0 }, .element = .n, .residue_idx = 2, .chain_key = 0, .res_num = 3, .order = 4 },
    };
    try conf.atoms.appendSlice(testing.allocator, &atoms);
    const b = burial(&conf, 0, .{ .x = 0, .y = 0, .z = 0 }, p);
    try testing.expectEqual(@as(f64, 2), b.count);
    const expected = 1.40 / std.math.pow(f64, 2.75, 4) + 2.64 / 256.0 + 1.06 / std.math.pow(f64, 16, 4);
    try testing.expectApproxEqAbs(expected, b.volume, 1e-15);
}
