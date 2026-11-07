//! Energy profiles and determinants of the empirical engine ("Distance and
//! angle factors", "Backbone hydrogen bonds", "Ions", "Local backbone
//! reorganisation term", "Side-chain pair interactions").

const std = @import("std");
const params_mod = @import("params.zig");
const sites = @import("sites.zig");
const desolvation = @import("desolvation.zig");

const Parameters = params_mod.Parameters;
const Group = sites.Group;
const Conformation = sites.Conformation;
const Vec = sites.Vec;

// ── Profiles ──────────────────────────────────────────────────────────────────

/// f(d): 1 below c₁, 0 above c₂, linear in between.
pub fn hbondDistanceFactor(d: f64, c1: f64, c2: f64) f64 {
    if (d < c1) return 1.0;
    if (d > c2) return 0.0;
    return 1.0 - (d - c1) / (c2 - c1);
}

/// E_hb = |Δ · f(d) · F_angle|.
pub fn hbondEnergy(d: f64, delta: f64, c1: f64, c2: f64, angle: f64) f64 {
    return @abs(delta * hbondDistanceFactor(d, c1, c2) * angle);
}

/// Cosine between the D→H direction and the H→A direction.
pub fn angleFactor(donor: Vec, hydrogen: Vec, acceptor: Vec) f64 {
    const u = hydrogen.sub(donor).normalize();
    const v = acceptor.sub(hydrogen).normalize();
    return u.dot(v);
}

/// Angle factor as used by the side-chain H-bond terms: the raw cosine in
/// parity mode (its sign is lost in E_hb), max(0, cos) with Q2.
pub fn sideChainAngle(cos: f64, params: Parameters) f64 {
    return if (params.fixes.q2) @max(0.0, cos) else cos;
}

/// E_coul = |244.12 / (ε · r′) · t|.
pub fn coulombEnergy(r: f64, w: f64, params: Parameters) f64 {
    const eps = params.eps_water - (params.eps_water - params.eps_protein) * w;
    const rp = @max(r, params.coulomb_min_distance);
    const t = std.math.clamp((rp - params.coulomb_cutoff_far) / (params.coulomb_min_distance - params.coulomb_cutoff_far), 0.0, 1.0);
    return @abs(params.coulomb_prefactor / (eps * rp) * t);
}

pub const ClosestPair = struct { a: u32, b: u32, d: f64 };

/// Minimum distance over the cross product of two atom lists.
pub fn closestPair(conf: *const Conformation, l1: []const u32, l2: []const u32) ?ClosestPair {
    var best: ?ClosestPair = null;
    for (l1) |a| for (l2) |b| {
        const d = conf.atoms.items[a].pos.distance(conf.atoms.items[b].pos);
        if (best == null or d < best.?.d) best = .{ .a = a, .b = b, .d = d };
    };
    return best;
}

// ── Backbone ──────────────────────────────────────────────────────────────────

/// Backbone H-bond determinant of titratable group `t` from backbone group `b`.
pub fn backboneDeterminant(conf: *const Conformation, b: *const Group, t: *const Group, params: Parameters) f64 {
    const bp = params_mod.backboneParams(b.gtype, t.gtype) orelse return 0;
    const cp = closestPair(conf, b.acid_list.slice(), t.acid_list.slice()) orelse return 0;
    if (cp.d >= bp.c2) return 0;
    const atoms = conf.atoms.items;
    var angle: f64 = 1.0;
    if (b.gtype == .bbn) {
        const bat = atoms[cp.a];
        angle = if (bat.isHydrogen()) angleFactor(atoms[bat.parent].pos, bat.pos, atoms[cp.b].pos) else 0.0;
    } else if (b.gtype == .bbc and t.gtype == .his) {
        const tat = atoms[cp.b];
        angle = if (tat.isHydrogen()) angleFactor(atoms[tat.parent].pos, tat.pos, atoms[cp.a].pos) else 0.0;
    }
    if (angle <= params.angle_min) return 0;
    return t.q * hbondEnergy(cp.d, bp.delta, bp.c1, bp.c2, angle);
}

// ── Ions ──────────────────────────────────────────────────────────────────────

/// Coulomb determinant of an ion on titratable group `t`.
pub fn ionDeterminant(t: *const Group, ion: *const Group, params: Parameters) f64 {
    const r = t.center.distance(ion.center);
    if (r >= params.ion_cutoff) return 0;
    const w = desolvation.pairWeight(t.count, ion.count, params);
    return -ion.q * coulombEnergy(r, w, params);
}

// ── Reorganisation ────────────────────────────────────────────────────────────

/// One backbone C=O contribution (before the burial weight) to an ASP/GLU
/// centre: 0.8 · min(1, 1 − (d − 3)/(6 − 3)) when d < 6 and cosθ > 0.001.
pub fn reorganisationTerm(center: Vec, c: Vec, o: Vec, params: Parameters) f64 {
    const d = center.distance(o);
    if (d >= params.reorg_far) return 0;
    const cos = o.sub(c).normalize().dot(center.sub(o).normalize());
    if (cos <= params.angle_min) return 0;
    const f = @min(1.0, 1.0 - (d - params.reorg_near) / (params.reorg_far - params.reorg_near));
    return params.reorg_prefactor * f;
}

// ── Side-chain pairs ──────────────────────────────────────────────────────────

/// Groups treated as buried pair: not (N₁ + N₂ ≤ 900 and (N₁ ≤ 400 or N₂ ≤ 400)).
pub fn isBuriedPair(n1: f64, n2: f64, params: Parameters) bool {
    return !(n1 + n2 <= params.buried_pair_sum and (n1 <= params.buried_pair_single or n2 <= params.buried_pair_single));
}

fn isPair(g1: *const Group, g2: *const Group, a: params_mod.GroupType, b: params_mod.GroupType) bool {
    return (g1.gtype == a and g2.gtype == b) or (g1.gtype == b and g2.gtype == a);
}

/// COO–ARG angle factor for ARG atom `arg_atom` and carboxyl oxygen `o`:
/// for a hydrogen, cos(N→H, H→O); for a nitrogen, the cosine of largest
/// magnitude over its bonded heavy neighbours X of cos(X→N, N→O) (the spec's
/// open point; amendment A9). 1 when the atom has no neighbour. With Q2 the
/// nitrogen uses the largest cosine and every value is clamped at 0.
fn argAngle(conf: *const Conformation, arg_atom: u32, o: u32, params: Parameters) f64 {
    const atoms = conf.atoms.items;
    const atom = atoms[arg_atom];
    if (atom.isHydrogen()) return sideChainAngle(angleFactor(atoms[atom.parent].pos, atom.pos, atoms[o].pos), params);
    var best: ?f64 = null;
    for (conf.neighbours(arg_atom)) |nb| {
        const c = angleFactor(atoms[nb].pos, atom.pos, atoms[o].pos);
        const prev = best orelse {
            best = c;
            continue;
        };
        if (if (params.fixes.q2) c > prev else @abs(c) > @abs(prev)) best = c;
    }
    return sideChainAngle(best orelse 1.0, params);
}

/// COO–ARG: sum over the two closest disjoint atom pairs.
fn carboxylArginine(conf: *const Conformation, coo: *const Group, arg: *const Group, params: Parameters) f64 {
    const cut = params_mod.sideChainCutoffs(.coo, .arg);
    var l1_buf: [8]u32 = undefined;
    var l2_buf: [8]u32 = undefined;
    const l1_src = coo.listAgainst(arg.gtype);
    const l2_src = arg.listAgainst(coo.gtype);
    @memcpy(l1_buf[0..l1_src.len], l1_src);
    @memcpy(l2_buf[0..l2_src.len], l2_src);
    var l1: []u32 = l1_buf[0..l1_src.len];
    var l2: []u32 = l2_buf[0..l2_src.len];
    var total: f64 = 0;
    var k: usize = 0;
    while (k < 2) : (k += 1) {
        const cp = closestPair(conf, l1, l2) orelse break;
        if (cp.d < cut.c2) {
            total += hbondEnergy(cp.d, params.hbond_max, cut.c1, cut.c2, argAngle(conf, cp.b, cp.a, params));
        }
        l1 = removeValue(l1, cp.a);
        l2 = removeValue(l2, cp.b);
    }
    return total;
}

fn removeValue(list: []u32, value: u32) []u32 {
    const i = std.mem.indexOfScalar(u32, list, value) orelse return list;
    list[i] = list[list.len - 1];
    return list[0 .. list.len - 1];
}

/// Side-chain H-bond magnitude h(G₁, G₂).
pub fn sideChainHbond(conf: *const Conformation, g1: *const Group, g2: *const Group, params: Parameters) f64 {
    const cut = params_mod.sideChainCutoffs(g1.gtype, g2.gtype);
    const cp = closestPair(conf, g1.listAgainst(g2.gtype), g2.listAgainst(g1.gtype)) orelse return 0;
    if (cp.d >= cut.c2) return 0;
    if (conf.withinBonds(g1.group_atom, g2.group_atom, params.bond_exclusion)) return 0;

    if (isPair(g1, g2, .coo, .arg)) {
        return if (g1.gtype == .coo) carboxylArginine(conf, g1, g2, params) else carboxylArginine(conf, g2, g1, params);
    }
    const buried = isBuriedPair(g1.count, g2.count, params);
    if ((isPair(g1, g2, .coo, .his) or isPair(g1, g2, .cys, .his)) and buried) {
        return if (g1.gtype == .coo or g2.gtype == .coo) params.coo_his_buried_value else params.cys_his_buried_value;
    }
    if (g1.gtype == .cys and g2.gtype == .cys and buried) return params.cys_cys_buried_value;

    const atoms = conf.atoms.items;
    var angle: f64 = 1.0;
    if (g2.gtype.isAngular()) {
        const a2 = atoms[cp.b];
        angle = if (a2.isHydrogen()) angleFactor(atoms[a2.parent].pos, a2.pos, atoms[cp.a].pos) else 0.0;
    } else if (g1.gtype.isAngular()) {
        const a1 = atoms[cp.a];
        angle = if (a1.isHydrogen()) angleFactor(atoms[a1.parent].pos, a1.pos, atoms[cp.b].pos) else 0.0;
    }
    angle = sideChainAngle(angle, params);
    if (g1.gtype == .coo and g2.gtype == .coo) {
        const w12 = desolvation.pairWeight(g1.count, g2.count, params);
        return hbondEnergy(cp.d, params.hbond_max, cut.c1, cut.c2, 1.0) * (1.0 + w12);
    }
    return hbondEnergy(cp.d, params.hbond_max, cut.c1, cut.c2, angle);
}

/// Coulomb magnitude c(G₁, G₂).
pub fn sideChainCoulomb(g1: *const Group, g2: *const Group, params: Parameters) f64 {
    if (!g1.isTitratable() or !g2.isTitratable()) return 0;
    const r = g1.center.distance(g2.center);
    if (r > params.pair_cutoff) return 0;
    if (g1.count + g2.count < params.coulomb_min_count) return 0;
    return coulombEnergy(r, desolvation.pairWeight(g1.count, g2.count, params), params);
}

pub const Shifts = struct {
    side_chain1: f64 = 0,
    side_chain2: f64 = 0,
    coulomb1: f64 = 0,
    coulomb2: f64 = 0,
};

/// Non-iterative assignment of h and c to G₁ (later group) and G₂.
pub fn nonIterative(g1: *const Group, g2: *const Group, h: f64, c: f64) Shifts {
    var s = Shifts{};
    if (h > 0) {
        if (g1.q == g2.q) {
            if (g1.model_pka < g2.model_pka) {
                s.side_chain1 = -h;
                s.side_chain2 = h;
            } else if (g2.model_pka < g1.model_pka) {
                s.side_chain1 = h;
                s.side_chain2 = -h;
            } else {
                s.side_chain1 = h;
                s.side_chain2 = -h;
            }
        } else {
            s.side_chain1 = g1.q * h;
            s.side_chain2 = g2.q * h;
        }
    }
    if (c > 0) {
        if (g1.q < 0 and g2.q < 0) {
            if (g1.model_pka > g2.model_pka) s.coulomb1 = c else s.coulomb2 = c;
        } else if (g1.q > 0 and g2.q > 0) {
            if (g1.model_pka < g2.model_pka) s.coulomb1 = -c else s.coulomb2 = -c;
        } else {
            s.coulomb1 = g1.q * c;
            s.coulomb2 = g2.q * c;
        }
    }
    return s;
}

const testing = std.testing;

test "H-bond distance profile" {
    try testing.expectEqual(@as(f64, 1), hbondDistanceFactor(1.9, 2.0, 3.0));
    try testing.expectApproxEqAbs(@as(f64, 0.25), hbondDistanceFactor(2.75, 2.0, 3.0), 1e-12);
    try testing.expectEqual(@as(f64, 0), hbondDistanceFactor(3.1, 2.0, 3.0));
    // Negative angle cosine still yields a positive magnitude (quirk Q2).
    try testing.expectApproxEqAbs(@as(f64, 0.85 * 0.25 * 0.5), hbondEnergy(2.75, 0.85, 2.0, 3.0, -0.5), 1e-12);
}

test "angle factor is the D-H / H-A cosine" {
    const d = Vec{ .x = 0, .y = 0, .z = 0 };
    const h = Vec{ .x = 1, .y = 0, .z = 0 };
    try testing.expectApproxEqAbs(@as(f64, 1), angleFactor(d, h, .{ .x = 3, .y = 0, .z = 0 }), 1e-12);
    try testing.expectApproxEqAbs(@as(f64, 0), angleFactor(d, h, .{ .x = 1, .y = 2, .z = 0 }), 1e-12);
    try testing.expectApproxEqAbs(@sqrt(0.5), angleFactor(d, h, .{ .x = 2, .y = 1, .z = 0 }), 1e-12);
}

test "Coulomb profile by hand" {
    const p = Parameters{};
    // r = 6, w = 0: ε = 160, t = (6 − 10)/(4 − 10) = 2/3 → 244.12/(960)·2/3
    try testing.expectApproxEqAbs(244.12 / 960.0 * 2.0 / 3.0, coulombEnergy(6, 0, p), 1e-12);
    // r = 3 → r' = 4, t = 1; w = 0.5: ε = 95 → 244.12/380
    try testing.expectApproxEqAbs(244.12 / 380.0, coulombEnergy(3, 0.5, p), 1e-12);
    try testing.expectEqual(@as(f64, 0), coulombEnergy(10.5, 1, p));
}

test "reorganisation term by hand" {
    const p = Parameters{};
    // C at origin, O on +x at 1.23, centre further along +x: cosθ = 1.
    const c = Vec{ .x = 0, .y = 0, .z = 0 };
    const o = Vec{ .x = 1.23, .y = 0, .z = 0 };
    // d = 4.5 → 0.8 · (1 − 1.5/3) = 0.4
    try testing.expectApproxEqAbs(@as(f64, 0.4), reorganisationTerm(.{ .x = 5.73, .y = 0, .z = 0 }, c, o, p), 1e-12);
    // d = 2.5 → capped at 0.8
    try testing.expectApproxEqAbs(@as(f64, 0.8), reorganisationTerm(.{ .x = 3.73, .y = 0, .z = 0 }, c, o, p), 1e-12);
    // centre behind the oxygen: cosθ < 0
    try testing.expectEqual(@as(f64, 0), reorganisationTerm(.{ .x = -1, .y = 0, .z = 0 }, c, o, p));
}

test "buried pair rule" {
    const p = Parameters{};
    // Sum <= 900 with one count <= 400: not buried.
    try testing.expect(!isBuriedPair(600, 300, p));
    try testing.expect(!isBuriedPair(100, 200, p));
    // Sum <= 900 but both counts > 400: buried.
    try testing.expect(isBuriedPair(450, 450, p));
    // Sum > 900: buried.
    try testing.expect(isBuriedPair(800, 101, p));
}

test "non-iterative assignment by charge and model pKa" {
    const asp = Group{ .gtype = .coo, .kind = .asp, .residue_idx = 0, .group_atom = 0, .center = .{ .x = 0, .y = 0, .z = 0 }, .q = -1, .model_pka = 3.8 };
    const tyr = Group{ .gtype = .tyr, .kind = .tyr, .residue_idx = 1, .group_atom = 1, .center = .{ .x = 0, .y = 0, .z = 0 }, .q = -1, .model_pka = 10.0 };
    const lys = Group{ .gtype = .lys, .kind = .lys, .residue_idx = 2, .group_atom = 2, .center = .{ .x = 0, .y = 0, .z = 0 }, .q = 1, .model_pka = 10.5 };
    const arg = Group{ .gtype = .arg, .kind = .arg, .residue_idx = 3, .group_atom = 3, .center = .{ .x = 0, .y = 0, .z = 0 }, .q = 1, .model_pka = 12.5 };
    const ser = Group{ .gtype = .roh, .residue_idx = 4, .group_atom = 4, .center = .{ .x = 0, .y = 0, .z = 0 } };
    // Two acids: the lower model pKa (ASP) loses the H-bond, TYR gains;
    // Coulomb goes to the higher model pKa.
    const s1 = nonIterative(&tyr, &asp, 0.3, 0.2);
    try testing.expectEqual(@as(f64, 0.3), s1.side_chain1);
    try testing.expectEqual(@as(f64, -0.3), s1.side_chain2);
    try testing.expectEqual(@as(f64, 0.2), s1.coulomb1);
    try testing.expectEqual(@as(f64, 0), s1.coulomb2);
    // Two bases: Coulomb −c to the lower model pKa (LYS).
    const s2 = nonIterative(&arg, &lys, 0, 0.4);
    try testing.expectEqual(@as(f64, -0.4), s2.coulomb2);
    try testing.expectEqual(@as(f64, 0), s2.coulomb1);
    // Acid–base: q_i · x for each.
    const s3 = nonIterative(&lys, &asp, 0.5, 0.1);
    try testing.expectEqual(@as(f64, 0.5), s3.side_chain1);
    try testing.expectEqual(@as(f64, -0.5), s3.side_chain2);
    try testing.expectEqual(@as(f64, 0.1), s3.coulomb1);
    try testing.expectEqual(@as(f64, -0.1), s3.coulomb2);
    // Neutral partner receives nothing.
    const s4 = nonIterative(&ser, &asp, 0.6, 0);
    try testing.expectEqual(@as(f64, 0), s4.side_chain1);
    try testing.expectEqual(@as(f64, -0.6), s4.side_chain2);
}

fn testAtom(x: f64, y: f64, element: sites.Element, res_num: i32) sites.Atom {
    return .{ .pos = .{ .x = x, .y = y, .z = 0 }, .element = element, .residue_idx = @intCast(res_num), .chain_key = 0, .res_num = res_num, .order = @intCast(res_num) };
}

test "backbone NH to carboxylate by hand" {
    const p = Parameters{};
    var conf = Conformation{ .allocator = testing.allocator };
    defer conf.deinit();
    // 0: backbone N, 1: its H along +x, 2: carboxyl O.
    try conf.atoms.append(testing.allocator, testAtom(0, 0, .n, 1));
    var h = testAtom(1.01, 0, .h, 1);
    h.parent = 0;
    try conf.atoms.append(testing.allocator, h);
    try conf.atoms.append(testing.allocator, testAtom(2.51, 0, .o, 5));
    var bbn = Group{ .gtype = .bbn, .residue_idx = 1, .group_atom = 0, .center = conf.atoms.items[0].pos };
    bbn.acid_list.append(0);
    bbn.acid_list.append(1);
    var asp = Group{ .gtype = .coo, .kind = .asp, .residue_idx = 5, .group_atom = 2, .center = conf.atoms.items[2].pos, .q = -1, .model_pka = 3.8 };
    asp.acid_list.append(2);
    // Linear N–H···O, d(H,O) = 1.5 < 2.0: f = 1, F = 1 → q · 0.85.
    try testing.expectApproxEqAbs(@as(f64, -0.85), backboneDeterminant(&conf, &bbn, &asp, p), 1e-12);
    // H→O perpendicular to N→H: F = 0 ≤ 0.001, no determinant.
    conf.atoms.items[2].pos = .{ .x = 1.01, .y = 2.5, .z = 0 };
    try testing.expectEqual(@as(f64, 0), backboneDeterminant(&conf, &bbn, &asp, p));
    // At 60° (cos = 0.5) and d = 2.5: f = 0.5 → −0.85 · 0.5 · 0.5.
    conf.atoms.items[2].pos = .{ .x = 1.01 + 1.25, .y = 2.5 * @sqrt(3.0) / 2.0, .z = 0 };
    try testing.expectApproxEqAbs(@as(f64, -0.2125), backboneDeterminant(&conf, &bbn, &asp, p), 1e-12);
    // A backbone C=O has no parameters for ASP.
    const bbc = Group{ .gtype = .bbc, .residue_idx = 1, .group_atom = 0, .center = conf.atoms.items[0].pos };
    try testing.expectEqual(@as(f64, 0), backboneDeterminant(&conf, &bbc, &asp, p));
}

test "side-chain H-bond COO–LYS by hand" {
    const p = Parameters{};
    var conf = Conformation{ .allocator = testing.allocator };
    defer conf.deinit();
    try conf.atoms.append(testing.allocator, testAtom(0, 0, .o, 1)); // ASP O
    try conf.atoms.append(testing.allocator, testAtom(3.0, 0, .n, 9)); // LYS NZ
    var asp = Group{ .gtype = .coo, .kind = .asp, .residue_idx = 1, .group_atom = 0, .center = conf.atoms.items[0].pos, .q = -1, .model_pka = 3.8 };
    asp.acid_list.append(0);
    asp.base_list.append(0);
    var lys = Group{ .gtype = .lys, .kind = .lys, .residue_idx = 9, .group_atom = 1, .center = conf.atoms.items[1].pos, .q = 1, .model_pka = 10.5 };
    lys.acid_list.append(1);
    lys.base_list.append(1);
    // COO–LYS cut-offs 2.85/3.85; d = 3.0 → f = 0.85; no angular group.
    const h = sideChainHbond(&conf, &lys, &asp, p);
    try testing.expectApproxEqAbs(@as(f64, 0.85 * 0.85), h, 1e-12);
    // Non-iterative acid–base: LYS +h, ASP −h.
    const s = nonIterative(&lys, &asp, h, 0);
    try testing.expectApproxEqAbs(h, s.side_chain1, 1e-12);
    try testing.expectApproxEqAbs(-h, s.side_chain2, 1e-12);
    // Beyond c₂ no interaction.
    conf.atoms.items[1].pos.x = 3.9;
    try testing.expectEqual(@as(f64, 0), sideChainHbond(&conf, &lys, &asp, p));
}

test "COO–COO H-bond is scaled by 1 + pair weight" {
    const p = Parameters{};
    var conf = Conformation{ .allocator = testing.allocator };
    defer conf.deinit();
    try conf.atoms.append(testing.allocator, testAtom(0, 0, .o, 1));
    try conf.atoms.append(testing.allocator, testAtom(3.0, 0, .o, 9));
    var a = Group{ .gtype = .coo, .kind = .asp, .residue_idx = 1, .group_atom = 0, .center = conf.atoms.items[0].pos, .q = -1, .model_pka = 3.8, .count = 400 };
    a.acid_list.append(0);
    a.base_list = a.acid_list;
    var b = Group{ .gtype = .coo, .kind = .glu, .residue_idx = 9, .group_atom = 1, .center = conf.atoms.items[1].pos, .q = -1, .model_pka = 4.5, .count = 440 };
    b.acid_list.append(1);
    b.base_list = b.acid_list;
    // c = 2.5/3.5, d = 3.0 → f = 0.5; w₁₂ = (840 − 560)/560 = 0.5 → 0.85·0.5·1.5.
    try testing.expectApproxEqAbs(@as(f64, 0.6375), sideChainHbond(&conf, &b, &a, p), 1e-12);
    // Coulomb: N₁ + N₂ ≥ 280, r = 3 → r′ = 4, t = 1, ε = 160 − 130·0.5 = 95.
    try testing.expectApproxEqAbs(244.12 / (95.0 * 4.0), sideChainCoulomb(&b, &a, p), 1e-12);
}

test "ion Coulomb determinant by hand" {
    const p = Parameters{};
    const asp = Group{ .gtype = .coo, .kind = .asp, .residue_idx = 1, .group_atom = 0, .center = .{ .x = 0, .y = 0, .z = 0 }, .q = -1, .model_pka = 3.8 };
    const zn = Group{ .gtype = .ion, .residue_idx = 2, .group_atom = 1, .center = .{ .x = 6, .y = 0, .z = 0 }, .q = 2 };
    // w = 0 (counts 0), r = 6: E = 244.12/(160·6)·(2/3); determinant −2·E.
    try testing.expectApproxEqAbs(-2.0 * 244.12 / 960.0 * 2.0 / 3.0, ionDeterminant(&asp, &zn, p), 1e-12);
    const far = Group{ .gtype = .ion, .residue_idx = 2, .group_atom = 1, .center = .{ .x = 10.5, .y = 0, .z = 0 }, .q = 2 };
    try testing.expectEqual(@as(f64, 0), ionDeterminant(&asp, &far, p));
}

test "Q2: a negative angle cosine gives no side-chain H-bond" {
    var conf = Conformation{ .allocator = testing.allocator };
    defer conf.deinit();
    // TRP NE1 at the origin with HE1 along +x; carboxyl O behind the H:
    // H→O = (−0.3, 1.6), |H→O| = 1.6279 < c₁ = 2.0 (COO–TRP), so f = 1 and
    // cos = −0.3 / 1.6279. N–O = 1.7505, so the H is the closest atom.
    try conf.atoms.append(testing.allocator, testAtom(0, 0, .n, 1));
    var h = testAtom(1.01, 0, .h, 1);
    h.parent = 0;
    try conf.atoms.append(testing.allocator, h);
    try conf.atoms.append(testing.allocator, testAtom(0.71, 1.6, .o, 5));
    var trp = Group{ .gtype = .trp, .residue_idx = 1, .group_atom = 0, .center = conf.atoms.items[0].pos };
    trp.acid_list.append(0);
    trp.acid_list.append(1);
    trp.base_list.append(0);
    var asp = Group{ .gtype = .coo, .kind = .asp, .residue_idx = 5, .group_atom = 2, .center = conf.atoms.items[2].pos, .q = -1, .model_pka = 3.8 };
    asp.acid_list.append(2);
    asp.base_list.append(2);
    const cos = 0.3 / @sqrt(0.3 * 0.3 + 1.6 * 1.6);
    // Parity: |0.85 · 1 · (−cos)|.
    try testing.expectApproxEqAbs(0.85 * cos, sideChainHbond(&conf, &asp, &trp, .{}), 1e-12);
    // Q2: F_angle = max(0, cos) = 0.
    const q2 = Parameters{ .fixes = .{ .q2 = true } };
    try testing.expectEqual(@as(f64, 0), sideChainHbond(&conf, &asp, &trp, q2));
    // Mirrored in front of the H (cos > 0): both modes agree.
    conf.atoms.items[2].pos.x = 1.31;
    try testing.expectApproxEqAbs(0.85 * cos, sideChainHbond(&conf, &asp, &trp, q2), 1e-12);
    try testing.expectApproxEqAbs(0.85 * cos, sideChainHbond(&conf, &asp, &trp, .{}), 1e-12);
}
