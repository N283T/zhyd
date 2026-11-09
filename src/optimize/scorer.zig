//! Dot-sphere bump/H-bond scoring for probe-based contact evaluation.
//! Reimplementation of the original reduce scoring algorithm — see
//! Word et al. (1999) and reduce's `probe.cpp` / `ScoreContact` for the
//! reference. No source code is copied; this is an independent Zig implementation.

const std = @import("std");
const math_mod = @import("../math.zig");
const Vec3 = math_mod.Vec3;
const element = @import("../element.zig");
const dot_sphere = @import("dot_sphere.zig");
const DotSphere = dot_sphere.DotSphere;

pub const ScoringParams = struct {
    gap_scale: f32 = 0.25, // Contact score gaussian width
    bump_weight: f32 = 10.0, // Clash penalty multiplier
    hb_weight: f32 = 4.0, // H-bond reward multiplier
    min_reg_hb_gap: f32 = 0.6, // Regular H-bond gap threshold (Å)
    min_charged_hb_gap: f32 = 0.8, // Charged H-bond gap threshold (Å)
    bad_bump_gap_cut: f32 = 0.4, // Bad bump classification threshold (Å)
    dot_density: f32 = 16.0, // Dots per Å²
    probe_radius: f32 = 0.25, // Probe sphere radius (reduce default)
};

pub const ScoreResult = struct {
    total: f32 = 0.0,
    contact_sub: f32 = 0.0,
    bump_sub: f32 = 0.0,
    hb_sub: f32 = 0.0,
    n_bad_bumps: u32 = 0,
};

/// H-bond overlap allowance for a pair, or null when the pair cannot H-bond.
/// One atom must donate and the other accept. When both are charged they must
/// have opposite charges, and the allowance is the charged one.
pub fn hbCutoff(flags_a: element.AtomFlags, flags_b: element.AtomFlags, params: ScoringParams) ?f32 {
    if (!((flags_a.donor and flags_b.acceptor) or (flags_a.acceptor and flags_b.donor))) return null;
    const a_charged = flags_a.positive or flags_a.negative;
    const b_charged = flags_b.positive or flags_b.negative;
    if (a_charged and b_charged) {
        const complementary = (flags_a.positive and flags_b.negative) or (flags_a.negative and flags_b.positive);
        if (!complementary) return null;
        return params.min_charged_hb_gap;
    }
    return params.min_reg_hb_gap;
}

/// Check if an overlap of `gap` (negative) is an H-bond within its allowance.
pub fn isHBond(flags_a: element.AtomFlags, flags_b: element.AtomFlags, gap: f32, params: ScoringParams) bool {
    const cut = hbCutoff(flags_a, flags_b, params) orelse return false;
    return (-gap) <= cut;
}

/// Score of one dot on the surface of the scored atom, `gap` from the surface
/// of the atom it is attributed to. Follows reduce's dot classification:
/// contacts score a gaussian of the gap, H-bonds reward the overlap, and
/// bumps penalize it. An H-bond overlap beyond its allowance is a bump of the
/// excess only.
fn dotScore(gap: f32, hb_cut: ?f32, params: ScoringParams) f32 {
    if (gap > 0) {
        const r = gap / params.gap_scale;
        return @exp(-r * r);
    }
    if (hb_cut) |c| {
        if (gap >= -c) return params.hb_weight * (-0.5 * gap);
        return -params.bump_weight * (-0.5 * (gap + c));
    }
    return -params.bump_weight * (-0.5 * gap);
}

/// Score a single atom against its environment by placing dots on its van
/// der Waals surface, as reduce does: each dot whose probe (at the probe
/// radius beyond the surface) touches a neighbor is attributed to the
/// neighbor with the smallest gap and scored with `dotScore`. The total is
/// divided by the dot density, so it is in square Angstroms.
///
/// The optimizer uses `pairDotScore`, which integrates the same dot scores
/// analytically for a single neighbor; this function is the discrete
/// reference it is tested against.
pub fn scoreAtom(
    atom_pos: Vec3(f32),
    atom_radius: f32,
    atom_flags: element.AtomFlags,
    neighbors: []const Vec3(f32),
    neighbor_radii: []const f32,
    neighbor_flags: []const element.AtomFlags,
    sphere: *const DotSphere,
    params: ScoringParams,
) ScoreResult {
    var result = ScoreResult{};
    const p = params.probe_radius;
    const probe_scale = (atom_radius + p) / atom_radius;

    for (sphere.points) |dot| {
        const q = atom_pos.add(dot);
        const probe_pos = atom_pos.add(dot.scale(probe_scale));

        var min_gap: f32 = std.math.floatMax(f32);
        var cause: ?usize = null;
        for (neighbors, 0..) |npos, i| {
            const reach = neighbor_radii[i] + p;
            const dp = probe_pos.sub(npos);
            if (dp.dot(dp) > reach * reach) continue;
            const gap = q.distance(npos) - neighbor_radii[i];
            if (gap < min_gap) {
                min_gap = gap;
                cause = i;
            }
        }
        const ci = cause orelse continue;
        const hb_cut = hbCutoff(atom_flags, neighbor_flags[ci], params);
        if (hb_cut == null and neighbor_flags[ci].hb_only_dummy) continue;

        const score = dotScore(min_gap, hb_cut, params) / sphere.density;
        result.total += score;
        if (min_gap > 0) {
            result.contact_sub += score;
        } else if (hb_cut != null and min_gap >= -hb_cut.?) {
            result.hb_sub += score;
        } else {
            result.bump_sub += score;
            if (-min_gap >= params.bad_bump_gap_cut) result.n_bad_bumps += 1;
        }
    }
    return result;
}

/// Error function (Abramowitz and Stegun 7.1.26, absolute error < 1.5e-7).
fn erf(x: f32) f32 {
    const sign: f32 = if (x < 0) -1 else 1;
    const ax = @abs(x);
    const t = 1.0 / (1.0 + 0.3275911 * ax);
    const poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    return sign * (1.0 - poly * @exp(-ax * ax));
}

/// Dot score of an atom of radius `ra` against a single neighbor of radius
/// `rb` at center distance `d`: the integral of `dotScore` over the dots of
/// `scoreAtom`, in closed form.
///
/// Dots are uniform on the sphere, so with u the cosine of a dot's angle to
/// the neighbor, the dot-to-neighbor distance is t = sqrt(ra^2 + d^2 - 2 ra d u)
/// and du = t dt / (ra d). A dot counts when its probe touches the neighbor
/// (u above a threshold), and its gap is g = t - rb. The area integral is
///   2 pi ra^2 * integral over u of dotScore = (2 pi ra / d) * integral of dotScore(g) (g + rb) dg.
/// The H-bond and bump pieces are linear in g and the contact piece is a
/// gaussian, so each piece integrates exactly.
///
/// With several neighbors reduce gives each dot to the closest one; summing
/// pair scores counts a dot once per neighbor it touches instead. That keeps
/// the score a sum over atom pairs, which exact clique search relies on.
pub fn pairDotScore(d: f32, ra: f32, rb: f32, hb_cut: ?f32, b_is_dummy: bool, params: ScoringParams) f32 {
    const p = params.probe_radius;
    if (d <= 1e-4) return 0;
    const rp = ra + p;
    const rbp = rb + p;
    if (d >= rp + rbp) return 0;
    // Dots with u above u_min have their probe within rb + p of the neighbor.
    const u_min = @max(-1.0, (rp * rp + d * d - rbp * rbp) / (2 * rp * d));
    if (u_min >= 1.0) return 0;
    const g_lo = @abs(d - ra) - rb;
    const g_hi = @sqrt(@max(0.0, ra * ra + d * d - 2 * ra * d * u_min)) - rb;

    var integral: f32 = 0;
    // Contact (g > 0). A dummy neighbor only counts as an H-bond partner.
    if (g_hi > 0 and !(b_is_dummy and hb_cut == null)) {
        const lo = @max(g_lo, 0.0);
        integral += gaussianPiece(g_hi, rb, params.gap_scale) - gaussianPiece(lo, rb, params.gap_scale);
    }
    // Overlap (g < 0).
    if (g_lo < 0) {
        const top = @min(g_hi, 0.0);
        if (hb_cut) |c| {
            // H-bond: hb_weight * -g/2 for -c <= g <= 0.
            const hb_lo = @max(g_lo, -c);
            if (top > hb_lo) integral += linearPiece(-0.5 * params.hb_weight, 0, hb_lo, top, rb);
            // Too close: bump of the excess, bump_weight * (g + c) / 2.
            const bump_hi = @min(top, -c);
            if (bump_hi > g_lo) integral += linearPiece(0.5 * params.bump_weight, 0.5 * params.bump_weight * c, g_lo, bump_hi, rb);
        } else if (!b_is_dummy) {
            integral += linearPiece(0.5 * params.bump_weight, 0, g_lo, top, rb);
        }
    }
    return 2 * std.math.pi * ra / d * integral;
}

/// Integral of (alpha g + beta)(g + rb) over [g1, g2].
fn linearPiece(alpha: f32, beta: f32, g1: f32, g2: f32, rb: f32) f32 {
    const F = struct {
        fn at(al: f32, be: f32, g: f32, r: f32) f32 {
            return al * g * g * g / 3 + (al * r + be) * g * g / 2 + be * r * g;
        }
    };
    return F.at(alpha, beta, g2, rb) - F.at(alpha, beta, g1, rb);
}

/// Antiderivative of exp(-(g/s)^2) (g + rb).
fn gaussianPiece(g: f32, rb: f32, s: f32) f32 {
    const x = g / s;
    return -0.5 * s * s * @exp(-x * x) + rb * s * (@sqrt(std.math.pi) / 2.0) * erf(x);
}

// ── Tests ─────────────────────────────────────────────────────────────────────

const testing = std.testing;

test "score no overlap gives positive contact" {
    // Two atoms at VDW contact distance should give positive score
    // Atom at origin with r=1.70, neighbor at distance 3.10 (slight contact)
    var sphere = try DotSphere.generate(testing.allocator, 1.70, 16.0);
    defer sphere.deinit();

    const atom_pos = Vec3(f32){ .x = 0.0, .y = 0.0, .z = 0.0 };
    const neighbor_pos = [_]Vec3(f32){.{ .x = 3.50, .y = 0.0, .z = 0.0 }}; // slight contact
    const neighbor_radii = [_]f32{1.40};
    const neighbor_flags = [_]element.AtomFlags{.{}};

    const result = scoreAtom(atom_pos, 1.70, .{}, &neighbor_pos, &neighbor_radii, &neighbor_flags, &sphere, .{});
    try testing.expect(result.total > 0);
    try testing.expect(result.bump_sub == 0);
}

test "score clash gives negative bump" {
    // Overlapping atoms should produce negative bump score
    var sphere = try DotSphere.generate(testing.allocator, 1.70, 16.0);
    defer sphere.deinit();

    const atom_pos = Vec3(f32){ .x = 0.0, .y = 0.0, .z = 0.0 };
    const neighbor_pos = [_]Vec3(f32){.{ .x = 2.0, .y = 0.0, .z = 0.0 }}; // heavy overlap
    const neighbor_radii = [_]f32{1.70};
    const neighbor_flags = [_]element.AtomFlags{.{}};

    const result = scoreAtom(atom_pos, 1.70, .{}, &neighbor_pos, &neighbor_radii, &neighbor_flags, &sphere, .{});
    try testing.expect(result.total < 0);
    try testing.expect(result.bump_sub < 0);
}

test "score H-bond overlap gives positive" {
    // Donor-acceptor overlap within H-bond threshold should be positive
    var sphere = try DotSphere.generate(testing.allocator, 1.05, 16.0);
    defer sphere.deinit();

    const atom_pos = Vec3(f32){ .x = 0.0, .y = 0.0, .z = 0.0 };
    const neighbor_pos = [_]Vec3(f32){.{ .x = 2.15, .y = 0.0, .z = 0.0 }}; // slight overlap, H-bond range
    const neighbor_radii = [_]f32{1.40};
    const donor_flags = element.AtomFlags{ .donor = true };
    const acc_flags = [_]element.AtomFlags{element.AtomFlags{ .acceptor = true }};

    const result = scoreAtom(atom_pos, 1.05, donor_flags, &neighbor_pos, &neighbor_radii, &acc_flags, &sphere, .{});
    try testing.expect(result.hb_sub > 0);
}

test "pairDotScore matches dot-sphere scoring of a single neighbor" {
    // Dense dots make the discrete sum close to the integral.
    const params = ScoringParams{};
    const Case = struct { ra: f32, rb: f32, a: element.AtomFlags, b: element.AtomFlags };
    const cases = [_]Case{
        .{ .ra = 1.05, .rb = 1.40, .a = .{ .donor = true }, .b = .{ .acceptor = true } }, // H...O H-bond
        .{ .ra = 1.05, .rb = 1.40, .a = .{ .donor = true, .positive = true }, .b = .{ .acceptor = true, .negative = true } }, // salt bridge
        .{ .ra = 1.22, .rb = 1.70, .a = .{}, .b = .{} }, // H...C contact/bump
        .{ .ra = 1.05, .rb = 1.55, .a = .{ .donor = true, .positive = true }, .b = .{ .donor = true, .positive = true } }, // like charges: no H-bond
    };
    for (cases) |c| {
        var sphere = try DotSphere.generate(testing.allocator, c.ra, 400.0);
        defer sphere.deinit();
        var d: f32 = 1.2;
        while (d < c.ra + c.rb + 0.6) : (d += 0.05) {
            const npos = [_]Vec3(f32){.{ .x = d, .y = 0, .z = 0 }};
            const discrete = scoreAtom(.{ .x = 0, .y = 0, .z = 0 }, c.ra, c.a, &npos, &.{c.rb}, &.{c.b}, &sphere, params).total;
            const analytic = pairDotScore(d, c.ra, c.rb, hbCutoff(c.a, c.b, params), false, params);
            const tol = 0.03 * @abs(discrete) + 0.03;
            if (@abs(discrete - analytic) > tol) {
                std.debug.print("ra={d} rb={d} d={d}: discrete={d} analytic={d}\n", .{ c.ra, c.rb, d, discrete, analytic });
                return error.TestUnexpectedResult;
            }
        }
    }
}

test "too-close H-bond only penalizes the excess overlap" {
    const params = ScoringParams{};
    const hb = hbCutoff(.{ .donor = true }, .{ .acceptor = true }, params);
    // H...O at 1.75 A: center overlap 0.70 A exceeds the 0.6 A allowance, but
    // most of the H surface still H-bonds, so the pair scores positive.
    try testing.expect(pairDotScore(1.75, 1.05, 1.40, hb, false, params) > 0);
    // Without H-bond flags the same overlap is a clash.
    try testing.expect(pairDotScore(1.75, 1.05, 1.40, null, false, params) < 0);
}

test "hbCutoff requires complementary charges" {
    const params = ScoringParams{};
    try testing.expectEqual(@as(?f32, params.min_reg_hb_gap), hbCutoff(.{ .donor = true }, .{ .acceptor = true }, params));
    try testing.expectEqual(@as(?f32, params.min_charged_hb_gap), hbCutoff(.{ .donor = true, .positive = true }, .{ .acceptor = true, .negative = true }, params));
    try testing.expectEqual(@as(?f32, null), hbCutoff(.{ .donor = true, .positive = true }, .{ .acceptor = true, .positive = true }, params));
    try testing.expectEqual(@as(?f32, null), hbCutoff(.{ .donor = true }, .{ .donor = true }, params));
}
