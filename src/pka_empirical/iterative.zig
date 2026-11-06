//! Iterative assignment of strongly coupled pair interactions
//! ("Iterative scheme").

const std = @import("std");
const params_mod = @import("params.zig");
const sites = @import("sites.zig");

const Parameters = params_mod.Parameters;
const Group = sites.Group;

pub const Pair = struct {
    /// Index of the later group (G₁) and the earlier group (G₂).
    g1: u32,
    g2: u32,
    h: f64,
    c: f64,
    /// Persistent corrections.
    a1: f64 = 0,
    a2: f64 = 0,
    /// Contributions of the latest iteration.
    side_chain1: f64 = 0,
    side_chain2: f64 = 0,
    coulomb1: f64 = 0,
    coulomb2: f64 = 0,
};

/// One iteration of one pair, using the previous-iteration values p1, p2.
pub fn step(pair: *Pair, q1: f64, q2: f64, p1: f64, p2: f64, params: Parameters) void {
    const h = pair.h;
    const c = pair.c;
    pair.side_chain1 = 0;
    pair.side_chain2 = 0;
    pair.coulomb1 = 0;
    pair.coulomb2 = 0;
    if (q1 < 0 and q2 < 0) {
        const d = c + 2.0 * h;
        const first_wins = p1 + pair.a1 > p2 + pair.a2;
        pair.a1 = 0;
        pair.a2 = 0;
        if (first_wins) {
            pair.side_chain1 = h;
            pair.coulomb1 = c;
            pair.side_chain2 = -h;
            pair.a1 = -d;
        } else {
            pair.side_chain2 = h;
            pair.coulomb2 = c;
            pair.side_chain1 = -h;
            pair.a2 = -d;
        }
    } else if (q1 > 0 and q2 > 0) {
        const d = -(c + 2.0 * h);
        const first_wins = p1 + pair.a1 + d < p2 + pair.a2 + d;
        pair.a1 = 0;
        pair.a2 = 0;
        if (first_wins) {
            pair.side_chain1 = -h;
            pair.coulomb1 = -c;
            pair.side_chain2 = h;
            pair.a1 = -d;
        } else {
            pair.side_chain2 = -h;
            pair.coulomb2 = -c;
            pair.side_chain1 = h;
            pair.a2 = -d;
        }
    } else {
        const k1 = p1 + pair.a1 + q1 * (c + h);
        const k2 = p2 + pair.a2 + q2 * (c + h);
        pair.a1 = 0;
        pair.a2 = 0;
        const apply = (q1 == -1 and k1 < k2) or (q1 == 1 and k1 > k2);
        if (apply) {
            if (c > params.min_contribution) {
                pair.coulomb1 = q1 * c;
                pair.coulomb2 = q2 * c;
                pair.a1 -= q1 * c;
                pair.a2 -= q2 * c;
            }
            if (h > params.min_contribution) {
                pair.side_chain1 = q1 * h;
                pair.side_chain2 = q2 * h;
                pair.a1 -= q1 * h;
                pair.a2 -= q2 * h;
            }
        }
    }
}

/// Run the scheme and add the final contributions above the threshold to the
/// groups' side-chain and Coulomb determinants.
pub fn solve(allocator: std.mem.Allocator, groups: []Group, pairs: []Pair, params: Parameters) !void {
    if (pairs.len == 0) return;
    const p0 = try allocator.alloc(f64, groups.len);
    defer allocator.free(p0);
    const p_old = try allocator.alloc(f64, groups.len);
    defer allocator.free(p_old);
    const p_new = try allocator.alloc(f64, groups.len);
    defer allocator.free(p_new);
    for (groups, 0..) |g, i| p0[i] = g.total();
    @memcpy(p_old, p0);

    var iteration: u32 = 0;
    while (iteration < params.max_iterations) : (iteration += 1) {
        @memcpy(p_new, p0);
        for (pairs) |*pair| {
            step(pair, groups[pair.g1].q, groups[pair.g2].q, p_old[pair.g1], p_old[pair.g2], params);
            p_new[pair.g1] += pair.side_chain1 + pair.coulomb1;
            p_new[pair.g2] += pair.side_chain2 + pair.coulomb2;
        }
        const converged = std.mem.eql(f64, p_new, p_old);
        @memcpy(p_old, p_new);
        if (converged) break;
    }

    for (pairs) |pair| {
        if (@abs(pair.side_chain1) > params.min_contribution) groups[pair.g1].side_chain += pair.side_chain1;
        if (@abs(pair.side_chain2) > params.min_contribution) groups[pair.g2].side_chain += pair.side_chain2;
        if (@abs(pair.coulomb1) > params.min_contribution) groups[pair.g1].coulomb += pair.coulomb1;
        if (@abs(pair.coulomb2) > params.min_contribution) groups[pair.g2].coulomb += pair.coulomb2;
    }
}

const testing = std.testing;

fn testGroup(kind: @import("../pka.zig").SiteKind, model: f64) Group {
    return .{
        .gtype = params_mod.groupTypeOf(kind),
        .kind = kind,
        .residue_idx = 0,
        .group_atom = 0,
        .center = .{ .x = 0, .y = 0, .z = 0 },
        .q = params_mod.charge(kind),
        .model_pka = model,
    };
}

test "two acids: the higher pKa wins the proton" {
    // ASP (3.8) earlier, GLU (4.5) later; h = 0.4, c = 0.3.
    var groups = [_]Group{ testGroup(.asp, 3.8), testGroup(.glu, 4.5) };
    var pairs = [_]Pair{.{ .g1 = 1, .g2 = 0, .h = 0.4, .c = 0.3 }};
    try solve(testing.allocator, &groups, &pairs, .{});
    // Iteration 1: 4.5 > 3.8 → GLU +h +c, ASP −h; a₁ = −(c + 2h) = −1.1.
    // Iteration 2: 5.2 − 1.1 = 4.1 > 3.4 → same result, p unchanged → stop.
    try testing.expectApproxEqAbs(@as(f64, 0.4), groups[1].side_chain, 1e-12);
    try testing.expectApproxEqAbs(@as(f64, 0.3), groups[1].coulomb, 1e-12);
    try testing.expectApproxEqAbs(@as(f64, -0.4), groups[0].side_chain, 1e-12);
    try testing.expectEqual(@as(f64, 0), groups[0].coulomb);
}

test "acid-base pair applies q·c and q·h when favourable" {
    // HIS (6.5, later, q = +1) with GLU (4.5, earlier, q = −1); h = 0.8, c = 0.5.
    var groups = [_]Group{ testGroup(.glu, 4.5), testGroup(.his, 6.5) };
    var pairs = [_]Pair{.{ .g1 = 1, .g2 = 0, .h = 0.8, .c = 0.5 }};
    try solve(testing.allocator, &groups, &pairs, .{});
    // k₁ = 6.5 + 1.3 = 7.8 > k₂ = 4.5 − 1.3 = 3.2 → applied every iteration;
    // a is reset and re-derived each time, so the values are stable.
    try testing.expectApproxEqAbs(@as(f64, 0.8), groups[1].side_chain, 1e-12);
    try testing.expectApproxEqAbs(@as(f64, 0.5), groups[1].coulomb, 1e-12);
    try testing.expectApproxEqAbs(@as(f64, -0.8), groups[0].side_chain, 1e-12);
    try testing.expectApproxEqAbs(@as(f64, -0.5), groups[0].coulomb, 1e-12);
}

test "three-group iterative case" {
    // Two carboxylates sharing a histidine: GLU A (4.5, index 0),
    // ASP B (3.8, index 1), HIS C (6.5, index 2).
    // Pairs in order (i ascending, then j): (1,0), (2,0), (2,1).
    var groups = [_]Group{ testGroup(.glu, 4.5), testGroup(.asp, 3.8), testGroup(.his, 6.5) };
    var pairs = [_]Pair{
        .{ .g1 = 1, .g2 = 0, .h = 0.5, .c = 0.2 },
        .{ .g1 = 2, .g2 = 0, .h = 0.6, .c = 0.3 },
        .{ .g1 = 2, .g2 = 1, .h = 0.0, .c = 0.4 },
    };
    try solve(testing.allocator, &groups, &pairs, .{});
    // Iteration 1 (p = 4.5, 3.8, 6.5):
    //  (ASP,GLU) acids: 3.8 > 4.5 false → GLU wins: GLU +0.5 +0.2, ASP −0.5; a_GLU = −1.2.
    //  (HIS,GLU): k₁ = 6.5 + 0.9 = 7.4, k₂ = 4.5 − 0.9 = 3.6 → HIS +0.3 +0.6, GLU −0.3 −0.6.
    //  (HIS,ASP): k₁ = 6.5 + 0.4, k₂ = 3.8 − 0.4 → HIS +0.4, ASP −0.4 (h = 0 skipped).
    //  p = GLU 4.5 + 0.7 − 0.9 = 4.3, ASP 3.8 − 0.9 = 2.9, HIS 6.5 + 1.3 = 7.8.
    // Iteration 2: (ASP,GLU): 2.9 + 0 > 4.3 − 1.2 = 3.1 false → GLU wins again.
    //  The acid–base pairs reset a and give the same result → p unchanged → stop.
    try testing.expectApproxEqAbs(@as(f64, 0.5 - 0.6), groups[0].side_chain, 1e-12);
    try testing.expectApproxEqAbs(@as(f64, 0.2 - 0.3), groups[0].coulomb, 1e-12);
    try testing.expectApproxEqAbs(@as(f64, -0.5), groups[1].side_chain, 1e-12);
    try testing.expectApproxEqAbs(@as(f64, -0.4), groups[1].coulomb, 1e-12);
    try testing.expectApproxEqAbs(@as(f64, 0.6), groups[2].side_chain, 1e-12);
    try testing.expectApproxEqAbs(@as(f64, 0.7), groups[2].coulomb, 1e-12);
}
