//! Q10: coupled microstates for the side-chain pairs of the iterative scheme
//! (docs/specs/zhyd-pka-empirical-engine.md, "Q10 coupled microstates").
//!
//! The iterative scheme gives each pair's interaction to one partner (or, with
//! Q9, splits an acid–base pair by its ion-pair population) and iterates
//! pKa values. Here the groups joined by iterative pairs are split into
//! connected components, and every protonation microstate of a component is
//! weighted by its free energy at each pH:
//!
//!   G(x, pH) / (ln 10 kT) = Σᵢ xᵢ (pH − p⁰ᵢ) + Σ_pairs E(x)
//!
//! with xᵢ = 1 when group i is protonated, p⁰ᵢ its pKa without the pairs, and
//! for a pair with H-bond h and Coulomb c (both ≥ 0):
//!   opposite charges, both charged:  −(c + h)
//!   like charges, both charged:      +c
//!   like charges, exactly one charged: −h
//! (an acid is charged when deprotonated, a base when protonated). A single
//! pair reproduces the iterative scheme's full-interaction shifts. A group's
//! pKa is the pH at which its protonated fraction is 0.5; its shift from p⁰
//! is added to the group's `coupling` term. Components larger than
//! `max_exact` fall back to the iterative scheme.

const std = @import("std");
const params_mod = @import("params.zig");
const sites = @import("sites.zig");
const iterative = @import("iterative.zig");

const Parameters = params_mod.Parameters;
const Group = sites.Group;
const Pair = iterative.Pair;

pub const max_exact: usize = 12;
const ph_min: f64 = -10.0;
const ph_max: f64 = 30.0;
const ph_step: f64 = 0.25;

pub const Stats = struct {
    exact_components: u32 = 0,
    fallback_components: u32 = 0,
};

pub fn solve(allocator: std.mem.Allocator, groups: []Group, pairs: []Pair, params: Parameters) !Stats {
    var stats: Stats = .{};
    if (pairs.len == 0) return stats;

    // Union-find over the groups that appear in pairs.
    const parent = try allocator.alloc(u32, groups.len);
    defer allocator.free(parent);
    for (parent, 0..) |*p, i| p.* = @intCast(i);
    for (pairs) |pair| unite(parent, pair.g1, pair.g2);

    const done = try allocator.alloc(bool, groups.len);
    defer allocator.free(done);
    @memset(done, false);

    var members: std.ArrayListUnmanaged(u32) = .empty;
    defer members.deinit(allocator);
    var comp_pairs: std.ArrayListUnmanaged(Pair) = .empty;
    defer comp_pairs.deinit(allocator);

    for (pairs) |first| {
        const root = find(parent, first.g1);
        if (done[root]) continue;
        done[root] = true;
        members.clearRetainingCapacity();
        comp_pairs.clearRetainingCapacity();
        for (groups, 0..) |_, gi| {
            if (find(parent, @intCast(gi)) == root and inPairs(pairs, @intCast(gi))) try members.append(allocator, @intCast(gi));
        }
        for (pairs) |pair| {
            if (find(parent, pair.g1) == root) try comp_pairs.append(allocator, pair);
        }
        if (members.items.len > max_exact) {
            stats.fallback_components += 1;
            try iterative.solve(allocator, groups, comp_pairs.items, params);
            continue;
        }
        stats.exact_components += 1;
        try solveExact(allocator, groups, members.items, comp_pairs.items);
    }
    return stats;
}

fn inPairs(pairs: []const Pair, g: u32) bool {
    for (pairs) |pair| {
        if (pair.g1 == g or pair.g2 == g) return true;
    }
    return false;
}

fn find(parent: []u32, x: u32) u32 {
    var r = x;
    while (parent[r] != r) r = parent[r];
    var y = x;
    while (parent[y] != r) {
        const next = parent[y];
        parent[y] = r;
        y = next;
    }
    return r;
}

fn unite(parent: []u32, a: u32, b: u32) void {
    const ra = find(parent, a);
    const rb = find(parent, b);
    if (ra != rb) parent[ra] = rb;
}

/// Pair energy (pKa units) of a microstate given the two groups' protonation.
pub fn pairEnergy(q1: f64, q2: f64, h: f64, c: f64, prot1: bool, prot2: bool) f64 {
    const charged1 = if (q1 < 0) !prot1 else prot1;
    const charged2 = if (q2 < 0) !prot2 else prot2;
    if (q1 * q2 < 0) return if (charged1 and charged2) -(c + h) else 0;
    if (charged1 and charged2) return c;
    if (charged1 != charged2) return -h;
    return 0;
}

fn solveExact(allocator: std.mem.Allocator, groups: []Group, members: []const u32, pairs: []const Pair) !void {
    const n = members.len;
    const n_states: usize = @as(usize, 1) << @intCast(n);
    // G(x, pH) = a[x] + pH · nprot[x].
    const a = try allocator.alloc(f64, n_states);
    defer allocator.free(a);
    const nprot = try allocator.alloc(f64, n_states);
    defer allocator.free(nprot);

    var p0: [max_exact]f64 = undefined;
    var q: [max_exact]f64 = undefined;
    for (members, 0..) |g, k| {
        p0[k] = groups[g].total();
        q[k] = groups[g].q;
    }
    var local: [max_exact * max_exact]struct { i: u8, j: u8, h: f64, c: f64 } = undefined;
    var n_local: usize = 0;
    for (pairs) |pair| {
        const i = std.mem.indexOfScalar(u32, members, pair.g1).?;
        const j = std.mem.indexOfScalar(u32, members, pair.g2).?;
        const h = if (pair.h > 0) pair.h else 0;
        const c = if (pair.c > 0) pair.c else 0;
        local[n_local] = .{ .i = @intCast(i), .j = @intCast(j), .h = h, .c = c };
        n_local += 1;
    }

    for (0..n_states) |x| {
        var sum: f64 = 0;
        var count: f64 = 0;
        for (0..n) |k| {
            if (x & (@as(usize, 1) << @intCast(k)) != 0) {
                sum -= p0[k];
                count += 1;
            }
        }
        for (local[0..n_local]) |lp| {
            const pi = x & (@as(usize, 1) << @intCast(lp.i)) != 0;
            const pj = x & (@as(usize, 1) << @intCast(lp.j)) != 0;
            sum += pairEnergy(q[lp.i], q[lp.j], lp.h, lp.c, pi, pj);
        }
        a[x] = sum;
        nprot[x] = count;
    }

    // A member's pKa is the last downward crossing of its protonated
    // fraction through 0.5: located on a coarse pH grid, then bisected.
    var pka: [max_exact]f64 = undefined;
    var prev: [max_exact]f64 = undefined;
    var crossing: [max_exact]?f64 = @splat(null);
    var ph = ph_min;
    var first = true;
    while (ph <= ph_max + 1e-9) : (ph += ph_step) {
        var frac: [max_exact]f64 = undefined;
        fractions(a, nprot, n, ph, &frac);
        for (0..n) |k| {
            if (!first and prev[k] >= 0.5 and frac[k] < 0.5) crossing[k] = ph - ph_step;
            prev[k] = frac[k];
        }
        first = false;
    }
    for (0..n) |k| {
        const lo_start = crossing[k] orelse {
            pka[k] = if (prev[k] >= 0.5) ph_max else ph_min;
            continue;
        };
        var lo = lo_start;
        var hi = lo_start + ph_step;
        for (0..30) |_| {
            const mid = 0.5 * (lo + hi);
            var frac: [max_exact]f64 = undefined;
            fractions(a, nprot, n, mid, &frac);
            if (frac[k] >= 0.5) lo = mid else hi = mid;
        }
        pka[k] = 0.5 * (lo + hi);
    }

    for (members, 0..) |g, k| groups[g].coupling += pka[k] - p0[k];
}

/// Protonated fraction of each of the `n` members at `ph`.
fn fractions(a: []const f64, nprot: []const f64, n: usize, ph: f64, out: *[max_exact]f64) void {
    var gmin: f64 = std.math.inf(f64);
    for (a, nprot) |ax, nx| gmin = @min(gmin, ax + ph * nx);
    var z: f64 = 0;
    var acc: [max_exact]f64 = @splat(0);
    for (a, nprot, 0..) |ax, nx, x| {
        const w = std.math.pow(f64, 10.0, -(ax + ph * nx - gmin));
        z += w;
        for (0..n) |k| {
            if (x & (@as(usize, 1) << @intCast(k)) != 0) acc[k] += w;
        }
    }
    for (0..n) |k| out[k] = acc[k] / z;
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

test "a well separated acid-base pair shifts both by the full interaction" {
    // ASP 3.8 and LYS 10.5, h + c = 1.0: the ion pair dominates.
    var groups = [_]Group{ testGroup(.asp, 3.8), testGroup(.lys, 10.5) };
    var pairs = [_]Pair{.{ .g1 = 1, .g2 = 0, .h = 0.6, .c = 0.4 }};
    const stats = try solve(testing.allocator, &groups, &pairs, .{});
    try testing.expectEqual(@as(u32, 1), stats.exact_components);
    try testing.expectApproxEqAbs(@as(f64, 2.8), groups[0].total(), 0.02);
    try testing.expectApproxEqAbs(@as(f64, 11.5), groups[1].total(), 0.02);
}

test "two acids: the higher pKa gains h + c, the other loses h" {
    var groups = [_]Group{ testGroup(.asp, 3.0), testGroup(.glu, 7.0) };
    var pairs = [_]Pair{.{ .g1 = 1, .g2 = 0, .h = 0.8, .c = 0.5 }};
    _ = try solve(testing.allocator, &groups, &pairs, .{});
    try testing.expectApproxEqAbs(@as(f64, 8.3), groups[1].total(), 0.02);
    try testing.expectApproxEqAbs(@as(f64, 2.2), groups[0].total(), 0.02);
}

test "an acid and a base with crossing pKa values share the interaction" {
    // HIS 6.5 and ASP 6.5 with W = 1: neither pair state dominates.
    var groups = [_]Group{ testGroup(.asp, 6.5), testGroup(.his, 6.5) };
    var pairs = [_]Pair{.{ .g1 = 1, .g2 = 0, .h = 0.5, .c = 0.5 }};
    _ = try solve(testing.allocator, &groups, &pairs, .{});
    try testing.expect(groups[1].total() > 6.5 and groups[1].total() < 7.5);
    try testing.expect(groups[0].total() < 6.5 and groups[0].total() > 5.5);
}

test "pair energies follow the charge rules" {
    try testing.expectEqual(@as(f64, -1.0), pairEnergy(-1, 1, 0.6, 0.4, false, true));
    try testing.expectEqual(@as(f64, 0.0), pairEnergy(-1, 1, 0.6, 0.4, true, true));
    try testing.expectEqual(@as(f64, 0.4), pairEnergy(-1, -1, 0.6, 0.4, false, false));
    try testing.expectEqual(@as(f64, -0.6), pairEnergy(-1, -1, 0.6, 0.4, true, false));
    try testing.expectEqual(@as(f64, 0.4), pairEnergy(1, 1, 0.6, 0.4, true, true));
}
