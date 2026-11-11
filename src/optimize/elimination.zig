//! Exact maximization of a sum of unary and pairwise score tables by
//! variable elimination (max-sum dynamic programming).
//!
//! A clique's score is a sum of one table per mover and one table per
//! interacting mover pair. Eliminating a mover replaces every table that
//! mentions it with one table over its remaining neighbors, holding the best
//! score over the eliminated mover's candidates. The cost depends on how
//! densely the movers interact (the size of the tables this creates), not on
//! the product of all candidate counts, so cliques far beyond brute force
//! are solved exactly. The mover to eliminate next is the one whose new table
//! is smallest.

const std = @import("std");
const Allocator = std.mem.Allocator;

/// A pair table: members `i < j`, stored row-major over (candidate of i,
/// candidate of j) at `offset` in `Problem.pair_table`.
pub const Pair = struct { i: u32, j: u32, offset: usize };

pub const Problem = struct {
    /// Candidate count per variable.
    card: []const u16,
    /// Variable v's terms are `unary[unary_offsets[v]..unary_offsets[v + 1]]`.
    unary_offsets: []const usize,
    unary: []const f32,
    pairs: []const Pair,
    pair_table: []const f32,
};

pub const Error = error{ TooLarge, OutOfMemory };

/// A table over a sorted set of variables, row-major (last variable fastest).
const Factor = struct {
    scope: []u32,
    table: []f32,
};

/// The choice made when a variable was eliminated: its best candidate for
/// each assignment of `scope`.
const Step = struct {
    variable: u32,
    scope: []u32,
    argmax: []u16,
};

/// Find the candidate per variable that maximizes the total score, written to
/// `out` (one entry per variable). Ties go to the lowest candidate index.
/// Returns error.TooLarge when a table would hold more than `table_limit`
/// entries.
pub fn solve(gpa: Allocator, problem: Problem, table_limit: usize, out: []u16) Error!void {
    var arena_state = std.heap.ArenaAllocator.init(gpa);
    defer arena_state.deinit();
    const arena = arena_state.allocator();

    const n = problem.card.len;
    std.debug.assert(out.len == n);

    var factors = std.ArrayListUnmanaged(Factor).empty;
    for (0..n) |v| {
        const scope = try arena.alloc(u32, 1);
        scope[0] = @intCast(v);
        const table = try arena.dupe(f32, problem.unary[problem.unary_offsets[v]..problem.unary_offsets[v + 1]]);
        try factors.append(arena, .{ .scope = scope, .table = table });
    }
    for (problem.pairs) |p| {
        const len = @as(usize, problem.card[p.i]) * problem.card[p.j];
        const table = problem.pair_table[p.offset..][0..len];
        // An all-zero pair table adds nothing and only couples the two movers.
        if (std.mem.allEqual(f32, table, 0)) continue;
        const scope = try arena.alloc(u32, 2);
        scope[0] = p.i;
        scope[1] = p.j;
        try factors.append(arena, .{ .scope = scope, .table = try arena.dupe(f32, table) });
    }

    const eliminated = try arena.alloc(bool, n);
    @memset(eliminated, false);
    var steps = std.ArrayListUnmanaged(Step).empty;
    const in_scope = try arena.alloc(bool, n);
    const assign = try arena.alloc(u16, n);

    for (0..n) |_| {
        // Pick the variable whose elimination creates the smallest table.
        var best_v: u32 = 0;
        var best_size: usize = std.math.maxInt(usize);
        for (0..n) |v| {
            if (eliminated[v]) continue;
            const size = newScopeSize(factors.items, @intCast(v), problem.card, in_scope);
            if (size < best_size) {
                best_size = size;
                best_v = @intCast(v);
            }
        }
        if (best_size > table_limit) return error.TooLarge;
        const v = best_v;

        // Collect the factors that mention v and the scope of the new one.
        @memset(in_scope, false);
        var touching = std.ArrayListUnmanaged(Factor).empty;
        var kept = std.ArrayListUnmanaged(Factor).empty;
        for (factors.items) |f| {
            if (std.mem.indexOfScalar(u32, f.scope, v) != null) {
                try touching.append(arena, f);
                for (f.scope) |u| in_scope[u] = true;
            } else {
                try kept.append(arena, f);
            }
        }
        in_scope[v] = false;
        var scope_list = std.ArrayListUnmanaged(u32).empty;
        for (0..n) |u| if (in_scope[u]) try scope_list.append(arena, @intCast(u));
        const scope = scope_list.items;

        const table = try arena.alloc(f32, best_size);
        const argmax = try arena.alloc(u16, best_size);

        // Enumerate assignments of the new scope (last variable fastest).
        for (scope) |u| assign[u] = 0;
        var idx: usize = 0;
        while (true) : (idx += 1) {
            var best_score: f32 = -std.math.inf(f32);
            var best_c: u16 = 0;
            for (0..problem.card[v]) |c| {
                assign[v] = @intCast(c);
                var total: f32 = 0;
                for (touching.items) |f| total += f.table[tableIndex(f.scope, assign, problem.card)];
                if (total > best_score) {
                    best_score = total;
                    best_c = @intCast(c);
                }
            }
            table[idx] = best_score;
            argmax[idx] = best_c;
            if (!increment(scope, assign, problem.card)) break;
        }

        try steps.append(arena, .{ .variable = v, .scope = scope, .argmax = argmax });
        eliminated[v] = true;
        if (scope.len > 0) try kept.append(arena, .{ .scope = scope, .table = table });
        factors = kept;
    }

    // Assign variables in reverse elimination order: each step's scope holds
    // only variables eliminated later, which are already assigned.
    var k = steps.items.len;
    while (k > 0) {
        k -= 1;
        const s = steps.items[k];
        assign[s.variable] = s.argmax[tableIndex(s.scope, assign, problem.card)];
    }
    @memcpy(out, assign);
}

/// Size of the table that eliminating `v` would create.
fn newScopeSize(factors: []const Factor, v: u32, card: []const u16, mark: []bool) usize {
    @memset(mark, false);
    for (factors) |f| {
        if (std.mem.indexOfScalar(u32, f.scope, v) == null) continue;
        for (f.scope) |u| mark[u] = true;
    }
    mark[v] = false;
    var size: usize = 1;
    for (mark, 0..) |m, u| {
        if (m) size = std.math.mul(usize, size, card[u]) catch return std.math.maxInt(usize);
    }
    return size;
}

/// Row-major index of the current assignment within a table over `scope`.
fn tableIndex(scope: []const u32, assign: []const u16, card: []const u16) usize {
    var idx: usize = 0;
    for (scope) |u| idx = idx * card[u] + assign[u];
    return idx;
}

/// Mixed-radix increment over `scope`, last variable fastest. False on wrap.
fn increment(scope: []const u32, assign: []u16, card: []const u16) bool {
    var k = scope.len;
    while (k > 0) {
        k -= 1;
        const u = scope[k];
        assign[u] += 1;
        if (assign[u] < card[u]) return true;
        assign[u] = 0;
    }
    return false;
}

// ── Tests ─────────────────────────────────────────────────────────────────────

const testing = std.testing;

fn bruteForce(problem: Problem, out: []u16) f32 {
    const n = problem.card.len;
    var cand: [12]u16 = @splat(0);
    var best: f32 = -std.math.inf(f32);
    while (true) {
        const sc = evaluate(problem, cand[0..n]);
        if (sc > best) {
            best = sc;
            @memcpy(out, cand[0..n]);
        }
        var k = n;
        const done = while (k > 0) {
            k -= 1;
            cand[k] += 1;
            if (cand[k] < problem.card[k]) break false;
            cand[k] = 0;
        } else true;
        if (done) break;
    }
    return best;
}

fn evaluate(problem: Problem, cand: []const u16) f32 {
    var sc: f32 = 0;
    for (cand, 0..) |c, v| sc += problem.unary[problem.unary_offsets[v] + c];
    for (problem.pairs) |p| sc += problem.pair_table[p.offset + @as(usize, cand[p.i]) * problem.card[p.j] + cand[p.j]];
    return sc;
}

test "variable elimination matches brute force on random problems" {
    var prng = std.Random.DefaultPrng.init(7);
    const rand = prng.random();
    for (0..50) |_| {
        const n = 3 + rand.uintLessThan(usize, 6);
        var card: [12]u16 = undefined;
        var unary_offsets: [13]usize = undefined;
        unary_offsets[0] = 0;
        for (0..n) |v| {
            card[v] = 1 + rand.uintLessThan(u16, 5);
            unary_offsets[v + 1] = unary_offsets[v] + card[v];
        }
        var unary: [64]f32 = undefined;
        for (unary[0..unary_offsets[n]]) |*x| x.* = rand.float(f32) * 4 - 2;

        var pairs: [66]Pair = undefined;
        var n_pairs: usize = 0;
        var pair_table: [66 * 25]f32 = undefined;
        var len: usize = 0;
        for (0..n) |i| for (i + 1..n) |j| {
            if (rand.float(f32) < 0.5) continue;
            pairs[n_pairs] = .{ .i = @intCast(i), .j = @intCast(j), .offset = len };
            n_pairs += 1;
            const sz = @as(usize, card[i]) * card[j];
            for (pair_table[len..][0..sz]) |*x| x.* = rand.float(f32) * 4 - 2;
            len += sz;
        };
        const problem = Problem{
            .card = card[0..n],
            .unary_offsets = unary_offsets[0 .. n + 1],
            .unary = unary[0..unary_offsets[n]],
            .pairs = pairs[0..n_pairs],
            .pair_table = pair_table[0..len],
        };

        var expected: [12]u16 = undefined;
        const best = bruteForce(problem, expected[0..n]);
        var got: [12]u16 = undefined;
        try solve(testing.allocator, problem, 1 << 20, got[0..n]);
        try testing.expectApproxEqAbs(best, evaluate(problem, got[0..n]), 1e-4);
    }
}

test "solve reports TooLarge when a table exceeds the limit" {
    // Three mutually coupled variables with 4 candidates: eliminating one
    // creates a 16-entry table.
    const card = [_]u16{ 4, 4, 4 };
    const unary_offsets = [_]usize{ 0, 4, 8, 12 };
    const unary: [12]f32 = @splat(0);
    const pairs = [_]Pair{ .{ .i = 0, .j = 1, .offset = 0 }, .{ .i = 0, .j = 2, .offset = 16 }, .{ .i = 1, .j = 2, .offset = 32 } };
    var pair_table: [48]f32 = undefined;
    for (&pair_table, 0..) |*x, k| x.* = @floatFromInt(k % 7);
    const problem = Problem{ .card = &card, .unary_offsets = &unary_offsets, .unary = &unary, .pairs = &pairs, .pair_table = &pair_table };
    var out: [3]u16 = undefined;
    try testing.expectError(error.TooLarge, solve(testing.allocator, problem, 8, &out));
    try solve(testing.allocator, problem, 16, &out);
}
