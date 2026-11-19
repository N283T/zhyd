//! Per-conformation pipeline and conformation averaging ("Pipeline").

const std = @import("std");
const model_mod = @import("../model.zig");
const pka = @import("../pka.zig");
const params_mod = @import("params.zig");
const sites = @import("sites.zig");
const desolvation = @import("desolvation.zig");
const interactions = @import("interactions.zig");
const iterative = @import("iterative.zig");

const Parameters = params_mod.Parameters;
const Group = sites.Group;
const Model = model_mod.Model;

/// Which determinant of a site a partner contributes to.
pub const Term = enum { side_chain, backbone, coulomb };

/// One partner's contribution to a site's pKa (PROPKA's "determinants").
pub const Determinant = struct {
    partner_residue: u32,
    partner_type: sites.GroupType,
    term: Term,
    value: f64,
};

/// A determinant of group `group` in one conformation.
pub const Recorded = struct {
    group: u32,
    determinant: Determinant,
};

const Recorder = std.ArrayListUnmanaged(Recorded);

fn record(allocator: std.mem.Allocator, recorder: ?*Recorder, groups: []const Group, group: usize, partner: usize, term: Term, value: f64) !void {
    const r = recorder orelse return;
    if (value == 0) return;
    try r.append(allocator, .{ .group = @intCast(group), .determinant = .{
        .partner_residue = groups[partner].residue_idx,
        .partner_type = groups[partner].gtype,
        .term = term,
        .value = value,
    } });
}

/// Run steps 3–8 of the pipeline on one conformation. With a recorder, every
/// partner contribution to a side-chain, backbone or Coulomb term is also
/// appended to it.
pub fn evaluate(allocator: std.mem.Allocator, conf: *sites.Conformation, params: Parameters, recorder: ?*Recorder) !void {
    const groups = conf.groups.items;
    desolvation.apply(conf, params);

    // Backbone H-bonds, ions, reorganisation.
    for (groups, 0..) |*t, ti| {
        if (!t.isTitratable()) continue;
        for (groups, 0..) |*other, oi| {
            switch (other.gtype) {
                .bbn, .bbc => {
                    const v = interactions.backboneDeterminant(conf, other, t, params);
                    t.backbone += v;
                    try record(allocator, recorder, groups, ti, oi, .backbone, v);
                },
                .ion => {
                    const v = interactions.ionDeterminant(t, other, params);
                    t.coulomb += v;
                    try record(allocator, recorder, groups, ti, oi, .coulomb, v);
                },
                else => {},
            }
        }
        if (t.kind == .asp or t.kind == .glu) {
            var sum: f64 = 0;
            for (groups) |*b| {
                if (b.gtype != .bbc) continue;
                sum += interactions.reorganisationTerm(t.center, conf.atoms.items[b.group_atom].pos, conf.atoms.items[b.acid_list.items[0]].pos, params);
            }
            t.reorganisation = sum * t.weight;
        }
    }

    markCovalentCoupling(conf);

    // Side-chain pairs: G₁ is the later group (groups are in input order).
    var pairs: std.ArrayListUnmanaged(iterative.Pair) = .empty;
    defer pairs.deinit(allocator);
    for (groups, 0..) |*g1, i| {
        if (!isPairGroup(g1.*)) continue;
        for (groups[0..i], 0..) |*g2, j| {
            if (!isPairGroup(g2.*)) continue;
            const class = params_mod.pairClass(g1.gtype, g2.gtype);
            if (class == .none) continue;
            if (g1.center.distance(g2.center) >= params.pair_cutoff) continue;
            const h = interactions.sideChainHbond(conf, g1, g2, params);
            const c = interactions.sideChainCoulomb(g1, g2, params);
            if (class == .non_iterative) {
                const s = interactions.nonIterative(g1, g2, h, c);
                g1.side_chain += s.side_chain1;
                g2.side_chain += s.side_chain2;
                g1.coulomb += s.coulomb1;
                g2.coulomb += s.coulomb2;
                try record(allocator, recorder, groups, i, j, .side_chain, s.side_chain1);
                try record(allocator, recorder, groups, j, i, .side_chain, s.side_chain2);
                try record(allocator, recorder, groups, i, j, .coulomb, s.coulomb1);
                try record(allocator, recorder, groups, j, i, .coulomb, s.coulomb2);
            } else if (h > 0 or c > 0) {
                try pairs.append(allocator, .{ .g1 = @intCast(i), .g2 = @intCast(j), .h = h, .c = c });
            }
        }
    }
    try iterative.solve(allocator, groups, pairs.items, params);
    if (recorder != null) {
        // The same threshold as `iterative.solve` uses when it adds them.
        for (pairs.items) |pair| {
            const parts = [_]struct { g: u32, p: u32, term: Term, v: f64 }{
                .{ .g = pair.g1, .p = pair.g2, .term = .side_chain, .v = pair.side_chain1 },
                .{ .g = pair.g2, .p = pair.g1, .term = .side_chain, .v = pair.side_chain2 },
                .{ .g = pair.g1, .p = pair.g2, .term = .coulomb, .v = pair.coulomb1 },
                .{ .g = pair.g2, .p = pair.g1, .term = .coulomb, .v = pair.coulomb2 },
            };
            for (parts) |part| {
                if (@abs(part.v) > params.min_contribution) try record(allocator, recorder, groups, part.g, part.p, part.term, part.v);
            }
        }
    }
}

/// Two titratable groups whose group atoms are within three bonds of each
/// other (in proteins: the N-terminus and the ASP/HIS/CYS side chain of the
/// same residue, or side chains joined by a clash shorter than the bond
/// distance): one of them does not titrate and takes no part in pair
/// interactions (amendment A1). Acid with base: the acid. Two bases: the
/// lower pKa so far (model, desolvation, reorganisation, backbone and ion
/// terms). Two acids: the higher pKa so far. Ties: the later group.
pub fn markCovalentCoupling(conf: *sites.Conformation) void {
    const groups = conf.groups.items;
    for (groups, 0..) |*g, i| {
        if (!g.isTitratable()) continue;
        for (groups[i + 1 ..]) |*other| {
            if (!other.isTitratable()) continue;
            // N-termini of one insertion-code family do not couple (A12).
            if (g.kind == .n_term and other.kind == .n_term) continue;
            if (!conf.withinBonds(g.group_atom, other.group_atom, sites.covalent_coupling_bonds)) continue;
            if (couplingDropsFirst(g.*, other.*)) {
                g.coupled_out = true;
                break;
            }
            other.coupled_out = true;
        }
    }
}

/// Whether `first` (earlier) rather than `second` stops titrating.
fn couplingDropsFirst(first: Group, second: Group) bool {
    if (first.q != second.q) return first.q < 0;
    if (first.q > 0) return first.total() < second.total();
    return first.total() > second.total();
}

fn isPairGroup(g: Group) bool {
    if (!g.gtype.isSideChain()) return false;
    return !g.disulfide and !g.coupled_out;
}

/// Reported result of one site, averaged over conformations.
pub const SiteSum = struct {
    residue_idx: u32,
    kind: pka.SiteKind,
    /// Conformations in which the site exists / titrates.
    n: u32 = 0,
    n_titrating: u32 = 0,
    disulfide: bool = false,
    /// Not titrating in the first conformation: reported as not titrating
    /// whatever the other conformations give (amendment A8).
    not_titrating_first: bool = false,
    predicted: f64 = 0,
    burial: f64 = 0,
    hydrogen_bond: f64 = 0,
    electrostatic: f64 = 0,
    weight: f64 = 0,
    center: sites.Vec = .{ .x = 0, .y = 0, .z = 0 },
    /// Conformation-mean determinants, `Result.determinants[det_start..][0..det_len]`
    /// (only from `runDetailed` with determinants requested).
    det_start: u32 = 0,
    det_len: u32 = 0,
};

/// Conformation average of one site identity: (chain, residue number, site
/// kind), insertion code ignored (amendment A10). Every residue that carried
/// the site in some conformation reports the same mean. With Q8 the identity
/// is (residue, site kind).
const Accumulator = struct {
    sum: SiteSum,
    key: Identity,
    residues: [16]u32 = undefined,
    n_residues: usize = 0,
    /// Determinant sums over conformations, one entry per (partner, term).
    dets: std.ArrayListUnmanaged(Determinant) = .empty,

    fn addDeterminant(self: *Accumulator, allocator: std.mem.Allocator, d: Determinant) !void {
        for (self.dets.items) |*e| {
            if (e.partner_residue == d.partner_residue and e.partner_type == d.partner_type and e.term == d.term) {
                e.value += d.value;
                return;
            }
        }
        try self.dets.append(allocator, d);
    }

    fn addResidue(self: *Accumulator, residue_idx: u32) void {
        if (std.mem.indexOfScalar(u32, self.residues[0..self.n_residues], residue_idx) != null) return;
        if (self.n_residues == self.residues.len) return;
        self.residues[self.n_residues] = residue_idx;
        self.n_residues += 1;
    }
};

/// Site identity used to match groups across and within conformations.
const Identity = struct {
    kind: pka.SiteKind,
    /// Parity: chain key and residue number. Q8: residue index (`b` unused).
    a: u32,
    b: i32,

    fn of(conf: *const sites.Conformation, g: Group, params: Parameters) Identity {
        if (params.fixes.q8) return .{ .kind = g.kind.?, .a = g.residue_idx, .b = 0 };
        const atom = conf.atoms.items[g.group_atom];
        return .{ .kind = g.kind.?, .a = atom.chain_key, .b = atom.res_num };
    }

    fn eql(x: Identity, y: Identity) bool {
        return x.kind == y.kind and x.a == y.a and x.b == y.b;
    }
};

fn findAcc(allocator: std.mem.Allocator, accs: *std.ArrayListUnmanaged(Accumulator), key: Identity, residue_idx: u32) !*Accumulator {
    for (accs.items) |*a| {
        if (a.key.eql(key)) return a;
    }
    try accs.append(allocator, .{ .sum = .{ .residue_idx = residue_idx, .kind = key.kind }, .key = key });
    return &accs.items[accs.items.len - 1];
}

/// Evaluate every conformation of `mdl` and average per site.
pub fn run(allocator: std.mem.Allocator, mdl: *const Model, params: Parameters) ![]SiteSum {
    const result = try runDetailed(allocator, mdl, params, false);
    allocator.free(result.determinants);
    return result.sums;
}

pub const Result = struct {
    sums: []SiteSum,
    /// Backing store of every site's determinants (`SiteSum.det_start`).
    determinants: []Determinant,
};

/// `run`, optionally also averaging each site's determinants over
/// conformations with the same denominator as its terms.
pub fn runDetailed(allocator: std.mem.Allocator, mdl: *const Model, params: Parameters, want_determinants: bool) !Result {
    var letters_buf: [64]u8 = undefined;
    const letters = sites.altlocLetters(mdl, &letters_buf);
    var accs: std.ArrayListUnmanaged(Accumulator) = .empty;
    defer {
        for (accs.items) |*a| a.dets.deinit(allocator);
        accs.deinit(allocator);
    }

    const n_conf = @max(letters.len, 1);
    for (0..n_conf) |ci| {
        const altloc: ?u8 = if (letters.len == 0) null else letters[ci];
        const first: ?u8 = if (letters.len == 0) null else letters[0];
        var conf = try sites.build(allocator, mdl, altloc, first, params);
        defer conf.deinit();
        var recorder: Recorder = .empty;
        defer recorder.deinit(allocator);
        try evaluate(allocator, &conf, params, if (want_determinants) &recorder else null);
        // Within one conformation the first group (input order) of each
        // identity supplies the value (amendment A12); every residue carrying
        // the identity reports it.
        var chosen: std.ArrayListUnmanaged(u32) = .empty;
        defer chosen.deinit(allocator);
        for (conf.groups.items, 0..) |g, gi| {
            if (g.kind == null) continue;
            const key = Identity.of(&conf, g, params);
            const slot = for (chosen.items) |*c| {
                if (Identity.of(&conf, conf.groups.items[c.*], params).eql(key)) break c;
            } else null;
            if (slot != null) {
                const acc = try findAcc(allocator, &accs, key, g.residue_idx);
                acc.addResidue(g.residue_idx);
            } else try chosen.append(allocator, @intCast(gi));
        }
        for (chosen.items) |gi| {
            const g = conf.groups.items[gi];
            const acc = try findAcc(allocator, &accs, Identity.of(&conf, g, params), g.residue_idx);
            acc.addResidue(g.residue_idx);
            addConformation(&acc.sum, g, ci == 0, params);
            if (g.disulfide or g.coupled_out) continue;
            for (recorder.items) |r| {
                if (r.group == gi) try acc.addDeterminant(allocator, r.determinant);
            }
        }
    }

    var out: std.ArrayListUnmanaged(SiteSum) = .empty;
    errdefer out.deinit(allocator);
    var store: std.ArrayListUnmanaged(Determinant) = .empty;
    errdefer store.deinit(allocator);
    for (accs.items) |*a| {
        var s = mean(a.sum, @intCast(n_conf), params);
        const inv = inverseDenominator(a.sum, @intCast(n_conf), params);
        s.det_start = @intCast(store.items.len);
        s.det_len = @intCast(a.dets.items.len);
        for (a.dets.items) |d| {
            var scaled = d;
            scaled.value *= inv;
            try store.append(allocator, scaled);
        }
        for (a.residues[0..a.n_residues]) |ri| {
            s.residue_idx = ri;
            try out.append(allocator, s);
        }
    }
    const sums = try out.toOwnedSlice(allocator);
    errdefer allocator.free(sums);
    return .{ .sums = sums, .determinants = try store.toOwnedSlice(allocator) };
}

/// 1 / the number of conformations a site's terms are averaged over.
fn inverseDenominator(sum: SiteSum, n_conf: u32, params: Parameters) f64 {
    const denominator: u32 = if (params.fixes.q7) sum.n_titrating else if (params.fixes.q8) sum.n else n_conf;
    return if (denominator == 0) 0.0 else 1.0 / @as(f64, @floatFromInt(denominator));
}

/// Add one conformation's group to its site sum. A group that does not
/// titrate (disulfide, covalently coupled) adds the fixed 99.99 with zero
/// terms in parity mode (amendment A8) and nothing with Q7.
pub fn addConformation(entry: *SiteSum, g: Group, first_conformation: bool, params: Parameters) void {
    entry.n += 1;
    if (g.disulfide or g.coupled_out) {
        entry.disulfide = entry.disulfide or g.disulfide;
        if (first_conformation) entry.not_titrating_first = true;
        if (!params.fixes.q7) entry.predicted += params_mod.not_titrating_pka;
        return;
    }
    entry.n_titrating += 1;
    entry.predicted += g.total();
    entry.burial += g.desolvation + g.reorganisation;
    entry.hydrogen_bond += g.side_chain + g.backbone;
    entry.electrostatic += g.coulomb;
    entry.weight += g.weight;
    entry.center = entry.center.add(g.center);
}

/// Conformation mean of a site sum. Parity: over all `n_conf`
/// conformations; a site absent from one contributes zero (amendment A10).
/// Q8: over the conformations in which the site exists. Q7: over those in
/// which it titrates (and a site not titrating in the first conformation is
/// still reported when it titrates in another).
pub fn mean(sum: SiteSum, n_conf: u32, params: Parameters) SiteSum {
    var s = sum;
    if (params.fixes.q7) s.not_titrating_first = false;
    const inv = inverseDenominator(sum, n_conf, params);
    s.predicted *= inv;
    s.burial *= inv;
    s.hydrogen_bond *= inv;
    s.electrostatic *= inv;
    s.weight *= inv;
    if (s.n_titrating > 0) s.center = s.center.scale(1.0 / @as(f64, @floatFromInt(s.n_titrating)));
    return s;
}

const testing = std.testing;

test "covalent coupling: acid yields to base, like charges by pKa" {
    const base = Group{ .gtype = .nterm, .kind = .n_term, .residue_idx = 0, .group_atom = 0, .center = .{ .x = 0, .y = 0, .z = 0 }, .q = 1, .model_pka = 8.0 };
    const acid = Group{ .gtype = .coo, .kind = .asp, .residue_idx = 0, .group_atom = 1, .center = .{ .x = 0, .y = 0, .z = 0 }, .q = -1, .model_pka = 3.8 };
    try testing.expect(!couplingDropsFirst(base, acid));
    try testing.expect(couplingDropsFirst(acid, base));
    // Two acids: the higher pKa so far stops titrating.
    const glu = Group{ .gtype = .coo, .kind = .glu, .residue_idx = 1, .group_atom = 2, .center = .{ .x = 0, .y = 0, .z = 0 }, .q = -1, .model_pka = 4.5 };
    try testing.expect(couplingDropsFirst(glu, acid));
    try testing.expect(!couplingDropsFirst(acid, glu));
    // Two bases: the lower pKa so far stops titrating.
    const his = Group{ .gtype = .his, .kind = .his, .residue_idx = 1, .group_atom = 3, .center = .{ .x = 0, .y = 0, .z = 0 }, .q = 1, .model_pka = 6.5 };
    try testing.expect(couplingDropsFirst(his, base));
}

test "N-terminal histidine is covalently coupled to the N-terminus" {
    const mmcif = @import("../mmcif.zig");
    var mdl = try mmcif.parseModel(testing.allocator, @embedFile("../test_data/his.cif"));
    defer mdl.deinit();
    for (mdl.residues.items) |*res| res.entity_type = .polymer;
    const sums = try run(testing.allocator, &mdl, .{});
    defer testing.allocator.free(sums);
    try testing.expectEqual(@as(usize, 2), sums.len);
    for (sums) |s| {
        if (s.kind == .his) {
            // N–CA–CB–CG: three bonds; the base (N-terminus) keeps titrating.
            try testing.expect(s.not_titrating_first);
            try testing.expectApproxEqAbs(params_mod.not_titrating_pka, s.predicted, 1e-9);
        } else {
            try testing.expectEqual(pka.SiteKind.n_term, s.kind);
            try testing.expect(!s.not_titrating_first);
        }
    }
}

test "Q7 and Q8: conformation means by hand" {
    const zero = sites.Vec{ .x = 0, .y = 0, .z = 0 };
    // ASP: model 3.80 + desolvation 0.40 = 4.20 in the first conformation;
    // covalently coupled (not titrating) in the second; absent in the third.
    var titrating = Group{ .gtype = .coo, .kind = .asp, .residue_idx = 0, .group_atom = 0, .center = zero, .q = -1, .model_pka = 3.8, .desolvation = 0.4, .weight = 0.5 };
    var coupled = titrating;
    coupled.coupled_out = true;
    titrating.center = .{ .x = 2, .y = 0, .z = 0 };
    const cases = [_]struct { fixes: params_mod.Fixes, predicted: f64, burial: f64 }{
        // Parity: (4.20 + 99.99 + 0) / 3.
        .{ .fixes = .none, .predicted = (4.2 + 99.99) / 3.0, .burial = 0.4 / 3.0 },
        // Q8: absent conformation dropped: (4.20 + 99.99) / 2.
        .{ .fixes = .{ .q8 = true }, .predicted = (4.2 + 99.99) / 2.0, .burial = 0.2 },
        // Q7: titrating conformation only.
        .{ .fixes = .{ .q7 = true }, .predicted = 4.2, .burial = 0.4 },
        .{ .fixes = .{ .q7 = true, .q8 = true }, .predicted = 4.2, .burial = 0.4 },
    };
    for (cases) |c| {
        const p = Parameters{ .fixes = c.fixes };
        var sum = SiteSum{ .residue_idx = 0, .kind = .asp };
        addConformation(&sum, titrating, true, p);
        addConformation(&sum, coupled, false, p);
        const m = mean(sum, 3, p);
        try testing.expectApproxEqAbs(c.predicted, m.predicted, 1e-12);
        try testing.expectApproxEqAbs(c.burial, m.burial, 1e-12);
        // The centre is always the mean over titrating conformations.
        try testing.expectApproxEqAbs(@as(f64, 2), m.center.x, 1e-12);
    }
    // Not titrating in the first conformation: parity reports the site as
    // not titrating; Q7 reports the mean of the titrating conformations.
    for ([_]params_mod.Fixes{ .none, .{ .q7 = true } }) |f| {
        const p = Parameters{ .fixes = f };
        var sum = SiteSum{ .residue_idx = 0, .kind = .asp };
        addConformation(&sum, coupled, true, p);
        addConformation(&sum, titrating, false, p);
        const m = mean(sum, 2, p);
        try testing.expectEqual(!f.q7, m.not_titrating_first);
        try testing.expectApproxEqAbs(@as(f64, if (f.q7) 4.2 else (4.2 + 99.99) / 2.0), m.predicted, 1e-12);
    }
}
