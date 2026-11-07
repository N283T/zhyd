//! Conformations, site detection and group geometry of the empirical engine
//! ("Input preparation", "Conformations", "Internal hydrogens",
//! "Group geometry").

const std = @import("std");
const model_mod = @import("../model.zig");
const math = @import("../math.zig");
const pka = @import("../pka.zig");
const topology = @import("../place/topology.zig");
const params_mod = @import("params.zig");
const hydrogens = @import("hydrogens.zig");

pub const Vec = math.Vec3(f64);
pub const GroupType = params_mod.GroupType;
const Parameters = params_mod.Parameters;
const Model = model_mod.Model;

pub const none: u32 = std.math.maxInt(u32);

pub const Element = enum { h, c, n, o, s, f, cl, p, other };

pub const Atom = struct {
    pos: Vec,
    element: Element,
    name: [4]u8 = .{ ' ', ' ', ' ', ' ' },
    name_len: u8 = 0,
    residue_idx: u32,
    /// Index of the residue's auth chain identifier (equal ids share a key).
    chain_key: u32,
    res_num: i32,
    /// Position in the input (model atom index); hydrogens get `none`.
    order: u32,
    /// For an internal hydrogen: the heavy atom it is bonded to.
    parent: u32 = none,
    /// Modelled atom absent from the input (Q3 C-terminal oxygen): an
    /// interaction atom only, never a burial atom.
    virtual: bool = false,

    pub fn nameSlice(self: *const Atom) []const u8 {
        return self.name[0..self.name_len];
    }

    pub fn isHydrogen(self: Atom) bool {
        return self.element == .h;
    }

    pub fn isHeavy(self: Atom) bool {
        return self.element != .h;
    }
};

/// Small fixed-capacity list of atom indices (largest: ARG acid list, 8).
pub const AtomList = struct {
    items: [8]u32 = undefined,
    len: u8 = 0,

    pub fn append(self: *AtomList, idx: u32) void {
        self.items[self.len] = idx;
        self.len += 1;
    }

    pub fn slice(self: *const AtomList) []const u32 {
        return self.items[0..self.len];
    }
};

pub const Group = struct {
    gtype: GroupType,
    /// Titratable site kind; null for ROH, AMD, TRP, backbone and ions.
    kind: ?pka.SiteKind = null,
    residue_idx: u32,
    group_atom: u32,
    center: Vec,
    acid_list: AtomList = .{},
    base_list: AtomList = .{},
    /// Charge of the charged form (ions: the ion charge; partners: 0).
    q: f64 = 0,
    model_pka: f64 = 0,
    disulfide: bool = false,
    /// Covalently coupled to another titratable group and chosen not to
    /// titrate (engine.markCovalentCoupling, amendment A1).
    coupled_out: bool = false,

    // Energy terms, filled by the engine.
    count: f64 = 0,
    weight: f64 = 0,
    desolvation: f64 = 0,
    reorganisation: f64 = 0,
    side_chain: f64 = 0,
    backbone: f64 = 0,
    coulomb: f64 = 0,

    pub fn isTitratable(self: Group) bool {
        return self.kind != null and !self.disulfide and !self.coupled_out;
    }

    /// Interaction atoms used against `other` (base list if `other` is a base).
    pub fn listAgainst(self: *const Group, other: GroupType) []const u32 {
        return if (other.isBase()) self.base_list.slice() else self.acid_list.slice();
    }

    pub fn total(self: Group) f64 {
        return self.model_pka + self.desolvation + self.reorganisation + self.side_chain + self.backbone + self.coulomb;
    }
};

pub const Conformation = struct {
    allocator: std.mem.Allocator,
    atoms: std.ArrayListUnmanaged(Atom) = .empty,
    groups: std.ArrayListUnmanaged(Group) = .empty,
    /// Heavy-atom bond graph (CSR): neighbours of atom i are
    /// `bond_targets[bond_offsets[i]..bond_offsets[i + 1]]`.
    bond_offsets: []u32 = &.{},
    bond_targets: []u32 = &.{},

    pub fn deinit(self: *Conformation) void {
        self.atoms.deinit(self.allocator);
        self.groups.deinit(self.allocator);
        self.allocator.free(self.bond_offsets);
        self.allocator.free(self.bond_targets);
    }

    pub fn neighbours(self: *const Conformation, idx: u32) []const u32 {
        if (idx + 1 >= self.bond_offsets.len) return &.{};
        return self.bond_targets[self.bond_offsets[idx]..self.bond_offsets[idx + 1]];
    }

    /// True if atoms `a` and `b` are connected by a path of at most
    /// `max_bonds` bonds. A virtual atom is bonded to its `parent` only.
    pub fn withinBonds(self: *const Conformation, a_in: u32, b_in: u32, max_bonds_in: u8) bool {
        var a = a_in;
        var b = b_in;
        var max_bonds = max_bonds_in;
        if (a == b) return true;
        for ([_]*u32{ &a, &b }) |end| {
            const atom = self.atoms.items[end.*];
            if (!atom.virtual) continue;
            if (max_bonds == 0) return false;
            end.* = atom.parent;
            max_bonds -= 1;
        }
        if (a == b) return true;
        var frontier: [1024]u32 = undefined;
        var visited: [1024]u32 = undefined;
        var n_visited: usize = 1;
        visited[0] = a;
        frontier[0] = a;
        var lo: usize = 0;
        var hi: usize = 1;
        var depth: u8 = 0;
        while (depth < max_bonds and lo < hi) : (depth += 1) {
            const end = hi;
            while (lo < end) : (lo += 1) {
                for (self.neighbours(frontier[lo])) |nb| {
                    if (nb == b) return true;
                    if (std.mem.indexOfScalar(u32, visited[0..n_visited], nb) != null) continue;
                    if (n_visited == visited.len or hi == frontier.len) return false;
                    visited[n_visited] = nb;
                    n_visited += 1;
                    frontier[hi] = nb;
                    hi += 1;
                }
            }
        }
        return false;
    }
};

// ── Residue classification ────────────────────────────────────────────────────

/// Residue names written as ATOM records in the oracle's input (standard
/// polymer residues); every other residue is a HETATM residue.
fn isStandardPolymerName(comp: []const u8) bool {
    const ids = std.StaticStringMap(void).initComptime(.{
        .{ "ALA", {} }, .{ "ARG", {} }, .{ "ASN", {} }, .{ "ASP", {} }, .{ "CYS", {} },
        .{ "GLN", {} }, .{ "GLU", {} }, .{ "GLY", {} }, .{ "HIS", {} }, .{ "ILE", {} },
        .{ "LEU", {} }, .{ "LYS", {} }, .{ "MET", {} }, .{ "PHE", {} }, .{ "PRO", {} },
        .{ "SER", {} }, .{ "THR", {} }, .{ "TRP", {} }, .{ "TYR", {} }, .{ "VAL", {} },
        .{ "UNK", {} }, .{ "DA", {} },  .{ "DC", {} },  .{ "DG", {} },  .{ "DT", {} },
        .{ "DI", {} },  .{ "A", {} },   .{ "C", {} },   .{ "G", {} },   .{ "U", {} },
        .{ "I", {} },
    });
    return ids.has(comp);
}

pub fn isAtomRecord(res: model_mod.Residue) bool {
    return res.entity_type == .polymer and isStandardPolymerName(res.compIdSlice());
}

/// Standard amino acids plus UNK (an ATOM-record residue the oracle treats
/// as an amino acid: termini and backbone groups; amendment A11).
fn isAminoAcidName(comp: []const u8) bool {
    return topology.getBonds(comp) != null or std.mem.eql(u8, comp, "UNK");
}

/// Residues dropped before anything else: zreduce waters plus the parity list.
pub fn isRemoved(res: model_mod.Residue) bool {
    if (res.entity_type == .water) return true;
    const comp = res.compIdSlice();
    if (std.mem.eql(u8, comp, "WAT") or std.mem.eql(u8, comp, "DOD")) return true;
    return params_mod.isRemovedResidue(comp);
}

fn elementOf(atom: model_mod.Atom) Element {
    return switch (atom.element_type) {
        .H, .Har, .Hpol, .Ha_p, .HOd => .h,
        .C, .Car, .C_eq_O => .c,
        .N, .Nacc => .n,
        .O => .o,
        .S => .s,
        .F => .f,
        .Cl => .cl,
        .P => .p,
        else => .other,
    };
}

pub fn volumeClass(atom: Atom) params_mod.VolumeClass {
    return switch (atom.element) {
        .c => if (std.mem.eql(u8, atom.nameSlice(), "CA") or std.mem.eql(u8, atom.nameSlice(), "C")) .c_backbone else .c_other,
        .n => .n,
        .o => .o,
        .s => .s,
        .f => .f,
        .cl => .cl,
        .p => .p,
        .h, .other => .other,
    };
}

/// Coordinates as the oracle sees them (PDB text, three decimals).
fn inputPosition(p: math.Vec3(f32)) Vec {
    return .{
        .x = @round(@as(f64, p.x) * 1000.0) / 1000.0,
        .y = @round(@as(f64, p.y) * 1000.0) / 1000.0,
        .z = @round(@as(f64, p.z) * 1000.0) / 1000.0,
    };
}

/// Distinct alternate-location letters among the atoms the engine reads.
/// An empty result means a single conformation.
pub fn altlocLetters(mdl: *const Model, buf: *[64]u8) []const u8 {
    var n: usize = 0;
    for (mdl.atoms.items) |atom| {
        if (atom.is_hydrogen or atom.altloc == ' ') continue;
        if (isRemoved(mdl.residues.items[atom.residue_idx])) continue;
        if (std.mem.indexOfScalar(u8, buf[0..n], atom.altloc) != null) continue;
        if (n == buf.len) break;
        buf[n] = atom.altloc;
        n += 1;
    }
    std.mem.sort(u8, buf[0..n], {}, std.sort.asc(u8));
    return buf[0..n];
}

// ── Conformation building ─────────────────────────────────────────────────────

/// Parity quirk (amendment A10): in every conformation after the first, of
/// the residues sharing a chain and residue number (differing in insertion
/// code) only the last one in input order is present; the others, atoms and
/// groups, are absent. Returns a per-residue flag.
fn shadowedResidues(allocator: std.mem.Allocator, mdl: *const Model, chain_keys: []const u32) ![]bool {
    const Key = struct { chain: u32, num: i32 };
    var last: std.AutoHashMapUnmanaged(Key, u32) = .empty;
    defer last.deinit(allocator);
    for (mdl.residues.items, 0..) |res, ri| {
        try last.put(allocator, .{ .chain = chain_keys[res.chain_idx], .num = res.auth_seq_id }, @intCast(ri));
    }
    const out = try allocator.alloc(bool, mdl.residues.items.len);
    for (mdl.residues.items, 0..) |res, ri| {
        out[ri] = last.get(.{ .chain = chain_keys[res.chain_idx], .num = res.auth_seq_id }).? != ri;
    }
    return out;
}

/// Whether model atom `ai` belongs to the conformation of `letter`: blank
/// altloc atoms always; otherwise the atom with that letter, or, when its
/// residue has no such alternate for this atom name, the first alternate
/// listed (proposed amendment A3).
fn inConformation(mdl: *const Model, ai: u32, letter: u8) bool {
    const atom = mdl.atoms.items[ai];
    if (atom.altloc == ' ' or atom.altloc == letter) return true;
    const res = mdl.residues.items[atom.residue_idx];
    var first: ?u32 = null;
    for (mdl.atoms.items[res.atom_start..res.atom_end], res.atom_start..) |other, oi| {
        if (other.altloc == ' ' or other.is_hydrogen) continue;
        if (!std.mem.eql(u8, other.nameSlice(), atom.nameSlice())) continue;
        if (other.altloc == letter) return false;
        if (first == null) first = @intCast(oi);
    }
    return first == ai;
}

const ResidueRange = struct { start: u32, end: u32 };

const Builder = struct {
    allocator: std.mem.Allocator,
    mdl: *const Model,
    params: Parameters,
    conf: *Conformation,
    /// Per model residue: range of conformation atoms.
    ranges: []ResidueRange,
    chain_keys: []const u32,

    fn find(self: *const Builder, residue_idx: u32, name: []const u8) ?u32 {
        const r = self.ranges[residue_idx];
        var i = r.start;
        while (i < r.end) : (i += 1) {
            if (std.mem.eql(u8, self.conf.atoms.items[i].nameSlice(), name)) return i;
        }
        return null;
    }

    fn pos(self: *const Builder, idx: u32) Vec {
        return self.conf.atoms.items[idx].pos;
    }

    fn addHydrogen(self: *Builder, parent: u32, p: Vec) !u32 {
        const heavy = self.conf.atoms.items[parent];
        const idx: u32 = @intCast(self.conf.atoms.items.len);
        try self.conf.atoms.append(self.allocator, .{
            .pos = p,
            .element = .h,
            .residue_idx = heavy.residue_idx,
            .chain_key = heavy.chain_key,
            .res_num = heavy.res_num,
            .order = none,
            .parent = parent,
        });
        return idx;
    }

    /// Bisector hydrogen on N. The two heavy neighbours are the first two
    /// bonded atoms in input order (proposed amendment A5); `a`, `b` are the
    /// topology neighbours used when N has fewer than two bonds.
    fn bisectorH(self: *Builder, n: u32, a_topo: u32, b_topo: u32) !u32 {
        const nbs = self.conf.neighbours(n);
        const a = if (nbs.len >= 2) nbs[0] else a_topo;
        const b = if (nbs.len >= 2) nbs[1] else b_topo;
        return self.addHydrogen(n, hydrogens.bisector(self.pos(n), self.pos(a), self.pos(b), self.params.nh_length));
    }

    /// Planar pair on N bonded to X; Y₁, Y₂ are X's other bonded neighbours
    /// (in input order).
    fn planarH(self: *Builder, n: u32, x: u32) !?[2]u32 {
        var others: [2]u32 = undefined;
        var n_others: usize = 0;
        for (self.conf.neighbours(x)) |nb| {
            if (nb == n or n_others == 2) continue;
            others[n_others] = nb;
            n_others += 1;
        }
        if (n_others == 0) return null;
        const p2: ?Vec = if (n_others == 2) self.pos(others[1]) else null;
        const hs = hydrogens.planarPair(self.pos(n), self.pos(x), self.pos(others[0]), p2, self.params.nh_length);
        return .{ try self.addHydrogen(n, hs[0]), try self.addHydrogen(n, hs[1]) };
    }

    fn addGroup(self: *Builder, group: Group) !void {
        try self.conf.groups.append(self.allocator, group);
    }

    fn mean(self: *const Builder, idxs: []const u32) Vec {
        var s = Vec{ .x = 0, .y = 0, .z = 0 };
        for (idxs) |i| s = s.add(self.pos(i));
        return s.scale(1.0 / @as(f64, @floatFromInt(idxs.len)));
    }

    fn site(kind: pka.SiteKind, residue_idx: u32, group_atom: u32, center: Vec) Group {
        return .{
            .gtype = params_mod.groupTypeOf(kind),
            .kind = kind,
            .residue_idx = residue_idx,
            .group_atom = group_atom,
            .center = center,
            .q = params_mod.charge(kind),
            .model_pka = params_mod.modelPka(kind),
        };
    }

    fn bonded(self: *const Builder, a: u32, b: u32) bool {
        return std.mem.indexOfScalar(u32, self.conf.neighbours(a), b) != null;
    }

    /// Carboxyl group: the named oxygens bonded to the carbon. Without a
    /// bonded oxygen the centre falls back to the group atom and the group
    /// has no interaction atoms (proposed amendment A2).
    fn carboxyl(self: *Builder, residue_idx: u32, kind: pka.SiteKind, c: []const u8, o1: []const u8, o2: []const u8) !void {
        const ci = self.find(residue_idx, c) orelse return;
        var oxygens: AtomList = .{};
        for ([_][]const u8{ o1, o2 }) |name| {
            const o = self.find(residue_idx, name) orelse continue;
            if (self.bonded(ci, o)) oxygens.append(o);
        }
        const center = if (oxygens.len > 0) self.mean(oxygens.slice()) else self.pos(ci);
        var g = site(kind, residue_idx, ci, center);
        g.acid_list = oxygens;
        g.base_list = oxygens;
        try self.addGroup(g);
    }

    /// Q3: C-terminus without OXT. The missing oxygen is modelled at
    /// `cterm_co_length` from C along −(û_CA + û_O) (in the C, CA, O plane,
    /// about 120° from both bonds). It is the group atom (bonded to C only,
    /// ordered as C), an interaction atom and part of the centre, but never
    /// a burial atom.
    fn carboxylCtermVirtual(self: *Builder, residue_idx: u32) !void {
        const c = self.find(residue_idx, "C") orelse return;
        const o = self.find(residue_idx, "O") orelse return;
        const ca = self.find(residue_idx, "CA") orelse return;
        if (!self.bonded(c, o)) return;
        const heavy = self.conf.atoms.items[c];
        const v: u32 = @intCast(self.conf.atoms.items.len);
        try self.conf.atoms.append(self.allocator, .{
            .pos = hydrogens.bisector(self.pos(c), self.pos(ca), self.pos(o), cterm_co_length),
            .element = .o,
            .name = .{ 'O', 'X', 'T', ' ' },
            .name_len = 3,
            .residue_idx = heavy.residue_idx,
            .chain_key = heavy.chain_key,
            .res_num = heavy.res_num,
            .order = heavy.order,
            .parent = c,
            .virtual = true,
        });
        var g = site(.c_term, residue_idx, v, self.mean(&.{ o, v }));
        g.acid_list.append(o);
        g.acid_list.append(v);
        g.base_list = g.acid_list;
        try self.addGroup(g);
    }

    fn carboxylCterm(self: *Builder, residue_idx: u32) !void {
        const oxt = self.find(residue_idx, "OXT") orelse return;
        const c = self.find(residue_idx, "C") orelse return;
        var oxygens: AtomList = .{};
        for ([_][]const u8{ "O", "OXT" }) |name| {
            const o = self.find(residue_idx, name) orelse continue;
            if (self.bonded(c, o)) oxygens.append(o);
        }
        const center = if (oxygens.len > 0) self.mean(oxygens.slice()) else self.pos(oxt);
        var g = site(.c_term, residue_idx, oxt, center);
        g.acid_list = oxygens;
        g.base_list = oxygens;
        try self.addGroup(g);
    }

    /// Site defined by one atom (CYS, TYR, LYS, N-terminus).
    fn singleAtom(self: *Builder, residue_idx: u32, kind: pka.SiteKind, name: []const u8) !void {
        const i = self.find(residue_idx, name) orelse return;
        var g = site(kind, residue_idx, i, self.pos(i));
        g.acid_list.append(i);
        g.base_list.append(i);
        try self.addGroup(g);
    }

    fn histidine(self: *Builder, residue_idx: u32) !void {
        const cg = self.find(residue_idx, "CG") orelse return;
        const nd1 = self.find(residue_idx, "ND1") orelse return;
        const cd2 = self.find(residue_idx, "CD2") orelse return;
        const ce1 = self.find(residue_idx, "CE1") orelse return;
        const ne2 = self.find(residue_idx, "NE2") orelse return;
        const hd1 = try self.bisectorH(nd1, cg, ce1);
        const he2 = try self.bisectorH(ne2, cd2, ce1);
        var g = site(.his, residue_idx, cg, self.mean(&.{ cg, nd1, cd2, ce1, ne2 }));
        for ([_]u32{ nd1, ne2, hd1, he2 }) |i| g.acid_list.append(i);
        g.base_list.append(nd1);
        g.base_list.append(ne2);
        try self.addGroup(g);
    }

    fn arginine(self: *Builder, residue_idx: u32) !void {
        const cd = self.find(residue_idx, "CD") orelse return;
        const ne = self.find(residue_idx, "NE") orelse return;
        const cz = self.find(residue_idx, "CZ") orelse return;
        const nh1 = self.find(residue_idx, "NH1") orelse return;
        const nh2 = self.find(residue_idx, "NH2") orelse return;
        const he = try self.bisectorH(ne, cd, cz);
        const hh1 = (try self.planarH(nh1, cz)) orelse return;
        const hh2 = (try self.planarH(nh2, cz)) orelse return;
        var g = site(.arg, residue_idx, cz, self.pos(cz));
        for ([_]u32{ ne, nh1, nh2, he, hh1[0], hh1[1], hh2[0], hh2[1] }) |i| g.acid_list.append(i);
        for ([_]u32{ ne, nh1, nh2 }) |i| g.base_list.append(i);
        try self.addGroup(g);
    }

    fn hydroxyl(self: *Builder, residue_idx: u32, name: []const u8) !void {
        const i = self.find(residue_idx, name) orelse return;
        var g = Group{ .gtype = .roh, .residue_idx = residue_idx, .group_atom = i, .center = self.pos(i) };
        g.acid_list.append(i);
        g.base_list.append(i);
        try self.addGroup(g);
    }

    /// Amide: hydrogens as a planar pair on N with X = the carbonyl carbon,
    /// Y₁ = X's carbon neighbour, Y₂ = the carbonyl oxygen.
    fn amide(self: *Builder, residue_idx: u32, c_prev: []const u8, c: []const u8, o: []const u8, n: []const u8) !void {
        const pi = self.find(residue_idx, c_prev) orelse return;
        const ci = self.find(residue_idx, c) orelse return;
        const oi = self.find(residue_idx, o) orelse return;
        const ni = self.find(residue_idx, n) orelse return;
        const hp = hydrogens.planarPair(self.pos(ni), self.pos(ci), self.pos(pi), self.pos(oi), self.params.nh_length);
        const h0 = try self.addHydrogen(ni, hp[0]);
        const h1 = try self.addHydrogen(ni, hp[1]);
        var g = Group{ .gtype = .amd, .residue_idx = residue_idx, .group_atom = ci, .center = self.mean(&.{ oi, ni }) };
        for ([_]u32{ ni, h0, h1 }) |i| g.acid_list.append(i);
        g.base_list.append(oi);
        try self.addGroup(g);
    }

    fn tryptophan(self: *Builder, residue_idx: u32) !void {
        const ne1 = self.find(residue_idx, "NE1") orelse return;
        const cd1 = self.find(residue_idx, "CD1") orelse return;
        const ce2 = self.find(residue_idx, "CE2") orelse return;
        const he1 = try self.bisectorH(ne1, cd1, ce2);
        var g = Group{ .gtype = .trp, .residue_idx = residue_idx, .group_atom = ne1, .center = self.pos(ne1) };
        g.acid_list.append(ne1);
        g.acid_list.append(he1);
        g.base_list.append(ne1);
        try self.addGroup(g);
    }

    /// `is_cterm_carboxyl`: the residue's C is a Q3 C-terminal carboxyl and
    /// forms no backbone C=O group.
    fn backbone(self: *Builder, residue_idx: u32, is_nterm: bool, is_cterm_carboxyl: bool) !void {
        const comp = self.mdl.residues.items[residue_idx].compIdSlice();
        if (!is_nterm and !std.mem.eql(u8, comp, "PRO")) bbn: {
            // Backbone N with exactly two heavy neighbours (previous C, CA).
            const n = self.find(residue_idx, "N") orelse break :bbn;
            const nbs = self.conf.neighbours(n);
            if (nbs.len != 2) break :bbn;
            const h = try self.bisectorH(n, nbs[0], nbs[1]);
            var g = Group{ .gtype = .bbn, .residue_idx = residue_idx, .group_atom = n, .center = self.pos(n) };
            g.acid_list.append(n);
            g.acid_list.append(h);
            g.base_list = g.acid_list;
            try self.addGroup(g);
        }
        bbc: {
            if (is_cterm_carboxyl) break :bbc;
            // Backbone C bonded to exactly one oxygen.
            const c = self.find(residue_idx, "C") orelse break :bbc;
            var oxygen: ?u32 = null;
            var n_oxygen: u32 = 0;
            for (self.conf.neighbours(c)) |nb| {
                if (self.conf.atoms.items[nb].element != .o) continue;
                n_oxygen += 1;
                oxygen = nb;
            }
            if (n_oxygen != 1) break :bbc;
            var g = Group{ .gtype = .bbc, .residue_idx = residue_idx, .group_atom = c, .center = self.pos(c) };
            g.acid_list.append(oxygen.?);
            g.base_list.append(oxygen.?);
            try self.addGroup(g);
        }
    }
};

/// Heavy atoms closer than or at this distance are bonded (proposed
/// amendment A2: the oracle bonds C–O at 2.000 Å but not at 2.005 Å).
pub const bond_max: f64 = 2.0;
/// C–O length of the modelled Q3 C-terminal oxygen (carboxylate C–O).
pub const cterm_co_length: f64 = 1.25;
/// Two sulfur atoms closer than this are bonded (disulfide).
pub const disulfide_bond_max: f64 = 2.5;

/// Build one conformation: the blank-altloc atoms plus those with `altloc`
/// (null: all atoms, for a structure without alternate locations).
pub fn build(allocator: std.mem.Allocator, mdl: *const Model, altloc: ?u8, first_altloc: ?u8, params: Parameters) !Conformation {
    const later_conformation = altloc != null and first_altloc != null and altloc.? != first_altloc.?;
    var conf = Conformation{ .allocator = allocator };
    errdefer conf.deinit();

    const ranges = try allocator.alloc(ResidueRange, mdl.residues.items.len);
    defer allocator.free(ranges);
    @memset(ranges, .{ .start = 0, .end = 0 });

    // Chain keys by auth chain id (the oracle aliases equal ids together).
    var chain_keys = try allocator.alloc(u32, mdl.chains.items.len);
    defer allocator.free(chain_keys);
    for (mdl.chains.items, 0..) |*chain, ci| {
        chain_keys[ci] = @intCast(ci);
        for (mdl.chains.items[0..ci], 0..) |*other, oi| {
            if (std.mem.eql(u8, chain.authSlice(), other.authSlice())) {
                chain_keys[ci] = chain_keys[oi];
                break;
            }
        }
    }

    const shadowed = try shadowedResidues(allocator, mdl, chain_keys);
    defer allocator.free(shadowed);

    for (mdl.atoms.items, 0..) |atom, ai| {
        if (atom.is_hydrogen) continue;
        if (ai >= mdl.original_atom_count and mdl.original_atom_count != 0) continue;
        const res = mdl.residues.items[atom.residue_idx];
        if (isRemoved(res)) continue;
        if (altloc) |letter| {
            if (!inConformation(mdl, @intCast(ai), letter)) continue;
        }
        if (later_conformation and shadowed[atom.residue_idx] and !params.fixes.q8) continue;
        const element = elementOf(atom);
        if (element == .h) continue;
        const idx: u32 = @intCast(conf.atoms.items.len);
        var a = Atom{
            .pos = inputPosition(atom.pos),
            .element = element,
            .residue_idx = atom.residue_idx,
            .chain_key = chain_keys[res.chain_idx],
            .res_num = res.auth_seq_id,
            .order = @intCast(ai),
        };
        const name = atom.nameSlice();
        @memcpy(a.name[0..name.len], name);
        a.name_len = @intCast(name.len);
        const r = &ranges[atom.residue_idx];
        if (r.end == r.start) r.start = idx;
        r.end = idx + 1;
        try conf.atoms.append(allocator, a);
    }
    try buildBonds(allocator, &conf);

    var b = Builder{ .allocator = allocator, .mdl = mdl, .params = params, .conf = &conf, .ranges = ranges, .chain_keys = chain_keys };

    // zreduce termini (Q3, Q4): the first residue of each zreduce chain is
    // an N-terminus, the last polymer residue a C-terminus; with
    // `break_termini` (--nterm aggressive) also the residues after / before
    // a chain break.
    const termini = try allocator.alloc(Termini, mdl.residues.items.len);
    defer allocator.free(termini);
    markTermini(mdl, termini, params.break_termini);

    // N-terminus: first ATOM-record residue after a chain start, a TER (the
    // oracle writes one after every polymer chain) or an OXT.
    var need_nterm = true;
    var nterm_family: ?struct { chain: u32, num: i32 } = null;
    var prev_chain: ?u32 = null;
    var prev_polymer = false;
    var prev_auth: ?u32 = null;
    for (mdl.residues.items, 0..) |res, ri| {
        const residue_idx: u32 = @intCast(ri);
        if (ranges[ri].end == ranges[ri].start) continue;
        if (prev_chain) |pc| {
            if (pc != res.chain_idx and prev_polymer) need_nterm = true;
        }
        if (prev_auth) |pa| {
            if (pa != chain_keys[res.chain_idx]) need_nterm = true;
        }
        prev_chain = res.chain_idx;
        prev_auth = chain_keys[res.chain_idx];
        prev_polymer = res.entity_type == .polymer;

        if (params.fixes.q5) {
            // Q5: a non-polymer residue with one heavy atom whose element is
            // a metal (or halogen) with a known ion charge.
            if (res.entity_type != .polymer and ranges[ri].end - ranges[ri].start == 1) {
                const i = ranges[ri].start;
                const element = mdl.atoms.items[conf.atoms.items[i].order].element_type;
                if (params_mod.ionChargeByElement(res.compIdSlice(), element)) |q| {
                    try b.addGroup(.{ .gtype = .ion, .residue_idx = residue_idx, .group_atom = i, .center = b.pos(i), .q = q });
                    continue;
                }
            }
        } else if (params_mod.ionCharge(res.compIdSlice())) |q| {
            var i = ranges[ri].start;
            while (i < ranges[ri].end) : (i += 1) {
                try b.addGroup(.{ .gtype = .ion, .residue_idx = residue_idx, .group_atom = i, .center = b.pos(i), .q = q });
            }
            continue;
        }
        if (!isAtomRecord(res)) continue;
        const comp = res.compIdSlice();
        if (!isAminoAcidName(comp)) continue;

        // Residues sharing the N-terminal residue's chain and number (other
        // insertion codes) are N-termini too (amendment A12).
        // Q8 drops the insertion-code family rule.
        var is_nterm = false;
        if (params.fixes.q4) {
            if (termini[ri].n and b.find(residue_idx, "N") != null) {
                is_nterm = true;
                try b.singleAtom(residue_idx, .n_term, "N");
            }
        } else if (need_nterm or (!params.fixes.q8 and nterm_family != null and nterm_family.?.chain == chain_keys[res.chain_idx] and nterm_family.?.num == res.auth_seq_id)) {
            if (need_nterm) nterm_family = .{ .chain = chain_keys[res.chain_idx], .num = res.auth_seq_id };
            need_nterm = false;
            if (b.find(residue_idx, "N") != null) {
                is_nterm = true;
                try b.singleAtom(residue_idx, .n_term, "N");
            }
        } else nterm_family = null;
        const has_oxt = b.find(residue_idx, "OXT") != null;
        const virtual_cterm = params.fixes.q3 and !has_oxt and termini[ri].c and
            b.find(residue_idx, "CA") != null and b.find(residue_idx, "O") != null;
        try b.backbone(residue_idx, is_nterm, virtual_cterm);

        if (std.mem.eql(u8, comp, "ASP")) {
            try b.carboxyl(residue_idx, .asp, "CG", "OD1", "OD2");
        } else if (std.mem.eql(u8, comp, "GLU")) {
            try b.carboxyl(residue_idx, .glu, "CD", "OE1", "OE2");
        } else if (std.mem.eql(u8, comp, "HIS")) {
            try b.histidine(residue_idx);
        } else if (std.mem.eql(u8, comp, "CYS")) {
            try b.singleAtom(residue_idx, .cys, "SG");
        } else if (std.mem.eql(u8, comp, "TYR")) {
            try b.singleAtom(residue_idx, .tyr, "OH");
        } else if (std.mem.eql(u8, comp, "LYS")) {
            try b.singleAtom(residue_idx, .lys, "NZ");
        } else if (std.mem.eql(u8, comp, "ARG")) {
            try b.arginine(residue_idx);
        } else if (std.mem.eql(u8, comp, "SER")) {
            try b.hydroxyl(residue_idx, "OG");
        } else if (std.mem.eql(u8, comp, "THR")) {
            try b.hydroxyl(residue_idx, "OG1");
        } else if (std.mem.eql(u8, comp, "ASN")) {
            try b.amide(residue_idx, "CB", "CG", "OD1", "ND2");
        } else if (std.mem.eql(u8, comp, "GLN")) {
            try b.amide(residue_idx, "CG", "CD", "OE1", "NE2");
        } else if (std.mem.eql(u8, comp, "TRP")) {
            try b.tryptophan(residue_idx);
        }
        if (has_oxt) {
            try b.carboxylCterm(residue_idx);
            need_nterm = true;
        } else if (virtual_cterm) try b.carboxylCtermVirtual(residue_idx);
    }

    // Disulfides: two CYS SG closer than the disulfide distance.
    for (conf.groups.items, 0..) |*g, gi| {
        if (g.gtype != .cys) continue;
        for (conf.groups.items, 0..) |other, oi| {
            if (oi == gi or other.gtype != .cys) continue;
            if (g.center.distance(other.center) < params.disulfide_distance) g.disulfide = true;
        }
    }

    // Model pKa from the parameter table (built-in or overridden).
    for (conf.groups.items) |*g| {
        if (g.kind) |k| g.model_pka = params.model_pka.get(k);
    }

    std.mem.sort(Group, conf.groups.items, &conf, groupOrderLess);
    return conf;
}

const Termini = struct { n: bool = false, c: bool = false };

/// zreduce termini: first residue of each chain and last polymer residue;
/// with `breaks`, also around `is_chain_break_before` (the same rule as the
/// hydrogen placer under --nterm aggressive).
fn markTermini(mdl: *const Model, out: []Termini, breaks: bool) void {
    @memset(out, .{});
    for (mdl.chains.items) |chain| {
        if (chain.residue_end <= chain.residue_start) continue;
        out[chain.residue_start].n = true;
        var last: ?u32 = null;
        var ri = chain.residue_start;
        while (ri < chain.residue_end) : (ri += 1) {
            const res = mdl.residues.items[ri];
            if (res.entity_type == .polymer) last = ri;
            if (breaks and ri > chain.residue_start and res.is_chain_break_before) {
                out[ri].n = true;
                out[ri - 1].c = true;
            }
        }
        if (last) |l| out[l].c = true;
    }
}

fn groupOrderLess(conf: *const Conformation, a: Group, b: Group) bool {
    return conf.atoms.items[a.group_atom].order < conf.atoms.items[b.group_atom].order;
}

/// Maximum bond separation of two titratable group atoms treated as
/// covalently coupled.
pub const covalent_coupling_bonds: u8 = 3;

const Edge = struct { a: u32, b: u32 };

/// Distance-based bond graph over the heavy atoms (uniform grid search).
fn buildBonds(allocator: std.mem.Allocator, conf: *Conformation) !void {
    const atoms = conf.atoms.items;
    const n = atoms.len;
    var edges: std.ArrayListUnmanaged(Edge) = .empty;
    defer edges.deinit(allocator);

    if (n > 0) {
        const cell = disulfide_bond_max;
        var lo = atoms[0].pos;
        var hi = atoms[0].pos;
        for (atoms) |a| {
            lo = .{ .x = @min(lo.x, a.pos.x), .y = @min(lo.y, a.pos.y), .z = @min(lo.z, a.pos.z) };
            hi = .{ .x = @max(hi.x, a.pos.x), .y = @max(hi.y, a.pos.y), .z = @max(hi.z, a.pos.z) };
        }
        const dims = [3]u64{
            @as(u64, @intFromFloat((hi.x - lo.x) / cell)) + 1,
            @as(u64, @intFromFloat((hi.y - lo.y) / cell)) + 1,
            @as(u64, @intFromFloat((hi.z - lo.z) / cell)) + 1,
        };
        const Keyed = struct { key: u64, idx: u32 };
        const keyed = try allocator.alloc(Keyed, n);
        defer allocator.free(keyed);
        for (atoms, 0..) |a, i| {
            const c = cellOf(a.pos, lo, cell);
            keyed[i] = .{ .key = c[0] + dims[0] * (c[1] + dims[1] * c[2]), .idx = @intCast(i) };
        }
        std.mem.sort(Keyed, keyed, {}, struct {
            fn less(_: void, x: Keyed, y: Keyed) bool {
                return x.key < y.key or (x.key == y.key and x.idx < y.idx);
            }
        }.less);
        for (atoms, 0..) |a, i| {
            const c = cellOf(a.pos, lo, cell);
            var dz: i64 = -1;
            while (dz <= 1) : (dz += 1) {
                var dy: i64 = -1;
                while (dy <= 1) : (dy += 1) {
                    var dx: i64 = -1;
                    while (dx <= 1) : (dx += 1) {
                        const cx = @as(i64, @intCast(c[0])) + dx;
                        const cy = @as(i64, @intCast(c[1])) + dy;
                        const cz = @as(i64, @intCast(c[2])) + dz;
                        if (cx < 0 or cy < 0 or cz < 0) continue;
                        if (cx >= dims[0] or cy >= dims[1] or cz >= dims[2]) continue;
                        const key = @as(u64, @intCast(cx)) + dims[0] * (@as(u64, @intCast(cy)) + dims[1] * @as(u64, @intCast(cz)));
                        var k = std.sort.lowerBound(Keyed, keyed, key, struct {
                            fn order(target: u64, x: Keyed) std.math.Order {
                                return std.math.order(target, x.key);
                            }
                        }.order);
                        while (k < n and keyed[k].key == key) : (k += 1) {
                            const j = keyed[k].idx;
                            if (j <= i) continue;
                            const b = atoms[j];
                            const d = a.pos.distance(b.pos);
                            const limit = if (a.element == .s and b.element == .s) disulfide_bond_max else bond_max;
                            const ok = if (a.element == .s and b.element == .s) d < limit else d <= limit;
                            if (ok) try edges.append(allocator, .{ .a = @intCast(i), .b = j });
                        }
                    }
                }
            }
        }
    }

    const offsets = try allocator.alloc(u32, n + 1);
    errdefer allocator.free(offsets);
    @memset(offsets, 0);
    for (edges.items) |e| {
        offsets[e.a + 1] += 1;
        offsets[e.b + 1] += 1;
    }
    for (1..n + 1) |i| offsets[i] += offsets[i - 1];
    const targets = try allocator.alloc(u32, offsets[n]);
    errdefer allocator.free(targets);
    const fill = try allocator.alloc(u32, n);
    defer allocator.free(fill);
    @memcpy(fill, offsets[0..n]);
    // Edges are generated in ascending (i, j) order; keep neighbours sorted
    // by input position.
    for (edges.items) |e| {
        targets[fill[e.a]] = e.b;
        fill[e.a] += 1;
    }
    for (edges.items) |e| {
        targets[fill[e.b]] = e.a;
        fill[e.b] += 1;
    }
    for (0..n) |i| std.mem.sort(u32, targets[offsets[i]..offsets[i + 1]], {}, std.sort.asc(u32));
    conf.bond_offsets = offsets;
    conf.bond_targets = targets;
}

fn cellOf(p: Vec, lo: Vec, cell: f64) [3]u64 {
    return .{
        @intFromFloat((p.x - lo.x) / cell),
        @intFromFloat((p.y - lo.y) / cell),
        @intFromFloat((p.z - lo.z) / cell),
    };
}

const testing = std.testing;
const mmcif = @import("../mmcif.zig");

fn parsePolymer(source: []const u8) !Model {
    const mdl = try mmcif.parseModel(testing.allocator, source);
    for (mdl.residues.items) |*res| res.entity_type = .polymer;
    return mdl;
}

fn findGroup(conf: *const Conformation, gtype: GroupType) ?*const Group {
    for (conf.groups.items) |*g| if (g.gtype == gtype) return g;
    return null;
}

test "histidine groups, internal hydrogens and centre" {
    var mdl = try parsePolymer(@embedFile("../test_data/his.cif"));
    defer mdl.deinit();
    var conf = try build(testing.allocator, &mdl, null, null, .{});
    defer conf.deinit();
    // First residue: N-terminus (no BBN), backbone C=O, HIS ring.
    try testing.expect(findGroup(&conf, .nterm) != null);
    try testing.expect(findGroup(&conf, .bbn) == null);
    try testing.expect(findGroup(&conf, .bbc) != null);
    const his = findGroup(&conf, .his).?;
    try testing.expectEqual(@as(u8, 4), his.acid_list.len);
    try testing.expectEqual(@as(u8, 2), his.base_list.len);
    // Centre: mean of CG ND1 CD2 CE1 NE2 (hand-computed from the fixture).
    try testing.expectApproxEqAbs(@as(f64, 0.08), his.center.x, 1e-9);
    try testing.expectApproxEqAbs(@as(f64, -3.36), his.center.y, 1e-9);
    try testing.expectApproxEqAbs(@as(f64, 0.96), his.center.z, 1e-9);
    // HD1 and HE2 sit 1.01 Å from their nitrogens.
    for (his.acid_list.slice()[2..]) |h| {
        const atom = conf.atoms.items[h];
        try testing.expect(atom.isHydrogen());
        try testing.expectApproxEqAbs(@as(f64, 1.01), atom.pos.distance(conf.atoms.items[atom.parent].pos), 1e-9);
    }
    // Groups are ordered by their group atom's input position: N (atom 0)
    // before the backbone C (atom 2) before CG (atom 5).
    try testing.expectEqual(GroupType.nterm, conf.groups.items[0].gtype);
    try testing.expectEqual(GroupType.bbc, conf.groups.items[1].gtype);
    try testing.expectEqual(GroupType.his, conf.groups.items[2].gtype);
}

/// One-residue ASP model; `od2_shift` moves OD2 along +x.
fn aspModel(od1_shift: f32, od2_shift: f32) !Model {
    var mdl = Model.init(testing.allocator);
    errdefer mdl.deinit();
    var chain = model_mod.Chain{ .residue_end = 1 };
    chain.setAuthAsymId("A");
    try mdl.chains.append(testing.allocator, chain);
    var res = model_mod.Residue{ .auth_seq_id = 1, .atom_end = 5, .entity_type = .polymer };
    res.setCompId("ASP");
    try mdl.residues.append(testing.allocator, res);
    const Spec = struct { name: []const u8, element: @import("../element.zig").AtomType, x: f32, y: f32 };
    const specs = [_]Spec{
        .{ .name = "CA", .element = .C, .x = -1.5, .y = 0 },
        .{ .name = "CB", .element = .C, .x = 0, .y = 0 },
        .{ .name = "CG", .element = .C, .x = 1.5, .y = 0 },
        .{ .name = "OD1", .element = .O, .x = 2.1 + od1_shift, .y = 1.1 },
        .{ .name = "OD2", .element = .O, .x = 2.1 + od2_shift, .y = -1.1 },
    };
    for (specs) |sp| {
        var atom = model_mod.Atom{ .pos = .{ .x = sp.x, .y = sp.y, .z = 0 }, .element_type = sp.element };
        atom.setName(sp.name);
        try mdl.atoms.append(testing.allocator, atom);
    }
    return mdl;
}

test "carboxyl uses only oxygens bonded to the carbon" {
    // CG–O = sqrt(0.6² + 1.1²) = 1.253 Å: both bonded, centre = (2.1, 0).
    {
        var mdl = try aspModel(0, 0);
        defer mdl.deinit();
        var conf = try build(testing.allocator, &mdl, null, null, .{});
        defer conf.deinit();
        const g = findGroup(&conf, .coo).?;
        try testing.expectEqual(@as(u8, 2), g.acid_list.len);
        try testing.expectApproxEqAbs(@as(f64, 2.1), g.center.x, 1e-6);
        try testing.expectApproxEqAbs(@as(f64, 0), g.center.y, 1e-6);
    }
    // OD2 moved 1.5 Å out: CG–OD2 = sqrt(2.1² + 1.1²) = 2.37 Å > 2.0, so the
    // group is OD1 alone.
    {
        var mdl = try aspModel(0, 1.5);
        defer mdl.deinit();
        var conf = try build(testing.allocator, &mdl, null, null, .{});
        defer conf.deinit();
        const g = findGroup(&conf, .coo).?;
        try testing.expectEqual(@as(u8, 1), g.acid_list.len);
        try testing.expectApproxEqAbs(@as(f64, 2.1), g.center.x, 1e-6);
        try testing.expectApproxEqAbs(@as(f64, 1.1), g.center.y, 1e-6);
    }
    // Both oxygens unbonded: the centre falls back to CG, no interaction atoms.
    {
        var mdl = try aspModel(1.5, 1.5);
        defer mdl.deinit();
        var conf = try build(testing.allocator, &mdl, null, null, .{});
        defer conf.deinit();
        const g = findGroup(&conf, .coo).?;
        try testing.expectEqual(@as(u8, 0), g.acid_list.len);
        try testing.expectApproxEqAbs(@as(f64, 1.5), g.center.x, 1e-6);
    }
}

test "conformations: blank atoms plus the letter, else the first alternate" {
    var mdl = try parsePolymer(@embedFile("../test_data/ala_altloc.cif"));
    defer mdl.deinit();
    // Rename CB B to CG: CB then only has an A alternate and CG only a B one.
    mdl.atoms.items[6].setName("CG");
    var buf: [64]u8 = undefined;
    try testing.expectEqualStrings("AB", altlocLetters(&mdl, &buf));
    var conf_b = try build(testing.allocator, &mdl, 'B', 'A', .{});
    defer conf_b.deinit();
    var names: [8][]const u8 = undefined;
    var n: usize = 0;
    for (conf_b.atoms.items) |*a| {
        if (a.isHydrogen()) continue;
        names[n] = a.nameSlice();
        n += 1;
    }
    // N, CA(B), C, O, CB(A, fallback), CG(B).
    try testing.expectEqual(@as(usize, 6), n);
    try testing.expectApproxEqAbs(@as(f64, 0.1), conf_b.atoms.items[1].pos.x, 1e-9);
    try testing.expectEqualStrings("CB", names[4]);
    try testing.expectEqualStrings("CG", names[5]);
}

const AtomSpec = struct {
    name: []const u8,
    element: @import("../element.zig").AtomType,
    x: f32,
    y: f32 = 0,
    z: f32 = 0,
    altloc: u8 = ' ',
};

const ResidueSpec = struct {
    comp: []const u8,
    num: i32,
    ins: u8 = ' ',
    entity: model_mod.residue.EntityType = .polymer,
    chain_break_before: bool = false,
    atoms: []const AtomSpec,
};

/// Model from residue specs; `chain_ends` lists the exclusive end residue
/// index of each chain (all chains have auth id "A" + index).
fn specModel(residues: []const ResidueSpec, chain_ends: []const u32) !Model {
    var mdl = Model.init(testing.allocator);
    errdefer mdl.deinit();
    var start: u32 = 0;
    for (chain_ends, 0..) |end, ci| {
        var chain = model_mod.Chain{ .residue_start = start, .residue_end = end };
        const ids = [_][]const u8{ "A", "B", "C" };
        chain.setAuthAsymId(ids[ci]);
        try mdl.chains.append(testing.allocator, chain);
        start = end;
    }
    for (residues, 0..) |spec, ri| {
        const chain_idx: u32 = for (chain_ends, 0..) |end, ci| {
            if (ri < end) break @intCast(ci);
        } else unreachable;
        var res = model_mod.Residue{
            .auth_seq_id = spec.num,
            .ins_code = spec.ins,
            .chain_idx = chain_idx,
            .entity_type = spec.entity,
            .is_chain_break_before = spec.chain_break_before,
            .atom_start = @intCast(mdl.atoms.items.len),
        };
        res.setCompId(spec.comp);
        for (spec.atoms) |a| {
            var atom = model_mod.Atom{ .pos = .{ .x = a.x, .y = a.y, .z = a.z }, .element_type = a.element, .altloc = a.altloc, .residue_idx = @intCast(ri) };
            atom.setName(a.name);
            try mdl.atoms.append(testing.allocator, atom);
        }
        res.atom_end = @intCast(mdl.atoms.items.len);
        try mdl.residues.append(testing.allocator, res);
    }
    return mdl;
}

/// Backbone of an amino acid at offset `x0` (N, CA, C, O; C–O 1.25 Å with
/// û_O = (0.6, 0.8) and û_CA = (−1, 0) at C).
fn backboneAtoms(x0: f32) [4]AtomSpec {
    return .{
        .{ .name = "N", .element = .N, .x = x0 - 2.9, .y = 0.5 },
        .{ .name = "CA", .element = .C, .x = x0 - 1.5 },
        .{ .name = "C", .element = .C, .x = x0 },
        .{ .name = "O", .element = .O, .x = x0 + 0.75, .y = 1.0 },
    };
}

fn countGroups(conf: *const Conformation, gtype: GroupType) usize {
    var n: usize = 0;
    for (conf.groups.items) |g| n += @intFromBool(g.gtype == gtype);
    return n;
}

test "Q3: chain-last residue without OXT is a C-terminus with a modelled oxygen" {
    const ala = backboneAtoms(0);
    var mdl = try specModel(&.{.{ .comp = "ALA", .num = 1, .atoms = &ala }}, &.{1});
    defer mdl.deinit();
    // Parity: no OXT, no C-terminus; the C=O is a backbone group.
    {
        var conf = try build(testing.allocator, &mdl, null, null, .{});
        defer conf.deinit();
        try testing.expect(findGroup(&conf, .coo) == null);
        try testing.expectEqual(@as(usize, 1), countGroups(&conf, .bbc));
    }
    var conf = try build(testing.allocator, &mdl, null, null, .{ .fixes = .{ .q3 = true } });
    defer conf.deinit();
    try testing.expectEqual(@as(usize, 0), countGroups(&conf, .bbc));
    const g = findGroup(&conf, .coo).?;
    try testing.expectEqual(pka.SiteKind.c_term, g.kind.?);
    // OXT = C + 1.25 · unit(−(û_CA + û_O)) = 1.25 · (0.4, −0.8)/√0.8
    //     = (0.559017, −1.118034); centre = mean(O, OXT).
    const oxt = conf.atoms.items[g.group_atom];
    try testing.expect(oxt.virtual);
    try testing.expectApproxEqAbs(@as(f64, 0.5590170), oxt.pos.x, 1e-6);
    try testing.expectApproxEqAbs(@as(f64, -1.1180340), oxt.pos.y, 1e-6);
    try testing.expectApproxEqAbs(@as(f64, 0.6545085), g.center.x, 1e-6);
    try testing.expectApproxEqAbs(@as(f64, -0.0590170), g.center.y, 1e-6);
    try testing.expectEqual(@as(u8, 2), g.acid_list.len);
    // The modelled oxygen is bonded to C only: OXT–C–CA is two bonds.
    const ca = for (conf.atoms.items, 0..) |a, i| {
        if (std.mem.eql(u8, a.nameSlice(), "CA")) break @as(u32, @intCast(i));
    } else unreachable;
    try testing.expect(conf.withinBonds(g.group_atom, ca, 2));
    try testing.expect(!conf.withinBonds(g.group_atom, ca, 1));
}

test "Q4: N-termini follow zreduce chains and chain breaks" {
    const ace = [_]AtomSpec{.{ .name = "C", .element = .C, .x = -20 }};
    const a1 = backboneAtoms(0);
    const a2 = backboneAtoms(20);
    const a3 = backboneAtoms(40);
    // Chain A: ACE cap, ALA 1, ALA 2. Chain B: ALA 10, ALA 12 after a break.
    var mdl = try specModel(&.{
        .{ .comp = "ACE", .num = 0, .atoms = &ace },
        .{ .comp = "ALA", .num = 1, .atoms = &a1 },
        .{ .comp = "ALA", .num = 2, .atoms = &a2 },
        .{ .comp = "ALA", .num = 10, .atoms = &a3 },
        .{ .comp = "ALA", .num = 12, .chain_break_before = true, .atoms = &backboneAtoms(60) },
    }, &.{ 3, 5 });
    defer mdl.deinit();
    const Case = struct { params: Parameters, nterm: []const u32 };
    const cases = [_]Case{
        // Parity: the first ATOM-record residue of each chain (ACE is not one).
        .{ .params = .{}, .nterm = &.{ 1, 3 } },
        // Q4: the capped chain has no N-terminus; breaks only when aggressive.
        .{ .params = .{ .fixes = .{ .q4 = true } }, .nterm = &.{3} },
        .{ .params = .{ .fixes = .{ .q4 = true }, .break_termini = true }, .nterm = &.{ 3, 4 } },
    };
    for (cases) |c| {
        var conf = try build(testing.allocator, &mdl, null, null, c.params);
        defer conf.deinit();
        var found: [4]u32 = undefined;
        var n: usize = 0;
        for (conf.groups.items) |g| {
            if (g.gtype != .nterm) continue;
            found[n] = g.residue_idx;
            n += 1;
        }
        try testing.expectEqualSlices(u32, c.nterm, found[0..n]);
    }
    // Q3 with break termini: the residue before the break is a C-terminus.
    var conf = try build(testing.allocator, &mdl, null, null, .{ .fixes = .{ .q3 = true }, .break_termini = true });
    defer conf.deinit();
    var cterm: [4]u32 = undefined;
    var n: usize = 0;
    for (conf.groups.items) |g| {
        if (g.gtype != .coo) continue;
        cterm[n] = g.residue_idx;
        n += 1;
    }
    try testing.expectEqualSlices(u32, &.{ 2, 3, 4 }, cterm[0..n]);
}

test "Q5: monatomic metal ions by element" {
    const ala = backboneAtoms(0);
    const zn = [_]AtomSpec{.{ .name = "ZN", .element = .Zn, .x = 10 }};
    const cluster = [_]AtomSpec{
        .{ .name = "FE1", .element = .Fe, .x = 20 },
        .{ .name = "S1", .element = .S, .x = 22 },
    };
    const ca = [_]AtomSpec{.{ .name = "CA", .element = .Ca, .x = 30 }};
    // "XZN" is not in the ion table; "CA" is (calcium, +2).
    var mdl = try specModel(&.{
        .{ .comp = "ALA", .num = 1, .atoms = &ala },
        .{ .comp = "XZN", .num = 101, .entity = .non_polymer, .atoms = &zn },
        .{ .comp = "FES", .num = 102, .entity = .non_polymer, .atoms = &cluster },
        .{ .comp = "CA", .num = 103, .entity = .non_polymer, .atoms = &ca },
    }, &.{4});
    defer mdl.deinit();
    for ([_]bool{ false, true }) |q5| {
        var conf = try build(testing.allocator, &mdl, null, null, .{ .fixes = .{ .q5 = q5 } });
        defer conf.deinit();
        var charges: [4]f64 = undefined;
        var n: usize = 0;
        for (conf.groups.items) |g| {
            if (g.gtype != .ion) continue;
            charges[n] = g.q;
            n += 1;
        }
        // Parity: CA by name only. Q5: XZN (+2 by element) and CA; the
        // two-atom FES cluster is not an ion.
        if (q5) try testing.expectEqualSlices(f64, &.{ 2, 2 }, charges[0..n]) else try testing.expectEqualSlices(f64, &.{2}, charges[0..n]);
    }
}

test "Q6: zreduce waters are removed in every mode" {
    const ala = backboneAtoms(0);
    const wat = [_]AtomSpec{.{ .name = "O", .element = .O, .x = 3 }};
    var mdl = try specModel(&.{
        .{ .comp = "ALA", .num = 1, .atoms = &ala },
        .{ .comp = "SOL", .num = 201, .entity = .water, .atoms = &wat },
    }, &.{2});
    defer mdl.deinit();
    for ([_]bool{ false, true }) |q6| {
        var conf = try build(testing.allocator, &mdl, null, null, .{ .fixes = .{ .q6 = q6 } });
        defer conf.deinit();
        try testing.expectEqual(@as(usize, 4), conf.atoms.items.len);
    }
}

test "Q8: insertion-code variants stay present in later conformations" {
    const asp = [_]AtomSpec{
        .{ .name = "CG", .element = .C, .x = 0 },
        .{ .name = "OD1", .element = .O, .x = 0.6, .y = 1.1 },
        .{ .name = "OD2", .element = .O, .x = 0.6, .y = -1.1 },
    };
    const ala = [_]AtomSpec{
        .{ .name = "CA", .element = .C, .x = 10 },
        .{ .name = "CB", .element = .C, .x = 11.5, .altloc = 'A' },
        .{ .name = "CB", .element = .C, .x = 11.5, .y = 0.5, .altloc = 'B' },
    };
    // ASP 5 and ALA 5A share chain and number; only ALA has alternates.
    var mdl = try specModel(&.{
        .{ .comp = "ASP", .num = 5, .atoms = &asp },
        .{ .comp = "ALA", .num = 5, .ins = 'A', .atoms = &ala },
    }, &.{2});
    defer mdl.deinit();
    for ([_]bool{ false, true }) |q8| {
        const p = Parameters{ .fixes = .{ .q8 = q8 } };
        var conf_a = try build(testing.allocator, &mdl, 'A', 'A', p);
        defer conf_a.deinit();
        try testing.expect(findGroup(&conf_a, .coo) != null);
        // Parity: in conformation B only the last residue of the (chain, 5)
        // family is present, so the ASP is absent.
        var conf_b = try build(testing.allocator, &mdl, 'B', 'A', p);
        defer conf_b.deinit();
        try testing.expectEqual(q8, findGroup(&conf_b, .coo) != null);
    }
}
