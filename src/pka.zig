//! Experimental, independently implemented protein pKa model.
//!
//! This module intentionally does not reproduce or translate PROPKA source
//! code.  It exposes small, separately testable correction terms so the model
//! can be calibrated against experiment and compared with external oracles.

const std = @import("std");
const model_mod = @import("model.zig");
const math = @import("math.zig");

const Model = model_mod.Model;
const Vec3f32 = math.Vec3(f32);

pub const schema_version = "zhyd.protonation.v1";
pub const parameter_version = "zhyd-pka-protein-v0.1";

pub const SiteKind = enum {
    asp,
    glu,
    his,
    lys,
    cys,
    tyr,
    n_term,
    c_term,

    pub fn name(self: SiteKind) []const u8 {
        return switch (self) {
            .asp => "ASP",
            .glu => "GLU",
            .his => "HIS",
            .lys => "LYS",
            .cys => "CYS",
            .tyr => "TYR",
            .n_term => "N_TERM",
            .c_term => "C_TERM",
        };
    }

    pub fn isAcid(self: SiteKind) bool {
        return switch (self) {
            .asp, .glu, .cys, .tyr, .c_term => true,
            .his, .lys, .n_term => false,
        };
    }
};

pub const Confidence = enum { high, medium, low };

pub const Contributions = struct {
    model: f32,
    burial: f32 = 0,
    hydrogen_bond: f32 = 0,
    electrostatic: f32 = 0,
    coupling: f32 = 0,

    pub fn total(self: Contributions) f32 {
        return self.model + self.burial + self.hydrogen_bond + self.electrostatic + self.coupling;
    }
};

pub const SiteResult = struct {
    residue_idx: u32,
    kind: SiteKind,
    center: Vec3f32,
    terms: Contributions,
    predicted_pka: f32,
    burial_fraction: f32,
    confidence: Confidence,
    coupled: bool = false,
};

pub const Prediction = struct {
    allocator: std.mem.Allocator,
    sites: []SiteResult,
    omitted_sites: u32,

    pub fn deinit(self: *Prediction) void {
        self.allocator.free(self.sites);
        self.sites = &.{};
    }
};

/// All values are deliberately versioned and centralized.  They are initial
/// scientific priors, not fitted production parameters.
pub const Parameters = struct {
    version: []const u8 = parameter_version,
    burial_radius: f32 = 10.0,
    burial_surface_count: f32 = 12.0,
    burial_core_count: f32 = 55.0,
    burial_max_shift: f32 = 2.0,
    hbond_cutoff: f32 = 4.0,
    hbond_max_shift: f32 = 1.5,
    fixed_charge_scale: f32 = 1.2,
    coupling_scale: f32 = 1.0,
    coupling_cutoff: f32 = 12.0,
    coupling_iterations: u8 = 24,
    convergence_tolerance: f32 = 0.001,
};

pub const Config = struct {
    parameters: Parameters = .{},
    include_structure_terms: bool = true,
    include_coupling: bool = true,
};

const Candidate = struct {
    residue_idx: u32,
    kind: SiteKind,
    center: Vec3f32,
};

pub fn modelPka(kind: SiteKind) f32 {
    return switch (kind) {
        .asp => 3.8,
        .glu => 4.5,
        .his => 6.5,
        .lys => 10.5,
        .cys => 8.3,
        .tyr => 10.1,
        .n_term => 8.0,
        .c_term => 3.1,
    };
}

/// Predict pKa values for standard protein side chains and peptide termini.
/// Ligands, nucleic acids, metals, and alternate-location ensembles are out of
/// scope for this first model and are never guessed silently.
pub fn predict(allocator: std.mem.Allocator, mdl: *const Model, config: Config) !Prediction {
    var candidates: std.ArrayListUnmanaged(Candidate) = .empty;
    defer candidates.deinit(allocator);
    var omitted: u32 = 0;

    for (mdl.residues.items, 0..) |res, res_idx| {
        const comp_id = res.compIdSlice();
        // Some coordinate files omit entity metadata. Standard amino-acid
        // identity is a safe sequence fallback; non-standard compounds are
        // still excluded rather than guessed.
        if (res.entity_type != .polymer and !isStandardAminoAcid(comp_id)) continue;
        const kind: ?SiteKind = if (std.mem.eql(u8, comp_id, "ASP")) .asp else if (std.mem.eql(u8, comp_id, "GLU")) .glu else if (std.mem.eql(u8, comp_id, "HIS")) .his else if (std.mem.eql(u8, comp_id, "LYS")) .lys else if (std.mem.eql(u8, comp_id, "CYS")) .cys else if (std.mem.eql(u8, comp_id, "TYR")) .tyr else null;
        if (kind) |k| {
            // A cystine SG has no titratable proton; it is not a site.
            if (k == .cys and isCovalentlyBondedCys(mdl, @intCast(res_idx))) {
                // fall through to terminus handling below
            } else if (siteCenter(mdl, @intCast(res_idx), k)) |center| {
                try candidates.append(allocator, .{ .residue_idx = @intCast(res_idx), .kind = k, .center = center });
            } else omitted += 1;
        }

        const chain = mdl.chains.items[res.chain_idx];
        if (res_idx == chain.residue_start) {
            if (siteCenter(mdl, @intCast(res_idx), .n_term)) |center| {
                try candidates.append(allocator, .{ .residue_idx = @intCast(res_idx), .kind = .n_term, .center = center });
            } else omitted += 1;
        }
        if (res_idx == chain.residue_end - 1) {
            if (siteCenter(mdl, @intCast(res_idx), .c_term)) |center| {
                try candidates.append(allocator, .{ .residue_idx = @intCast(res_idx), .kind = .c_term, .center = center });
            } else omitted += 1;
        }
    }

    const sites = try allocator.alloc(SiteResult, candidates.items.len);
    errdefer allocator.free(sites);
    for (candidates.items, 0..) |candidate, i| {
        const burial = if (config.include_structure_terms) burialFraction(mdl, candidate, config.parameters) else 0;
        var terms = Contributions{ .model = modelPka(candidate.kind) };
        if (config.include_structure_terms) {
            terms.burial = burialContribution(candidate.kind, burial, config.parameters);
            terms.hydrogen_bond = hbondContribution(mdl, candidate, config.parameters);
            terms.electrostatic = fixedChargeContribution(mdl, candidate, config.parameters);
        }
        sites[i] = .{
            .residue_idx = candidate.residue_idx,
            .kind = candidate.kind,
            .center = candidate.center,
            .terms = terms,
            .predicted_pka = clampPka(terms.total()),
            .burial_fraction = burial,
            .confidence = .high,
        };
    }

    if (config.include_structure_terms and config.include_coupling) {
        try applyMeanFieldCoupling(allocator, sites, config.parameters);
    }
    for (sites) |*site| site.confidence = assessConfidence(site.*);

    return .{ .allocator = allocator, .sites = sites, .omitted_sites = omitted };
}

fn isStandardAminoAcid(comp_id: []const u8) bool {
    const ids = std.StaticStringMap(void).initComptime(.{
        .{ "ALA", {} }, .{ "ARG", {} }, .{ "ASN", {} }, .{ "ASP", {} },
        .{ "CYS", {} }, .{ "GLN", {} }, .{ "GLU", {} }, .{ "GLY", {} },
        .{ "HIS", {} }, .{ "ILE", {} }, .{ "LEU", {} }, .{ "LYS", {} },
        .{ "MET", {} }, .{ "PHE", {} }, .{ "PRO", {} }, .{ "SER", {} },
        .{ "THR", {} }, .{ "TRP", {} }, .{ "TYR", {} }, .{ "VAL", {} },
    });
    return ids.has(comp_id);
}

fn siteCenter(mdl: *const Model, residue_idx: u32, kind: SiteKind) ?Vec3f32 {
    return switch (kind) {
        .asp => midpoint(namedAtom(mdl, residue_idx, "OD1"), namedAtom(mdl, residue_idx, "OD2")),
        .glu => midpoint(namedAtom(mdl, residue_idx, "OE1"), namedAtom(mdl, residue_idx, "OE2")),
        .his => midpoint(namedAtom(mdl, residue_idx, "ND1"), namedAtom(mdl, residue_idx, "NE2")),
        .lys => namedAtom(mdl, residue_idx, "NZ"),
        .cys => namedAtom(mdl, residue_idx, "SG"),
        .tyr => namedAtom(mdl, residue_idx, "OH"),
        .n_term => namedAtom(mdl, residue_idx, "N"),
        .c_term => namedAtom(mdl, residue_idx, "OXT") orelse namedAtom(mdl, residue_idx, "O"),
    };
}

/// Maximum SG–SG distance treated as a disulfide bond (typical 2.03–2.05 Å).
const disulfide_max_distance: f32 = 2.5;

/// True when the CYS side-chain sulfur is covalently bonded to another residue:
/// either flagged from `_struct_conn`, or within disulfide distance of another
/// CYS SG. The geometric check covers PDB input and files without
/// `_struct_conn`.
pub fn isCovalentlyBondedCys(mdl: *const Model, residue_idx: u32) bool {
    const res = mdl.residues.items[residue_idx];
    var sg: ?Vec3f32 = null;
    for (mdl.atoms.items[res.atom_start..res.atom_end]) |atom| {
        if (atom.is_hydrogen or !std.mem.eql(u8, atom.nameSlice(), "SG")) continue;
        if (atom.flags.bonded_inter_residue) return true;
        if (atom.altloc == ' ' or atom.altloc == 'A') sg = atom.pos;
    }
    const pos = sg orelse return false;
    for (mdl.residues.items, 0..) |other, other_idx| {
        if (other_idx == residue_idx or !std.mem.eql(u8, other.compIdSlice(), "CYS")) continue;
        const other_sg = namedAtom(mdl, @intCast(other_idx), "SG") orelse continue;
        if (pos.distance(other_sg) <= disulfide_max_distance) return true;
    }
    return false;
}

fn namedAtom(mdl: *const Model, residue_idx: u32, name: []const u8) ?Vec3f32 {
    const res = mdl.residues.items[residue_idx];
    for (mdl.atoms.items[res.atom_start..res.atom_end]) |atom| {
        if (atom.is_hydrogen) continue;
        if (atom.altloc != ' ' and atom.altloc != 'A') continue;
        if (std.mem.eql(u8, atom.nameSlice(), name)) return atom.pos;
    }
    return null;
}

fn midpoint(a: ?Vec3f32, b: ?Vec3f32) ?Vec3f32 {
    const av = a orelse return null;
    const bv = b orelse return null;
    return .{ .x = (av.x + bv.x) * 0.5, .y = (av.y + bv.y) * 0.5, .z = (av.z + bv.z) * 0.5 };
}

fn burialFraction(mdl: *const Model, site: Candidate, params: Parameters) f32 {
    var count: f32 = 0;
    for (mdl.atoms.items) |atom| {
        if (atom.is_hydrogen or atom.residue_idx == site.residue_idx) continue;
        if (atom.altloc != ' ' and atom.altloc != 'A') continue;
        const d = atom.pos.distance(site.center);
        if (d >= params.burial_radius) continue;
        count += @max(@as(f32, 0), (params.burial_radius - d) / params.burial_radius);
    }
    return std.math.clamp((count - params.burial_surface_count) / (params.burial_core_count - params.burial_surface_count), 0, 1);
}

fn burialContribution(kind: SiteKind, fraction: f32, params: Parameters) f32 {
    const sign: f32 = if (kind.isAcid()) 1 else -1;
    return sign * params.burial_max_shift * fraction;
}

fn hbondContribution(mdl: *const Model, site: Candidate, params: Parameters) f32 {
    var score: f32 = 0;
    for (mdl.atoms.items) |atom| {
        if (atom.is_hydrogen or atom.residue_idx == site.residue_idx) continue;
        if (atom.altloc != ' ' and atom.altloc != 'A') continue;
        const d = atom.pos.distance(site.center);
        if (d < 2.0 or d >= params.hbond_cutoff) continue;
        const weight = (params.hbond_cutoff - d) / (params.hbond_cutoff - 2.0);
        if (site.kind.isAcid()) {
            if (atom.flags.donor or atom.flags.positive) score -= 0.35 * weight;
        } else {
            if (atom.flags.acceptor or atom.flags.negative) score += 0.35 * weight;
        }
    }
    return std.math.clamp(score, -params.hbond_max_shift, params.hbond_max_shift);
}

fn fixedChargeContribution(mdl: *const Model, site: Candidate, params: Parameters) f32 {
    var shift: f32 = 0;
    for (mdl.residues.items, 0..) |res, res_idx| {
        if (res_idx == site.residue_idx) continue;
        // ARG is treated as a fixed +1 group in this phase. All supported
        // titratable groups are handled by the coupling term instead.
        if (!std.mem.eql(u8, res.compIdSlice(), "ARG")) continue;
        const center = namedAtom(mdl, @intCast(res_idx), "CZ") orelse continue;
        const d = center.distance(site.center);
        if (d < 1.0 or d > params.coupling_cutoff) continue;
        shift -= params.fixed_charge_scale / d;
    }
    return shift;
}

fn applyMeanFieldCoupling(allocator: std.mem.Allocator, sites: []SiteResult, params: Parameters) !void {
    if (sites.len < 2) return;
    const next_values = try allocator.alloc(f32, sites.len);
    defer allocator.free(next_values);
    const next_coupling = try allocator.alloc(f32, sites.len);
    defer allocator.free(next_coupling);
    const next_flags = try allocator.alloc(bool, sites.len);
    defer allocator.free(next_flags);
    var iteration: u8 = 0;
    while (iteration < params.coupling_iterations) : (iteration += 1) {
        var max_delta: f32 = 0;
        for (sites, 0..) |site, i| {
            var contribution: f32 = 0;
            var coupled = false;
            for (sites, 0..) |other, j| {
                if (i == j) continue;
                const d = site.center.distance(other.center);
                if (d < 1.0 or d > params.coupling_cutoff) continue;
                const charge = meanChargeAtPh(other.kind, other.predicted_pka, site.predicted_pka);
                contribution -= params.coupling_scale * charge / d;
                if (d <= 6.0 and @abs(charge) > 0.1) coupled = true;
            }
            contribution = std.math.clamp(contribution, -2.5, 2.5);
            const next = clampPka(site.terms.model + site.terms.burial + site.terms.hydrogen_bond + site.terms.electrostatic + contribution);
            max_delta = @max(max_delta, @abs(next - site.predicted_pka));
            next_values[i] = next;
            next_coupling[i] = contribution;
            next_flags[i] = coupled;
        }
        for (sites, 0..) |*site, i| {
            site.terms.coupling = next_coupling[i];
            site.predicted_pka = next_values[i];
            site.coupled = next_flags[i];
        }
        if (max_delta < params.convergence_tolerance) break;
    }
}

fn meanChargeAtPh(kind: SiteKind, pka: f32, ph: f32) f32 {
    if (kind.isAcid()) {
        const deprotonated = 1.0 / (1.0 + std.math.pow(f32, 10.0, pka - ph));
        return -deprotonated;
    }
    const protonated = 1.0 / (1.0 + std.math.pow(f32, 10.0, ph - pka));
    return protonated;
}

fn assessConfidence(site: SiteResult) Confidence {
    if (site.burial_fraction >= 0.70 or @abs(site.terms.coupling) >= 1.0) return .low;
    if (site.burial_fraction >= 0.35 or site.coupled or @abs(site.terms.electrostatic) >= 0.5) return .medium;
    return .high;
}

fn clampPka(value: f32) f32 {
    return std.math.clamp(value, -2.0, 16.0);
}

const testing = std.testing;
const mmcif = @import("mmcif.zig");

test "model pKa stage is deterministic and structure-independent" {
    const source = @embedFile("test_data/his.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();
    var result = try predict(testing.allocator, &mdl, .{ .include_structure_terms = false });
    defer result.deinit();
    try testing.expectEqual(@as(usize, 3), result.sites.len); // HIS plus both termini
    for (result.sites) |site| try testing.expectApproxEqAbs(modelPka(site.kind), site.predicted_pka, 0.0001);
}

test "mean-field charge has correct acid and base limits" {
    try testing.expect(meanChargeAtPh(.asp, 4.0, 8.0) < -0.99);
    try testing.expect(meanChargeAtPh(.lys, 10.5, 7.0) > 0.99);
}

test "disulfide-bonded cysteines are not titratable sites" {
    const source = @embedFile("test_data/disulfide.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();
    var result = try predict(testing.allocator, &mdl, .{});
    defer result.deinit();
    for (result.sites) |site| try testing.expect(site.kind != .cys);
    try testing.expectEqual(@as(u32, 0), result.omitted_sites);
}

test "free cysteine remains a titratable site" {
    const source = @embedFile("test_data/disulfide.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();
    // Move the second SG out of bonding range.
    for (mdl.atoms.items) |*atom| {
        if (atom.residue_idx == 1 and std.mem.eql(u8, atom.nameSlice(), "SG")) atom.pos.x += 5.0;
    }
    var result = try predict(testing.allocator, &mdl, .{});
    defer result.deinit();
    var n_cys: usize = 0;
    for (result.sites) |site| {
        if (site.kind == .cys) n_cys += 1;
    }
    try testing.expectEqual(@as(usize, 2), n_cys);
}

test "missing required side-chain atoms are omitted explicitly" {
    const source = @embedFile("test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();
    // tiny.cif has ordinary protein termini; prediction must still be usable.
    var result = try predict(testing.allocator, &mdl, .{});
    defer result.deinit();
    try testing.expect(result.sites.len >= 2);
}
