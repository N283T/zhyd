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
pub const parameter_version = "zhyd-pka-protein-v0.2";
pub const engine_name = "zhyd-builtin";

pub const SiteKind = enum {
    asp,
    glu,
    his,
    lys,
    /// Only produced by the PROPKA engine; the built-in model treats ARG as a
    /// fixed +1 charge.
    arg,
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
            .arg => "ARG",
            .cys => "CYS",
            .tyr => "TYR",
            .n_term => "N_TERM",
            .c_term => "C_TERM",
        };
    }

    pub fn isAcid(self: SiteKind) bool {
        return switch (self) {
            .asp, .glu, .cys, .tyr, .c_term => true,
            .his, .lys, .arg, .n_term => false,
        };
    }
};

pub const Confidence = enum { high, medium, low };

/// Why a titratable residue was deliberately not given a pKa.
pub const ExclusionReason = enum {
    /// Cystine sulfur: no titratable proton.
    disulfide,
    /// Side-chain ligand of a metal ion; the metal-bound state is outside
    /// this model and must not be guessed from protein-only terms.
    metal_coordinated,
    /// SG covalently bonded to a non-cysteine residue or ligand (e.g. a
    /// covalent inhibitor or modification), or flagged covalent by
    /// `_struct_conn` without a resolvable partner.
    covalently_bonded,
};

pub const ExcludedSite = struct {
    residue_idx: u32,
    kind: SiteKind,
    reason: ExclusionReason,
};

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
    excluded: []ExcludedSite = &.{},
    omitted_sites: u32,
    /// Version string of the parameter set that produced these values.
    parameter_version: []const u8 = parameter_version,
    /// Which pKa engine produced these values.
    engine: []const u8 = engine_name,

    pub fn deinit(self: *Prediction) void {
        self.allocator.free(self.sites);
        self.allocator.free(self.excluded);
        self.sites = &.{};
        self.excluded = &.{};
    }
};

/// All values are deliberately versioned and centralized.
///
/// v0.2: `burial_max_shift`, `hbond_strength` and `hbond_max_shift` were fitted
/// by coordinate descent on the PKAD-R training split only (train MAE; see
/// scripts/pka/params/zhyd-pka-protein-v0.2.json). The fit removed the burial
/// shift: with the current neighbor-count descriptor it did not improve
/// training error. Burial is still computed for confidence and strata. This
/// is a known weakness for deeply buried groups (SNase stress set).
/// v0.1 priors were burial 2.0, hbond strength 0.35, hbond max 1.5.
pub const Parameters = struct {
    version: []const u8 = parameter_version,
    burial_radius: f32 = 10.0,
    burial_surface_count: f32 = 12.0,
    burial_core_count: f32 = 55.0,
    burial_max_shift: f32 = 0.0,
    hbond_cutoff: f32 = 4.0,
    hbond_max_shift: f32 = 0.75,
    hbond_strength: f32 = 0.5,
    fixed_charge_scale: f32 = 1.2,
    coupling_scale: f32 = 1.0,
    coupling_cutoff: f32 = 12.0,
    coupling_iterations: u8 = 24,
    convergence_tolerance: f32 = 0.001,
};

/// Parse a parameter set from JSON. Missing fields keep their defaults, unknown
/// fields are rejected so typos cannot silently fall back to defaults. A custom
/// set must carry its own `version` so reports never mislabel it.
pub fn parseParameters(allocator: std.mem.Allocator, json: []const u8) !std.json.Parsed(Parameters) {
    const parsed = try std.json.parseFromSlice(Parameters, allocator, json, .{ .ignore_unknown_fields = false, .allocate = .alloc_always });
    errdefer parsed.deinit();
    if (std.mem.eql(u8, parsed.value.version, parameter_version)) {
        // Only the built-in defaults may use the built-in version label.
        if (!sameValues(parsed.value, .{})) return error.ParameterVersionRequired;
    }
    return parsed;
}

/// Compare every numeric field, ignoring the version label.
fn sameValues(a: Parameters, b: Parameters) bool {
    inline for (std.meta.fields(Parameters)) |field| {
        if (comptime std.mem.eql(u8, field.name, "version")) continue;
        if (@field(a, field.name) != @field(b, field.name)) return false;
    }
    return true;
}

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
        .arg => 12.5,
        .cys => 8.3,
        .tyr => 10.1,
        .n_term => 8.0,
        .c_term => 3.1,
    };
}

/// Predict pKa values for standard protein side chains and peptide termini.
/// Ligands, nucleic acids, metals, and alternate-location ensembles are out of
/// scope for this first model and are never guessed silently: disulfide and
/// metal-coordinated residues are listed in `excluded` with a reason.
pub fn predict(allocator: std.mem.Allocator, mdl: *const Model, config: Config) !Prediction {
    var candidates: std.ArrayListUnmanaged(Candidate) = .empty;
    defer candidates.deinit(allocator);
    var excluded: std.ArrayListUnmanaged(ExcludedSite) = .empty;
    errdefer excluded.deinit(allocator);
    var omitted: u32 = 0;
    const metals = try metalPositions(allocator, mdl);
    defer allocator.free(metals);

    for (mdl.residues.items, 0..) |res, res_idx| {
        const comp_id = res.compIdSlice();
        // Some coordinate files omit entity metadata. Standard amino-acid
        // identity is a safe sequence fallback; non-standard compounds are
        // still excluded rather than guessed.
        if (res.entity_type != .polymer and !isStandardAminoAcid(comp_id)) continue;
        const kind: ?SiteKind = if (std.mem.eql(u8, comp_id, "ASP")) .asp else if (std.mem.eql(u8, comp_id, "GLU")) .glu else if (std.mem.eql(u8, comp_id, "HIS")) .his else if (std.mem.eql(u8, comp_id, "LYS")) .lys else if (std.mem.eql(u8, comp_id, "CYS")) .cys else if (std.mem.eql(u8, comp_id, "TYR")) .tyr else null;
        if (kind) |k| {
            if (exclusionReason(mdl, @intCast(res_idx), k, metals)) |reason| {
                try excluded.append(allocator, .{ .residue_idx = @intCast(res_idx), .kind = k, .reason = reason });
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

    const excluded_slice = try excluded.toOwnedSlice(allocator);
    return .{ .allocator = allocator, .sites = sites, .excluded = excluded_slice, .omitted_sites = omitted, .parameter_version = config.parameters.version };
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
        .arg => namedAtom(mdl, residue_idx, "CZ"),
        .cys => namedAtom(mdl, residue_idx, "SG"),
        .tyr => namedAtom(mdl, residue_idx, "OH"),
        .n_term => namedAtom(mdl, residue_idx, "N"),
        .c_term => namedAtom(mdl, residue_idx, "OXT") orelse namedAtom(mdl, residue_idx, "O"),
    };
}

/// Maximum ligand-atom–metal distance treated as coordination. Typical
/// protein metal-ligand bonds are 1.9–2.6 Å; 2.8 Å tolerates coordinate error.
const metal_coordination_max_distance: f32 = 2.8;

pub fn metalPositions(allocator: std.mem.Allocator, mdl: *const Model) ![]Vec3f32 {
    var list: std.ArrayListUnmanaged(Vec3f32) = .empty;
    errdefer list.deinit(allocator);
    for (mdl.atoms.items) |atom| {
        if (atom.element_type.info().flags.metallic or atom.flags.metallic) try list.append(allocator, atom.pos);
    }
    return list.toOwnedSlice(allocator);
}

/// Side-chain atoms that can coordinate a metal for each site kind.
fn ligandAtomNames(kind: SiteKind) []const []const u8 {
    return switch (kind) {
        .asp => &.{ "OD1", "OD2" },
        .glu => &.{ "OE1", "OE2" },
        .his => &.{ "ND1", "NE2" },
        .lys => &.{"NZ"},
        .arg => &.{ "NE", "NH1", "NH2" },
        .cys => &.{"SG"},
        .tyr => &.{"OH"},
        .n_term, .c_term => &.{},
    };
}

pub fn isMetalCoordinated(mdl: *const Model, residue_idx: u32, kind: SiteKind, metals: []const Vec3f32) bool {
    for (ligandAtomNames(kind)) |name| {
        const pos = namedAtom(mdl, residue_idx, name) orelse continue;
        for (metals) |metal| {
            if (pos.distance(metal) <= metal_coordination_max_distance) return true;
        }
    }
    return false;
}

pub fn exclusionReason(mdl: *const Model, residue_idx: u32, kind: SiteKind, metals: []const Vec3f32) ?ExclusionReason {
    // Metal first: `_struct_conn` flags metalc bonds too, and a metal ligand
    // must not be mistaken for a disulfide.
    if (isMetalCoordinated(mdl, residue_idx, kind, metals)) return .metal_coordinated;
    if (kind == .cys) return cysteineBond(mdl, residue_idx);
    return null;
}

/// Maximum SG–heavy-atom distance treated as a covalent bond to a
/// non-cysteine partner (S–C ≈ 1.8 Å).
const covalent_max_distance: f32 = 2.1;

/// Classify a non-metal covalent bond on a CYS SG by its nearest partner.
/// Callers check metal coordination first.
pub fn cysteineBond(mdl: *const Model, residue_idx: u32) ?ExclusionReason {
    const res = mdl.residues.items[residue_idx];
    var sg: ?Vec3f32 = null;
    var flagged = false;
    for (mdl.atoms.items[res.atom_start..res.atom_end]) |atom| {
        if (atom.is_hydrogen or !std.mem.eql(u8, atom.nameSlice(), "SG")) continue;
        if (atom.flags.bonded_inter_residue) flagged = true;
        if (sg == null and (atom.altloc == ' ' or atom.altloc == 'A')) sg = atom.pos;
    }
    const pos = sg orelse return if (flagged) .covalently_bonded else null;
    var covalent_partner = false;
    for (mdl.atoms.items) |atom| {
        if (atom.is_hydrogen or atom.residue_idx == residue_idx) continue;
        if (atom.altloc != ' ' and atom.altloc != 'A') continue;
        if (atom.element_type.info().flags.metallic or atom.flags.metallic) continue;
        const d = pos.distance(atom.pos);
        if (d > disulfide_max_distance) continue;
        const other = mdl.residues.items[atom.residue_idx];
        if (std.mem.eql(u8, atom.nameSlice(), "SG") and std.mem.eql(u8, other.compIdSlice(), "CYS")) return .disulfide;
        if (d <= covalent_max_distance and !isSolvent(mdl, atom.residue_idx)) covalent_partner = true;
    }
    if (covalent_partner) return .covalently_bonded;
    // A flag without a partner in this model (e.g. a symmetry-mate
    // disulfide) is still covalent; do not guess which kind.
    return if (flagged) .covalently_bonded else null;
}

/// Maximum SG–SG distance treated as a disulfide bond (typical 2.03–2.05 Å).
const disulfide_max_distance: f32 = 2.5;

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

/// Water is the reference environment of the model pKa values, so explicit
/// solvent must not count as burial or as a protein hydrogen-bond partner.
fn isSolvent(mdl: *const Model, residue_idx: u32) bool {
    const res = mdl.residues.items[residue_idx];
    if (res.entity_type == .water) return true;
    const comp = res.compIdSlice();
    return std.mem.eql(u8, comp, "HOH") or std.mem.eql(u8, comp, "WAT") or
        std.mem.eql(u8, comp, "DOD") or std.mem.eql(u8, comp, "H2O");
}

fn burialFraction(mdl: *const Model, site: Candidate, params: Parameters) f32 {
    var count: f32 = 0;
    for (mdl.atoms.items) |atom| {
        if (atom.is_hydrogen or atom.residue_idx == site.residue_idx) continue;
        if (isSolvent(mdl, atom.residue_idx)) continue;
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
        if (isSolvent(mdl, atom.residue_idx)) continue;
        if (atom.altloc != ' ' and atom.altloc != 'A') continue;
        const d = atom.pos.distance(site.center);
        if (d < 2.0 or d >= params.hbond_cutoff) continue;
        const weight = (params.hbond_cutoff - d) / (params.hbond_cutoff - 2.0);
        if (site.kind.isAcid()) {
            if (atom.flags.donor or atom.flags.positive) score -= params.hbond_strength * weight;
        } else {
            if (atom.flags.acceptor or atom.flags.negative) score += params.hbond_strength * weight;
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

pub fn assessConfidence(site: SiteResult) Confidence {
    if (site.burial_fraction >= 0.70 or @abs(site.terms.coupling) >= 1.0) return .low;
    if (site.burial_fraction >= 0.35 or site.coupled or @abs(site.terms.electrostatic) >= 0.5) return .medium;
    return .high;
}

/// Values outside this range are not credible for protein residues.
pub fn isPhysicalPka(value: f32) bool {
    return value >= -5.0 and value <= 20.0;
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
    try testing.expectEqual(@as(usize, 2), result.excluded.len);
    for (result.excluded) |site| try testing.expectEqual(ExclusionReason.disulfide, site.reason);
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

test "parameter files keep defaults and require a custom version" {
    var parsed = try parseParameters(testing.allocator, "{\"version\": \"custom-1\", \"hbond_strength\": 0.9}");
    defer parsed.deinit();
    try testing.expectEqualStrings("custom-1", parsed.value.version);
    try testing.expectApproxEqAbs(@as(f32, 0.9), parsed.value.hbond_strength, 1e-6);
    const defaults = Parameters{};
    try testing.expectApproxEqAbs(defaults.burial_radius, parsed.value.burial_radius, 1e-6);
    try testing.expectError(error.ParameterVersionRequired, parseParameters(testing.allocator, "{\"hbond_strength\": 0.9}"));

    // The parsed version must outlive the caller's JSON buffer.
    const buffer = try testing.allocator.dupe(u8, "{\"version\": \"heap-1\"}");
    var owned = try parseParameters(testing.allocator, buffer);
    defer owned.deinit();
    @memset(buffer, 0);
    testing.allocator.free(buffer);
    try testing.expectEqualStrings("heap-1", owned.value.version);
    try testing.expectError(error.UnknownField, parseParameters(testing.allocator, "{\"version\": \"x\", \"hbond_strenght\": 0.5}"));
}

test "metal-coordinated cysteine is excluded with a reason" {
    const source = @embedFile("test_data/disulfide.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();
    // Break the disulfide, then place a zinc 2.3 Å from the first SG.
    var first_sg: Vec3f32 = undefined;
    for (mdl.atoms.items) |*atom| {
        if (!std.mem.eql(u8, atom.nameSlice(), "SG")) continue;
        if (atom.residue_idx == 0) first_sg = atom.pos else atom.pos.x += 5.0;
    }
    var zn = mdl.atoms.items[0];
    zn.pos = .{ .x = first_sg.x, .y = first_sg.y + 2.3, .z = first_sg.z };
    zn.element_type = .Zn;
    zn.residue_idx = 1;
    try mdl.atoms.append(testing.allocator, zn);
    // struct_conn also flags metalc bonds; the flag must not turn a metal
    // ligand into a "disulfide".
    for (mdl.atoms.items) |*atom| {
        if (atom.residue_idx == 0 and std.mem.eql(u8, atom.nameSlice(), "SG")) atom.flags.bonded_inter_residue = true;
    }
    var result = try predict(testing.allocator, &mdl, .{});
    defer result.deinit();
    try testing.expectEqual(@as(usize, 1), result.excluded.len);
    try testing.expectEqual(ExclusionReason.metal_coordinated, result.excluded[0].reason);
    try testing.expectEqual(@as(u32, 0), result.excluded[0].residue_idx);
}

test "cysteine bonded to a non-cysteine heavy atom is excluded as covalently bonded" {
    const source = @embedFile("test_data/disulfide.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();
    var first_sg: Vec3f32 = undefined;
    for (mdl.atoms.items) |*atom| {
        if (!std.mem.eql(u8, atom.nameSlice(), "SG")) continue;
        if (atom.residue_idx == 0) first_sg = atom.pos else atom.pos.x += 5.0;
    }
    // A ligand-like carbon 1.8 Å from the first SG (S–C single bond).
    var c = mdl.atoms.items[0];
    c.pos = .{ .x = first_sg.x, .y = first_sg.y + 1.8, .z = first_sg.z };
    c.element_type = .C;
    c.setName("C1");
    c.residue_idx = 1;
    try mdl.atoms.append(testing.allocator, c);
    var result = try predict(testing.allocator, &mdl, .{});
    defer result.deinit();
    try testing.expectEqual(@as(usize, 1), result.excluded.len);
    try testing.expectEqual(ExclusionReason.covalently_bonded, result.excluded[0].reason);
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
