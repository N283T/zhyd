//! Versioned parameter set of the empirical pKa engine.
//!
//! Values restate the published PROPKA 3 empirical model as described in
//! docs/specs/zhyd-pka-empirical-engine.md (section numbers in comments refer
//! to that document). Nothing here is derived from PROPKA source code.

const std = @import("std");
const pka = @import("../pka.zig");
const AtomType = @import("../element.zig").AtomType;

pub const parity_version = "zhyd-empirical-v1-parity";
/// All adopted corrections (Q1–Q8) on the parity parameter table.
pub const zhyd_v1_version = "zhyd-empirical-v1";
/// All adopted corrections on the recalibrated table (`calibrated`).
pub const zhyd_version = "zhyd-empirical-v1.3";

/// zhyd-mode corrections of the parity quirks ("Known parity quirks"), one
/// switch each. All off reproduces parity mode exactly.
pub const Fixes = packed struct(u8) {
    /// Same-residue exclusion compares the full residue identity.
    q1: bool = false,
    /// F_angle = max(0, cos).
    q2: bool = false,
    /// The chain-last amino-acid residue is a C-terminus without OXT.
    q3: bool = false,
    /// N-termini follow zhyd chains and chain-break detection.
    q4: bool = false,
    /// Monatomic metal ions by element.
    q5: bool = false,
    /// Waters removed by zhyd's classification (already part of the
    /// parity input, so no further effect).
    q6: bool = false,
    /// Conformation mean over the conformations where the site titrates.
    q7: bool = false,
    /// Insertion codes: full residue identity in conformation matching, no
    /// last-of-family shadowing, absent conformations do not count.
    q8: bool = false,

    pub const none: Fixes = .{};
    pub const names = [_][]const u8{ "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Q8" };

    pub fn bits(self: Fixes) u8 {
        return @bitCast(self);
    }

    pub fn fromBits(b: u8) Fixes {
        return @bitCast(b);
    }

    pub fn isNone(self: Fixes) bool {
        return self.bits() == 0;
    }

    /// Parse a comma-separated list such as "Q2,Q7" (case-insensitive) or
    /// "none".
    pub fn parse(list: []const u8) error{InvalidFix}!Fixes {
        if (std.ascii.eqlIgnoreCase(list, "none")) return .none;
        var b: u8 = 0;
        var it = std.mem.splitScalar(u8, list, ',');
        while (it.next()) |raw| {
            const item = std.mem.trim(u8, raw, " ");
            const i = for (names, 0..) |n, i| {
                if (std.ascii.eqlIgnoreCase(n, item)) break i;
            } else return error.InvalidFix;
            b |= @as(u8, 1) << @intCast(i);
        }
        return fromBits(b);
    }
};

/// `base` followed by "+Qn" for every enabled correction ("+none" for the
/// empty set), for every possible correction set (index: `Fixes.bits()`).
fn customVersions(comptime base: []const u8) [256][]const u8 {
    @setEvalBranchQuota(100_000);
    var table: [256][]const u8 = undefined;
    for (0..256) |i| {
        var v: []const u8 = base;
        for (Fixes.names, 0..) |n, k| {
            if ((i >> k) & 1 == 1) v = v ++ "+" ++ n;
        }
        table[i] = if (i == 0) base ++ "+none" else v;
    }
    return table;
}

const v1_versions = customVersions(zhyd_v1_version);
const calibrated_versions = customVersions(zhyd_version);

/// Parameter version of a correction set on the parity table
/// (`calibrated` false) or the recalibrated table. Parity table: parity for
/// no correction, "zhyd-empirical-v1" for the adopted set, otherwise
/// "zhyd-empirical-v1" with the enabled corrections appended (e.g.
/// "zhyd-empirical-v1+Q2+Q7"). Recalibrated table: "zhyd-empirical-v1.3"
/// for the adopted set, otherwise "zhyd-empirical-v1.3+..." ("+none" for
/// no correction).
pub fn versionFor(fixes: Fixes, adopted: Fixes, is_calibrated: bool) []const u8 {
    if (is_calibrated) {
        if (fixes.bits() == adopted.bits()) return zhyd_version;
        return calibrated_versions[fixes.bits()];
    }
    if (fixes.isNone()) return parity_version;
    if (fixes.bits() == adopted.bits()) return zhyd_v1_version;
    return v1_versions[fixes.bits()];
}

/// Per-site-kind parameter values.
pub const KindTable = std.enums.EnumArray(pka.SiteKind, f64);

/// Model pKa values of the parity table ("Model values and charges").
pub const default_model_pka = KindTable.init(.{
    .asp = 3.80,
    .glu = 4.50,
    .his = 6.50,
    .lys = 10.50,
    .arg = 12.50,
    .cys = 9.00,
    .tyr = 10.00,
    .n_term = 8.00,
    .c_term = 3.20,
});

/// Recalibrated zhyd table (`zhyd-empirical-v1.3`): CYS model pKa and
/// desolvation scale and HIS desolvation scale, fitted in order CYS, TYR, HIS
/// from the parity table on PKAD-R merged with the external set (all
/// non-SNase rows) and kept because they lower the nested cluster CV MAE
/// (spec "Recalibration (r7)"; fit record
/// scripts/pka/params/zhyd-empirical-v1.3.json). TYR keeps the parity values:
/// its refit did not beat the default under nested CV. Every other value is
/// the parity table's.
pub const calibrated: Parameters = blk: {
    var p = Parameters{ .version = zhyd_version };
    p.model_pka.set(.cys, 7.5);
    p.desolvation_scale.set(.cys, 0.45);
    p.desolvation_scale.set(.his, 0.6);
    break :blk p;
};

/// Per-kind values in a parameter override file; null keeps the base value.
pub const KindOverrides = struct {
    ASP: ?f64 = null,
    GLU: ?f64 = null,
    HIS: ?f64 = null,
    LYS: ?f64 = null,
    ARG: ?f64 = null,
    CYS: ?f64 = null,
    TYR: ?f64 = null,
    N_TERM: ?f64 = null,
    C_TERM: ?f64 = null,

    fn apply(self: KindOverrides, table: *KindTable) void {
        inline for (std.meta.fields(KindOverrides)) |field| {
            if (@field(self, field.name)) |v| {
                const kind = comptime for (std.enums.values(pka.SiteKind)) |k| {
                    if (std.mem.eql(u8, k.name(), field.name)) break k;
                } else unreachable;
                table.set(kind, v);
            }
        }
    }

    fn valid(self: KindOverrides, min: f64, max: f64) bool {
        inline for (std.meta.fields(KindOverrides)) |field| {
            if (@field(self, field.name)) |v| {
                if (!std.math.isFinite(v) or v < min or v > max) return false;
            }
        }
        return true;
    }
};

/// Parameter override file for zhyd mode (`zhyd pka --parameters`):
/// a required version label plus per-kind model pKa values and desolvation
/// scale factors. Unknown keys are rejected.
pub const Overrides = struct {
    version: []const u8,
    model_pka: KindOverrides = .{},
    desolvation_scale: KindOverrides = .{},

    pub fn applyTo(self: Overrides, base: Parameters) Parameters {
        var p = base;
        p.version = self.version;
        self.model_pka.apply(&p.model_pka);
        self.desolvation_scale.apply(&p.desolvation_scale);
        return p;
    }
};

pub const OverrideError = error{ ParameterVersionRequired, ParameterOutOfRange };

/// Parse and validate an override file. The version must not be one of the
/// built-in labels; model pKa values must lie in −5..20 and scale factors in
/// 0..10.
pub fn parseOverrides(allocator: std.mem.Allocator, json: []const u8) !std.json.Parsed(Overrides) {
    const parsed = try std.json.parseFromSlice(Overrides, allocator, json, .{ .ignore_unknown_fields = false, .allocate = .alloc_always });
    errdefer parsed.deinit();
    const v = parsed.value.version;
    if (v.len == 0 or std.mem.startsWith(u8, v, zhyd_v1_version)) return error.ParameterVersionRequired;
    if (!parsed.value.model_pka.valid(-5.0, 20.0) or !parsed.value.desolvation_scale.valid(0.0, 10.0)) return error.ParameterOutOfRange;
    return parsed;
}

/// Group types of the model ("Group geometry" table).
pub const GroupType = enum {
    coo,
    his,
    cys,
    tyr,
    lys,
    arg,
    nterm,
    roh,
    amd,
    trp,
    bbn,
    bbc,
    ion,

    /// Bases are ARG, LYS, HIS and the N-terminus.
    pub fn isBase(self: GroupType) bool {
        return switch (self) {
            .arg, .lys, .his, .nterm => true,
            else => false,
        };
    }

    /// Groups whose H-bond strength depends on a donor hydrogen direction.
    pub fn isAngular(self: GroupType) bool {
        return switch (self) {
            .his, .arg, .amd, .trp => true,
            else => false,
        };
    }

    pub fn isSideChain(self: GroupType) bool {
        return switch (self) {
            .bbn, .bbc, .ion => false,
            else => true,
        };
    }
};

pub fn groupTypeOf(kind: pka.SiteKind) GroupType {
    return switch (kind) {
        .asp, .glu, .c_term => .coo,
        .his => .his,
        .cys => .cys,
        .tyr => .tyr,
        .lys => .lys,
        .arg => .arg,
        .n_term => .nterm,
    };
}

/// Interaction class of a side-chain pair ("Side-chain pair interactions").
pub const PairClass = enum { none, iterative, non_iterative };

pub const Parameters = struct {
    version: []const u8 = parity_version,
    /// zhyd-mode corrections; none in parity mode.
    fixes: Fixes = .{},
    /// Residues after a chain break are termini (zhyd `--nterm
    /// aggressive`); used by Q3 and Q4 only.
    break_termini: bool = false,

    // Per-site-kind model pKa and desolvation scale factor (multiplies the
    // desolvation shift).
    model_pka: KindTable = default_model_pka,
    desolvation_scale: KindTable = KindTable.initFill(1.0),

    // Internal hydrogens
    nh_length: f64 = 1.01,

    // Desolvation
    desolv_count_radius: f64 = 15.0,
    desolv_volume_radius: f64 = 20.0,
    desolv_min_distance: f64 = 2.75,
    buried_count_min: f64 = 280.0,
    buried_count_max: f64 = 560.0,
    surface_scale: f64 = 0.25,
    desolv_prefactor: f64 = -13.0,
    desolv_volume_offset: f64 = 0.0,

    // Coulomb
    eps_water: f64 = 160.0,
    eps_protein: f64 = 30.0,
    coulomb_min_distance: f64 = 4.0,
    coulomb_cutoff_far: f64 = 10.0,
    coulomb_prefactor: f64 = 244.12,
    coulomb_min_count: f64 = 280.0,

    // H-bonds
    hbond_max: f64 = 0.85,
    angle_min: f64 = 0.001,

    // Reorganisation (ASP/GLU)
    reorg_prefactor: f64 = 0.8,
    reorg_near: f64 = 3.0,
    reorg_far: f64 = 6.0,

    // Pair cut-offs
    pair_cutoff: f64 = 10.0,
    ion_cutoff: f64 = 10.0,
    bond_exclusion: u8 = 4,

    // Buried pair ("not (N1 + N2 <= 900 and (N1 <= 400 or N2 <= 400))")
    buried_pair_sum: f64 = 900.0,
    buried_pair_single: f64 = 400.0,
    coo_his_buried_value: f64 = 1.60,
    cys_his_buried_value: f64 = 1.60,
    cys_cys_buried_value: f64 = 3.60,

    // Iterative scheme
    max_iterations: u32 = 10,
    min_contribution: f64 = 0.005,

    // Disulfide SG–SG distance
    disulfide_distance: f64 = 2.5,
};

/// Value of a site that does not titrate in a conformation (disulfide or
/// covalently coupled); conformation averages include it (amendment A8).
pub const not_titrating_pka: f64 = 99.99;

/// Model pKa of the parity table.
pub fn modelPka(kind: pka.SiteKind) f64 {
    return default_model_pka.get(kind);
}

/// Charge of the charged form; the sign of every shift.
pub fn charge(kind: pka.SiteKind) f64 {
    return if (kind.isAcid()) -1.0 else 1.0;
}

/// Relative atomic volume used by the desolvation sum.
pub const VolumeClass = enum { c_backbone, c_other, n, o, s, f, cl, p, other };

pub fn relativeVolume(class: VolumeClass) f64 {
    return switch (class) {
        .c_backbone => 1.40,
        .c_other => 2.64,
        .n => 1.06,
        .o => 1.00,
        .s => 1.66,
        .f => 0.90,
        .cl => 1.53,
        .p => 1.66,
        .other => 1.00,
    };
}

/// Ion table ("Ions"): residue name -> charge.
pub fn ionCharge(comp_id: []const u8) ?f64 {
    const map = std.StaticStringMap(f64).initComptime(.{
        .{ "MG", 2.0 },  .{ "CA", 2.0 },   .{ "ZN", 2.0 }, .{ "NA", 1.0 },
        .{ "CL", -1.0 }, .{ "MN", 2.0 },   .{ "K", 1.0 },  .{ "CO", 2.0 },
        .{ "FE2", 2.0 }, .{ "CD", 2.0 },   .{ "FE", 3.0 }, .{ "SR", 2.0 },
        .{ "CU", 2.0 },  .{ "IOD", -1.0 }, .{ "HG", 2.0 }, .{ "BR", -1.0 },
        .{ "NI", 2.0 },  .{ "1P", 1.0 },   .{ "2P", 2.0 }, .{ "1N", -1.0 },
        .{ "2N", -2.0 },
    });
    return map.get(comp_id);
}

/// Residues removed before anything else (parity input preparation).
pub fn isRemovedResidue(comp_id: []const u8) bool {
    const map = std.StaticStringMap(void).initComptime(.{
        .{ "HOH", {} }, .{ "H2O", {} }, .{ "SO4", {} }, .{ "PO4", {} },
        .{ "PEG", {} }, .{ "EPE", {} }, .{ "TRS", {} },
    });
    return map.has(comp_id);
}

/// Index of a side-chain group type in the interaction matrix.
fn matrixIndex(t: GroupType) ?usize {
    return switch (t) {
        .cys => 0,
        .nterm => 1,
        .his => 2,
        .lys => 3,
        .amd => 4,
        .coo => 5,
        .arg => 6,
        .trp => 7,
        .roh => 8,
        .tyr => 9,
        .bbn, .bbc, .ion => null,
    };
}

// Rows/columns: CYS N+ HIS LYS AMD COO ARG TRP ROH TYR
const I = PairClass.iterative;
const N = PairClass.non_iterative;
const X = PairClass.none;
const matrix = [10][10]PairClass{
    .{ I, N, I, N, N, I, N, N, N, N }, // CYS
    .{ N, I, N, N, X, N, N, X, X, I }, // N+
    .{ I, N, I, N, N, I, N, X, X, I }, // HIS
    .{ N, N, N, I, X, N, N, X, X, I }, // LYS
    .{ N, X, N, X, X, N, X, X, X, N }, // AMD
    .{ I, N, I, N, N, I, N, N, N, N }, // COO
    .{ N, N, N, N, X, N, I, X, X, N }, // ARG
    .{ N, X, X, X, X, N, X, X, X, N }, // TRP
    .{ N, X, X, X, X, N, X, X, X, N }, // ROH
    .{ N, I, I, I, N, N, N, N, N, I }, // TYR
};

pub fn pairClass(a: GroupType, b: GroupType) PairClass {
    const i = matrixIndex(a) orelse return .none;
    const j = matrixIndex(b) orelse return .none;
    return matrix[i][j];
}

pub const Cutoffs = struct { c1: f64, c2: f64 };

/// Side-chain H-bond cut-offs by pair type (symmetric; default 3.0/4.0).
pub fn sideChainCutoffs(a: GroupType, b: GroupType) Cutoffs {
    const Entry = struct { a: GroupType, b: GroupType, c: Cutoffs };
    const table = [_]Entry{
        .{ .a = .coo, .b = .coo, .c = .{ .c1 = 2.50, .c2 = 3.50 } },
        .{ .a = .coo, .b = .arg, .c = .{ .c1 = 1.85, .c2 = 2.85 } },
        .{ .a = .coo, .b = .lys, .c = .{ .c1 = 2.85, .c2 = 3.85 } },
        .{ .a = .coo, .b = .his, .c = .{ .c1 = 2.00, .c2 = 3.00 } },
        .{ .a = .coo, .b = .amd, .c = .{ .c1 = 2.00, .c2 = 3.00 } },
        .{ .a = .coo, .b = .trp, .c = .{ .c1 = 2.00, .c2 = 3.00 } },
        .{ .a = .coo, .b = .roh, .c = .{ .c1 = 2.65, .c2 = 3.65 } },
        .{ .a = .coo, .b = .tyr, .c = .{ .c1 = 2.65, .c2 = 3.65 } },
        .{ .a = .coo, .b = .nterm, .c = .{ .c1 = 2.85, .c2 = 3.85 } },
        .{ .a = .arg, .b = .cys, .c = .{ .c1 = 2.50, .c2 = 4.00 } },
        .{ .a = .arg, .b = .tyr, .c = .{ .c1 = 2.50, .c2 = 4.00 } },
        .{ .a = .his, .b = .amd, .c = .{ .c1 = 2.00, .c2 = 3.00 } },
        .{ .a = .his, .b = .tyr, .c = .{ .c1 = 2.00, .c2 = 3.00 } },
        .{ .a = .cys, .b = .cys, .c = .{ .c1 = 3.00, .c2 = 5.00 } },
        .{ .a = .cys, .b = .trp, .c = .{ .c1 = 2.50, .c2 = 3.50 } },
        .{ .a = .cys, .b = .roh, .c = .{ .c1 = 3.50, .c2 = 4.50 } },
        .{ .a = .cys, .b = .amd, .c = .{ .c1 = 2.50, .c2 = 3.50 } },
        .{ .a = .cys, .b = .tyr, .c = .{ .c1 = 3.50, .c2 = 4.50 } },
        .{ .a = .cys, .b = .nterm, .c = .{ .c1 = 3.00, .c2 = 4.50 } },
        .{ .a = .tyr, .b = .tyr, .c = .{ .c1 = 3.50, .c2 = 4.50 } },
        .{ .a = .tyr, .b = .nterm, .c = .{ .c1 = 3.00, .c2 = 4.50 } },
        .{ .a = .tyr, .b = .amd, .c = .{ .c1 = 2.50, .c2 = 3.50 } },
        .{ .a = .tyr, .b = .trp, .c = .{ .c1 = 2.50, .c2 = 3.50 } },
        .{ .a = .tyr, .b = .roh, .c = .{ .c1 = 3.50, .c2 = 4.50 } },
    };
    for (table) |e| {
        if ((e.a == a and e.b == b) or (e.a == b and e.b == a)) return e.c;
    }
    return .{ .c1 = 3.0, .c2 = 4.0 };
}

pub const BackboneParams = struct { delta: f64, c1: f64, c2: f64 };

/// Backbone H-bond parameters by (backbone type, site type).
pub fn backboneParams(backbone: GroupType, site: GroupType) ?BackboneParams {
    return switch (backbone) {
        .bbn => switch (site) {
            .coo => .{ .delta = 0.85, .c1 = 2.00, .c2 = 3.00 },
            .cys => .{ .delta = 0.85, .c1 = 3.00, .c2 = 4.00 },
            .tyr => .{ .delta = 0.85, .c1 = 2.20, .c2 = 3.20 },
            else => null,
        },
        .bbc => switch (site) {
            .his => .{ .delta = 0.85, .c1 = 2.00, .c2 = 3.00 },
            else => null,
        },
        else => null,
    };
}

const testing = std.testing;

test "interaction matrix is symmetric" {
    const types = [_]GroupType{ .cys, .nterm, .his, .lys, .amd, .coo, .arg, .trp, .roh, .tyr };
    for (types) |a| for (types) |b| try testing.expectEqual(pairClass(a, b), pairClass(b, a));
    try testing.expectEqual(PairClass.iterative, pairClass(.coo, .his));
    try testing.expectEqual(PairClass.non_iterative, pairClass(.coo, .arg));
    try testing.expectEqual(PairClass.none, pairClass(.roh, .roh));
    try testing.expectEqual(PairClass.none, pairClass(.bbn, .coo));
}

test "cut-off table is symmetric with a 3/4 default" {
    try testing.expectEqual(@as(f64, 1.85), sideChainCutoffs(.arg, .coo).c1);
    try testing.expectEqual(@as(f64, 2.85), sideChainCutoffs(.coo, .arg).c2);
    try testing.expectEqual(@as(f64, 3.0), sideChainCutoffs(.lys, .lys).c1);
    try testing.expectEqual(@as(f64, 4.0), sideChainCutoffs(.his, .his).c2);
}

test "correction sets parse and name their version" {
    const f = try Fixes.parse("Q2,q7");
    try testing.expect(f.q2 and f.q7 and !f.q1);
    try testing.expectEqual(@as(u8, 0b0100_0010), f.bits());
    try testing.expect((try Fixes.parse("none")).isNone());
    try testing.expectError(error.InvalidFix, Fixes.parse("Q9"));
    try testing.expectError(error.InvalidFix, Fixes.parse(""));
    try testing.expectEqualStrings("zhyd-empirical-v1+Q2+Q7", versionFor(f, .{ .q1 = true }, false));
    try testing.expectEqualStrings(zhyd_v1_version, versionFor(f, f, false));
    try testing.expectEqualStrings(parity_version, versionFor(.none, f, false));
    try testing.expectEqualStrings(zhyd_version, versionFor(f, f, true));
    try testing.expectEqualStrings("zhyd-empirical-v1.3+Q2+Q7", versionFor(f, .{ .q1 = true }, true));
    try testing.expectEqualStrings("zhyd-empirical-v1.3+none", versionFor(.none, f, true));
}

test "parameter overrides: explicit keys, required version" {
    const json =
        \\{"version": "cys-fit", "model_pka": {"CYS": 8.2}, "desolvation_scale": {"CYS": 0.5}}
    ;
    const parsed = try parseOverrides(testing.allocator, json);
    defer parsed.deinit();
    const p = parsed.value.applyTo(.{});
    try testing.expectEqualStrings("cys-fit", p.version);
    try testing.expectEqual(@as(f64, 8.2), p.model_pka.get(.cys));
    try testing.expectEqual(@as(f64, 0.5), p.desolvation_scale.get(.cys));
    // Other kinds keep the base values.
    try testing.expectEqual(@as(f64, 10.0), p.model_pka.get(.tyr));
    try testing.expectEqual(@as(f64, 1.0), p.desolvation_scale.get(.asp));
    // Unknown keys, a missing or built-in version, and out-of-range values fail.
    try testing.expectError(error.UnknownField, parseOverrides(testing.allocator, "{\"version\": \"x\", \"model_pka\": {\"SER\": 1}}"));
    try testing.expectError(error.UnknownField, parseOverrides(testing.allocator, "{\"version\": \"x\", \"eps_water\": 80}"));
    try testing.expectError(error.MissingField, parseOverrides(testing.allocator, "{\"model_pka\": {}}"));
    try testing.expectError(error.ParameterVersionRequired, parseOverrides(testing.allocator, "{\"version\": \"zhyd-empirical-v1.3\"}"));
    try testing.expectError(error.ParameterOutOfRange, parseOverrides(testing.allocator, "{\"version\": \"x\", \"desolvation_scale\": {\"CYS\": -1}}"));
}

/// Default charge of a monatomic metal ion by element (Q5); null for
/// elements that are not treated as free ions.
pub fn metalCharge(element: AtomType) ?f64 {
    return switch (element) {
        .Li, .Na, .K, .Rb, .Cs, .Ag => 1.0,
        .Mg, .Ca, .Sr, .Ba, .Mn, .Co, .Ni, .Cu, .Zn, .Cd, .Hg, .Pb => 2.0,
        .Fe => 3.0,
        else => null,
    };
}

/// Charge of a monatomic ion residue under Q5: the parity ion table, a few
/// CCD ids that name another oxidation state, else the element default.
pub fn ionChargeByElement(comp_id: []const u8, element: AtomType) ?f64 {
    if (ionCharge(comp_id)) |q| return q;
    const named = std.StaticStringMap(f64).initComptime(.{
        .{ "CU1", 1.0 }, .{ "MN3", 3.0 }, .{ "CO3", 3.0 }, .{ "3CO", 3.0 }, .{ "FE2", 2.0 },
    });
    if (named.get(comp_id)) |q| return q;
    return metalCharge(element);
}

test "ion charges by element" {
    try testing.expectEqual(@as(?f64, 2.0), ionChargeByElement("ZN2", .Zn));
    try testing.expectEqual(@as(?f64, 1.0), ionChargeByElement("CU1", .Cu));
    try testing.expectEqual(@as(?f64, 3.0), ionChargeByElement("FE", .Fe));
    try testing.expectEqual(@as(?f64, -1.0), ionChargeByElement("CL", .Cl));
    try testing.expectEqual(@as(?f64, null), ionChargeByElement("XYZ", .Mo));
}

test "ion table" {
    try testing.expectEqual(@as(?f64, 2.0), ionCharge("ZN"));
    try testing.expectEqual(@as(?f64, -2.0), ionCharge("2N"));
    try testing.expectEqual(@as(?f64, null), ionCharge("SO4"));
}

test "recalibrated table changes only CYS and HIS" {
    for (std.enums.values(pka.SiteKind)) |k| {
        try testing.expectEqual(k == .cys, calibrated.model_pka.get(k) != default_model_pka.get(k));
        try testing.expectEqual(k == .cys or k == .his, calibrated.desolvation_scale.get(k) != 1.0);
    }
    try testing.expectEqual(@as(f64, 7.5), calibrated.model_pka.get(.cys));
    try testing.expectEqual(@as(f64, 0.6), calibrated.desolvation_scale.get(.his));
    // TYR keeps its desolvation: a zero scale drops buried tyrosines to pKa
    // 6-7 and lets them pair with neighbouring histidines (spec r7).
    try testing.expectEqual(@as(f64, 1.0), calibrated.desolvation_scale.get(.tyr));
}
