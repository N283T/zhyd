//! zhyd empirical pKa engine (MIT): parity mode and zhyd mode.
//!
//! An independent implementation of the PROPKA 3 empirical model as written
//! in docs/specs/zhyd-pka-empirical-engine.md. It was implemented from that
//! specification, the published papers and black-box oracle reports only; it
//! contains no code from PROPKA or zpropka. zhyd mode applies the corrections
//! Q1–Q8 of the specification's "Known parity quirks", each behind its own
//! switch (`params.Fixes`).
//!
//! Disulfide, metal-coordinated and covalently bonded cysteines are listed in
//! `excluded` exactly as for the other engines.

const std = @import("std");
const model_mod = @import("model.zig");
const pka = @import("pka.zig");

pub const params = @import("pka_empirical/params.zig");
pub const hydrogens = @import("pka_empirical/hydrogens.zig");
pub const sites = @import("pka_empirical/sites.zig");
pub const desolvation = @import("pka_empirical/desolvation.zig");
pub const interactions = @import("pka_empirical/interactions.zig");
pub const iterative = @import("pka_empirical/iterative.zig");
pub const engine = @import("pka_empirical/engine.zig");

pub const parity_engine_name = "zhyd-empirical (parity)";
pub const zhyd_engine_name = "zhyd-empirical (zhyd)";
pub const parameter_version = params.parity_version;
pub const Fixes = params.Fixes;

/// Corrections enabled by zhyd mode: Q1–Q8, adopted on the PKAD-R train
/// split under the r4 protocol (see the specification's "zhyd mode
/// evaluation protocol" and "zhyd mode ablation").
pub const zhyd_fixes: Fixes = .{ .q1 = true, .q2 = true, .q3 = true, .q4 = true, .q5 = true, .q6 = true, .q7 = true, .q8 = true };

pub const Mode = enum {
    /// Reproduce the PROPKA 3 model, quirks included.
    parity,
    /// Adopted corrections (`zhyd_fixes`) on the recalibrated parameter
    /// table (`params.calibrated`).
    zhyd,

    pub fn defaultFixes(self: Mode) Fixes {
        return switch (self) {
            .parity => .none,
            .zhyd => zhyd_fixes,
        };
    }
};

pub const Options = struct {
    mode: Mode = .parity,
    /// Overrides the mode's correction set (ablation).
    fixes: ?Fixes = null,
    /// zhyd N-terminal policy; under `aggressive` chain breaks create
    /// termini for Q3 and Q4.
    nterm_mode: NtermMode = .auto,
    /// Per-kind model pKa / desolvation scale overrides (`--parameters`,
    /// zhyd mode only), applied on top of the mode's parameter table.
    overrides: ?params.Overrides = null,

    pub fn effectiveFixes(self: Options) Fixes {
        return self.fixes orelse self.mode.defaultFixes();
    }

    /// Full parameter set: the mode's table (parity table in parity mode,
    /// recalibrated table in zhyd mode), the correction set and overrides.
    pub fn parameters(self: Options) params.Parameters {
        const fixes = self.effectiveFixes();
        var p: params.Parameters = if (self.mode == .zhyd) params.calibrated else .{};
        p.fixes = fixes;
        p.break_termini = self.nterm_mode.treatsBreakAsNterm();
        p.version = params.versionFor(fixes, zhyd_fixes, self.mode == .zhyd);
        if (self.overrides) |o| p = o.applyTo(p);
        return p;
    }

    /// Engine name: parity only for the unmodified parity model.
    pub fn engineName(self: Options) []const u8 {
        const parity = self.mode == .parity and self.effectiveFixes().isNone() and self.overrides == null;
        return if (parity) parity_engine_name else zhyd_engine_name;
    }
};

pub fn isEngineName(name: []const u8) bool {
    return std.mem.eql(u8, name, parity_engine_name) or std.mem.eql(u8, name, zhyd_engine_name);
}

const Model = model_mod.Model;
const NtermMode = @import("place/terminal.zig").NtermMode;

pub fn predict(allocator: std.mem.Allocator, mdl: *const Model, options: Options) !pka.Prediction {
    const parameters = options.parameters();
    const fixes = parameters.fixes;
    const sums = try engine.run(allocator, mdl, parameters);
    defer allocator.free(sums);

    const metals = try pka.metalPositions(allocator, mdl);
    defer allocator.free(metals);

    var out_sites: std.ArrayListUnmanaged(pka.SiteResult) = .empty;
    errdefer out_sites.deinit(allocator);
    var excluded: std.ArrayListUnmanaged(pka.ExcludedSite) = .empty;
    errdefer excluded.deinit(allocator);
    var omitted: u32 = 0;

    for (sums) |s| {
        if (pka.exclusionReason(mdl, s.residue_idx, s.kind, metals)) |reason| {
            try excluded.append(allocator, .{ .residue_idx = s.residue_idx, .kind = s.kind, .reason = reason });
            continue;
        }
        // Not titrating in every conformation (the oracle reports 99.99). Q7:
        // reported when it titrates in any conformation.
        if (s.n == 0 or (fixes.q7 and s.n_titrating == 0) or s.not_titrating_first or s.predicted >= params.not_titrating_pka - 0.005) {
            omitted += 1;
            continue;
        }
        var site = pka.SiteResult{
            .residue_idx = s.residue_idx,
            .kind = s.kind,
            .center = s.center.cast(f32),
            .terms = .{
                .model = @floatCast(parameters.model_pka.get(s.kind)),
                .burial = @floatCast(s.burial),
                .hydrogen_bond = @floatCast(s.hydrogen_bond),
                .electrostatic = @floatCast(s.electrostatic),
            },
            .predicted_pka = @floatCast(s.predicted),
            .burial_fraction = @floatCast(s.weight),
            .confidence = .high,
        };
        site.confidence = pka.assessConfidence(site);
        if (!pka.isPhysicalPka(site.predicted_pka)) site.confidence = .low;
        try out_sites.append(allocator, site);
    }

    return .{
        .allocator = allocator,
        .sites = try out_sites.toOwnedSlice(allocator),
        .excluded = try excluded.toOwnedSlice(allocator),
        .omitted_sites = omitted,
        .parameter_version = parameters.version,
        .engine = options.engineName(),
    };
}

const testing = std.testing;

test {
    testing.refAllDecls(@This());
}

test "prediction excludes disulfide cysteines and labels the engine" {
    const mmcif = @import("mmcif.zig");
    var mdl = try mmcif.parseModel(testing.allocator, @embedFile("test_data/disulfide.cif"));
    defer mdl.deinit();
    for (mdl.residues.items) |*res| res.entity_type = .polymer;
    var result = try predict(testing.allocator, &mdl, .{});
    defer result.deinit();
    try testing.expectEqualStrings(parity_engine_name, result.engine);
    try testing.expectEqualStrings(parameter_version, result.parameter_version);
    for (result.sites) |site| try testing.expect(site.kind != .cys);
    try testing.expectEqual(@as(usize, 2), result.excluded.len);
    for (result.excluded) |site| try testing.expectEqual(pka.ExclusionReason.disulfide, site.reason);
}

test "prediction terms add up to the predicted value" {
    const mmcif = @import("mmcif.zig");
    var mdl = try mmcif.parseModel(testing.allocator, @embedFile("test_data/his.cif"));
    defer mdl.deinit();
    for (mdl.residues.items) |*res| res.entity_type = .polymer;
    var result = try predict(testing.allocator, &mdl, .{});
    defer result.deinit();
    // The N-terminal HIS is covalently coupled and omitted; the N-terminus
    // remains.
    try testing.expectEqual(@as(usize, 1), result.sites.len);
    try testing.expectEqual(@as(u32, 1), result.omitted_sites);
    const site = result.sites[0];
    try testing.expectEqual(pka.SiteKind.n_term, site.kind);
    try testing.expectApproxEqAbs(site.predicted_pka, site.terms.total(), 1e-4);
}

test "zhyd mode names its engine and correction set" {
    const mmcif = @import("mmcif.zig");
    var mdl = try mmcif.parseModel(testing.allocator, @embedFile("test_data/his.cif"));
    defer mdl.deinit();
    for (mdl.residues.items) |*res| res.entity_type = .polymer;
    const Case = struct { options: Options, engine: []const u8, version: []const u8 };
    const cases = [_]Case{
        .{ .options = .{}, .engine = parity_engine_name, .version = params.parity_version },
        .{ .options = .{ .mode = .zhyd }, .engine = zhyd_engine_name, .version = params.zhyd_version },
        .{ .options = .{ .mode = .zhyd, .fixes = .none }, .engine = zhyd_engine_name, .version = "zhyd-empirical-v1.3+none" },
        .{ .options = .{ .fixes = zhyd_fixes }, .engine = zhyd_engine_name, .version = params.zhyd_v1_version },
        .{ .options = .{ .fixes = .{ .q2 = true, .q7 = true } }, .engine = zhyd_engine_name, .version = "zhyd-empirical-v1+Q2+Q7" },
        .{ .options = .{ .mode = .zhyd, .overrides = .{ .version = "my-fit" } }, .engine = zhyd_engine_name, .version = "my-fit" },
    };
    for (cases) |c| {
        var result = try predict(testing.allocator, &mdl, c.options);
        defer result.deinit();
        try testing.expectEqualStrings(c.engine, result.engine);
        try testing.expectEqualStrings(c.version, result.parameter_version);
    }
}

test "overrides change the model pKa and scale the desolvation shift" {
    const mmcif = @import("mmcif.zig");
    var mdl = try mmcif.parseModel(testing.allocator, @embedFile("test_data/his.cif"));
    defer mdl.deinit();
    for (mdl.residues.items) |*res| res.entity_type = .polymer;
    var base = try predict(testing.allocator, &mdl, .{ .mode = .zhyd });
    defer base.deinit();
    // N-terminus: model 8.00 → 7.00; desolvation × 0.5. The site has no
    // pair partner here (the HIS is covalently coupled out), so the total
    // moves by exactly −1 + desolvation / 2.
    var o = params.Overrides{ .version = "test" };
    o.model_pka.N_TERM = 7.0;
    o.desolvation_scale.N_TERM = 0.5;
    var fit = try predict(testing.allocator, &mdl, .{ .mode = .zhyd, .overrides = o });
    defer fit.deinit();
    try testing.expectEqual(pka.SiteKind.n_term, base.sites[0].kind);
    const b = base.sites[0];
    const f = fit.sites[0];
    try testing.expectApproxEqAbs(@as(f32, 7.0), f.terms.model, 1e-6);
    try testing.expectApproxEqAbs(b.terms.burial * 0.5, f.terms.burial, 1e-5);
    try testing.expectApproxEqAbs(b.predicted_pka - 1.0 - b.terms.burial * 0.5, f.predicted_pka, 1e-5);
}
