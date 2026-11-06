//! zhyd empirical pKa engine (MIT), parity mode.
//!
//! An independent implementation of the PROPKA 3 empirical model as written
//! in docs/specs/zhyd-pka-empirical-engine.md. It was implemented from that
//! specification, the published papers and black-box oracle reports only; it
//! contains no code from PROPKA or zpropka.
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

pub const engine_name = "zhyd-empirical (parity)";
pub const parameter_version = params.parity_version;

pub const Mode = enum {
    /// Reproduce the PROPKA 3 model, quirks included.
    parity,
};

pub const Options = struct {
    mode: Mode = .parity,
    parameters: params.Parameters = .{},
};

const Model = model_mod.Model;

pub fn predict(allocator: std.mem.Allocator, mdl: *const Model, options: Options) !pka.Prediction {
    const sums = try engine.run(allocator, mdl, options.parameters);
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
        // Not titrating in every conformation (the oracle reports 99.99).
        if (s.n == 0 or s.not_titrating_first or s.predicted >= params.not_titrating_pka - 0.005) {
            omitted += 1;
            continue;
        }
        var site = pka.SiteResult{
            .residue_idx = s.residue_idx,
            .kind = s.kind,
            .center = s.center.cast(f32),
            .terms = .{
                .model = @floatCast(params.modelPka(s.kind)),
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
        .parameter_version = options.parameters.version,
        .engine = engine_name,
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
    try testing.expectEqualStrings(engine_name, result.engine);
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
