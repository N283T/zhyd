//! pH-driven protonation states (docs/specs/zhyd-protonation-pka.md,
//! Phase 4): turn a pKa prediction into protonation overrides for placement.
//!
//! The pKa decides only the macrostate (charge). Tautomers and orientations
//! stay with the optimizer: a pKa-derived HIP keeps its ring flip, and a
//! protonated Asp/Glu leaves the proton's oxygen to the carboxyl flip mover.
//! A site is protonated when its predicted pKa is at or above the pH, the
//! same rule as the `pka` report's recommended state. Only sites whose state
//! differs from the default get an entry; TYR, ARG and the termini keep their
//! default states.

const std = @import("std");
const model_mod = @import("model.zig");
const pka = @import("pka.zig");
const protonation = @import("place/protonation.zig");

const Model = model_mod.Model;
const ProtonationOverrides = protonation.ProtonationOverrides;

/// The non-default state a predicted pKa implies at `ph`, if any.
pub fn stateAt(kind: pka.SiteKind, predicted_pka: f32, ph: f32) ?protonation.ResidueState {
    if (!pka.isPhysicalPka(predicted_pka)) return null;
    const protonated = predicted_pka >= ph;
    return switch (kind) {
        .his => if (protonated) .{ .his = .hip } else null,
        .asp => if (protonated) .{ .asp = .protonated } else null,
        .glu => if (protonated) .{ .glu = .protonated } else null,
        .lys => if (protonated) null else .{ .lys = .neutral },
        .cys => if (protonated) null else .{ .cys = .thiolate },
        .tyr, .arg, .n_term, .c_term => null,
    };
}

/// Overrides for one model: the pKa-derived states first, then copies of the
/// user's entries, so a user entry wins over a prediction for the same residue.
pub fn fromPrediction(
    allocator: std.mem.Allocator,
    mdl: *const Model,
    prediction: *const pka.Prediction,
    ph: f32,
    user: ?*const ProtonationOverrides,
) !ProtonationOverrides {
    var entries: std.ArrayListUnmanaged(protonation.Entry) = .empty;
    errdefer {
        for (entries.items) |entry| allocator.free(entry.selector.chain_id);
        entries.deinit(allocator);
    }

    for (prediction.sites) |site| {
        const state = stateAt(site.kind, site.predicted_pka, ph) orelse continue;
        const res = mdl.residues.items[site.residue_idx];
        const chain = mdl.chains.items[res.chain_idx];
        var entry = protonation.Entry{
            .selector = .{
                .chain_id = try allocator.dupe(u8, chain.authSlice()),
                .auth_seq_id = res.auth_seq_id,
                .ins_code = res.ins_code,
            },
            .state = state,
            .source = .pka,
        };
        entry.comp_id.set(res.compIdSlice());
        entries.append(allocator, entry) catch |err| {
            allocator.free(entry.selector.chain_id);
            return err;
        };
    }

    if (user) |ov| {
        for (ov.entries) |user_entry| {
            var entry = user_entry;
            entry.selector.chain_id = try allocator.dupe(u8, user_entry.selector.chain_id);
            entries.append(allocator, entry) catch |err| {
                allocator.free(entry.selector.chain_id);
                return err;
            };
        }
    }

    return .{ .allocator = allocator, .entries = try entries.toOwnedSlice(allocator) };
}

const testing = std.testing;

test "states follow the predicted pKa relative to the pH" {
    try testing.expectEqual(protonation.ResidueState{ .his = .hip }, stateAt(.his, 7.2, 7.0).?);
    try testing.expectEqual(@as(?protonation.ResidueState, null), stateAt(.his, 6.4, 7.0));
    try testing.expectEqual(protonation.ResidueState{ .asp = .protonated }, stateAt(.asp, 7.5, 7.0).?);
    try testing.expectEqual(@as(?protonation.ResidueState, null), stateAt(.glu, 4.4, 7.0));
    try testing.expectEqual(protonation.ResidueState{ .lys = .neutral }, stateAt(.lys, 6.0, 7.0).?);
    try testing.expectEqual(protonation.ResidueState{ .cys = .thiolate }, stateAt(.cys, 6.5, 7.0).?);
    try testing.expectEqual(@as(?protonation.ResidueState, null), stateAt(.cys, 9.0, 7.0));
    try testing.expectEqual(@as(?protonation.ResidueState, null), stateAt(.tyr, 6.0, 7.0));
    try testing.expectEqual(@as(?protonation.ResidueState, null), stateAt(.his, 30.0, 7.0));
}

test "user overrides win over pKa-derived states" {
    const mmcif = @import("mmcif.zig");
    var mdl = try mmcif.parseModel(testing.allocator, @embedFile("test_data/his.cif"));
    defer mdl.deinit();
    var his_idx: ?u32 = null;
    for (mdl.residues.items, 0..) |res, i| {
        if (std.mem.eql(u8, res.compIdSlice(), "HIS")) his_idx = @intCast(i);
    }
    const sites = [_]pka.SiteResult{.{
        .residue_idx = his_idx.?,
        .kind = .his,
        .center = .{ .x = 0, .y = 0, .z = 0 },
        .terms = .{ .model = 6.5 },
        .predicted_pka = 7.5,
        .burial_fraction = 0,
        .confidence = .high,
    }};
    const prediction = pka.Prediction{
        .allocator = testing.allocator,
        .sites = @constCast(&sites),
        .excluded = &.{},
        .omitted_sites = 0,
        .parameter_version = "test",
        .engine = "test",
    };

    var predicted = try fromPrediction(testing.allocator, &mdl, &prediction, 7.0, null);
    defer predicted.deinit();
    const entry = predicted.findEntry(&mdl, his_idx.?).?;
    try testing.expectEqual(protonation.Source.pka, entry.source);
    try testing.expectEqual(protonation.HisState.hip, entry.state.his);

    const res = mdl.residues.items[his_idx.?];
    const line = try std.fmt.allocPrint(testing.allocator, "{s}:{d} HIS HID\n", .{ mdl.chains.items[res.chain_idx].authSlice(), res.auth_seq_id });
    defer testing.allocator.free(line);
    var user = try protonation.parseString(testing.allocator, line);
    defer user.deinit();
    var merged = try fromPrediction(testing.allocator, &mdl, &prediction, 7.0, &user);
    defer merged.deinit();
    const winner = merged.findEntry(&mdl, his_idx.?).?;
    try testing.expectEqual(protonation.Source.user, winner.source);
    try testing.expectEqual(protonation.HisState.hid, winner.state.his);
}
