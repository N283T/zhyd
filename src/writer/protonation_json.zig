//! Versioned machine-readable protonation and pKa report.

const std = @import("std");
const model_mod = @import("../model.zig");
const pka = @import("../pka.zig");
const pka_empirical = @import("../pka_empirical.zig");
const protonation = @import("../place/protonation.zig");
const terminal = @import("../place/terminal.zig");
const math = @import("../math.zig");

pub const ReportModel = struct {
    model_num: u32,
    model: *const model_mod.Model,
    prediction: *const pka.Prediction,
};

pub const StateReportModel = struct {
    model_num: u32,
    model: *const model_mod.Model,
};

pub fn writeStateReport(
    writer: anytype,
    input: []const u8,
    coordinate_file: ?[]const u8,
    models: []const StateReportModel,
    overrides: ?*const protonation.ProtonationOverrides,
    nterm_mode: terminal.NtermMode,
) !void {
    try writer.writeAll("{\n  \"schema_version\": \"");
    try writer.writeAll(pka.schema_version);
    try writer.writeAll("\",\n  \"report_kind\": \"applied_protonation_state\",\n  \"status\": \"coordinate_assignment\",\n  \"input\": \"");
    try writeJsonString(writer, input);
    try writer.writeAll("\",\n  \"coordinate_file\": ");
    if (coordinate_file) |path| {
        try writer.writeAll("\"");
        try writeJsonString(writer, path);
        try writer.writeAll("\"");
    } else try writer.writeAll("null");
    try writer.writeAll(",\n  \"parameter_version\": null,\n  \"models\": [\n");
    for (models, 0..) |entry, model_i| {
        try writer.print("    {{\"model\": {d}, \"sites\": [\n", .{entry.model_num});
        const metals = try pka.metalPositions(entry.model.allocator, entry.model);
        defer entry.model.allocator.free(metals);
        var first = true;
        for (entry.model.residues.items, 0..) |res, res_idx| {
            const comp = res.compIdSlice();
            const kind: ?pka.SiteKind = if (std.mem.eql(u8, comp, "ASP")) .asp else if (std.mem.eql(u8, comp, "GLU")) .glu else if (std.mem.eql(u8, comp, "HIS")) .his else if (std.mem.eql(u8, comp, "LYS")) .lys else if (std.mem.eql(u8, comp, "CYS")) .cys else if (std.mem.eql(u8, comp, "TYR")) .tyr else null;
            if (kind) |k| {
                if (!first) try writer.writeAll(",\n");
                first = false;
                try writeAppliedSite(writer, entry.model, @intCast(res_idx), k, overrides, nterm_mode, metals);
            }
            const chain = entry.model.chains.items[res.chain_idx];
            const is_nterm = res_idx == chain.residue_start or
                (res.is_chain_break_before and nterm_mode.treatsBreakAsNterm());
            const is_cterm = res_idx == chain.residue_end - 1 or
                (res_idx + 1 < entry.model.residues.items.len and entry.model.residues.items[res_idx + 1].is_chain_break_before);
            if (is_nterm) {
                if (!first) try writer.writeAll(",\n");
                first = false;
                try writeAppliedSite(writer, entry.model, @intCast(res_idx), .n_term, overrides, nterm_mode, metals);
            }
            if (is_cterm) {
                if (!first) try writer.writeAll(",\n");
                first = false;
                try writeAppliedSite(writer, entry.model, @intCast(res_idx), .c_term, overrides, nterm_mode, metals);
            }
        }
        try writer.writeAll("\n    ]}");
        if (model_i + 1 < models.len) try writer.writeAll(",");
        try writer.writeAll("\n");
    }
    try writer.writeAll("  ]\n}\n");
}

const AppliedState = struct {
    macrostate: []const u8,
    tautomer: ?[]const u8 = null,
    charge: i8,
    source: []const u8,
    warning: ?[]const u8 = null,
};

fn writeAppliedSite(
    writer: anytype,
    mdl: *const model_mod.Model,
    residue_idx: u32,
    kind: pka.SiteKind,
    overrides: ?*const protonation.ProtonationOverrides,
    nterm_mode: terminal.NtermMode,
    metals: []const math.Vec3(f32),
) !void {
    const res = mdl.residues.items[residue_idx];
    const chain = mdl.chains.items[res.chain_idx];
    var state = appliedState(mdl, residue_idx, kind, overrides, nterm_mode, metals);
    // The applied state comes from protein-only rules; flag metal ligands so
    // downstream charge models do not trust it blindly.
    if (state.warning == null and pka.isMetalCoordinated(mdl, residue_idx, kind, metals)) {
        state.warning = "metal_coordinated_state_not_modeled";
    }
    try writer.writeAll("      {\"id\": {\"chain_label\": \"");
    try writeJsonString(writer, chain.labelSlice());
    try writer.writeAll("\", \"chain_auth\": \"");
    try writeJsonString(writer, chain.authSlice());
    try writer.writeAll("\", \"comp_id\": \"");
    try writeJsonString(writer, res.compIdSlice());
    try writer.print("\", \"seq_id\": {d}, \"auth_seq_id\": {d}, \"ins_code\": ", .{ res.seq_id, res.auth_seq_id });
    if (res.ins_code == ' ') try writer.writeAll("null") else try writer.print("\"{c}\"", .{res.ins_code});
    try writer.writeAll(", \"site\": \"");
    try writer.writeAll(kind.name());
    try writer.writeAll("\"}, \"state\": {\"macrostate\": \"");
    try writer.writeAll(state.macrostate);
    try writer.writeAll("\", \"tautomer\": ");
    if (state.tautomer) |value| {
        try writer.print("\"{s}\"", .{value});
    } else try writer.writeAll("null");
    try writer.print(", \"formal_charge\": {d}, \"source\": \"{s}\", \"applies_to_coordinates\": true}}, \"pka\": null, \"confidence\": null, \"warnings\": [", .{ state.charge, state.source });
    if (state.warning) |warning| try writer.print("\"{s}\"", .{warning});
    try writer.writeAll("]}");
}

fn appliedState(
    mdl: *const model_mod.Model,
    residue_idx: u32,
    kind: pka.SiteKind,
    overrides: ?*const protonation.ProtonationOverrides,
    nterm_mode: terminal.NtermMode,
    metals: []const math.Vec3(f32),
) AppliedState {
    const override_state = if (overrides) |value| value.find(mdl, residue_idx) else null;

    // A covalently bonded SG (cystine or covalent adduct) has no titratable
    // proton; a thiol/thiolate override cannot apply and is reported rather
    // than silently honored. Metal ligands are handled by the caller.
    if (kind == .cys and !pka.isMetalCoordinated(mdl, residue_idx, kind, metals)) {
        if (pka.cysteineBond(mdl, residue_idx)) |bond| {
            const overridden = if (override_state) |value| value == .cys else false;
            return .{
                .macrostate = if (bond == .disulfide) "disulfide_bonded" else "covalently_bonded",
                .charge = 0,
                .source = "covalent_bond",
                .warning = if (overridden) "override_ignored_covalent_bond" else null,
            };
        }
    }

    if (override_state) |value| switch (value) {
        .his => |his| if (kind == .his and his != .auto) return switch (his) {
            .hid => .{ .macrostate = "neutral", .tautomer = "HID", .charge = 0, .source = "user_override" },
            .hie => .{ .macrostate = "neutral", .tautomer = "HIE", .charge = 0, .source = "user_override" },
            .hip => .{ .macrostate = "protonated_cation", .tautomer = "HIP", .charge = 1, .source = "user_override" },
            .auto => unreachable,
        },
        .asp => |acid| if (kind == .asp) return acidAppliedState(acid, "OD1", "OD2"),
        .glu => |acid| if (kind == .glu) return acidAppliedState(acid, "OE1", "OE2"),
        .lys => |lys| if (kind == .lys) return if (lys == .charged)
            .{ .macrostate = "protonated_cation", .charge = 1, .source = "user_override" }
        else
            .{ .macrostate = "neutral", .charge = 0, .source = "user_override" },
        .cys => |cys| if (kind == .cys) return if (cys == .thiol)
            .{ .macrostate = "protonated_neutral", .charge = 0, .source = "user_override" }
        else
            .{ .macrostate = "deprotonated_anion", .charge = -1, .source = "user_override" },
    };

    return switch (kind) {
        .asp, .glu => .{ .macrostate = "deprotonated_anion", .charge = -1, .source = "default" },
        .lys, .arg => .{ .macrostate = "protonated_cation", .charge = 1, .source = "default" },
        .cys, .tyr => .{ .macrostate = "protonated_neutral", .charge = 0, .source = "default" },
        .his => histidineAppliedState(mdl, residue_idx),
        .n_term => if (nterm_mode == .neutral and !std.mem.eql(u8, mdl.residues.items[residue_idx].compIdSlice(), "PRO"))
            .{ .macrostate = "neutral", .charge = 0, .source = "user_override" }
        else
            .{ .macrostate = "protonated_cation", .charge = 1, .source = if (nterm_mode == .auto) "default" else "user_override" },
        .c_term => .{ .macrostate = "deprotonated_anion", .charge = -1, .source = "default" },
    };
}

fn acidAppliedState(state: protonation.AcidState, atom1: []const u8, atom2: []const u8) AppliedState {
    return switch (state) {
        .deprotonated => .{ .macrostate = "deprotonated_anion", .charge = -1, .source = "user_override" },
        .atom1 => .{ .macrostate = "protonated_neutral", .tautomer = atom1, .charge = 0, .source = "user_override" },
        .atom2 => .{ .macrostate = "protonated_neutral", .tautomer = atom2, .charge = 0, .source = "user_override" },
    };
}

fn histidineAppliedState(mdl: *const model_mod.Model, residue_idx: u32) AppliedState {
    const hd1 = hydrogenPresent(mdl, residue_idx, "HD1");
    const he2 = hydrogenPresent(mdl, residue_idx, "HE2");
    if (hd1 and he2) return .{ .macrostate = "protonated_cation", .tautomer = "HIP", .charge = 1, .source = "geometry_selected" };
    if (hd1) return .{ .macrostate = "neutral", .tautomer = "HID", .charge = 0, .source = "geometry_selected" };
    if (he2) return .{ .macrostate = "neutral", .tautomer = "HIE", .charge = 0, .source = "geometry_selected" };
    return .{ .macrostate = "neutral", .tautomer = "UNKNOWN", .charge = 0, .source = "geometry_selected", .warning = "histidine_tautomer_not_observable_from_hydrogens" };
}

fn hydrogenPresent(mdl: *const model_mod.Model, residue_idx: u32, name: []const u8) bool {
    for (mdl.atoms.items) |atom| {
        if (atom.residue_idx != residue_idx or !atom.is_hydrogen) continue;
        if (!std.mem.eql(u8, atom.nameSlice(), name)) continue;
        if (atom.pos.x > 999 and atom.pos.y > 999 and atom.pos.z > 999) continue;
        return true;
    }
    return false;
}

pub fn writePredictionReport(writer: anytype, input: []const u8, models: []const ReportModel) !void {
    try writer.writeAll("{\n  \"schema_version\": \"");
    try writeJsonString(writer, pka.schema_version);
    try writer.writeAll("\",\n  \"report_kind\": \"pka_prediction\",\n");
    const engine = if (models.len > 0) models[0].prediction.engine else pka.engine_name;
    const builtin_engine = std.mem.eql(u8, engine, pka.engine_name);
    const empirical_engine = std.mem.eql(u8, engine, pka_empirical.parity_engine_name);
    const zhyd_engine = std.mem.eql(u8, engine, pka_empirical.zhyd_engine_name);
    try writer.writeAll(if (builtin_engine or zhyd_engine)
        "  \"status\": \"experimental_unvalidated\",\n"
    else if (empirical_engine)
        "  \"status\": \"experimental_parity_with_propka3_model\",\n"
    else
        "  \"status\": \"external_reference_method\",\n");
    try writer.writeAll("  \"input\": \"");
    try writeJsonString(writer, input);
    try writer.writeAll("\",\n  \"parameter_version\": \"");
    const version = if (models.len > 0) models[0].prediction.parameter_version else pka.parameter_version;
    try writeJsonString(writer, version);
    try writer.writeAll("\",\n  \"engine\": \"");
    try writeJsonString(writer, engine);
    try writer.writeAll("\",\n  \"reference_ph_for_recommended_state\": 7.0,\n");
    try writer.writeAll(if (builtin_engine)
        "  \"warnings\": [\"not_for_clinical_or_production_charge_assignment\", \"buried_sites_poorly_modeled\"],\n"
    else if (empirical_engine)
        "  \"warnings\": [\"not_for_clinical_or_production_charge_assignment\", \"experimental_engine_reproduces_propka3_model_quirks\"],\n"
    else if (zhyd_engine)
        "  \"warnings\": [\"not_for_clinical_or_production_charge_assignment\", \"experimental_engine_corrects_propka3_model_quirks\"],\n"
    else
        "  \"warnings\": [\"not_for_clinical_or_production_charge_assignment\"],\n");
    try writer.writeAll("  \"models\": [\n");

    for (models, 0..) |entry, model_i| {
        try writer.print("    {{\"model\": {d}, \"omitted_sites\": {d}, \"sites\": [\n", .{ entry.model_num, entry.prediction.omitted_sites });
        for (entry.prediction.sites, 0..) |site, site_i| {
            const state = recommendedState(site.kind, site.predicted_pka, 7.0);
            try writer.writeAll("      {");
            try writeSiteId(writer, entry.model, site.residue_idx, site.kind);
            try writer.writeAll(", \"state\": {\"macrostate\": \"");
            try writer.writeAll(state.macrostate);
            try writer.writeAll("\", \"tautomer\": ");
            if (state.tautomer) |tautomer| {
                try writer.writeAll("\"");
                try writer.writeAll(tautomer);
                try writer.writeAll("\"");
            } else try writer.writeAll("null");
            try writer.print(", \"formal_charge\": {d}, \"source\": \"pka_predicted\", \"applies_to_coordinates\": false}}, ", .{state.formal_charge});
            try writer.print("\"pka\": {{\"predicted\": {d:.4}, \"terms\": {{\"model\": {d:.4}, \"burial\": {d:.4}, \"hydrogen_bond\": {d:.4}, \"electrostatic\": {d:.4}, \"coupling\": {d:.4}}}, \"burial_fraction\": {d:.4}}}, ", .{
                site.predicted_pka,
                site.terms.model,
                site.terms.burial,
                site.terms.hydrogen_bond,
                site.terms.electrostatic,
                site.terms.coupling,
                site.burial_fraction,
            });
            try writer.writeAll("\"confidence\": \"");
            try writer.writeAll(@tagName(site.confidence));
            try writer.writeAll("\", \"warnings\": [");
            var wrote_warning = false;
            if (site.coupled) {
                try writer.writeAll("\"coupled_group_mean_field_approximation\"");
                wrote_warning = true;
            }
            if (site.kind == .his and state.formal_charge == 0) {
                if (wrote_warning) try writer.writeAll(", ");
                try writer.writeAll("\"neutral_histidine_tautomer_requires_geometry_selection\"");
                wrote_warning = true;
            }
            if (site.confidence == .low) {
                if (wrote_warning) try writer.writeAll(", ");
                try writer.writeAll("\"low_confidence_environment\"");
                wrote_warning = true;
            }
            if (!pka.isPhysicalPka(site.predicted_pka)) {
                if (wrote_warning) try writer.writeAll(", ");
                try writer.writeAll("\"pka_outside_physical_range\"");
            }
            try writer.writeAll("]}");
            if (site_i + 1 < entry.prediction.sites.len) try writer.writeAll(",");
            try writer.writeAll("\n");
        }
        try writer.writeAll("    ], \"excluded_sites\": [");
        for (entry.prediction.excluded, 0..) |site, i| {
            if (i > 0) try writer.writeAll(",");
            try writer.writeAll("\n      {");
            try writeSiteId(writer, entry.model, site.residue_idx, site.kind);
            try writer.print(", \"reason\": \"{s}\"}}", .{@tagName(site.reason)});
        }
        if (entry.prediction.excluded.len > 0) try writer.writeAll("\n    ");
        try writer.writeAll("]}");
        if (model_i + 1 < models.len) try writer.writeAll(",");
        try writer.writeAll("\n");
    }
    try writer.writeAll("  ]\n}\n");
}

fn writeSiteId(writer: anytype, mdl: *const model_mod.Model, residue_idx: u32, kind: pka.SiteKind) !void {
    const res = mdl.residues.items[residue_idx];
    const chain = mdl.chains.items[res.chain_idx];
    try writer.writeAll("\"id\": {\"chain_label\": \"");
    try writeJsonString(writer, chain.labelSlice());
    try writer.writeAll("\", \"chain_auth\": \"");
    try writeJsonString(writer, chain.authSlice());
    try writer.writeAll("\", \"comp_id\": \"");
    try writeJsonString(writer, res.compIdSlice());
    try writer.print("\", \"seq_id\": {d}, \"auth_seq_id\": {d}, \"ins_code\": ", .{ res.seq_id, res.auth_seq_id });
    if (res.ins_code == ' ') {
        try writer.writeAll("null");
    } else {
        try writer.print("\"{c}\"", .{res.ins_code});
    }
    try writer.writeAll(", \"site\": \"");
    try writer.writeAll(kind.name());
    try writer.writeAll("\"}");
}

const RecommendedState = struct {
    macrostate: []const u8,
    tautomer: ?[]const u8,
    formal_charge: i8,
};

fn recommendedState(kind: pka.SiteKind, predicted_pka: f32, ph: f32) RecommendedState {
    const protonated = predicted_pka >= ph;
    if (kind.isAcid()) {
        return if (protonated)
            .{ .macrostate = "protonated_neutral", .tautomer = if (kind == .asp or kind == .glu) "unspecified_carboxyl_oxygen" else null, .formal_charge = 0 }
        else
            .{ .macrostate = "deprotonated_anion", .tautomer = null, .formal_charge = -1 };
    }
    if (protonated) {
        return .{ .macrostate = "protonated_cation", .tautomer = if (kind == .his) "HIP" else null, .formal_charge = 1 };
    }
    return .{ .macrostate = "neutral", .tautomer = if (kind == .his) "HID_OR_HIE" else null, .formal_charge = 0 };
}

fn writeJsonString(writer: anytype, value: []const u8) !void {
    for (value) |c| switch (c) {
        '"' => try writer.writeAll("\\\""),
        '\\' => try writer.writeAll("\\\\"),
        '\n' => try writer.writeAll("\\n"),
        '\r' => try writer.writeAll("\\r"),
        '\t' => try writer.writeAll("\\t"),
        else => if (c < 0x20) try writer.print("\\u{x:0>4}", .{c}) else try writer.writeByte(c),
    };
}

const testing = std.testing;

test "prediction report escapes input and carries stable schema" {
    var out: std.Io.Writer.Allocating = .init(testing.allocator);
    defer out.deinit();
    try writePredictionReport(&out.writer, "a\"b.cif", &.{});
    const json = out.writer.buffered();
    try testing.expect(std.mem.find(u8, json, "zhyd.protonation.v1") != null);
    try testing.expect(std.mem.find(u8, json, "a\\\"b.cif") != null);
}

test "prediction report distinguishes empirical parity and zhyd modes" {
    const mmcif = @import("../mmcif.zig");
    var mdl = try mmcif.parseModel(testing.allocator, @embedFile("../test_data/his.cif"));
    defer mdl.deinit();
    for (mdl.residues.items) |*res| res.entity_type = .polymer;
    const Case = struct { mode: pka_empirical.Mode, status: []const u8, version: []const u8 };
    for ([_]Case{
        .{ .mode = .parity, .status = "experimental_parity_with_propka3_model", .version = "\"zhyd-empirical-v1-parity\"" },
        .{ .mode = .zhyd, .status = "experimental_unvalidated", .version = "\"zhyd-empirical-v1\"" },
    }) |c| {
        var prediction = try pka_empirical.predict(testing.allocator, &mdl, .{ .mode = c.mode });
        defer prediction.deinit();
        var out: std.Io.Writer.Allocating = .init(testing.allocator);
        defer out.deinit();
        try writePredictionReport(&out.writer, "his.cif", &.{.{ .model_num = 1, .model = &mdl, .prediction = &prediction }});
        const json = out.writer.buffered();
        try testing.expect(std.mem.find(u8, json, c.status) != null);
        try testing.expect(std.mem.find(u8, json, c.version) != null);
    }
}

test "applied state report records geometry-selected histidine tautomer" {
    var mdl = model_mod.Model.init(testing.allocator);
    defer mdl.deinit();
    var chain = model_mod.Chain{};
    chain.setLabelAsymId("A");
    chain.setAuthAsymId("A");
    chain.residue_end = 1;
    try mdl.chains.append(testing.allocator, chain);
    var res = model_mod.Residue{ .chain_idx = 0, .seq_id = 1, .auth_seq_id = 1, .atom_end = 1, .entity_type = .polymer };
    res.setCompId("HIS");
    try mdl.residues.append(testing.allocator, res);
    var hd1 = model_mod.Atom{ .pos = .{ .x = 1, .y = 2, .z = 3 }, .residue_idx = 0, .is_hydrogen = true };
    hd1.setName("HD1");
    try mdl.atoms.append(testing.allocator, hd1);
    const reports = [_]StateReportModel{.{ .model_num = 1, .model = &mdl }};
    var out: std.Io.Writer.Allocating = .init(testing.allocator);
    defer out.deinit();
    try writeStateReport(&out.writer, "his.cif", "his-out.cif", &reports, null, .auto);
    const json = out.writer.buffered();
    try testing.expect(std.mem.find(u8, json, "\"tautomer\": \"HID\"") != null);
    try testing.expect(std.mem.find(u8, json, "\"source\": \"geometry_selected\"") != null);
}

test "applied state report marks disulfide cysteines as covalently bonded" {
    const mmcif = @import("../mmcif.zig");
    var mdl = try mmcif.parseModel(testing.allocator, @embedFile("../test_data/disulfide.cif"));
    defer mdl.deinit();
    const reports = [_]StateReportModel{.{ .model_num = 1, .model = &mdl }};
    var out: std.Io.Writer.Allocating = .init(testing.allocator);
    defer out.deinit();
    try writeStateReport(&out.writer, "disulfide.cif", null, &reports, null, .auto);
    const json = out.writer.buffered();
    try testing.expectEqual(@as(usize, 2), std.mem.count(u8, json, "\"macrostate\": \"disulfide_bonded\""));
    try testing.expect(std.mem.find(u8, json, "\"protonated_neutral\"") == null);
}
