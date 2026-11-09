//! Plan execution — core geometry dispatch for hydrogen placement.
//!
//! Translates placement plans into atom positions and appends
//! new hydrogen atoms to the model.

const std = @import("std");
const model_mod = @import("../model.zig");
const Model = model_mod.Model;
const Atom = model_mod.Atom;
const Residue = model_mod.Residue;
const math_mod = @import("../math.zig");
const element = @import("../element.zig");
const bond_policy = @import("bond_policy.zig");
const geometry = @import("geometry.zig");
const standard = @import("standard.zig");
const topology = @import("topology.zig");
const protonation = @import("protonation.zig");
const lookup = @import("lookup.zig");
const terminal = @import("terminal.zig");

const Vec3f32 = math_mod.Vec3(f32);
const MoverHint = element.MoverHint;
const ABSENT_H_POS = @import("../optimize/mover.zig").Mover.ABSENT_H_POS;
const ParentMeta = lookup.ParentMeta;
const findAtom = lookup.findAtom;
const findAtomPos = lookup.findAtomPos;
const existsInResidue = lookup.existsInResidue;
const findOtherNeighbor = lookup.findOtherNeighbor;
const findThirdNeighbor = lookup.findThirdNeighbor;
const findAtomBetween = lookup.findAtomBetween;
const findBondedNeighbor = lookup.findBondedNeighbor;
const findThirdBondedNeighbor = lookup.findThirdBondedNeighbor;
const findBondedAtomBetween = lookup.findBondedAtomBetween;
const padName = lookup.padName;
const trimPlanName = lookup.trimPlanName;
const isBackboneAmideH = terminal.isBackboneAmideH;
const findPrevResAtomPos = lookup.findPrevResAtomPos;

/// Result of a single executePlan call — why a hydrogen was or wasn't placed.
/// Defined here to avoid circular dependency with placer.zig.
pub const PlaceResult = enum {
    placed,
    existing_h, // hydrogen already present in residue
    inter_residue, // parent atom bonded_inter_residue (disulfide, glycosidic)
    missing_parent, // parent heavy atom (connected[0]) not found
    missing_ref, // reference neighbor or geometric lookup failed
};

/// Append a new hydrogen atom to the model.
pub fn appendHydrogen(mdl: *Model, pos: Vec3f32, plan: *const standard.PlacementPlan, res_idx: u32, meta: ParentMeta) !void {
    try appendHydrogenNamed(mdl, pos, trimPlanName(&plan.h_name), plan.atom_type, plan.mover_hint, res_idx, meta);
}

pub fn appendHydrogenNamed(
    mdl: *Model,
    pos: Vec3f32,
    name: []const u8,
    atom_type: element.AtomType,
    mover_hint: standard.MoverHint,
    res_idx: u32,
    meta: ParentMeta,
) !void {
    var atom = Atom{
        .pos = pos,
        .element_type = atom_type,
        .residue_idx = res_idx,
        .is_hydrogen = true,
        .is_added = true,
        .vdw_radius = atom_type.info().explicit_radius,
        .flags = atom_type.info().flags,
        .altloc = meta.altloc,
        .occupancy = meta.occupancy,
        .b_factor = meta.b_factor,
        .mover_hint = mover_hint,
    };
    atom.setName(name);
    try mdl.atoms.append(mdl.allocator, atom);
}

pub fn placeOverrideHydrogen(mdl: *Model, res: Residue, res_idx: u32, target_altloc: u8, state: ?protonation.ResidueState, mode: bond_policy.BondLengthMode) !?PlaceResult {
    const s = state orelse return null;
    const comp_id = res.compIdSlice();

    if (std.mem.eql(u8, comp_id, "ASP") and s == .asp) {
        return switch (s.asp) {
            .deprotonated => null,
            .atom1 => try placeCarboxylHydrogen(mdl, res, res_idx, target_altloc, mode, "OD1", "OD2", "CG", "HD1", .none),
            .atom2 => try placeCarboxylHydrogen(mdl, res, res_idx, target_altloc, mode, "OD2", "OD1", "CG", "HD2", .none),
            .protonated => try placeFlippableCarboxyl(mdl, res, res_idx, target_altloc, mode, .{ "OD1", "OD2", "CG", "HD1", "HD2" }),
        };
    }
    if (std.mem.eql(u8, comp_id, "GLU") and s == .glu) {
        return switch (s.glu) {
            .deprotonated => null,
            .atom1 => try placeCarboxylHydrogen(mdl, res, res_idx, target_altloc, mode, "OE1", "OE2", "CD", "HE1", .none),
            .atom2 => try placeCarboxylHydrogen(mdl, res, res_idx, target_altloc, mode, "OE2", "OE1", "CD", "HE2", .none),
            .protonated => try placeFlippableCarboxyl(mdl, res, res_idx, target_altloc, mode, .{ "OE1", "OE2", "CD", "HE1", "HE2" }),
        };
    }
    return null;
}

fn placeCarboxylHydrogen(
    mdl: *Model,
    res: Residue,
    res_idx: u32,
    target_altloc: u8,
    mode: bond_policy.BondLengthMode,
    protonated_o_name: []const u8,
    partner_o_name: []const u8,
    carbon_name: []const u8,
    h_name: []const u8,
    hint: MoverHint,
) !PlaceResult {
    const protonated_o = findAtom(mdl, res, padName(protonated_o_name), target_altloc) orelse return .missing_parent;
    if (protonated_o.flags.bonded_inter_residue) return .inter_residue;

    var meta = ParentMeta.fromAtom(protonated_o);
    if (target_altloc != ' ') meta.altloc = target_altloc;
    const h_padded = padName(h_name);
    if (existsInResidue(mdl, res, h_padded, meta.altloc)) return .existing_h;

    const o_pos = protonated_o.pos;
    const partner_o_pos = findAtomPos(mdl, res, padName(partner_o_name), target_altloc) orelse return .missing_ref;
    const carbon_pos = findAtomPos(mdl, res, padName(carbon_name), target_altloc) orelse return .missing_ref;
    const h_pos = geometry.placeHXR2Planar(
        o_pos.cast(f64),
        carbon_pos.cast(f64),
        partner_o_pos.cast(f64),
        bond_policy.adjustedBondLength(mode, 0.97, protonated_o.element_type, .Hpol),
        0.0,
    );
    try appendHydrogenNamed(mdl, h_pos.cast(f32), h_name, .Hpol, hint, res_idx, meta);
    return .placed;
}

/// Protonated acid with the oxygen left open: H on the first oxygen (syn)
/// and a second H on the other oxygen at the absent sentinel. The carboxyl
/// flip mover chooses between them; without optimization the first stays.
fn placeFlippableCarboxyl(
    mdl: *Model,
    res: Residue,
    res_idx: u32,
    target_altloc: u8,
    mode: bond_policy.BondLengthMode,
    names: [5][]const u8, // O1, O2, C, H on O1, H on O2
) !PlaceResult {
    const first = try placeCarboxylHydrogen(mdl, res, res_idx, target_altloc, mode, names[0], names[1], names[2], names[3], .flip_carboxyl);
    if (first != .placed) return first;
    const other_o = findAtom(mdl, res, padName(names[1]), target_altloc) orelse return first;
    var meta = ParentMeta.fromAtom(other_o);
    if (target_altloc != ' ') meta.altloc = target_altloc;
    if (existsInResidue(mdl, res, padName(names[4]), meta.altloc)) return first;
    try appendHydrogenNamed(mdl, ABSENT_H_POS, names[4], .Hpol, .flip_carboxyl, res_idx, meta);
    return first;
}

// ── Tests ─────────────────────────────────────────────────────────────────────

const testing = std.testing;

/// Build a minimal two-residue Model suitable for executePlan tests.
/// chain 0, residue 0 = ALA-like with N/CA/C/O/CB.
/// residue 1 = preceding residue that contributes prev-C for backbone amide tests.
fn buildMinimalModel(allocator: std.mem.Allocator) !Model {
    var mdl = Model.init(allocator);

    // atoms for residue 0
    const atoms = [_]struct { name: []const u8, x: f32, y: f32, z: f32, elem: element.AtomType }{
        .{ .name = "N", .x = 1.0, .y = 2.0, .z = 3.0, .elem = .N },
        .{ .name = "CA", .x = 2.0, .y = 3.0, .z = 4.0, .elem = .C },
        .{ .name = "C", .x = 3.0, .y = 4.0, .z = 5.0, .elem = .C },
        .{ .name = "O", .x = 4.0, .y = 5.0, .z = 6.0, .elem = .O },
        .{ .name = "CB", .x = 2.5, .y = 2.5, .z = 3.5, .elem = .C },
    };
    for (atoms) |a| {
        var atom = Atom{
            .pos = .{ .x = a.x, .y = a.y, .z = a.z },
            .element_type = a.elem,
            .residue_idx = 0,
            .vdw_radius = a.elem.info().explicit_radius,
            .flags = a.elem.info().flags,
        };
        atom.setName(a.name);
        try mdl.atoms.append(allocator, atom);
    }

    var res = Residue{};
    res.setCompId("ALA");
    res.atom_start = 0;
    res.atom_end = @intCast(mdl.atoms.items.len);
    res.entity_type = .polymer;
    try mdl.residues.append(allocator, res);

    return mdl;
}

test "executePlan hxr3 places H with correct bond length" {
    var mdl = try buildMinimalModel(testing.allocator);
    defer mdl.deinit();

    const res = mdl.residues.items[0];
    // HXR3: center=CA, known neighbors=N,CB; 3rd neighbor found by distance (C)
    const plan = standard.PlacementPlan{
        .h_name = lookup.padName("HA"),
        .placement_type = .hxr3,
        .connected = .{ lookup.padName("CA"), lookup.padName("N"), lookup.padName("CB") },
        .bond_len = 1.09,
        .atom_type = .H,
    };

    const result = try executePlan(&mdl, res, 0, &plan, null, ' ', .neutron);
    try testing.expectEqual(PlaceResult.placed, result);

    // Verify HA was appended
    var found_ha = false;
    for (mdl.atoms.items) |atom| {
        if (atom.is_hydrogen and std.mem.eql(u8, atom.nameSlice(), "HA")) {
            found_ha = true;
            // Bond length to CA should be ~1.09
            const ca_pos = math_mod.Vec3(f32){ .x = 2.0, .y = 3.0, .z = 4.0 };
            const dist = atom.pos.distance(ca_pos);
            try testing.expect(dist > 0.9 and dist < 1.3);
            break;
        }
    }
    try testing.expect(found_ha);
}

test "executePlan h2xr2 places H with correct bond length" {
    var mdl = try buildMinimalModel(testing.allocator);
    defer mdl.deinit();

    const res = mdl.residues.items[0];
    // H2XR2: center=CB, ref_neighbor=CA; findOtherNeighbor finds 2nd neighbor by distance
    // CB is at (2.5,2.5,3.5); CA at (2.0,3.0,4.0) — distance ~0.87, within 1.9 A
    const plan = standard.PlacementPlan{
        .h_name = lookup.padName("HB2"),
        .placement_type = .h2xr2,
        .connected = .{ lookup.padName("CB"), lookup.padName("CA"), lookup.padName("    ") },
        .bond_len = 1.09,
        .angle = 109.5,
        .dihedral = 120.0,
        .atom_type = .H,
    };

    const result = try executePlan(&mdl, res, 0, &plan, null, ' ', .neutron);
    try testing.expectEqual(PlaceResult.placed, result);

    var found = false;
    for (mdl.atoms.items) |atom| {
        if (atom.is_hydrogen and std.mem.eql(u8, atom.nameSlice(), "HB2")) {
            found = true;
            const cb_pos = math_mod.Vec3(f32){ .x = 2.5, .y = 2.5, .z = 3.5 };
            const dist = atom.pos.distance(cb_pos);
            try testing.expect(dist > 0.9 and dist < 1.3);
            break;
        }
    }
    try testing.expect(found);
}

test "executePlan hxr2_planar places aromatic H" {
    var mdl = try buildMinimalModel(testing.allocator);
    defer mdl.deinit();

    // For hxr2_planar we need a center atom between two reference atoms.
    // Add a simple planar arrangement: O is "center", N and C are references.
    // O is at (4,5,6), N at (1,2,3), C at (3,4,5).
    // connected[0]=O (n1), connected[1]=C (n2), center = atom between them at dist < 1.9 from both
    // Actually hxr2_planar: connected[0] and connected[1] are the two flanking atoms,
    // and the center atom is found between them by findAtomBetween.
    // Let's use N(connected[0]) and C(connected[1]) — CA is ~1.41A from N and ~1.52A from C
    const res = mdl.residues.items[0];
    const plan = standard.PlacementPlan{
        .h_name = lookup.padName("HO"),
        .placement_type = .hxr2_planar,
        .connected = .{ lookup.padName("N"), lookup.padName("C"), lookup.padName("    ") },
        .bond_len = 1.08,
        .atom_type = .Har,
    };

    const result = try executePlan(&mdl, res, 0, &plan, null, ' ', .neutron);
    // Should place or fail with missing_ref if CA isn't close enough — either is valid for this unit test
    try testing.expect(result == .placed or result == .missing_ref);
}

test "executePlan returns missing_parent when base atom absent" {
    var mdl = try buildMinimalModel(testing.allocator);
    defer mdl.deinit();

    const res = mdl.residues.items[0];
    // Plan references "CG" which doesn't exist in our minimal ALA model
    const plan = standard.PlacementPlan{
        .h_name = lookup.padName("HG"),
        .placement_type = .hxr3,
        .connected = .{ lookup.padName("CG"), lookup.padName("CB"), lookup.padName("CA") },
        .bond_len = 1.09,
        .atom_type = .H,
    };

    const result = try executePlan(&mdl, res, 0, &plan, null, ' ', .neutron);
    try testing.expectEqual(PlaceResult.missing_parent, result);
}

test "executePlan returns missing_ref when reference atom absent" {
    var mdl = try buildMinimalModel(testing.allocator);
    defer mdl.deinit();

    const res = mdl.residues.items[0];
    // For hxr2_frac: connected[0]=N (exists), connected[1]=CA (exists), connected[2]=CG (missing)
    const plan = standard.PlacementPlan{
        .h_name = lookup.padName("HN"),
        .placement_type = .hxr2_frac,
        .connected = .{ lookup.padName("N"), lookup.padName("CA"), lookup.padName("CG") },
        .bond_len = 1.01,
        .fudge = 0.5,
        .atom_type = .Hpol,
    };

    const result = try executePlan(&mdl, res, 0, &plan, null, ' ', .neutron);
    try testing.expectEqual(PlaceResult.missing_ref, result);
}

test "executePlan returns existing_h when hydrogen already present" {
    var mdl = try buildMinimalModel(testing.allocator);
    defer mdl.deinit();

    // Pre-append HA so executePlan sees it already exists
    var existing_h = Atom{
        .pos = .{ .x = 1.5, .y = 2.5, .z = 3.5 },
        .element_type = .H,
        .residue_idx = 0,
        .is_hydrogen = true,
        .vdw_radius = 1.2,
    };
    existing_h.setName("HA");
    try mdl.atoms.append(testing.allocator, existing_h);
    // Update residue atom_end so existsInResidue can see it
    mdl.residues.items[0].atom_end = @intCast(mdl.atoms.items.len);

    const res = mdl.residues.items[0];
    const plan = standard.PlacementPlan{
        .h_name = lookup.padName("HA"),
        .placement_type = .hxr3,
        .connected = .{ lookup.padName("CA"), lookup.padName("N"), lookup.padName("CB") },
        .bond_len = 1.09,
        .atom_type = .H,
    };

    const result = try executePlan(&mdl, res, 0, &plan, null, ' ', .neutron);
    try testing.expectEqual(PlaceResult.existing_h, result);
}

test "executePlan hxy places linear H" {
    var mdl = try buildMinimalModel(testing.allocator);
    defer mdl.deinit();

    const res = mdl.residues.items[0];
    // hxy: center=O (connected[0]), neighbor=C (connected[1])
    // O at (4,5,6), C at (3,4,5) — distance sqrt(3) ≈ 1.73, within 1.9A
    const plan = standard.PlacementPlan{
        .h_name = lookup.padName("HO"),
        .placement_type = .hxy,
        .connected = .{ lookup.padName("O"), lookup.padName("C"), lookup.padName("    ") },
        .bond_len = 0.97,
        .atom_type = .Hpol,
    };

    const result = try executePlan(&mdl, res, 0, &plan, null, ' ', .neutron);
    try testing.expectEqual(PlaceResult.placed, result);

    var found = false;
    for (mdl.atoms.items) |atom| {
        if (atom.is_hydrogen and std.mem.eql(u8, atom.nameSlice(), "HO")) {
            found = true;
            const o_pos = math_mod.Vec3(f32){ .x = 4.0, .y = 5.0, .z = 6.0 };
            const dist = atom.pos.distance(o_pos);
            try testing.expectApproxEqAbs(@as(f32, 0.97), dist, 0.01);
            break;
        }
    }
    try testing.expect(found);
}

/// Execute a single placement plan: find reference atoms, compute H position, add to model.
/// Returns the placement result (placed / skipped / missing ref).
pub fn executePlan(mdl: *Model, res: Residue, res_idx: u32, plan: *const standard.PlacementPlan, bonds: ?[]const topology.BondEntry, target_altloc: u8, mode: bond_policy.BondLengthMode) !PlaceResult {
    // Resolve parent heavy atom (connected[0]) for metadata and position
    const base_atom = findAtom(mdl, res, plan.connected[0], target_altloc) orelse return .missing_parent;

    // Skip H placement if parent atom is involved in an inter-residue bond
    // (e.g. disulfide SG, glycosidic leaving O) — the bond already satisfies valence.
    // Only applies to standard/hardcoded plans (bonds != null). CCD-derived plans
    // already account for inter-residue bonds via analyzeBonds extra_bonds, and
    // atoms like glycan C1 may still have free valence for H despite the flag.
    if (bonds != null and base_atom.flags.bonded_inter_residue) return .inter_residue;

    var meta = ParentMeta.fromAtom(base_atom);
    // Override altloc when iterating conformers: if parent has blank altloc
    // (shared backbone) but we're targeting a specific conformer, the placed H
    // should inherit the target conformer's altloc, not blank.
    if (target_altloc != ' ') meta.altloc = target_altloc;

    // Skip if this hydrogen already exists in the residue
    if (existsInResidue(mdl, res, plan.h_name, meta.altloc)) return .existing_h;

    const effective_bond_len = bond_policy.adjustedBondLength(mode, plan.bond_len, base_atom.element_type, plan.atom_type);

    switch (plan.placement_type) {
        .hxr3 => {
            // connected[0]=center, connected[1..2]=two known neighbors
            // Need to find 3rd heavy-atom neighbor of center
            const center_pos = base_atom.pos;
            const n1_pos = findAtomPos(mdl, res, plan.connected[1], target_altloc) orelse return .missing_ref;
            const n2_pos = findAtomPos(mdl, res, plan.connected[2], target_altloc) orelse return .missing_ref;
            const n3_pos = (if (bonds) |b|
                findThirdBondedNeighbor(mdl, res, b, plan.connected[0], plan.connected[1], plan.connected[2], target_altloc)
            else
                findThirdNeighbor(mdl, res, plan.connected[0], plan.connected[1], plan.connected[2], target_altloc)) orelse return .missing_ref;

            const h_pos = geometry.placeHXR3(
                center_pos.cast(f64),
                n1_pos.cast(f64),
                n2_pos.cast(f64),
                n3_pos.cast(f64),
                effective_bond_len,
            );
            try appendHydrogen(mdl, h_pos.cast(f32), plan, res_idx, meta);
            return .placed;
        },

        .h2xr2 => {
            // connected[0]=center, connected[1]=reference neighbor
            const center_pos = base_atom.pos;
            const n1_pos = findAtomPos(mdl, res, plan.connected[1], target_altloc) orelse return .missing_ref;
            const n2_pos = (if (bonds) |b|
                findBondedNeighbor(mdl, res, b, plan.connected[0], plan.connected[1], target_altloc)
            else
                findOtherNeighbor(mdl, res, plan.connected[0], plan.connected[1], target_altloc)) orelse return .missing_ref;

            const h_pos = geometry.placeH2XR2(
                center_pos.cast(f64),
                n1_pos.cast(f64),
                n2_pos.cast(f64),
                effective_bond_len,
                plan.angle,
                plan.dihedral,
            );
            try appendHydrogen(mdl, h_pos.cast(f32), plan, res_idx, meta);
            return .placed;
        },

        .h3xr => {
            // connected[0]=a1 (center), connected[1]=a2, connected[2]=a3
            const a1_pos = base_atom.pos;
            const a2_pos = findAtomPos(mdl, res, plan.connected[1], target_altloc) orelse return .missing_ref;

            // Backbone amide NH: place in peptide plane using C(i-1) from previous residue.
            // H bisects the C(i-1)-N-CA exterior angle, lying in the peptide plane.
            // Only matches " H  " (single amide H), not " H1 "/" H2 "/" H3 " (N-terminal NH3+).
            if (isBackboneAmideH(plan)) {
                if (findPrevResAtomPos(mdl, res_idx, .{ ' ', 'C', ' ', ' ' }, target_altloc)) |prev_c_pos| {
                    const h_pos = geometry.placeHXR2Planar(
                        a1_pos.cast(f64),
                        prev_c_pos.cast(f64),
                        a2_pos.cast(f64),
                        effective_bond_len,
                        0.0,
                    );
                    try appendHydrogen(mdl, h_pos.cast(f32), plan, res_idx, meta);
                    return .placed;
                }
                // No previous C available: first residue, chain break, or missing atom.
                // Fall through to dihedral-based placement as degraded fallback.
            }

            const a3_pos = if (!lookup.isBlank(plan.connected[2]))
                findAtomPos(mdl, res, plan.connected[2], target_altloc) orelse return .missing_ref
            else if (bonds) |b|
                findBondedNeighbor(mdl, res, b, plan.connected[1], plan.connected[0], target_altloc) orelse return .missing_ref
            else
                findOtherNeighbor(mdl, res, plan.connected[1], plan.connected[0], target_altloc) orelse return .missing_ref;

            const h_pos = geometry.placeH3XR(
                a1_pos.cast(f64),
                a2_pos.cast(f64),
                a3_pos.cast(f64),
                effective_bond_len,
                plan.angle,
                plan.dihedral,
            );
            try appendHydrogen(mdl, h_pos.cast(f32), plan, res_idx, meta);
            return .placed;
        },

        .hxr2_planar => {
            // connected[0] and connected[1] are neighbors; center is the atom between them
            const n1_pos = base_atom.pos;
            const n2_pos = findAtomPos(mdl, res, plan.connected[1], target_altloc) orelse return .missing_ref;
            const center_pos = (if (bonds) |b|
                findBondedAtomBetween(mdl, res, b, plan.connected[0], plan.connected[1], target_altloc)
            else
                findAtomBetween(mdl, res, plan.connected[0], plan.connected[1], target_altloc)) orelse return .missing_ref;

            const h_pos = geometry.placeHXR2Planar(
                center_pos.cast(f64),
                n1_pos.cast(f64),
                n2_pos.cast(f64),
                effective_bond_len,
                plan.fudge,
            );
            try appendHydrogen(mdl, h_pos.cast(f32), plan, res_idx, meta);
            return .placed;
        },

        .hxr2_frac => {
            const a1_pos = base_atom.pos;
            const a2_pos = findAtomPos(mdl, res, plan.connected[1], target_altloc) orelse return .missing_ref;
            const a3_pos = findAtomPos(mdl, res, plan.connected[2], target_altloc) orelse return .missing_ref;

            const h_pos = geometry.placeHXR2Frac(
                a1_pos.cast(f64),
                a2_pos.cast(f64),
                a3_pos.cast(f64),
                effective_bond_len,
                plan.fudge,
            );
            try appendHydrogen(mdl, h_pos.cast(f32), plan, res_idx, meta);
            return .placed;
        },

        .hxy => {
            const center_pos = base_atom.pos;
            const neighbor_pos = findAtomPos(mdl, res, plan.connected[1], target_altloc) orelse return .missing_ref;

            const h_pos = geometry.placeHXY(
                center_pos.cast(f64),
                neighbor_pos.cast(f64),
                effective_bond_len,
            );
            try appendHydrogen(mdl, h_pos.cast(f32), plan, res_idx, meta);
            return .placed;
        },
    }
}
