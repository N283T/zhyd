//! CCD-derived hydrogen placement plan generation.
//!
//! Derives PlacementPlan from CCD component topology. For each hydrogen atom
//! in the component that is NOT already present in the model, determine the
//! placement type from hybridization (inferred from bond orders).
//! Used as fallback when no hardcoded plan exists (standard, nucleotide, modified).

const std = @import("std");
const ccd = @import("../ccd.zig");
const standard = @import("standard.zig");
const PlacementPlan = standard.PlacementPlan;
const PlacementType = standard.PlacementType;
const MoverHint = standard.MoverHint;
const element_mod = @import("../element.zig");
const math_mod = @import("../math.zig");
const Vec3 = math_mod.Vec3;

// ---------------------------------------------------------------------------
// Bond analysis types
// ---------------------------------------------------------------------------

const Hybridization = enum { sp3, sp2, sp, unknown };

const BondInfo = struct {
    total_bonds: u8,
    heavy_neighbor_count: u8,
    h_neighbor_count: u8,
    has_double: bool,
    has_triple: bool,
    has_aromatic: bool,
    hybridization: Hybridization,
};

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Derive placement plans for non-standard residue hydrogens from CCD component definition.
/// Only generates plans for H atoms not already present in existing_atom_names.
pub fn derivePlans(
    allocator: std.mem.Allocator,
    component: *const ccd.Component,
    existing_atom_names: []const [4]u8,
    inter_residue_atoms: []const [4]u8,
) ![]PlacementPlan {
    var plans = std.ArrayListUnmanaged(PlacementPlan).empty;
    errdefer plans.deinit(allocator);

    for (component.atoms, 0..) |atom, atom_idx| {
        // Skip non-hydrogen atoms
        if (atom.element_symbol[0] != 'H') continue;

        // Skip if already present in model
        if (nameExists(atom.name, existing_atom_names)) continue;

        // Find the heavy atom this H is bonded to
        const heavy_idx = findBondedHeavyAtom(component, @intCast(atom_idx)) orelse continue;
        const heavy_atom = component.atoms[heavy_idx];

        // Skip leaving H atoms whose parent heavy atom is absent from the model
        // or involved in an inter-residue bond. The leaving_atom_flag marks atoms
        // that *may* leave during polymerization (e.g. H2 on backbone N for peptide
        // bond, HO1 on glycosidic O1, HXT on C-terminal OXT). We skip when:
        //  - the parent heavy atom is not in the structure (consumed by polymer bond), or
        //  - the parent has an explicit inter-residue bond.
        // Non-leaving H and leaving H with intact parent are still placed (e.g.
        // HO2-HO6 on sugars where the OH is intact).
        if (atom.leaving and (!nameExists(heavy_atom.name, existing_atom_names) or
            nameExists(heavy_atom.name, inter_residue_atoms))) continue;

        // Count inter-residue bonds on this heavy atom
        var extra: u8 = 0;
        for (inter_residue_atoms) |bonded_name| {
            if (std.mem.eql(u8, &heavy_atom.name, &bonded_name)) extra += 1;
        }
        // Count bonds on heavy atom and determine hybridization.
        // Use model-aware neighbor count: CCD neighbors that don't exist in the model
        // (e.g. leaving atoms like O1 in glycosidic bonds) are excluded so that the
        // placement type matches the actual number of available reference atoms.
        const missing_neighbors = countMissingNeighbors(component, @intCast(heavy_idx), existing_atom_names);
        const bond_info = analyzeBonds(component, @intCast(heavy_idx), extra, missing_neighbors);

        // Determine placement type from hybridization
        const plan = deriveSinglePlan(component, atom, @intCast(heavy_idx), heavy_atom, bond_info, existing_atom_names) orelse continue;
        try plans.append(allocator, plan);
    }

    const result = try allocator.dupe(PlacementPlan, plans.items);
    plans.deinit(allocator);
    return result;
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Check if atom name is in existing names list.
fn nameExists(name: [4]u8, list: []const [4]u8) bool {
    for (list) |existing| {
        if (std.mem.eql(u8, &name, &existing)) return true;
    }
    return false;
}

/// Find the non-H atom bonded to this H atom.
fn findBondedHeavyAtom(component: *const ccd.Component, h_idx: u16) ?u16 {
    for (component.bonds) |bond| {
        if (bond.atom_idx_1 == h_idx) {
            const other = bond.atom_idx_2;
            if (component.atoms[other].element_symbol[0] != 'H') return other;
        }
        if (bond.atom_idx_2 == h_idx) {
            const other = bond.atom_idx_1;
            if (component.atoms[other].element_symbol[0] != 'H') return other;
        }
    }
    return null;
}

/// Maximum valence for common elements; returns 0 for unknown/metals (skip check).
fn maxValence(element_symbol: [2]u8) u8 {
    const e0 = element_symbol[0];
    const e1 = element_symbol[1];
    if (e0 == 'H' and e1 == ' ') return 1;
    if (e0 == 'C' and e1 == ' ') return 4;
    if (e0 == 'C' and e1 == 'l') return 1;
    if (e0 == 'N' and e1 == ' ') return 4;
    if (e0 == 'O' and e1 == ' ') return 2;
    if (e0 == 'S' and e1 == ' ') return 6;
    if (e0 == 'S' and e1 == 'e') return 2;
    if (e0 == 'P' and e1 == ' ') return 5;
    if (e0 == 'F' and e1 == ' ') return 1;
    if (e0 == 'B' and e1 == 'r') return 1;
    if (e0 == 'B' and e1 == ' ') return 3;
    if (e0 == 'I' and e1 == ' ') return 1;
    return 0;
}

/// Bond order as a valence contribution.
fn bondValence(order: ccd.BondOrder) u8 {
    return switch (order) {
        .double => 2,
        .triple => 3,
        else => 1,
    };
}

/// Count heavy neighbors in CCD that are NOT present in the model atom list.
/// Used to adjust placement type when leaving atoms (e.g. O1 in glycans) are absent.
fn countMissingNeighbors(component: *const ccd.Component, atom_idx: u16, existing_atom_names: []const [4]u8) u8 {
    var missing: u8 = 0;
    for (component.bonds) |bond| {
        const neighbor_idx: ?u16 = if (bond.atom_idx_1 == atom_idx)
            bond.atom_idx_2
        else if (bond.atom_idx_2 == atom_idx)
            bond.atom_idx_1
        else
            null;

        if (neighbor_idx) |ni| {
            const neighbor = component.atoms[ni];
            if (neighbor.element_symbol[0] != 'H' and !nameExists(neighbor.name, existing_atom_names)) {
                missing += 1;
            }
        }
    }
    return missing;
}

/// Analyze bonds on the given atom index.
/// extra_bonds: number of inter-residue single bonds not present in the CCD template.
/// missing_neighbors: number of CCD heavy neighbors not present in the model.
fn analyzeBonds(component: *const ccd.Component, atom_idx: u16, extra_bonds: u8, missing_neighbors: u8) BondInfo {
    var info = BondInfo{
        .total_bonds = 0,
        .heavy_neighbor_count = 0,
        .h_neighbor_count = 0,
        .has_double = false,
        .has_triple = false,
        .has_aromatic = false,
        .hybridization = .unknown,
    };

    var template_valence: u8 = 0;
    for (component.bonds) |bond| {
        const neighbor_idx: ?u16 = if (bond.atom_idx_1 == atom_idx)
            bond.atom_idx_2
        else if (bond.atom_idx_2 == atom_idx)
            bond.atom_idx_1
        else
            null;

        if (neighbor_idx) |ni| {
            info.total_bonds += 1;
            template_valence += bondValence(bond.order);
            if (component.atoms[ni].element_symbol[0] == 'H') {
                info.h_neighbor_count += 1;
            } else {
                info.heavy_neighbor_count += 1;
            }
            switch (bond.order) {
                .double => info.has_double = true,
                .triple => info.has_triple = true,
                .aromatic => info.has_aromatic = true,
                else => {},
            }
        }
    }

    // Account for inter-residue bonds not in the CCD template
    info.total_bonds += extra_bonds;
    info.heavy_neighbor_count += extra_bonds;

    // Subtract CCD neighbors absent from the model (e.g. leaving atoms)
    if (missing_neighbors > 0) {
        info.heavy_neighbor_count -|= missing_neighbors;
        info.total_bonds -|= missing_neighbors;
    }

    // Note: total_valence is intentionally NOT adjusted for missing neighbors.
    // It reflects the CCD template's full valence for hybridization determination,
    // since the missing atom's bond (typically single) doesn't change the electronic
    // hybridization of the center atom. The neighbor count adjustment above only
    // affects which placement geometry is selected (hxr3 vs h2xr2 etc.).
    const total_valence = template_valence + extra_bonds;
    const max_val = maxValence(component.atoms[atom_idx].element_symbol);

    // Determine hybridization, with valence overflow correction
    if (max_val > 0 and total_valence > max_val) {
        // Inter-residue bond makes template bond orders impossible.
        // Demote hybridization: sp→sp2, sp2→sp3, sp3→sp3.
        if (info.has_triple or (info.has_double and countDoubleBonds(component, atom_idx) >= 2)) {
            info.hybridization = .sp2;
        } else {
            info.hybridization = .sp3;
        }
    } else {
        // Normal hybridization
        if (info.has_triple or (info.has_double and countDoubleBonds(component, atom_idx) >= 2)) {
            info.hybridization = .sp;
        } else if (info.has_double or info.has_aromatic) {
            info.hybridization = .sp2;
        } else {
            info.hybridization = .sp3;
        }
    }

    return info;
}

/// Count the number of double bonds on an atom.
fn countDoubleBonds(component: *const ccd.Component, atom_idx: u16) u8 {
    var count: u8 = 0;
    for (component.bonds) |bond| {
        if ((bond.atom_idx_1 == atom_idx or bond.atom_idx_2 == atom_idx) and bond.order == .double) {
            count += 1;
        }
    }
    return count;
}

/// Select placement type based on hybridization and neighbor configuration.
fn deriveSinglePlan(
    component: *const ccd.Component,
    h_atom: ccd.CompAtom,
    heavy_idx: u16,
    heavy_atom: ccd.CompAtom,
    bond_info: BondInfo,
    existing_atom_names: []const [4]u8,
) ?PlacementPlan {
    // Determine bond length from element type
    const bond_len: f32 = switch (heavy_atom.element_symbol[0]) {
        'O' => 0.97,
        'N' => 1.02,
        'S' => 1.33,
        else => 1.10, // C-H default
    };

    // Determine atom type for the H
    const atom_type: element_mod.AtomType = switch (heavy_atom.element_symbol[0]) {
        'N', 'O', 'S' => .Hpol,
        else => if (bond_info.has_aromatic) .Har else .H,
    };

    // Determine mover hint
    const mover_hint: MoverHint = switch (heavy_atom.element_symbol[0]) {
        'O' => .rotate,
        'S' => .rotate,
        else => .none,
    };

    // Plan type selection uses fallback: when inter-residue extras inflate
    // heavy_neighbor_count beyond what findHeavyNeighborNamesFiltered can
    // satisfy from the CCD template, fall through to a lower-arity plan.
    switch (bond_info.hybridization) {
        .sp3 => {
            if (bond_info.heavy_neighbor_count >= 3) {
                // HXR3: tetrahedral with 3 heavy neighbors.
                // Only use hxr3 if CCD template has 3+ heavy neighbors in the model,
                // because executePlan needs to geometrically find the 3rd neighbor
                // within the same residue. When the 3rd is an inter-residue partner
                // (not in CCD), it won't be found — fall through to h2xr2 instead.
                if (findHeavyNeighborNamesFiltered(component, heavy_idx, 3, existing_atom_names)) |refs| {
                    return PlacementPlan{
                        .h_name = h_atom.name,
                        .placement_type = .hxr3,
                        .connected = .{ heavy_atom.name, refs[0], refs[1] },
                        .bond_len = bond_len,
                        .atom_type = atom_type,
                    };
                }
                // Fallback: not enough CCD refs for hxr3, try h2xr2
            }
            if (bond_info.heavy_neighbor_count >= 2) {
                // H2XR2: 2 heavy neighbors, dihedral-controlled
                if (findHeavyNeighborNamesFiltered(component, heavy_idx, 2, existing_atom_names)) |refs| {
                    // Placed as h3xr with the second neighbor as explicit dihedral
                    // reference, so execution does not have to rediscover it.
                    // Two H on a tetrahedral center sit at +-120 deg from it.
                    const g = dihedralGeometry(component, h_atom, heavy_idx, heavy_atom, refs[0], refs[1], 109.5, &.{ 120.0, -120.0 });
                    return PlacementPlan{
                        .h_name = h_atom.name,
                        .placement_type = .h3xr,
                        .connected = .{ heavy_atom.name, refs[0], refs[1] },
                        .bond_len = bond_len,
                        .angle = g.angle,
                        .dihedral = g.dihedral,
                        .atom_type = atom_type,
                    };
                }
                // Fallback: not enough CCD refs for h2xr2, try h3xr
            }
            if (bond_info.heavy_neighbor_count >= 1) {
                // H3XR: dihedral-controlled (e.g., methyl, or backbone H with inter-residue bond)
                const refs = findHeavyNeighborNamesFiltered(component, heavy_idx, 1, existing_atom_names) orelse return null;
                const second_ref = findSecondReferenceFiltered(component, heavy_idx, refs[0], existing_atom_names) orelse return null;
                // Up to three H staggered around the X-R bond.
                const g = dihedralGeometry(component, h_atom, heavy_idx, heavy_atom, refs[0], second_ref, 109.5, &.{ 180.0, 60.0, -60.0 });
                return PlacementPlan{
                    .h_name = h_atom.name,
                    .placement_type = .h3xr,
                    .connected = .{ heavy_atom.name, refs[0], second_ref },
                    .bond_len = bond_len,
                    .angle = g.angle,
                    .dihedral = g.dihedral,
                    .atom_type = atom_type,
                    .mover_hint = if (bond_info.h_neighbor_count >= 3) .rotate_methyl else mover_hint,
                };
            }
            return null;
        },
        .sp2 => {
            if (bond_info.heavy_neighbor_count >= 2) {
                // Planar placement
                if (findHeavyNeighborNamesFiltered(component, heavy_idx, 2, existing_atom_names)) |refs| {
                    return PlacementPlan{
                        .h_name = h_atom.name,
                        .placement_type = .hxr2_planar,
                        .connected = .{ refs[0], refs[1], blank },
                        .bond_len = bond_len,
                        .atom_type = atom_type,
                    };
                }
                // Fallback: not enough CCD refs, try dihedral placement
            }
            if (bond_info.heavy_neighbor_count >= 1) {
                // Only one heavy neighbor on sp2 — use dihedral placement
                const refs = findHeavyNeighborNamesFiltered(component, heavy_idx, 1, existing_atom_names) orelse return null;
                const second_ref = findSecondReferenceFiltered(component, heavy_idx, refs[0], existing_atom_names) orelse return null;
                // In the plane of X, R and R2 (e.g. both H of a terminal =CH2).
                const g = dihedralGeometry(component, h_atom, heavy_idx, heavy_atom, refs[0], second_ref, 120.0, &.{ 180.0, 0.0 });
                return PlacementPlan{
                    .h_name = h_atom.name,
                    .placement_type = .h3xr,
                    .connected = .{ heavy_atom.name, refs[0], second_ref },
                    .bond_len = bond_len,
                    .angle = g.angle,
                    .dihedral = g.dihedral,
                    .atom_type = atom_type,
                };
            }
            return null;
        },
        .sp => {
            // Linear placement (HXY)
            if (bond_info.heavy_neighbor_count >= 1) {
                const refs = findHeavyNeighborNamesFiltered(component, heavy_idx, 1, existing_atom_names) orelse return null;
                return PlacementPlan{
                    .h_name = h_atom.name,
                    .placement_type = .hxy,
                    .connected = .{ heavy_atom.name, refs[0], blank },
                    .bond_len = bond_len,
                    .atom_type = atom_type,
                };
            }
            return null;
        },
        .unknown => return null,
    }
}

const blank: [4]u8 = .{ ' ', ' ', ' ', ' ' };

/// Find heavy neighbor names, optionally filtering to only atoms present in model.
fn findHeavyNeighborNamesFiltered(component: *const ccd.Component, atom_idx: u16, comptime max: u8, existing_atom_names: ?[]const [4]u8) ?[max][4]u8 {
    var result: [max][4]u8 = undefined;
    var count: u8 = 0;

    for (component.bonds) |bond| {
        const neighbor_idx: ?u16 = if (bond.atom_idx_1 == atom_idx)
            bond.atom_idx_2
        else if (bond.atom_idx_2 == atom_idx)
            bond.atom_idx_1
        else
            null;

        if (neighbor_idx) |ni| {
            if (component.atoms[ni].element_symbol[0] != 'H') {
                // Skip atoms not present in model (e.g. leaving atoms)
                if (existing_atom_names) |names| {
                    if (!nameExists(component.atoms[ni].name, names)) continue;
                }
                if (count < max) {
                    result[count] = component.atoms[ni].name;
                    count += 1;
                    if (count == max) return result;
                }
            }
        }
    }

    if (count >= max) return result;
    return null;
}

/// Find a dihedral reference atom (neighbor of first_ref that isn't the heavy atom).
/// Optionally filters to atoms present in the model.
fn findSecondReferenceFiltered(component: *const ccd.Component, heavy_idx: u16, first_ref_name: [4]u8, existing_atom_names: ?[]const [4]u8) ?[4]u8 {
    // Find index of first_ref
    var first_ref_idx: ?u16 = null;
    for (component.atoms, 0..) |a, i| {
        if (std.mem.eql(u8, &a.name, &first_ref_name)) {
            first_ref_idx = @intCast(i);
            break;
        }
    }
    const ref_idx = first_ref_idx orelse return null;

    // Find a neighbor of ref_idx that is not heavy_idx and not H
    for (component.bonds) |bond| {
        const neighbor_idx: ?u16 = if (bond.atom_idx_1 == ref_idx)
            bond.atom_idx_2
        else if (bond.atom_idx_2 == ref_idx)
            bond.atom_idx_1
        else
            null;

        if (neighbor_idx) |ni| {
            if (ni != heavy_idx and component.atoms[ni].element_symbol[0] != 'H') {
                if (existing_atom_names) |names| {
                    if (!nameExists(component.atoms[ni].name, names)) continue;
                }
                return component.atoms[ni].name;
            }
        }
    }
    return null;
}

/// Bond angle and dihedral for an `h3xr` plan placing `h_atom` on `heavy`
/// (connected = heavy, r, r2).
///
/// When the component has ideal coordinates for all four atoms they are
/// measured there, in the convention of `geometry.placeH3XR`: angle is
/// H-X-R and dihedral is minus the H-X-R-R2 torsion. Otherwise the H gets
/// `fallback_dihedrals[k]`, k being its rank among the H on the same heavy
/// atom, so that several H on one atom never share a position.
fn dihedralGeometry(
    component: *const ccd.Component,
    h_atom: ccd.CompAtom,
    heavy_idx: u16,
    heavy: ccd.CompAtom,
    r: [4]u8,
    r2: [4]u8,
    fallback_angle: f32,
    fallback_dihedrals: []const f32,
) struct { angle: f32, dihedral: f32 } {
    if (idealPos(h_atom)) |h| if (idealPos(heavy)) |x| {
        if (atomByName(component, r)) |ra| if (idealPos(ra)) |rp| {
            if (atomByName(component, r2)) |r2a| if (idealPos(r2a)) |r2p| {
                return .{
                    .angle = @floatCast(math_mod.angle(f64, h, x, rp)),
                    .dihedral = @floatCast(-math_mod.dihedralAngle(f64, h, x, rp, r2p)),
                };
            };
        };
    };
    const k = hRank(component, heavy_idx, h_atom.name);
    return .{ .angle = fallback_angle, .dihedral = fallback_dihedrals[k % fallback_dihedrals.len] };
}

/// Ideal position of a component atom, or null when the component has none
/// (all coordinates zero: missing in the CCD entry, or not provided by the
/// source of the component).
fn idealPos(a: ccd.CompAtom) ?Vec3(f64) {
    if (a.ideal_x == 0 and a.ideal_y == 0 and a.ideal_z == 0) return null;
    return .{ .x = a.ideal_x, .y = a.ideal_y, .z = a.ideal_z };
}

fn atomByName(component: *const ccd.Component, name: [4]u8) ?ccd.CompAtom {
    for (component.atoms) |a| {
        if (std.mem.eql(u8, &a.name, &name)) return a;
    }
    return null;
}

/// Rank of the H named `h_name` among the H bonded to `heavy_idx`, in bond order.
fn hRank(component: *const ccd.Component, heavy_idx: u16, h_name: [4]u8) usize {
    var k: usize = 0;
    for (component.bonds) |bond| {
        const ni: ?u16 = if (bond.atom_idx_1 == heavy_idx) bond.atom_idx_2 else if (bond.atom_idx_2 == heavy_idx) bond.atom_idx_1 else null;
        const idx = ni orelse continue;
        const a = component.atoms[idx];
        if (a.element_symbol[0] != 'H') continue;
        if (std.mem.eql(u8, &a.name, &h_name)) return k;
        k += 1;
    }
    return k;
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

const testing = std.testing;

test "nameExists finds matching names" {
    const list = [_][4]u8{
        .{ 'N', ' ', ' ', ' ' },
        .{ 'C', 'A', ' ', ' ' },
    };
    try testing.expect(nameExists(.{ 'N', ' ', ' ', ' ' }, &list));
    try testing.expect(nameExists(.{ 'C', 'A', ' ', ' ' }, &list));
    try testing.expect(!nameExists(.{ 'H', ' ', ' ', ' ' }, &list));
}

test "derivePlans on empty component returns empty" {
    const comp = ccd.Component{
        .comp_id = "TST",
        .comp_type = "non-polymer",
        .atoms = &.{},
        .bonds = &.{},
    };
    const plans = try derivePlans(testing.allocator, &comp, &.{}, &.{});
    defer testing.allocator.free(plans);
    try testing.expectEqual(@as(usize, 0), plans.len);
}

test "derivePlans skips existing H atoms" {
    // Build a simple component: C bonded to H
    var atoms = [_]ccd.CompAtom{
        .{ .name = .{ 'C', ' ', ' ', ' ' }, .name_len = 1, .element_symbol = .{ 'C', ' ' } },
        .{ .name = .{ 'H', '1', ' ', ' ' }, .name_len = 2, .element_symbol = .{ 'H', ' ' } },
    };
    var bonds = [_]ccd.CompBond{
        .{ .atom_idx_1 = 0, .atom_idx_2 = 1, .order = .single },
    };
    const comp = ccd.Component{
        .comp_id = "TST",
        .comp_type = "non-polymer",
        .atoms = &atoms,
        .bonds = &bonds,
    };

    // H1 already exists
    const existing = [_][4]u8{.{ 'H', '1', ' ', ' ' }};
    const plans = try derivePlans(testing.allocator, &comp, &existing, &.{});
    defer testing.allocator.free(plans);
    try testing.expectEqual(@as(usize, 0), plans.len);
}

test "analyzeBonds detects sp2" {
    var atoms = [_]ccd.CompAtom{
        .{ .name = .{ 'C', '1', ' ', ' ' }, .name_len = 2, .element_symbol = .{ 'C', ' ' } },
        .{ .name = .{ 'C', '2', ' ', ' ' }, .name_len = 2, .element_symbol = .{ 'C', ' ' } },
        .{ .name = .{ 'O', ' ', ' ', ' ' }, .name_len = 1, .element_symbol = .{ 'O', ' ' } },
        .{ .name = .{ 'H', ' ', ' ', ' ' }, .name_len = 1, .element_symbol = .{ 'H', ' ' } },
    };
    var bonds = [_]ccd.CompBond{
        .{ .atom_idx_1 = 0, .atom_idx_2 = 1, .order = .single },
        .{ .atom_idx_1 = 0, .atom_idx_2 = 2, .order = .double },
        .{ .atom_idx_1 = 0, .atom_idx_2 = 3, .order = .single },
    };
    const comp = ccd.Component{
        .comp_id = "TST",
        .comp_type = "non-polymer",
        .atoms = &atoms,
        .bonds = &bonds,
    };
    const info = analyzeBonds(&comp, 0, 0, 0);
    try testing.expectEqual(Hybridization.sp2, info.hybridization);
    try testing.expectEqual(@as(u8, 2), info.heavy_neighbor_count);
    try testing.expectEqual(@as(u8, 1), info.h_neighbor_count);
}

test "analyzeBonds detects sp3" {
    var atoms = [_]ccd.CompAtom{
        .{ .name = .{ 'C', ' ', ' ', ' ' }, .name_len = 1, .element_symbol = .{ 'C', ' ' } },
        .{ .name = .{ 'N', ' ', ' ', ' ' }, .name_len = 1, .element_symbol = .{ 'N', ' ' } },
        .{ .name = .{ 'H', '1', ' ', ' ' }, .name_len = 2, .element_symbol = .{ 'H', ' ' } },
        .{ .name = .{ 'H', '2', ' ', ' ' }, .name_len = 2, .element_symbol = .{ 'H', ' ' } },
        .{ .name = .{ 'H', '3', ' ', ' ' }, .name_len = 2, .element_symbol = .{ 'H', ' ' } },
    };
    var bonds = [_]ccd.CompBond{
        .{ .atom_idx_1 = 0, .atom_idx_2 = 1, .order = .single },
        .{ .atom_idx_1 = 0, .atom_idx_2 = 2, .order = .single },
        .{ .atom_idx_1 = 0, .atom_idx_2 = 3, .order = .single },
        .{ .atom_idx_1 = 0, .atom_idx_2 = 4, .order = .single },
    };
    const comp = ccd.Component{
        .comp_id = "TST",
        .comp_type = "non-polymer",
        .atoms = &atoms,
        .bonds = &bonds,
    };
    const info = analyzeBonds(&comp, 0, 0, 0);
    try testing.expectEqual(Hybridization.sp3, info.hybridization);
    try testing.expectEqual(@as(u8, 1), info.heavy_neighbor_count);
    try testing.expectEqual(@as(u8, 3), info.h_neighbor_count);
}

test "findBondedHeavyAtom returns correct index" {
    var atoms = [_]ccd.CompAtom{
        .{ .name = .{ 'C', ' ', ' ', ' ' }, .name_len = 1, .element_symbol = .{ 'C', ' ' } },
        .{ .name = .{ 'H', ' ', ' ', ' ' }, .name_len = 1, .element_symbol = .{ 'H', ' ' } },
    };
    var bonds = [_]ccd.CompBond{
        .{ .atom_idx_1 = 0, .atom_idx_2 = 1, .order = .single },
    };
    const comp = ccd.Component{
        .comp_id = "TST",
        .comp_type = "non-polymer",
        .atoms = &atoms,
        .bonds = &bonds,
    };
    const result = findBondedHeavyAtom(&comp, 1);
    try testing.expectEqual(@as(?u16, 0), result);
}

test "analyzeBonds demotes sp to sp2 on valence overflow" {
    // C1 has: C2 (triple) → template valence = 3, sp
    // With 2 extra bonds → valence = 5 > 4 → demote sp→sp2
    var atoms = [_]ccd.CompAtom{
        .{ .name = .{ 'C', '1', ' ', ' ' }, .name_len = 2, .element_symbol = .{ 'C', ' ' } },
        .{ .name = .{ 'C', '2', ' ', ' ' }, .name_len = 2, .element_symbol = .{ 'C', ' ' } },
    };
    var bonds = [_]ccd.CompBond{
        .{ .atom_idx_1 = 0, .atom_idx_2 = 1, .order = .triple },
    };
    const comp = ccd.Component{
        .comp_id = "TST",
        .comp_type = "non-polymer",
        .atoms = &atoms,
        .bonds = &bonds,
    };

    const info_normal = analyzeBonds(&comp, 0, 0, 0);
    try testing.expectEqual(Hybridization.sp, info_normal.hybridization);

    const info_overflow = analyzeBonds(&comp, 0, 2, 0);
    try testing.expectEqual(Hybridization.sp2, info_overflow.hybridization);
}

test "analyzeBonds demotes N sp2 to sp3 on valence overflow" {
    // N has: C (double), C2 (single) → template valence = 3, sp2
    // maxValence for N is now 4 (covers NH4+/quaternary N), so 1 extra bond
    // brings valence to 4 which does NOT overflow. 2 extra bonds → valence = 5 > 4
    // → demote sp2→sp3.
    var atoms = [_]ccd.CompAtom{
        .{ .name = .{ 'N', ' ', ' ', ' ' }, .name_len = 1, .element_symbol = .{ 'N', ' ' } },
        .{ .name = .{ 'C', '1', ' ', ' ' }, .name_len = 2, .element_symbol = .{ 'C', ' ' } },
        .{ .name = .{ 'C', '2', ' ', ' ' }, .name_len = 2, .element_symbol = .{ 'C', ' ' } },
    };
    var bonds = [_]ccd.CompBond{
        .{ .atom_idx_1 = 0, .atom_idx_2 = 1, .order = .double },
        .{ .atom_idx_1 = 0, .atom_idx_2 = 2, .order = .single },
    };
    const comp = ccd.Component{
        .comp_id = "TST",
        .comp_type = "non-polymer",
        .atoms = &atoms,
        .bonds = &bonds,
    };

    const info_normal = analyzeBonds(&comp, 0, 0, 0);
    try testing.expectEqual(Hybridization.sp2, info_normal.hybridization);

    // 1 extra bond: valence = 4, not > maxValence(N)=4, no demotion
    const info_no_overflow = analyzeBonds(&comp, 0, 1, 0);
    try testing.expectEqual(Hybridization.sp2, info_no_overflow.hybridization);

    // 2 extra bonds: valence = 5 > maxValence(N)=4, demote sp2→sp3
    const info_overflow = analyzeBonds(&comp, 0, 2, 0);
    try testing.expectEqual(Hybridization.sp3, info_overflow.hybridization);
}

test "analyzeBonds demotes sp2 to sp3 on valence overflow" {
    // C1 has: C2 (single), O (double), H (single) → template valence = 4
    // With 1 extra inter-residue bond → total valence = 5 > 4 → demote sp2→sp3
    var atoms = [_]ccd.CompAtom{
        .{ .name = .{ 'C', '1', ' ', ' ' }, .name_len = 2, .element_symbol = .{ 'C', ' ' } },
        .{ .name = .{ 'C', '2', ' ', ' ' }, .name_len = 2, .element_symbol = .{ 'C', ' ' } },
        .{ .name = .{ 'O', ' ', ' ', ' ' }, .name_len = 1, .element_symbol = .{ 'O', ' ' } },
        .{ .name = .{ 'H', ' ', ' ', ' ' }, .name_len = 1, .element_symbol = .{ 'H', ' ' } },
    };
    var bonds = [_]ccd.CompBond{
        .{ .atom_idx_1 = 0, .atom_idx_2 = 1, .order = .single },
        .{ .atom_idx_1 = 0, .atom_idx_2 = 2, .order = .double },
        .{ .atom_idx_1 = 0, .atom_idx_2 = 3, .order = .single },
    };
    const comp = ccd.Component{
        .comp_id = "TST",
        .comp_type = "non-polymer",
        .atoms = &atoms,
        .bonds = &bonds,
    };

    // Without extra bonds: sp2 (has double bond, valence = 4 which equals max for C)
    const info_normal = analyzeBonds(&comp, 0, 0, 0);
    try testing.expectEqual(Hybridization.sp2, info_normal.hybridization);

    // With 1 extra inter-residue bond: valence = 5 > 4 → demote to sp3
    const info_overflow = analyzeBonds(&comp, 0, 1, 0);
    try testing.expectEqual(Hybridization.sp3, info_overflow.hybridization);
    try testing.expectEqual(@as(u8, 3), info_overflow.heavy_neighbor_count); // 2 template + 1 extra
}

test "dihedrals from ideal coordinates reproduce the ideal H positions" {
    const geometry = @import("geometry.zig");
    // Ethanol-like: C1 methyl (3 H) bonded to C2 (2 H) bonded to O1.
    // Ideal coordinates from a staggered geometry.
    var atoms = [_]ccd.CompAtom{
        .{ .name = .{ 'C', '1', ' ', ' ' }, .name_len = 2, .element_symbol = .{ 'C', ' ' }, .ideal_x = -1.2, .ideal_y = 0.1, .ideal_z = 0.3 },
        .{ .name = .{ 'C', '2', ' ', ' ' }, .name_len = 2, .element_symbol = .{ 'C', ' ' }, .ideal_x = 0.29, .ideal_y = -0.1, .ideal_z = 0.05 },
        .{ .name = .{ 'O', '1', ' ', ' ' }, .name_len = 2, .element_symbol = .{ 'O', ' ' }, .ideal_x = 0.95, .ideal_y = 1.12, .ideal_z = 0.2 },
        .{ .name = .{ 'H', '1', '1', ' ' }, .name_len = 3, .element_symbol = .{ 'H', ' ' }, .ideal_x = -1.67, .ideal_y = -0.84, .ideal_z = 0.05 },
        .{ .name = .{ 'H', '1', '2', ' ' }, .name_len = 3, .element_symbol = .{ 'H', ' ' }, .ideal_x = -1.6, .ideal_y = 0.86, .ideal_z = -0.37 },
        .{ .name = .{ 'H', '1', '3', ' ' }, .name_len = 3, .element_symbol = .{ 'H', ' ' }, .ideal_x = -1.38, .ideal_y = 0.37, .ideal_z = 1.33 },
        .{ .name = .{ 'H', '2', '1', ' ' }, .name_len = 3, .element_symbol = .{ 'H', ' ' }, .ideal_x = 0.47, .ideal_y = -0.4, .ideal_z = -0.98 },
        .{ .name = .{ 'H', '2', '2', ' ' }, .name_len = 3, .element_symbol = .{ 'H', ' ' }, .ideal_x = 0.68, .ideal_y = -0.87, .ideal_z = 0.71 },
    };
    var bonds = [_]ccd.CompBond{
        .{ .atom_idx_1 = 0, .atom_idx_2 = 1, .order = .single },
        .{ .atom_idx_1 = 1, .atom_idx_2 = 2, .order = .single },
        .{ .atom_idx_1 = 0, .atom_idx_2 = 3, .order = .single },
        .{ .atom_idx_1 = 0, .atom_idx_2 = 4, .order = .single },
        .{ .atom_idx_1 = 0, .atom_idx_2 = 5, .order = .single },
        .{ .atom_idx_1 = 1, .atom_idx_2 = 6, .order = .single },
        .{ .atom_idx_1 = 1, .atom_idx_2 = 7, .order = .single },
    };
    const comp = ccd.Component{ .comp_id = "ETX", .comp_type = "non-polymer", .atoms = &atoms, .bonds = &bonds };
    const existing = [_][4]u8{ atoms[0].name, atoms[1].name, atoms[2].name };
    const plans = try derivePlans(testing.allocator, &comp, &existing, &.{});
    defer testing.allocator.free(plans);
    try testing.expectEqual(@as(usize, 5), plans.len);

    for (plans) |plan| {
        try testing.expectEqual(PlacementType.h3xr, plan.placement_type);
        const pos = struct {
            fn of(c: ccd.Component, name: [4]u8) Vec3(f64) {
                const a = atomByName(&c, name).?;
                return .{ .x = a.ideal_x, .y = a.ideal_y, .z = a.ideal_z };
            }
        };
        const h_ideal = pos.of(comp, plan.h_name);
        const x = pos.of(comp, plan.connected[0]);
        const placed = geometry.placeH3XR(x, pos.of(comp, plan.connected[1]), pos.of(comp, plan.connected[2]), h_ideal.distance(x), plan.angle, plan.dihedral);
        try testing.expect(placed.distance(h_ideal) < 1e-3);
    }
}
