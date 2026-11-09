//! Mover generator: scans placed H atoms and creates Mover instances.
//!
//! Inspects `mover_hint` on each added hydrogen atom, resolves axis/center
//! atoms from placement plans, and delegates to rotator constructors.
//! Supports single-H rotators, NH3 rotators, NH2 rotators,
//! amide flippers (ASN/GLN), and histidine ring flippers. Methyl groups get
//! no mover: they stay staggered as placed. The scoring has no torsion term,
//! so letting contacts rotate them only pushed them toward eclipsed, and
//! their thousands of movers merged most of a protein into one clique.

const std = @import("std");
const log = std.log.scoped(.mover_gen);
const model_mod = @import("../model.zig");
const Model = model_mod.Model;
const Atom = model_mod.Atom;
const Residue = model_mod.Residue;
const Vec3 = @import("../math.zig").Vec3(f32);
const mover_mod = @import("mover.zig");
const Mover = mover_mod.Mover;
const rotator = @import("rotator.zig");
const flipper = @import("flipper.zig");
const standard = @import("../place/standard.zig");
const nucleotide = @import("../place/nucleotide.zig");
const modified = @import("../place/modified.zig");
const protonation = @import("../place/protonation.zig");
const bond_policy = @import("../place/bond_policy.zig");
const ccd_mod = @import("../ccd.zig");
const ComponentDict = ccd_mod.ComponentDict;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Trim leading/trailing spaces from a 4-byte padded atom name.
/// Takes a pointer to avoid returning a slice into a stack-local copy.
fn trimName(name: *const [4]u8) []const u8 {
    var start: usize = 0;
    while (start < 4 and name[start] == ' ') start += 1;
    var end: usize = 4;
    while (end > start and name[end - 1] == ' ') end -= 1;
    return name[start..end];
}

/// Find a PlacementPlan for a given H atom name within a comp_id.
fn findPlanForH(comp_id: []const u8, h_name: []const u8) ?*const standard.PlacementPlan {
    const plans = standard.getPlans(comp_id) orelse
        nucleotide.getPlans(comp_id) orelse
        modified.getPlans(comp_id) orelse
        return null;
    for (plans) |*plan| {
        if (std.mem.eql(u8, trimName(&plan.h_name), h_name)) return plan;
    }
    return null;
}

/// Find atom index for a given atom name and altloc within a residue.
///
/// Performs two bounded scans instead of a full O(N) scan over all atoms:
///   - Pass 1: atoms[res.atom_start..res.atom_end] for original atoms (no
///             residue_idx check needed since the range is per-residue).
///   - Pass 2: atoms[original_atom_count..] for added H atoms, filtered by
///             residue_idx (added atoms live after all original atoms and carry
///             a residue_idx back-pointer).
///
/// Prefers an exact altloc match and falls back to blank altloc for shared atoms.
fn findAtomIdx(
    atoms: []const Atom,
    res: Residue,
    residue_idx: u32,
    name: []const u8,
    target_altloc: u8,
    original_atom_count: u32,
) ?u32 {
    var blank_match: ?u32 = null;

    // Pass 1: original atoms for this residue (tight range, no residue_idx filter needed)
    for (atoms[res.atom_start..res.atom_end], res.atom_start..) |*a, i| {
        if (!std.mem.eql(u8, a.nameSlice(), name)) continue;
        if (a.altloc == target_altloc) return @intCast(i);
        if (a.altloc == ' ' and blank_match == null) blank_match = @intCast(i);
    }

    // Pass 2: added H atoms appended after all original atoms
    for (atoms[original_atom_count..], original_atom_count..) |*a, i| {
        if (a.residue_idx != residue_idx) continue;
        if (!std.mem.eql(u8, a.nameSlice(), name)) continue;
        if (a.altloc == target_altloc) return @intCast(i);
        if (a.altloc == ' ' and blank_match == null) blank_match = @intCast(i);
    }

    return blank_match;
}

/// Resolve center and axis atom names for a rotatable H using CCD topology.
/// For a single-H rotator: center = heavy atom bonded to H, axis = heavy atom bonded to center.
fn resolveCcdRotatorAtoms(
    ccd_dict: ?*const ComponentDict,
    comp_id: []const u8,
    h_name: []const u8,
    center_out: *[4]u8,
    axis_out: *[4]u8,
) bool {
    const dict = ccd_dict orelse return false;
    const component = dict.get(comp_id) orelse return false;

    // Find H atom index in CCD component
    var h_ccd_idx: ?u16 = null;
    for (component.atoms, 0..) |a, i| {
        const a_name = trimName(&a.name);
        if (std.mem.eql(u8, a_name, h_name)) {
            h_ccd_idx = @intCast(i);
            break;
        }
    }
    const h_idx = h_ccd_idx orelse return false;

    // Find the heavy atom bonded to this H (center)
    var center_idx: ?u16 = null;
    for (component.bonds) |bond| {
        const other: u16 = if (bond.atom_idx_1 == h_idx)
            bond.atom_idx_2
        else if (bond.atom_idx_2 == h_idx)
            bond.atom_idx_1
        else
            continue;
        if (component.atoms[other].element_symbol[0] != 'H') {
            center_idx = other;
            break;
        }
    }
    const ci = center_idx orelse return false;
    center_out.* = component.atoms[ci].name;

    // Find the first heavy atom bonded to center (axis), excluding H atoms
    for (component.bonds) |bond| {
        const other: u16 = if (bond.atom_idx_1 == ci)
            bond.atom_idx_2
        else if (bond.atom_idx_2 == ci)
            bond.atom_idx_1
        else
            continue;
        if (other == h_idx) continue;
        if (component.atoms[other].element_symbol[0] != 'H') {
            axis_out.* = component.atoms[other].name;
            return true;
        }
    }
    return false;
}

/// Resolve center and axis atom names for N-terminal H atoms (H1, H2, H3).
/// N-terminal NH3+ and PRO NH2+ groups rotate around the N-CA bond:
///   center = N, axis = CA.
/// Returns true if h_name matches a known N-terminal H name.
const RotatorAtoms = struct { center: u32, axis: u32 };

/// Longest H–heavy bond accepted when the rotator center is found by distance.
const max_h_bond = 1.25;
/// Longest heavy–heavy bond accepted for the rotation axis.
const max_heavy_bond = 2.0;

/// Resolve a single-H rotator's center and axis from coordinates: the center
/// is the heavy atom of the residue nearest the H (within a bond length) and
/// the axis is the heavy atom nearest the center. Used when no plan or
/// dictionary names them, e.g. the 3'-terminal HO3' or H placed by
/// distance-based bond inference.
fn resolveGeometricRotatorAtoms(atoms: []const Atom, res: Residue, h: *const Atom) ?RotatorAtoms {
    const center = nearestHeavy(atoms, res, h.pos, h.altloc, null, max_h_bond) orelse return null;
    const axis = nearestHeavy(atoms, res, atoms[center].pos, h.altloc, center, max_heavy_bond) orelse return null;
    return .{ .center = center, .axis = axis };
}

fn nearestHeavy(atoms: []const Atom, res: Residue, pos: Vec3, altloc: u8, exclude: ?u32, max_dist: f32) ?u32 {
    var best: ?u32 = null;
    var best_d2: f32 = max_dist * max_dist;
    for (atoms[res.atom_start..res.atom_end], res.atom_start..) |*a, i| {
        if (a.is_hydrogen) continue;
        if (exclude) |e| if (e == i) continue;
        if (a.altloc != ' ' and altloc != ' ' and a.altloc != altloc) continue;
        const v = a.pos.sub(pos);
        const d2 = v.dot(v);
        if (d2 <= best_d2) {
            best_d2 = d2;
            best = @intCast(i);
        }
    }
    return best;
}

fn resolveNtermRotatorAtoms(h_name: []const u8, center_out: *[4]u8, axis_out: *[4]u8) bool {
    const is_nterm_h = std.mem.eql(u8, h_name, "H1") or
        std.mem.eql(u8, h_name, "H2") or
        std.mem.eql(u8, h_name, "H3");
    if (!is_nterm_h) return false;
    center_out.* = .{ ' ', 'N', ' ', ' ' };
    axis_out.* = .{ ' ', 'C', 'A', ' ' };
    return true;
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

pub const MoverGenResult = struct {
    movers: []Mover,
    n_skipped: u32,
};

/// Scan all placed H atoms in the model and create Mover instances for
/// optimization. Caller owns the returned movers slice and must call
/// `deinit()` on each Mover, then free the slice with the same allocator.
/// inline_dict is checked first per component, then ccd_dict as fallback.
pub fn generateMovers(
    allocator: std.mem.Allocator,
    mdl: *const Model,
    no_flip: bool,
    ccd_dict: ?*const ComponentDict,
    inline_dict: ?*const ComponentDict,
    overrides: ?*const protonation.ProtonationOverrides,
    mode: bond_policy.BondLengthMode,
) !MoverGenResult {
    var movers: std.ArrayListUnmanaged(Mover) = .empty;
    errdefer {
        for (movers.items) |*m| m.deinit();
        movers.deinit(allocator);
    }

    var n_skipped: u32 = 0;

    // Track which group rotators we have already created to avoid duplicates.
    // Key: (residue_idx, center_name) packed into u64.
    var seen_groups: std.AutoHashMapUnmanaged(GroupKey, void) = .empty;
    defer seen_groups.deinit(allocator);

    const atoms = mdl.atoms.items;

    for (atoms, 0..) |*atom, idx| {
        if (!atom.is_added) continue;
        if (atom.mover_hint == .none) continue;

        const residue_idx = atom.residue_idx;
        const res = mdl.residues.items[residue_idx];
        const comp_id = res.compIdSlice();
        const h_name = atom.nameSlice();
        const target_altloc = atom.altloc;
        const residue_state = if (overrides) |ov| ov.find(mdl, residue_idx) else null;

        // Per-component dict fallback: inline_dict first, then ccd_dict.
        const effective_dict: ?*const ComponentDict = blk: {
            if (inline_dict) |d| {
                if (d.get(comp_id) != null) break :blk inline_dict;
            }
            break :blk ccd_dict;
        };

        switch (atom.mover_hint) {
            .rotate => {
                // Single H rotator (OH, SH, etc.). The standard plan or the
                // dictionary names the center and axis; otherwise they come
                // from coordinates (3'-terminal HO3', distance-inferred H).
                // H names alone are not used: ligand H named H1/H2/H3 are
                // not N-terminal amine H.
                var center_name_buf: [4]u8 = undefined;
                var axis_name_buf: [4]u8 = undefined;
                const named = if (findPlanForH(comp_id, h_name)) |plan| blk: {
                    center_name_buf = plan.connected[0];
                    axis_name_buf = plan.connected[1];
                    break :blk true;
                } else resolveCcdRotatorAtoms(effective_dict, comp_id, h_name, &center_name_buf, &axis_name_buf);

                var center_idx: u32 = undefined;
                var axis_idx: u32 = undefined;
                if (named) {
                    const center_name = trimName(&center_name_buf);
                    const axis_name = trimName(&axis_name_buf);
                    center_idx = findAtomIdx(atoms, res, residue_idx, center_name, target_altloc, mdl.original_atom_count) orelse {
                        log.warn("center atom '{s}' not found for H '{s}' in {s} (res {d}), skipping", .{ center_name, h_name, comp_id, residue_idx });
                        n_skipped += 1;
                        continue;
                    };
                    axis_idx = findAtomIdx(atoms, res, residue_idx, axis_name, target_altloc, mdl.original_atom_count) orelse {
                        log.warn("axis atom '{s}' not found for H '{s}' in {s} (res {d}), skipping", .{ axis_name, h_name, comp_id, residue_idx });
                        n_skipped += 1;
                        continue;
                    };
                } else if (resolveGeometricRotatorAtoms(atoms, res, atom)) |ra| {
                    center_idx = ra.center;
                    axis_idx = ra.axis;
                } else {
                    log.warn("no bonded heavy atoms found for H '{s}' in {s} (res {d}), skipping rotator", .{ h_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                }
                const h_idx: u32 = @intCast(idx);

                const mover = try rotator.createSingleHRotator(
                    allocator,
                    atoms,
                    h_idx,
                    center_idx,
                    axis_idx,
                    residue_idx,
                );
                try movers.append(allocator, mover);
            },
            .rotate_nh2 => {
                // Group rotator for NH2 (non-PRO neutral N-term or PRO NH2+).
                // 2 H atoms share the same N center and axis CA. Unlike
                // .rotate_nh3, there are only 2 H to bundle and we use 12
                // orientations at 30° (no 3-fold symmetry).
                var center_name_raw: [4]u8 = undefined;
                var axis_name_buf: [4]u8 = undefined;
                var center_name: []const u8 = undefined;
                var axis_name: []const u8 = undefined;

                if (resolveNtermRotatorAtoms(h_name, &center_name_raw, &axis_name_buf)) {
                    center_name = trimName(&center_name_raw);
                    axis_name = trimName(&axis_name_buf);
                } else {
                    log.warn("no N-terminal resolver for NH2 H '{s}' in {s} (res {d}), skipping rotator", .{ h_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                }

                // Dedup key — one mover per (residue, N center, altloc)
                const group_key = GroupKey{ .residue_idx = residue_idx, .center_name = normalizeNameKey(center_name_raw), .altloc = target_altloc };
                const gop = try seen_groups.getOrPut(allocator, group_key);
                if (gop.found_existing) continue;

                const center_idx = findAtomIdx(atoms, res, residue_idx, center_name, target_altloc, mdl.original_atom_count) orelse {
                    log.warn("center atom '{s}' not found for NH2 group in {s} (res {d}), skipping", .{ center_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                };
                const axis_idx = findAtomIdx(atoms, res, residue_idx, axis_name, target_altloc, mdl.original_atom_count) orelse {
                    log.warn("axis atom '{s}' not found for NH2 group in {s} (res {d}), skipping", .{ axis_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                };

                // Collect exactly the 2 NH2 H atoms for this residue+altloc.
                var h_indices: [2]u32 = undefined;
                var h_count: u32 = 0;
                for (atoms[mdl.original_atom_count..], mdl.original_atom_count..) |*a, ai| {
                    if (a.residue_idx != residue_idx) continue;
                    if (!a.is_added) continue;
                    if (a.mover_hint != .rotate_nh2) continue;
                    if (a.altloc != target_altloc) continue;
                    if (h_count < 2) {
                        h_indices[h_count] = @intCast(ai);
                    }
                    h_count += 1;
                }
                if (h_count != 2) {
                    log.warn("expected 2 H for NH2 group in {s} (res {d}), found {d}, skipping", .{ comp_id, residue_idx, h_count });
                    n_skipped += 1;
                    continue;
                }

                const mover = try rotator.createNH2Rotator(
                    allocator,
                    atoms,
                    h_indices,
                    center_idx,
                    axis_idx,
                    residue_idx,
                );
                try movers.append(allocator, mover);
            },
            .rotate_methyl => continue,
            .rotate_nh3 => {
                if (residue_state) |state| {
                    if (state == .lys and state.lys == .neutral) continue;
                }
                // Group rotator: 3 H atoms share same center.
                // Try standard plan first, fall back to CCD topology.
                var center_name_raw: [4]u8 = undefined;
                var axis_name_buf: [4]u8 = undefined;
                var center_name: []const u8 = undefined;
                var axis_name: []const u8 = undefined;

                if (findPlanForH(comp_id, h_name)) |plan| {
                    center_name_raw = plan.connected[0];
                    axis_name_buf = plan.connected[1];
                    center_name = trimName(&center_name_raw);
                    axis_name = trimName(&axis_name_buf);
                } else if (resolveCcdRotatorAtoms(effective_dict, comp_id, h_name, &center_name_raw, &axis_name_buf)) {
                    center_name = trimName(&center_name_raw);
                    axis_name = trimName(&axis_name_buf);
                } else if (resolveNtermRotatorAtoms(h_name, &center_name_raw, &axis_name_buf)) {
                    center_name = trimName(&center_name_raw);
                    axis_name = trimName(&axis_name_buf);
                } else {
                    log.warn("no plan for group H '{s}' in {s} (res {d}), skipping rotator", .{ h_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                }

                // Dedup key (normalize padding so CCD 'N   ' and PDB ' N  ' match)
                const group_key = GroupKey{ .residue_idx = residue_idx, .center_name = normalizeNameKey(center_name_raw), .altloc = target_altloc };
                const gop = try seen_groups.getOrPut(allocator, group_key);
                if (gop.found_existing) continue;

                const center_idx = findAtomIdx(atoms, res, residue_idx, center_name, target_altloc, mdl.original_atom_count) orelse {
                    log.warn("center atom '{s}' not found for group in {s} (res {d}), skipping", .{ center_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                };
                const axis_idx = findAtomIdx(atoms, res, residue_idx, axis_name, target_altloc, mdl.original_atom_count) orelse {
                    log.warn("axis atom '{s}' not found for group in {s} (res {d}), skipping", .{ axis_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                };

                // Find all 3 H atoms with same (residue_idx, center, hint).
                // Use both standard plan lookup and CCD fallback to match center.
                // Scan only added atoms (atoms[original_atom_count..]) for efficiency.
                var h_indices: [3]u32 = undefined;
                var h_count: u32 = 0;
                for (atoms[mdl.original_atom_count..], mdl.original_atom_count..) |*a, ai| {
                    if (a.residue_idx == residue_idx and a.is_added and a.mover_hint == atom.mover_hint and a.altloc == target_altloc) {
                        // Resolve this H's center via standard plan, CCD, or N-terminal fallback
                        var a_center: [4]u8 = undefined;
                        var a_axis: [4]u8 = undefined;
                        const has_center = if (findPlanForH(comp_id, a.nameSlice())) |ap| blk: {
                            a_center = ap.connected[0];
                            break :blk true;
                        } else if (resolveCcdRotatorAtoms(effective_dict, comp_id, a.nameSlice(), &a_center, &a_axis))
                            true
                        else
                            resolveNtermRotatorAtoms(a.nameSlice(), &a_center, &a_axis);

                        if (has_center and std.mem.eql(u8, trimName(&a_center), center_name)) {
                            if (h_count < 3) {
                                h_indices[h_count] = @intCast(ai);
                                h_count += 1;
                            }
                        }
                    }
                }
                if (h_count != 3) {
                    log.warn("expected 3 H for group at '{s}' in {s} (res {d}), found {d}, skipping", .{ center_name, comp_id, residue_idx, h_count });
                    n_skipped += 1;
                    continue;
                }

                const mover = try rotator.createNH3Rotator(
                    allocator,
                    atoms,
                    h_indices,
                    center_idx,
                    axis_idx,
                    residue_idx,
                );
                try movers.append(allocator, mover);
            },
            .flip_amide => {
                if (no_flip) continue;

                // Deduplicate: one amide flipper per residue
                const group_key = GroupKey{ .residue_idx = residue_idx, .center_name = .{ 'A', 'M', 'D', ' ' }, .altloc = target_altloc };
                const gop = try seen_groups.getOrPut(allocator, group_key);
                if (gop.found_existing) continue;

                const plan = findPlanForH(comp_id, h_name) orelse {
                    log.warn("no plan for flip H '{s}' in {s} (res {d}), skipping", .{ h_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                };

                // N atom = connected[0] (ND2 for ASN, NE2 for GLN)
                const n_name = trimName(&plan.connected[0]);
                // C atom = connected[1] (CG for ASN, CD for GLN)
                const c_name = trimName(&plan.connected[1]);
                // O atom: OD1 for ASN, OE1 for GLN
                const o_name: []const u8 = if (std.mem.eql(u8, comp_id, "ASN")) "OD1" else "OE1";

                const n_at_idx = findAtomIdx(atoms, res, residue_idx, n_name, target_altloc, mdl.original_atom_count) orelse {
                    log.warn("N atom '{s}' not found for amide flip in {s} (res {d})", .{ n_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                };
                const o_at_idx = findAtomIdx(atoms, res, residue_idx, o_name, target_altloc, mdl.original_atom_count) orelse {
                    log.warn("O atom '{s}' not found for amide flip in {s} (res {d})", .{ o_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                };
                const c_at_idx = findAtomIdx(atoms, res, residue_idx, c_name, target_altloc, mdl.original_atom_count) orelse {
                    log.warn("C atom '{s}' not found for amide flip in {s} (res {d})", .{ c_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                };

                // Find the two H atoms for this amide group
                var h_indices: [2]u32 = undefined;
                var h_count: u32 = 0;
                const all_plans = standard.getPlans(comp_id) orelse continue;
                for (all_plans) |*p| {
                    if (p.mover_hint != .flip_amide) continue;
                    const p_h_name = trimName(&p.h_name);
                    if (findAtomIdx(atoms, res, residue_idx, p_h_name, target_altloc, mdl.original_atom_count)) |hi| {
                        if (h_count < 2) {
                            h_indices[h_count] = hi;
                            h_count += 1;
                        }
                    }
                }
                if (h_count != 2) {
                    log.warn("expected 2 H for amide flip in {s} (res {d}), found {d}", .{ comp_id, residue_idx, h_count });
                    n_skipped += 1;
                    continue;
                }

                const m = try flipper.createAmideFlipper(
                    allocator,
                    atoms,
                    o_at_idx,
                    n_at_idx,
                    h_indices[0],
                    h_indices[1],
                    c_at_idx,
                    residue_idx,
                    mode,
                );
                try movers.append(allocator, m);
            },
            .flip_his => {
                if (residue_state) |state| {
                    if (state == .his and state.his != .auto) continue;
                }
                if (no_flip) continue;

                // Deduplicate: one His flipper per residue
                const group_key = GroupKey{ .residue_idx = residue_idx, .center_name = .{ 'H', 'I', 'S', ' ' }, .altloc = target_altloc };
                const gop = try seen_groups.getOrPut(allocator, group_key);
                if (gop.found_existing) continue;

                // Find ring heavy atoms by name
                const nd1_idx = findAtomIdx(atoms, res, residue_idx, "ND1", target_altloc, mdl.original_atom_count) orelse {
                    log.warn("ND1 not found for His flip (res {d})", .{residue_idx});
                    n_skipped += 1;
                    continue;
                };
                const cd2_idx = findAtomIdx(atoms, res, residue_idx, "CD2", target_altloc, mdl.original_atom_count) orelse {
                    log.warn("CD2 not found for His flip (res {d})", .{residue_idx});
                    n_skipped += 1;
                    continue;
                };
                const ce1_idx = findAtomIdx(atoms, res, residue_idx, "CE1", target_altloc, mdl.original_atom_count) orelse {
                    log.warn("CE1 not found for His flip (res {d})", .{residue_idx});
                    n_skipped += 1;
                    continue;
                };
                const ne2_idx = findAtomIdx(atoms, res, residue_idx, "NE2", target_altloc, mdl.original_atom_count) orelse {
                    log.warn("NE2 not found for His flip (res {d})", .{residue_idx});
                    n_skipped += 1;
                    continue;
                };

                // H atoms are optional
                const hd1_idx = findAtomIdx(atoms, res, residue_idx, "HD1", target_altloc, mdl.original_atom_count);
                const he2_idx = findAtomIdx(atoms, res, residue_idx, "HE2", target_altloc, mdl.original_atom_count);

                const m = try flipper.createHisFlipper(
                    allocator,
                    atoms,
                    nd1_idx,
                    cd2_idx,
                    ce1_idx,
                    ne2_idx,
                    hd1_idx,
                    he2_idx,
                    residue_idx,
                    mode,
                );
                try movers.append(allocator, m);
            },
            .none => continue,
        }
    }

    return .{
        .movers = try movers.toOwnedSlice(allocator),
        .n_skipped = n_skipped,
    };
}

const GroupKey = struct {
    residue_idx: u32,
    center_name: [4]u8,
    altloc: u8,
};

/// Normalize a 4-byte atom name for use as a hash key.
/// Trims spaces and left-aligns so that ' N  ' and 'N   ' produce the same key.
fn normalizeNameKey(raw: [4]u8) [4]u8 {
    const trimmed = trimName(&raw);
    var out: [4]u8 = .{ ' ', ' ', ' ', ' ' };
    for (trimmed, 0..) |c, i| out[i] = c;
    return out;
}

// ── Tests ─────────────────────────────────────────────────────────────────────

const testing = std.testing;
const mmcif = @import("../mmcif.zig");
const placer = @import("../place/placer.zig");

test "trimName trims spaces" {
    try testing.expectEqualStrings("CB", trimName(&.{ ' ', 'C', 'B', ' ' }));
    try testing.expectEqualStrings("CA", trimName(&.{ ' ', 'C', 'A', ' ' }));
    try testing.expectEqualStrings("N", trimName(&.{ ' ', 'N', ' ', ' ' }));
    try testing.expectEqualStrings("HB1", trimName(&.{ 'H', 'B', '1', ' ' }));
}

test "findPlanForH returns correct plan for ALA HB1" {
    const plan = findPlanForH("ALA", "HB1");
    try testing.expect(plan != null);
    try testing.expectEqualStrings("HB1", trimName(&plan.?.h_name));
    try testing.expectEqual(standard.MoverHint.rotate_methyl, plan.?.mover_hint);
}

test "findPlanForH returns null for unknown atom" {
    const plan = findPlanForH("ALA", "HX9");
    try testing.expect(plan == null);
}

test "findPlanForH finds nucleotide plans" {
    // Aromatic C-H on adenine
    const h8 = findPlanForH("DA", "H8");
    try testing.expect(h8 != null);
    try testing.expectEqualStrings("H8", trimName(&h8.?.h_name));

    // OH rotator on RNA ribose
    const ho2 = findPlanForH("A", "HO2'");
    try testing.expect(ho2 != null);
    try testing.expectEqual(standard.MoverHint.rotate, ho2.?.mover_hint);

    // Thymine methyl
    const h71 = findPlanForH("DT", "H71");
    try testing.expect(h71 != null);
    try testing.expectEqual(standard.MoverHint.rotate_methyl, h71.?.mover_hint);

    // Standard AA still works
    try testing.expect(findPlanForH("ALA", "HB1") != null);

    // Unknown returns null
    try testing.expect(findPlanForH("DA", "HX9") == null);
}

test "resolveGeometricRotatorAtoms finds the bonded O and its heavy neighbor, not N" {
    // A ligand hydroxyl H named H1, as in FAD: the H name must not send the
    // lookup to an N-terminal N. Atom 0 (P) is the O's neighbor, atom 2 (N)
    // is farther from the O than P, atom 3 is the H.
    var atoms = [_]Atom{
        .{ .pos = .{ .x = 0, .y = 0, .z = 0 }, .is_hydrogen = false },
        .{ .pos = .{ .x = 1.6, .y = 0, .z = 0 }, .is_hydrogen = false },
        .{ .pos = .{ .x = 4.0, .y = 0, .z = 0 }, .is_hydrogen = false },
        .{ .pos = .{ .x = 1.9, .y = 0.92, .z = 0 }, .is_hydrogen = true },
    };
    atoms[0].setName("P");
    atoms[1].setName("O1A");
    atoms[2].setName("N");
    atoms[3].setName("H1");
    const res = Residue{ .atom_start = 0, .atom_end = 3 };
    const ra = resolveGeometricRotatorAtoms(&atoms, res, &atoms[3]).?;
    try testing.expectEqual(@as(u32, 1), ra.center);
    try testing.expectEqual(@as(u32, 0), ra.axis);

    // An H with no heavy atom within bonding distance resolves to nothing.
    var far = atoms[3];
    far.pos = .{ .x = 10, .y = 10, .z = 10 };
    try testing.expect(resolveGeometricRotatorAtoms(&atoms, res, &far) == null);
}

test "findAtomIdx prefers exact altloc and falls back to blank" {
    const source = @embedFile("../test_data/ala_altloc.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    const res0 = mdl.residues.items[0];
    const orig = mdl.original_atom_count;

    try testing.expectEqual(findAtomIdx(mdl.atoms.items, res0, 0, "CA", 'A', orig).?, 1);
    try testing.expectEqual(findAtomIdx(mdl.atoms.items, res0, 0, "CA", 'B', orig).?, 2);
    try testing.expectEqual(findAtomIdx(mdl.atoms.items, res0, 0, "N", 'A', orig).?, 0);
    try testing.expectEqual(findAtomIdx(mdl.atoms.items, res0, 0, "N", 'B', orig).?, 0);
}

test "generateMovers creates no mover for ALA methyl H" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    placer.applyChemistry(&mdl);
    _ = try placer.addHydrogens(&mdl, null, null);

    const gen_result = try generateMovers(testing.allocator, &mdl, false, null, null, null, .neutron);
    const movers = gen_result.movers;
    defer {
        for (0..movers.len) |i| @constCast(&movers[i]).deinit();
        testing.allocator.free(movers);
    }

    // The CB methyl stays staggered as placed: no mover controls its H.
    try expectNoMethylMover(&mdl, movers);
}

fn expectNoMethylMover(mdl: *const Model, movers: []const Mover) !void {
    var n_methyl_h: u32 = 0;
    for (mdl.atoms.items) |a| {
        if (a.is_added and a.mover_hint == .rotate_methyl) n_methyl_h += 1;
    }
    try testing.expect(n_methyl_h > 0);
    for (movers) |m| {
        for (m.atom_indices) |ai| {
            try testing.expect(mdl.atoms.items[ai].mover_hint != .rotate_methyl);
        }
    }
}

test "generateMovers creates no methyl mover for any altloc conformer" {
    const source = @embedFile("../test_data/ala_altloc.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    placer.applyChemistry(&mdl);
    _ = try placer.addHydrogens(&mdl, null, null);

    const gen_result = try generateMovers(testing.allocator, &mdl, false, null, null, null, .neutron);
    const movers = gen_result.movers;
    defer {
        for (0..movers.len) |i| @constCast(&movers[i]).deinit();
        testing.allocator.free(movers);
    }

    try expectNoMethylMover(&mdl, movers);
}

test "generateMovers total count for ALA" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    placer.applyChemistry(&mdl);
    _ = try placer.addHydrogens(&mdl, null, null);

    const gen_result = try generateMovers(testing.allocator, &mdl, false, null, null, null, .neutron);
    const movers = gen_result.movers;
    defer {
        for (0..movers.len) |i| @constCast(&movers[i]).deinit();
        testing.allocator.free(movers);
    }

    // N-terminal ALA has only the NH3+ rotator (H1/H2/H3): its CB methyl gets
    // no mover and it has no OH/SH groups.
    try testing.expectEqual(@as(usize, 1), movers.len);
    try testing.expectEqual(mover_mod.MoverKind.nh3_rotator, movers[0].kind);
}

test "optimizer pipeline runs without error on ALA" {
    const opt_mod = @import("optimizer.zig");

    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    placer.applyChemistry(&mdl);
    _ = try placer.addHydrogens(&mdl, null, null);

    const gen_result = try generateMovers(testing.allocator, &mdl, false, null, null, null, .neutron);
    const movers = gen_result.movers;
    defer {
        for (0..movers.len) |i| @constCast(&movers[i]).deinit();
        testing.allocator.free(movers);
    }

    try testing.expect(movers.len > 0);

    // Run optimizer - should complete without error
    const opt_config = opt_mod.OptConfig{};
    _ = try opt_mod.optimize(testing.allocator, movers, &mdl, opt_config);
}

test "NH3 rotator controls 3 atoms" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    placer.applyChemistry(&mdl);
    _ = try placer.addHydrogens(&mdl, null, null);

    const gen_result = try generateMovers(testing.allocator, &mdl, false, null, null, null, .neutron);
    const movers = gen_result.movers;
    defer {
        for (0..movers.len) |i| @constCast(&movers[i]).deinit();
        testing.allocator.free(movers);
    }

    for (movers) |m| {
        if (m.kind == .nh3_rotator) {
            try testing.expectEqual(@as(usize, 3), m.atom_indices.len);
            try testing.expectEqual(@as(usize, 4), m.orientations.len);
        }
    }
}

test "generateMovers creates amide flipper for ASN" {
    const source = @embedFile("../test_data/asn.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    placer.applyChemistry(&mdl);
    _ = try placer.addHydrogens(&mdl, null, null);

    const gen_result = try generateMovers(testing.allocator, &mdl, false, null, null, null, .neutron);
    const movers = gen_result.movers;
    defer {
        for (0..movers.len) |i| @constCast(&movers[i]).deinit();
        testing.allocator.free(movers);
    }

    var amide_count: u32 = 0;
    for (movers) |m| {
        if (m.kind == .amide_flip) amide_count += 1;
    }
    try testing.expectEqual(@as(u32, 1), amide_count);
}

test "no_flip suppresses flip movers" {
    const source = @embedFile("../test_data/asn.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    placer.applyChemistry(&mdl);
    _ = try placer.addHydrogens(&mdl, null, null);

    const gen_result = try generateMovers(testing.allocator, &mdl, true, null, null, null, .neutron);
    const movers = gen_result.movers;
    defer {
        for (0..movers.len) |i| @constCast(&movers[i]).deinit();
        testing.allocator.free(movers);
    }

    for (movers) |m| {
        try testing.expect(m.kind != .amide_flip);
        try testing.expect(m.kind != .his_flip);
    }
}

test "generateMovers creates His flipper" {
    const source = @embedFile("../test_data/his.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    placer.applyChemistry(&mdl);
    _ = try placer.addHydrogens(&mdl, null, null);

    const gen_result = try generateMovers(testing.allocator, &mdl, false, null, null, null, .neutron);
    const movers = gen_result.movers;
    defer {
        for (0..movers.len) |i| @constCast(&movers[i]).deinit();
        testing.allocator.free(movers);
    }

    var his_count: u32 = 0;
    for (movers) |m| {
        if (m.kind == .his_flip) {
            his_count += 1;
            try testing.expectEqual(@as(usize, 4), m.orientations.len);
        }
    }
    try testing.expectEqual(@as(u32, 1), his_count);
}

test "generateMovers skips fixed His override" {
    const source = @embedFile("../test_data/his.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();
    mdl.residues.items[0].auth_seq_id = 1;

    var overrides = try protonation.parseString(testing.allocator,
        \\A:1 HIS HID
    );
    defer overrides.deinit();

    placer.applyChemistryWithConfig(&mdl, .{ .protonation = &overrides });
    _ = try placer.addHydrogensWithConfig(&mdl, null, null, .{ .protonation = &overrides });

    const gen_result = try generateMovers(testing.allocator, &mdl, false, null, null, &overrides, .neutron);
    const movers = gen_result.movers;
    defer {
        for (0..movers.len) |i| @constCast(&movers[i]).deinit();
        testing.allocator.free(movers);
    }

    for (movers) |m| {
        try testing.expect(m.kind != .his_flip);
    }
}

// ── NH2 group rotator (issue #253) ───────────────────────────────────────────

test "generateMovers creates a single nh2_rotator for N-terminal PRO" {
    // Inline minimal PRO N-terminal fixture. N-terminal PRO should yield
    // one nh2_rotator covering H2 and H3 together, not two independent
    // single_h_rotators.
    const source =
        \\data_PRO_NH2
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\_atom_site.auth_asym_id
        \\ATOM 1 N N   PRO A 1 1.458 0.000 0.000 1.00 10.0 . A
        \\ATOM 2 C CA  PRO A 1 0.000 0.000 0.000 1.00 10.0 . A
        \\ATOM 3 C C   PRO A 1 -0.523 1.413 0.000 1.00 10.0 . A
        \\ATOM 4 O O   PRO A 1 -1.742 1.560 0.000 1.00 10.0 . A
        \\ATOM 5 C CB  PRO A 1 -0.700 -1.100 1.800 1.00 10.0 . A
        \\ATOM 6 C CG  PRO A 1 0.200 -1.900 2.500 1.00 10.0 . A
        \\ATOM 7 C CD  PRO A 1 1.500 -1.200 2.200 1.00 10.0 . A
        \\#
    ;
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    placer.applyChemistry(&mdl);
    _ = try placer.addHydrogens(&mdl, null, null);

    const gen_result = try generateMovers(testing.allocator, &mdl, false, null, null, null, .neutron);
    const movers = gen_result.movers;
    defer {
        for (0..movers.len) |i| @constCast(&movers[i]).deinit();
        testing.allocator.free(movers);
    }

    var nh2_count: u32 = 0;
    var single_h_count: u32 = 0;
    for (movers) |m| {
        if (m.kind == .nh2_rotator) nh2_count += 1;
        if (m.kind == .single_h_rotator) single_h_count += 1;
    }
    // Exactly one nh2_rotator covering both H, no single-H rotators for
    // the backbone amine (would indicate the old per-H pattern).
    try testing.expectEqual(@as(u32, 1), nh2_count);
    try testing.expectEqual(@as(u32, 0), single_h_count);

    // The nh2_rotator controls exactly 2 H atoms at 12 orientations.
    for (movers) |m| {
        if (m.kind != .nh2_rotator) continue;
        try testing.expectEqual(@as(usize, 2), m.atom_indices.len);
        try testing.expectEqual(@as(usize, 12), m.orientations.len);
        try testing.expect(m.center_idx != null);
        try testing.expect(m.axis_idx != null);
    }
}

test "generateMovers creates nh2_rotator for non-PRO neutral N-term" {
    const source = @embedFile("../test_data/ala_stretched.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    placer.applyChemistryWithConfig(&mdl, .{ .nterm_mode = .neutral });
    _ = try placer.addHydrogensWithConfig(&mdl, null, null, .{ .nterm_mode = .neutral });

    const gen_result = try generateMovers(testing.allocator, &mdl, false, null, null, null, .neutron);
    const movers = gen_result.movers;
    defer {
        for (0..movers.len) |i| @constCast(&movers[i]).deinit();
        testing.allocator.free(movers);
    }

    var nh2_count: u32 = 0;
    for (movers) |m| {
        if (m.kind == .nh2_rotator) nh2_count += 1;
    }
    try testing.expectEqual(@as(u32, 1), nh2_count);
}

test "nh2_rotator preserves HNH angle across all orientations" {
    // Critical regression test for the pre-existing bug fixed by this change:
    // two independent single-H rotators would let H2 and H3 rotate freely
    // around the same N-CA axis and destroy the fixed HNH angle. A group
    // rotator must preserve the relative positions of both H atoms.
    const source = @embedFile("../test_data/ala_stretched.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    placer.applyChemistryWithConfig(&mdl, .{ .nterm_mode = .neutral });
    _ = try placer.addHydrogensWithConfig(&mdl, null, null, .{ .nterm_mode = .neutral });

    const gen_result = try generateMovers(testing.allocator, &mdl, false, null, null, null, .neutron);
    const movers = gen_result.movers;
    defer {
        for (0..movers.len) |i| @constCast(&movers[i]).deinit();
        testing.allocator.free(movers);
    }

    // Grab the nh2_rotator.
    var target: ?*const mover_mod.Mover = null;
    for (movers) |*m| {
        if (m.kind == .nh2_rotator) {
            target = m;
            break;
        }
    }
    try testing.expect(target != null);
    const m = target.?;

    // Reference HNH distance from the initial orientation.
    const ref_h2 = m.orientations[0].positions[0];
    const ref_h3 = m.orientations[0].positions[1];
    const ref_d2 = ref_h2.distance(ref_h3);

    // Every orientation must preserve the same H2-H3 distance (rigid body).
    for (m.orientations) |orient| {
        const d2 = orient.positions[0].distance(orient.positions[1]);
        try testing.expectApproxEqAbs(ref_d2, d2, 0.005);
    }
}
