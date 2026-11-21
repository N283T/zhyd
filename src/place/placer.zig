//! Unified hydrogen placement entry point.
//!
//! Adds hydrogens to a Model using standard plans for known residues (20 AA)
//! and CCD-derived plans for non-standard residues.  When no dictionary entry
//! exists, falls back to distance-based bond inference.

const std = @import("std");
const model_mod = @import("../model.zig");
const Model = model_mod.Model;
const Atom = model_mod.Atom;
const Residue = model_mod.Residue;
const ccd_mod = @import("../ccd.zig");
const ComponentDict = ccd_mod.ComponentDict;
const standard = @import("standard.zig");
const nucleotide = @import("nucleotide.zig");
const modified = @import("modified.zig");
const ccd_derive = @import("ccd_derive.zig");
const distance_derive = @import("distance_derive.zig");
const geometry = @import("geometry.zig");
const bond_policy = @import("bond_policy.zig");
const math_mod = @import("../math.zig");
const element = @import("../element.zig");

const Vec3f32 = math_mod.Vec3(f32);
const topology = @import("topology.zig");
const chemistry = @import("chemistry.zig");
const protonation = @import("protonation.zig");
const neighbor_mod = @import("../model/neighbor.zig");
const CellList = neighbor_mod.CellList;

pub const lookup = @import("lookup.zig");
pub const terminal = @import("terminal.zig");
pub const execute = @import("execute.zig");
pub const water_mod = @import("water.zig");

const findAtom = lookup.findAtom;
const findAtomPos = lookup.findAtomPos;
const existsInResidue = lookup.existsInResidue;
const padName = lookup.padName;

/// Result of a single executePlan call — why a hydrogen was or wasn't placed.
/// Re-exported from execute.zig; defined there to avoid circular imports.
pub const PlaceResult = execute.PlaceResult;

pub const PlacementResult = struct {
    n_placed: u32 = 0,
    n_skipped_existing: u32 = 0,
    n_skipped_inter_residue: u32 = 0,
    n_skipped_missing_ref: u32 = 0,
    n_skipped_quality_filter: u32 = 0,
    n_residues: u32 = 0,
    n_distance_derived: u32 = 0,

    pub fn totalSkipped(self: PlacementResult) u32 {
        return self.n_skipped_existing + self.n_skipped_inter_residue + self.n_skipped_missing_ref + self.n_skipped_quality_filter;
    }

    fn tally(self: *PlacementResult, r: PlaceResult) void {
        switch (r) {
            .placed => self.n_placed += 1,
            .existing_h => self.n_skipped_existing += 1,
            .inter_residue => self.n_skipped_inter_residue += 1,
            .missing_parent, .missing_ref => self.n_skipped_missing_ref += 1,
        }
    }
};

pub const WaterConfig = water_mod.WaterConfig;

pub const PlacementConfig = struct {
    water: WaterConfig = .{},
    protonation: ?*const protonation.ProtonationOverrides = null,
    bond_policy: bond_policy.BondPolicy = .{},
    nterm_mode: terminal.NtermMode = .auto,
    sdf_dict: ?*const ComponentDict = null,
};

const AltlocSet = struct { locs: [10]u8, count: usize };

/// Collect distinct non-blank altloc values from a residue's atoms.
fn collectAltlocs(mdl: *const Model, res: Residue) AltlocSet {
    var result: AltlocSet = .{ .locs = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, .count = 0 };
    const atoms = mdl.atoms.items[res.atom_start..res.atom_end];
    for (atoms) |a| {
        if (a.altloc == ' ') continue;
        var found = false;
        for (result.locs[0..result.count]) |existing| {
            if (existing == a.altloc) {
                found = true;
                break;
            }
        }
        if (!found) {
            if (result.count >= result.locs.len) {
                std.log.warn("residue has more than {d} altlocs, excess ignored", .{result.locs.len});
                break;
            }
            result.locs[result.count] = a.altloc;
            result.count += 1;
        }
    }
    return result;
}

/// Add hydrogens to the model.
/// Uses hardcoded plans (standard AA, nucleotides, modified), CCD-derived plans as fallback.
/// For CCD-derived placement, inline_dict is checked first, then ccd_dict (per-component fallback).
/// New atoms are appended to the end of model.atoms. Each new atom carries its residue_idx.
pub fn addHydrogens(
    mdl: *Model,
    ccd_dict: ?*const ComponentDict,
    inline_dict: ?*const ComponentDict,
) !PlacementResult {
    return addHydrogensWithConfig(mdl, ccd_dict, inline_dict, .{});
}

/// Add hydrogens to the model with placement options.
pub fn addHydrogensWithConfig(
    mdl: *Model,
    ccd_dict: ?*const ComponentDict,
    inline_dict: ?*const ComponentDict,
    config: PlacementConfig,
) !PlacementResult {
    var result = PlacementResult{};

    // Build a CellList over all non-hydrogen atoms so that water neighbor
    // queries (max radius 4.0 Å) can run in O(1) instead of O(N).
    // We collect positions and a mapping from CellList index → atom index.
    var water_cell_ctx: ?water_mod.WaterCellCtx = null;
    defer if (water_cell_ctx) |*ctx| {
        ctx.cell_list.deinit();
        mdl.allocator.free(ctx.positions);
        mdl.allocator.free(ctx.atom_indices);
    };

    if (config.water.enabled) blk: {
        var positions = std.ArrayListUnmanaged(Vec3f32).empty;
        var atom_idx_map = std.ArrayListUnmanaged(u32).empty;
        errdefer positions.deinit(mdl.allocator);
        errdefer atom_idx_map.deinit(mdl.allocator);

        for (mdl.atoms.items, 0..) |atom, i| {
            if (atom.is_hydrogen) continue;
            try positions.append(mdl.allocator, atom.pos);
            try atom_idx_map.append(mdl.allocator, @intCast(i));
        }

        const pos_slice = try positions.toOwnedSlice(mdl.allocator);
        errdefer mdl.allocator.free(pos_slice);
        const idx_slice = try atom_idx_map.toOwnedSlice(mdl.allocator);
        errdefer mdl.allocator.free(idx_slice);

        const cl = CellList.init(mdl.allocator, pos_slice, 5.0) catch |err| switch (err) {
            error.GridTooLarge => {
                // Fall back to linear scan: free arrays and leave water_cell_ctx null.
                mdl.allocator.free(pos_slice);
                mdl.allocator.free(idx_slice);
                break :blk;
            },
            else => return err,
        };
        water_cell_ctx = water_mod.WaterCellCtx{
            .cell_list = cl,
            .positions = pos_slice,
            .atom_indices = idx_slice,
        };
    }

    const n_residues = mdl.residues.items.len;
    for (0..n_residues) |res_idx| {
        const res = mdl.residues.items[res_idx];
        const comp_id = res.compIdSlice();
        const protonation_state = if (config.protonation) |overrides| overrides.find(mdl, res_idx) else null;

        if (config.water.enabled and water_mod.isWaterResidue(res, comp_id)) {
            const altlocs = collectAltlocs(mdl, res);
            const wr = try water_mod.placeWaterHydrogens(
                mdl,
                res,
                @intCast(res_idx),
                config.water,
                config.bond_policy.mode,
                water_cell_ctx,
                altlocs.count,
                &altlocs.locs,
            );
            result.n_placed += wr.n_placed;
            result.n_skipped_existing += wr.n_skipped_existing;
            result.n_skipped_inter_residue += wr.n_skipped_inter_residue;
            result.n_skipped_missing_ref += wr.n_skipped_missing_ref;
            result.n_skipped_quality_filter += wr.n_skipped_quality_filter;
            result.n_residues += 1;
            continue;
        }

        // Detect terminal residues.
        // The first residue in a chain is always a real N-terminus (gets NH3+/NH2+),
        // even if is_chain_break_before is set due to unobserved leading residues.
        // Mid-chain residues after a gap are terminal for bond topology (→
        // `is_nterm_for_bonding` drives CCD leaving-atom logic and the peptide-bond
        // inter-residue list below), but they keep a single backbone amide H by
        // default. Only `--nterm aggressive` promotes them to full N-termini
        // for amine placement and charge annotation (→ `treat_as_nterm`).
        const chain = mdl.chains.items[res.chain_idx];
        const is_real_nterm = (res_idx == chain.residue_start);
        const is_break_nterm = res.is_chain_break_before and (res_idx != chain.residue_start);
        const is_nterm_for_bonding = is_real_nterm or is_break_nterm;
        const treat_as_nterm = is_real_nterm or
            (is_break_nterm and config.nterm_mode.treatsBreakAsNterm());
        const is_cterm = (res_idx == chain.residue_end - 1) or
            (res_idx + 1 < n_residues and mdl.residues.items[res_idx + 1].is_chain_break_before);

        if (standard.getPlans(comp_id) orelse
            nucleotide.getPlans(comp_id) orelse
            modified.getPlans(comp_id)) |plans|
        {
            const bonds = topology.getBonds(comp_id); // null for non-standard
            const altlocs = collectAltlocs(mdl, res);

            const targets: []const u8 = if (altlocs.count == 0)
                &[_]u8{' '}
            else
                altlocs.locs[0..altlocs.count];

            // Detect whether plans include a backbone amide H (peptide residues).
            // Nucleotides and residues like PCA (cyclized N) have no backbone H.
            const has_backbone = for (plans) |plan| {
                if (terminal.isBackboneH(&plan)) break true;
            } else false;

            for (targets) |alt| {
                for (plans) |plan| {
                    if (shouldSkipPlanForProtonation(comp_id, &plan, protonation_state)) continue;
                    // Skip the single backbone amide H only when this residue
                    // will receive an N-terminal amine (NH3+/NH2+/NH2). In the
                    // default `auto` mode this matches `is_real_nterm`; in
                    // `aggressive` mode chain-break residues are also covered.
                    if (treat_as_nterm and terminal.isBackboneH(&plan)) continue;

                    result.tally(try execute.executePlan(mdl, res, @intCast(res_idx), &plan, bonds, alt, config.bond_policy.mode));
                }

                // Place the N-terminal amine hydrogens according to the
                // configured mode. PRO keeps the NH2+ geometry under every
                // mode — a true neutral secondary amine (single H) would
                // need a separate placement path that does not exist today
                // (see NtermMode.placesPositiveCharge for the matching
                // chemistry-flag logic).
                //
                // PRO needs an explicit gate: PRO's plans contain no backbone
                // amide H entry (its N is already bonded to CA and CD), so
                // `has_backbone` is false — without the `is_pro` bypass, the
                // NH2Pro call below would be dead code and N-terminal PRO
                // would receive no backbone H at all. This bypass wires up
                // placement that was always intended but previously skipped.
                const is_pro = std.mem.eql(u8, comp_id, "PRO");
                if (treat_as_nterm and (has_backbone or is_pro)) {
                    const nterm = if (is_pro)
                        try terminal.placeNtermNH2Pro(mdl, res, @intCast(res_idx), alt, config.bond_policy.mode)
                    else switch (config.nterm_mode) {
                        .auto, .aggressive => try terminal.placeNtermNH3(mdl, res, @intCast(res_idx), alt, config.bond_policy.mode),
                        .neutral => try terminal.placeNtermNH2Neutral(mdl, res, @intCast(res_idx), alt, config.bond_policy.mode),
                    };
                    result.n_placed += nterm.placed;
                    result.n_skipped_existing += nterm.skipped;
                    result.n_skipped_missing_ref += nterm.missing_ref;
                }

                // 3' terminal nucleotide: place HO3' (leaving atom, absent mid-chain).
                if (is_cterm and nucleotide.getPlans(comp_id) != null) {
                    const oh = try terminal.place3primeOH(mdl, res, @intCast(res_idx), alt, config.bond_policy.mode);
                    result.n_placed += oh.placed;
                    result.n_skipped_existing += oh.skipped;
                    result.n_skipped_missing_ref += oh.missing_ref;
                }

                if (try execute.placeOverrideHydrogen(mdl, res, @intCast(res_idx), alt, protonation_state, config.bond_policy.mode)) |place_result| {
                    result.tally(place_result);
                }
            }

            result.n_residues += 1;
        } else {
            // Try inline dict first, then SDF dict, then external CCD (per-component fallback)
            const component = if (inline_dict) |d| d.get(comp_id) else null;
            const sdf_component = if (config.sdf_dict) |d| d.get(comp_id) else null;
            const effective_component = component orelse sdf_component orelse if (ccd_dict) |d| d.get(comp_id) else null;
            if (effective_component) |comp| {
                const existing = try collectAtomNames(mdl.allocator, mdl, res);
                defer mdl.allocator.free(existing);
                const explicit_ir = try collectInterResidueAtomNames(mdl.allocator, mdl, res);
                defer mdl.allocator.free(explicit_ir);

                // For amino-acid-like residues (have N, CA, C), add implicit peptide
                // bond atoms to inter-residue list. Peptide bonds are not in struct_conn
                // but backbone N (non-N-term) and C (non-C-term) participate in them.
                // This ensures leaving H on backbone N (e.g. H2 on 2MR) are skipped.
                const ir_atoms = try addImplicitPeptideBondAtoms(
                    mdl.allocator,
                    mdl,
                    res,
                    explicit_ir,
                    is_nterm_for_bonding,
                    is_cterm,
                );
                defer mdl.allocator.free(ir_atoms);
                const plans = try ccd_derive.derivePlans(mdl.allocator, &comp, existing, ir_atoms);
                defer mdl.allocator.free(plans);

                // Topology-aware H-bond flags (charged acid O, lone-pair N)
                // on top of the element defaults set by applyChemistry.
                for (mdl.atoms.items[res.atom_start..res.atom_end]) |*atom| {
                    if (atom.is_hydrogen) continue;
                    if (ccd_derive.heavyAtomFlags(&comp, atom.nameRaw(), existing)) |extra| {
                        atom.flags = element.mergeFlags(atom.flags, extra);
                    }
                }

                // For non-5'-terminal nucleotides, skip H on phosphate O (HOP2/HOP3).
                // Phosphodiester bonds are implicit polymer bonds (not in struct_conn),
                // so phosphate H should only be placed at the 5' terminus.
                const skip_phosphate_h = !is_nterm_for_bonding and terminal.hasPhosphorusAtom(mdl, res);

                for (plans) |plan| {
                    if (skip_phosphate_h and terminal.isPhosphateH(plan.h_name)) {
                        result.n_skipped_inter_residue += 1;
                        continue;
                    }
                    result.tally(try execute.executePlan(mdl, res, @intCast(res_idx), &plan, null, ' ', config.bond_policy.mode));
                }

                // 3' terminal nucleotide (CCD path): place HO3' only if CCD
                // plans did not already include it (avoids duplicate placement).
                if (is_cterm and terminal.isNucleotideResidue(mdl, res)) {
                    const ho3_name = padName("HO3'");
                    const already_in_plans = for (plans) |plan| {
                        if (std.mem.eql(u8, &plan.h_name, &ho3_name)) break true;
                    } else false;
                    if (!already_in_plans) {
                        const oh = try terminal.place3primeOH(mdl, res, @intCast(res_idx), ' ', config.bond_policy.mode);
                        result.n_placed += oh.placed;
                        result.n_skipped_existing += oh.skipped;
                        result.n_skipped_missing_ref += oh.missing_ref;
                    }
                }

                result.n_residues += 1;
            } else {
                // Last resort: infer bonds from coordinates for unknown residues
                if (try distance_derive.deriveComponentFromCoordinates(mdl.allocator, mdl, res)) |derived| {
                    defer {
                        mdl.allocator.free(derived.atoms);
                        mdl.allocator.free(derived.bonds);
                        mdl.allocator.free(derived.comp_id);
                        mdl.allocator.free(derived.comp_type);
                    }

                    const existing = try collectAtomNames(mdl.allocator, mdl, res);
                    defer mdl.allocator.free(existing);

                    const plans = try ccd_derive.derivePlans(mdl.allocator, &derived, existing, &.{});
                    defer mdl.allocator.free(plans);

                    for (plans) |plan| {
                        result.tally(try execute.executePlan(mdl, res, @intCast(res_idx), &plan, null, ' ', config.bond_policy.mode));
                    }

                    result.n_residues += 1;
                    result.n_distance_derived += 1;
                }
            }
        }
    }

    return result;
}

/// Apply residue/atom-specific chemical annotations to heavy atoms.
/// Updates element_type, flags, and vdw_radius on standard-residue heavy atoms.
/// Also applies terminal charge annotations (N-terminal positive, C-terminal negative).
/// Should be called once after parsing, before hydrogen placement.
pub fn applyChemistry(mdl: *Model) void {
    applyChemistryWithConfig(mdl, .{});
}

pub fn applyChemistryWithConfig(mdl: *Model, config: PlacementConfig) void {
    const n_residues = mdl.residues.items.len;
    for (0..n_residues) |res_idx| {
        const res = mdl.residues.items[res_idx];
        const comp_id = res.compIdSlice();
        const protonation_state = if (config.protonation) |overrides| overrides.find(mdl, res_idx) else null;

        // Detect terminal residues.
        // The positive-charge flag on backbone N is driven by the same rule
        // used for N-terminal H placement (see `addHydrogensWithConfig`):
        // residues receive the POS_D flag iff they get an NH3+/NH2+ placement.
        // PRO is the sole exception under `neutral` mode — it keeps the NH2+
        // geometry (secondary amine at physiological pH) and therefore also
        // keeps the positive charge flag. The two call sites must stay in sync
        // or else chemistry annotation will drift from placed hydrogen atoms.
        const chain = mdl.chains.items[res.chain_idx];
        const is_real_nterm = (res_idx == chain.residue_start);
        const is_break_nterm = res.is_chain_break_before and (res_idx != chain.residue_start);
        const treat_as_nterm = is_real_nterm or
            (is_break_nterm and config.nterm_mode.treatsBreakAsNterm());
        const is_pro = std.mem.eql(u8, comp_id, "PRO");
        const is_positive_nterm = treat_as_nterm and
            config.nterm_mode.placesPositiveCharge(is_pro);
        const is_cterm = (res_idx == chain.residue_end - 1) or
            (res_idx + 1 < n_residues and mdl.residues.items[res_idx + 1].is_chain_break_before);

        const atoms = mdl.atoms.items[res.atom_start..res.atom_end];
        for (atoms) |*atom| {
            if (atom.is_hydrogen) continue;

            // Apply standard residue annotations (replace flags, but preserve bonded_inter_residue)
            const has_std_ann = if (chemistry.getAnnotationWithOverride(comp_id, atom.nameRaw(), protonation_state) orelse
                chemistry.getNucleotideAnnotation(comp_id, atom.nameRaw())) |ann|
            blk: {
                atom.element_type = ann.atom_type;
                const keep_bonded = atom.flags.bonded_inter_residue;
                atom.flags = ann.flags;
                atom.flags.bonded_inter_residue = keep_bonded;
                atom.vdw_radius = ann.atom_type.info().explicit_radius;
                break :blk true;
            } else blk: {
                // No residue-specific annotation (ligands, water, nucleic acids,
                // unannotated amino acid atoms): fall back to the element-type
                // defaults so that e.g. ligand and water oxygens are scored as
                // acceptors. Parsed atoms otherwise keep empty flags and the
                // generic 1.70 Å radius.
                const info = atom.element_type.info();
                const keep_bonded = atom.flags.bonded_inter_residue;
                atom.flags = element.mergeFlags(atom.flags, info.flags);
                atom.flags.bonded_inter_residue = keep_bonded;
                if (atom.element_type != .unknown) atom.vdw_radius = info.explicit_radius;
                break :blk false;
            };

            // Apply terminal annotations (merge flags via OR)
            // Only set element_type/vdw_radius for atoms without standard annotation (e.g. OXT)
            if (chemistry.getTerminalAnnotation(atom.nameRaw(), is_positive_nterm, is_cterm)) |term_ann| {
                atom.flags = element.mergeFlags(atom.flags, term_ann.flags);
                if (!has_std_ann) {
                    atom.element_type = term_ann.atom_type;
                    atom.vdw_radius = term_ann.atom_type.info().explicit_radius;
                }
            }
        }
    }
}

fn shouldSkipPlanForProtonation(comp_id: []const u8, plan: *const standard.PlacementPlan, state: ?protonation.ResidueState) bool {
    const s = state orelse return false;

    if (std.mem.eql(u8, comp_id, "HIS") and s == .his) {
        if (s.his == .hid and std.mem.eql(u8, &plan.h_name, " HE2")) return true;
        if (s.his == .hie and std.mem.eql(u8, &plan.h_name, " HD1")) return true;
        return false;
    }
    if (std.mem.eql(u8, comp_id, "LYS") and s == .lys and s.lys == .neutral) {
        return std.mem.eql(u8, &plan.h_name, " HZ3");
    }
    if (std.mem.eql(u8, comp_id, "CYS") and s == .cys and s.cys == .thiolate) {
        return std.mem.eql(u8, &plan.h_name, " HG ");
    }
    return false;
}

/// Find the raw [4]u8 name of an atom in the model (using nameMatch for lookup).
/// Returns the actual stored name bytes, which may differ from PDB-padded format.
fn findAtomName(mdl: *const Model, res: Residue, name: [4]u8, target_altloc: u8) ?[4]u8 {
    const a = findAtom(mdl, res, name, target_altloc) orelse return null;
    return a.nameRaw();
}

/// Extend inter-residue atom names with implicit peptide bond atoms.
/// For amino-acid-like residues (having N, CA, C), backbone N participates in
/// a peptide bond to the previous residue (unless N-terminal), and backbone C
/// participates in a peptide bond to the next residue (unless C-terminal).
/// These implicit bonds are not recorded in struct_conn but affect leaving atom logic.
fn addImplicitPeptideBondAtoms(
    allocator: std.mem.Allocator,
    mdl: *const Model,
    res: Residue,
    explicit_ir: []const [4]u8,
    is_nterm: bool,
    is_cterm: bool,
) ![]const [4]u8 {
    // Check if residue has backbone atoms (amino-acid-like)
    const has_n = findAtomPos(mdl, res, .{ ' ', 'N', ' ', ' ' }, ' ') != null;
    const has_ca = findAtomPos(mdl, res, .{ ' ', 'C', 'A', ' ' }, ' ') != null;
    const has_c = findAtomPos(mdl, res, .{ ' ', 'C', ' ', ' ' }, ' ') != null;

    if (!(has_n and has_ca and has_c)) return try allocator.dupe([4]u8, explicit_ir);

    var list = std.ArrayListUnmanaged([4]u8).empty;
    defer list.deinit(allocator);
    try list.appendSlice(allocator, explicit_ir);

    // Non-N-terminal: backbone N has implicit peptide bond to previous C.
    // Use model atom name format (setName style, not PDB-padded).
    if (!is_nterm) try list.append(allocator, findAtomName(mdl, res, .{ ' ', 'N', ' ', ' ' }, ' ') orelse .{ 'N', ' ', ' ', ' ' });
    // Non-C-terminal: backbone C has implicit peptide bond to next N
    if (!is_cterm) try list.append(allocator, findAtomName(mdl, res, .{ ' ', 'C', ' ', ' ' }, ' ') orelse .{ 'C', ' ', ' ', ' ' });

    return try allocator.dupe([4]u8, list.items);
}

/// Collect existing atom names in a residue as [4]u8 arrays.
/// Caller must free the returned slice with the provided allocator.
fn collectAtomNames(allocator: std.mem.Allocator, mdl: *const Model, res: Residue) ![][4]u8 {
    const atoms = mdl.atoms.items[res.atom_start..res.atom_end];
    const result = try allocator.alloc([4]u8, atoms.len);
    for (atoms, 0..) |a, i| {
        result[i] = a.nameRaw();
    }
    return result;
}

/// Collect names of heavy atoms in this residue that carry an inter-residue bond.
/// Iterates mdl.bonds so that an atom involved in multiple inter-residue bonds
/// produces multiple entries (one per bond), giving correct valence accounting.
fn collectInterResidueAtomNames(allocator: std.mem.Allocator, mdl: *const Model, res: Residue) ![]const [4]u8 {
    const bond_mod = @import("../model/bond.zig");
    var names = std.ArrayListUnmanaged([4]u8).empty;
    defer names.deinit(allocator);
    for (mdl.bonds.items) |bond| {
        if (bond.source == bond_mod.BondSource.struct_conn or bond.source == bond_mod.BondSource.branch_link) {
            if (bond.atom_1 >= res.atom_start and bond.atom_1 < res.atom_end) {
                const atom = mdl.atoms.items[bond.atom_1];
                if (!atom.is_hydrogen) try names.append(allocator, atom.nameRaw());
            }
            if (bond.atom_2 >= res.atom_start and bond.atom_2 < res.atom_end) {
                const atom = mdl.atoms.items[bond.atom_2];
                if (!atom.is_hydrogen) try names.append(allocator, atom.nameRaw());
            }
        }
    }
    return try allocator.dupe([4]u8, names.items);
}
