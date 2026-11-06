//! PROPKA 3.5.1 engine for `zreduce pka`, via the optional zpropka package.
//!
//! zpropka is a separate LGPL-2.1 package; it is only compiled in with
//! `zig build -Dpropka=true`. Without it, `available` is false and `predict`
//! returns `error.PropkaNotBuilt`, so the default build stays MIT-only.
//!
//! The model is written as PDB text (heavy atoms only; PROPKA re-protonates),
//! with chain IDs aliased to single characters, and PROPKA groups are mapped
//! back to residues. Disulfide and metal-coordinated residues are moved to
//! `excluded` exactly as for the built-in engine.

const std = @import("std");
const build_options = @import("build_options");
const model_mod = @import("model.zig");
const pka = @import("pka.zig");
const format = @import("writer/format.zig");
const pdb_writer = @import("writer/pdb_writer.zig");

const zpropka = if (build_options.propka) @import("zpropka") else struct {};

pub const available = build_options.propka;
pub const engine_version = "PROPKA 3.5.1 (zpropka)";
pub const parameter_version = "propka-3.5.1-default";

const Model = model_mod.Model;

/// Residue key used to map PROPKA groups back: chain alias, number, icode.
const Key = struct { chain: u8, res_num: i32, icode: u8 };

const Encoded = struct {
    pdb: []u8,
    /// (chain alias, auth_seq_id, ins_code) -> residue index
    residues: std.AutoHashMapUnmanaged(Key, u32),

    fn deinit(self: *Encoded, allocator: std.mem.Allocator) void {
        allocator.free(self.pdb);
        self.residues.deinit(allocator);
    }
};

const alias_pool = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";

/// Encode heavy atoms as PDB records. Auth chain IDs that are a single,
/// unique character are kept; others get an unused alias.
fn encodePdb(allocator: std.mem.Allocator, mdl: *const Model, options: Options) !Encoded {
    var aliases: std.StringHashMapUnmanaged(u8) = .empty;
    defer aliases.deinit(allocator);
    var used = [_]bool{false} ** 256;
    // Keep single-character auth IDs first so they map to themselves.
    for (mdl.chains.items) |*chain| {
        const auth = chain.authSlice();
        if (auth.len == 1 and !used[auth[0]] and std.mem.indexOfScalar(u8, alias_pool, auth[0]) != null) {
            if (aliases.get(auth) == null) {
                try aliases.put(allocator, auth, auth[0]);
                used[auth[0]] = true;
            }
        }
    }
    for (mdl.chains.items) |*chain| {
        const auth = chain.authSlice();
        if (aliases.get(auth) != null) continue;
        const alias = for (alias_pool) |c| {
            if (!used[c]) break c;
        } else return error.TooManyChains;
        used[alias] = true;
        try aliases.put(allocator, auth, alias);
    }

    var out: std.Io.Writer.Allocating = .init(allocator);
    errdefer out.deinit();
    var residues: std.AutoHashMapUnmanaged(Key, u32) = .empty;
    errdefer residues.deinit(allocator);
    const w = &out.writer;
    var serial: u32 = 0;
    // PROPKA resets its N-terminus detection at TER (and after OXT), so a TER
    // must follow every polymer chain or the next chain loses its N+ group.
    var prev_chain: ?u32 = null;
    var prev_polymer = false;
    for (mdl.atoms.items) |atom| {
        if (atom.is_hydrogen) continue;
        const res = mdl.residues.items[atom.residue_idx];
        if (!options.keep_waters and isWater(res)) continue;
        if (res.auth_seq_id < -999 or res.auth_seq_id > 9999) return error.ResidueNumberOutOfRange;
        if (prev_chain) |pc| {
            if (pc != res.chain_idx and prev_polymer) try w.writeAll("TER   \n");
        }
        prev_chain = res.chain_idx;
        prev_polymer = res.entity_type == .polymer;
        const alias = aliases.get(mdl.chains.items[res.chain_idx].authSlice()).?;
        try residues.put(allocator, .{ .chain = alias, .res_num = res.auth_seq_id, .icode = res.ins_code }, atom.residue_idx);
        serial = serial % 99999 + 1;
        // PDB convention (and PROPKA's ligand handling): only standard
        // residues are ATOM; modified polymer residues are HETATM.
        const record = if (res.entity_type == .polymer and isStandardResidue(res.compIdSlice())) "ATOM  " else "HETATM";
        const sym = format.elementSymbol(atom.element_type);
        const name4 = pdb_writer.formatAtomNamePdb(atom.nameSlice(), sym);
        const comp = res.compIdSlice();
        // Zig prints '+' for padded signed integers, so format the residue
        // number as unsigned text first.
        var num_buf: [8]u8 = undefined;
        const num = if (res.auth_seq_id < 0)
            try std.fmt.bufPrint(&num_buf, "-{d}", .{@abs(res.auth_seq_id)})
        else
            try std.fmt.bufPrint(&num_buf, "{d}", .{@as(u32, @intCast(res.auth_seq_id))});
        try w.print("{s}{d:>5} {s}{c}{s:>3} {c}{s:>4}{c}   {d:>8.3}{d:>8.3}{d:>8.3}{d:>6.2}{d:>6.2}          {s:>2}\n", .{
            record,        serial,       &name4,     atom.altloc, comp[0..@min(comp.len, 3)], alias,
            num,           res.ins_code, atom.pos.x, atom.pos.y,  atom.pos.z,                 atom.occupancy,
            atom.b_factor, sym,
        });
    }
    if (prev_polymer) try w.writeAll("TER   \n");
    try w.writeAll("END\n");
    return .{ .pdb = try out.toOwnedSlice(), .residues = residues };
}

fn isStandardResidue(comp: []const u8) bool {
    const ids = std.StaticStringMap(void).initComptime(.{
        .{ "ALA", {} }, .{ "ARG", {} }, .{ "ASN", {} }, .{ "ASP", {} }, .{ "CYS", {} },
        .{ "GLN", {} }, .{ "GLU", {} }, .{ "GLY", {} }, .{ "HIS", {} }, .{ "ILE", {} },
        .{ "LEU", {} }, .{ "LYS", {} }, .{ "MET", {} }, .{ "PHE", {} }, .{ "PRO", {} },
        .{ "SER", {} }, .{ "THR", {} }, .{ "TRP", {} }, .{ "TYR", {} }, .{ "VAL", {} },
        .{ "UNK", {} }, .{ "DA", {} },  .{ "DC", {} },  .{ "DG", {} },  .{ "DT", {} },
        .{ "DI", {} },  .{ "A", {} },   .{ "C", {} },   .{ "G", {} },   .{ "U", {} },
        .{ "I", {} },
    });
    return ids.has(comp);
}

fn isWater(res: model_mod.Residue) bool {
    if (res.entity_type == .water) return true;
    const comp = res.compIdSlice();
    return std.mem.eql(u8, comp, "HOH") or std.mem.eql(u8, comp, "WAT") or std.mem.eql(u8, comp, "DOD");
}

pub const Options = struct {
    /// PROPKA has no special handling for water: crystal waters are read as
    /// HETATM and enter its burial and interaction terms. Waters are dropped
    /// by default, as PDB2PQR does before calling PROPKA (on PKAD-R this made
    /// no measurable difference: test MAE 0.677 vs 0.678). Keep them to
    /// reproduce PROPKA run directly on a deposited file.
    keep_waters: bool = false,
};

fn siteKind(residue_type: []const u8) ?pka.SiteKind {
    const map = std.StaticStringMap(pka.SiteKind).initComptime(.{
        .{ "ASP", .asp },   .{ "GLU", .glu }, .{ "HIS", .his }, .{ "LYS", .lys },
        .{ "ARG", .arg },   .{ "CYS", .cys }, .{ "TYR", .tyr }, .{ "N+", .n_term },
        .{ "C-", .c_term },
    });
    return map.get(residue_type);
}

/// Run PROPKA on `mdl` and return a prediction in the same form as the
/// built-in engine. Ligand groups are not reported (counted as omitted).
pub fn predict(allocator: std.mem.Allocator, mdl: *const Model, options: Options) !pka.Prediction {
    if (!available) return error.PropkaNotBuilt;

    var encoded = try encodePdb(allocator, mdl, options);
    defer encoded.deinit(allocator);
    var result = try zpropka.api.predictPdb(allocator, encoded.pdb, "zreduce", .{});
    defer result.deinit();

    const metals = try pka.metalPositions(allocator, mdl);
    defer allocator.free(metals);

    var sites: std.ArrayListUnmanaged(pka.SiteResult) = .empty;
    errdefer sites.deinit(allocator);
    var excluded: std.ArrayListUnmanaged(pka.ExcludedSite) = .empty;
    errdefer excluded.deinit(allocator);
    var omitted: u32 = 0;

    for (result.groups) |g| {
        const kind = siteKind(g.residue_type) orelse {
            omitted += 1; // ligand groups: outside this report's scope
            continue;
        };
        if (g.ligand) {
            omitted += 1;
            continue;
        }
        const icode: u8 = if (g.icode.len == 1) g.icode[0] else ' ';
        const chain: u8 = if (g.chain_id.len == 1) g.chain_id[0] else ' ';
        const res_num = std.math.cast(i32, g.res_num) orelse {
            omitted += 1;
            continue;
        };
        const residue_idx = encoded.residues.get(.{ .chain = chain, .res_num = res_num, .icode = icode }) orelse {
            omitted += 1;
            continue;
        };
        if (pka.exclusionReason(mdl, residue_idx, kind, metals)) |reason| {
            try excluded.append(allocator, .{ .residue_idx = residue_idx, .kind = kind, .reason = reason });
            continue;
        }
        if (g.pka >= 99.0) {
            // PROPKA's "not titrating" marker for a covalently coupled group
            // that our exclusion rules did not recognise.
            omitted += 1;
            continue;
        }
        var site = pka.SiteResult{
            .residue_idx = residue_idx,
            .kind = kind,
            .center = .{ .x = 0, .y = 0, .z = 0 },
            .terms = .{
                .model = @floatCast(g.model_pka),
                .burial = @floatCast(g.desolvation),
                .hydrogen_bond = @floatCast(g.sidechain_hbond + g.backbone_hbond),
                .electrostatic = @floatCast(g.coulomb),
            },
            .predicted_pka = @floatCast(g.pka),
            .burial_fraction = @floatCast(std.math.clamp(g.buried, 0.0, 1.0)),
            .confidence = .high,
            .coupled = g.coupled,
        };
        site.confidence = pka.assessConfidence(site);
        // PROPKA can return unphysical values (e.g. > 50 for some coupled
        // cysteines); keep them, flagged, rather than silently clipping.
        if (!pka.isPhysicalPka(site.predicted_pka)) site.confidence = .low;
        try sites.append(allocator, site);
    }

    return .{
        .allocator = allocator,
        .sites = try sites.toOwnedSlice(allocator),
        .excluded = try excluded.toOwnedSlice(allocator),
        .omitted_sites = omitted,
        .parameter_version = parameter_version,
        .engine = engine_version,
    };
}

test "PROPKA engine reports availability consistently" {
    if (available) return;
    var mdl = Model.init(std.testing.allocator);
    defer mdl.deinit();
    try std.testing.expectError(error.PropkaNotBuilt, predict(std.testing.allocator, &mdl, .{}));
}
