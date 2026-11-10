//! Residue/atom-specific chemical annotations for the 20 standard amino acids.
//!
//! Atom names follow PDB convention (4-character, space-padded).
//! Provides atom type and flag overrides beyond the generic element-based defaults.

const std = @import("std");
const element = @import("../element.zig");
const protonation = @import("protonation.zig");

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

pub const ChemAnnotation = struct {
    atom_type: element.AtomType,
    flags: element.AtomFlags,
};

// ---------------------------------------------------------------------------
// Flag constants (shorthand for common combinations)
// ---------------------------------------------------------------------------

const A = element.AtomFlags{ .acceptor = true };
const D = element.AtomFlags{ .donor = true };
const DA = element.AtomFlags{ .donor = true, .acceptor = true };
const ARA = element.AtomFlags{ .aromatic = true, .acceptor = true };
const ARD = element.AtomFlags{ .aromatic = true, .donor = true };
const ARDA = element.AtomFlags{ .aromatic = true, .donor = true, .acceptor = true };
const ARP = element.AtomFlags{ .aromatic = true, .donor = true, .positive = true };
const NEG_A = element.AtomFlags{ .negative = true, .acceptor = true };
const POS_D = element.AtomFlags{ .positive = true, .donor = true };
const NONE = element.AtomFlags{};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

const n = @import("lookup.zig").padName;

const AnnotEntry = struct {
    name: [4]u8,
    ann: ChemAnnotation,
};

fn a(comptime atom_name: []const u8, atom_type: element.AtomType, flags: element.AtomFlags) AnnotEntry {
    return .{ .name = n(atom_name), .ann = .{ .atom_type = atom_type, .flags = flags } };
}

/// Returns start and end indices of the non-space content in a 4-char name.
fn trimBounds(name: *const [4]u8) struct { start: usize, end: usize } {
    var start: usize = 0;
    while (start < 4 and name[start] == ' ') start += 1;
    var end: usize = 4;
    while (end > start and name[end - 1] == ' ') end -= 1;
    return .{ .start = start, .end = end };
}

/// Check whether two 4-char names match after trimming leading/trailing spaces.
fn nameEql(x: *const [4]u8, y: *const [4]u8) bool {
    const xb = trimBounds(x);
    const yb = trimBounds(y);
    const x_len = xb.end - xb.start;
    const y_len = yb.end - yb.start;
    if (x_len != y_len) return false;
    for (0..x_len) |i| {
        if (x[xb.start + i] != y[yb.start + i]) return false;
    }
    return true;
}

fn lookupName(entries: []const AnnotEntry, atom_name: [4]u8) ?ChemAnnotation {
    for (entries) |*entry| {
        if (nameEql(&entry.name, &atom_name)) return entry.ann;
    }
    return null;
}

// ---------------------------------------------------------------------------
// Backbone annotations (shared by all standard residues)
// ---------------------------------------------------------------------------

const backbone_annotations = [_]AnnotEntry{
    a(" C  ", .C_eq_O, NONE),
    a(" O  ", .O, A),
    a(" N  ", .N, D),
};

// ---------------------------------------------------------------------------
// Side-chain annotation tables per residue
// ---------------------------------------------------------------------------

const asp_sc = [_]AnnotEntry{
    a("OD1 ", .O, NEG_A),
    a("OD2 ", .O, NEG_A),
    a(" CG ", .C_eq_O, NONE),
};

const glu_sc = [_]AnnotEntry{
    a("OE1 ", .O, NEG_A),
    a("OE2 ", .O, NEG_A),
    a(" CD ", .C_eq_O, NONE),
};

const asn_sc = [_]AnnotEntry{
    a("OD1 ", .O, A),
    a(" CG ", .C_eq_O, NONE),
    a("ND2 ", .N, D),
};

const gln_sc = [_]AnnotEntry{
    a("OE1 ", .O, A),
    a(" CD ", .C_eq_O, NONE),
    a("NE2 ", .N, D),
};

const lys_sc = [_]AnnotEntry{
    a(" NZ ", .N, POS_D),
};

const arg_sc = [_]AnnotEntry{
    a(" NE ", .N, POS_D),
    a("NH1 ", .N, POS_D),
    a("NH2 ", .N, POS_D),
    a(" CZ ", .C_eq_O, NONE),
};

const ser_sc = [_]AnnotEntry{
    a(" OG ", .O, DA),
};

const thr_sc = [_]AnnotEntry{
    a("OG1 ", .O, DA),
};

const cys_sc = [_]AnnotEntry{
    a(" SG ", .S, A),
};

const met_sc = [_]AnnotEntry{
    a(" SD ", .S, A),
};

const tyr_sc = [_]AnnotEntry{
    a(" OH ", .O, DA),
    a(" CG ", .Car, ARA),
    a("CD1 ", .Car, ARA),
    a("CD2 ", .Car, ARA),
    a("CE1 ", .Car, ARA),
    a("CE2 ", .Car, ARA),
    a(" CZ ", .Car, ARA),
};

const phe_sc = [_]AnnotEntry{
    a(" CG ", .Car, ARA),
    a("CD1 ", .Car, ARA),
    a("CD2 ", .Car, ARA),
    a("CE1 ", .Car, ARA),
    a("CE2 ", .Car, ARA),
    a(" CZ ", .Car, ARA),
};

const trp_sc = [_]AnnotEntry{
    a(" CG ", .Car, ARA),
    a("CD1 ", .Car, ARA),
    a("CD2 ", .Car, ARA),
    a("CE2 ", .Car, ARA),
    a("CE3 ", .Car, ARA),
    a("CZ2 ", .Car, ARA),
    a("CZ3 ", .Car, ARA),
    a("CH2 ", .Car, ARA),
    a("NE1 ", .N, ARD),
};

const his_sc = [_]AnnotEntry{
    a(" CG ", .Car, ARA),
    a("CD2 ", .Car, ARA),
    a("CE1 ", .Car, ARA),
    a("ND1 ", .Nacc, ARDA), // donor+acceptor: protonation-dependent
    a("NE2 ", .Nacc, ARDA), // donor+acceptor: protonation-dependent
};

// ---------------------------------------------------------------------------
// Per-residue side-chain dispatch
// ---------------------------------------------------------------------------

fn getSideChainSlice(comp_id: []const u8) ?[]const AnnotEntry {
    if (std.mem.eql(u8, comp_id, "ASP")) return &asp_sc;
    if (std.mem.eql(u8, comp_id, "GLU")) return &glu_sc;
    if (std.mem.eql(u8, comp_id, "ASN")) return &asn_sc;
    if (std.mem.eql(u8, comp_id, "GLN")) return &gln_sc;
    if (std.mem.eql(u8, comp_id, "LYS")) return &lys_sc;
    if (std.mem.eql(u8, comp_id, "ARG")) return &arg_sc;
    if (std.mem.eql(u8, comp_id, "SER")) return &ser_sc;
    if (std.mem.eql(u8, comp_id, "THR")) return &thr_sc;
    if (std.mem.eql(u8, comp_id, "CYS")) return &cys_sc;
    if (std.mem.eql(u8, comp_id, "MET")) return &met_sc;
    if (std.mem.eql(u8, comp_id, "TYR")) return &tyr_sc;
    if (std.mem.eql(u8, comp_id, "PHE")) return &phe_sc;
    if (std.mem.eql(u8, comp_id, "TRP")) return &trp_sc;
    if (std.mem.eql(u8, comp_id, "HIS")) return &his_sc;
    // ALA, GLY, VAL, LEU, ILE, PRO — no side-chain annotations
    return null;
}

/// Returns true if comp_id is one of the 20 standard amino acids.
fn isStandardAminoAcid(comp_id: []const u8) bool {
    const standard_map = std.StaticStringMap(void).initComptime(.{
        .{ "ALA", {} }, .{ "GLY", {} }, .{ "VAL", {} }, .{ "LEU", {} },
        .{ "ILE", {} }, .{ "PRO", {} }, .{ "PHE", {} }, .{ "TYR", {} },
        .{ "TRP", {} }, .{ "SER", {} }, .{ "THR", {} }, .{ "CYS", {} },
        .{ "MET", {} }, .{ "ASP", {} }, .{ "GLU", {} }, .{ "ASN", {} },
        .{ "GLN", {} }, .{ "LYS", {} }, .{ "ARG", {} }, .{ "HIS", {} },
    });
    return standard_map.has(comp_id);
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Returns chemical annotation for a given residue and atom name.
///
/// Returns null if:
/// - comp_id is not one of the 20 standard amino acids
/// - atom_name has no specific annotation (e.g. ALA CB — use generic element type)
pub fn getAnnotation(comp_id: []const u8, atom_name: [4]u8) ?ChemAnnotation {
    return getAnnotationWithOverride(comp_id, atom_name, null);
}

pub fn getAnnotationWithOverride(comp_id: []const u8, atom_name: [4]u8, state: ?protonation.ResidueState) ?ChemAnnotation {
    if (!isStandardAminoAcid(comp_id)) return null;

    if (state) |s| {
        if (getOverrideAnnotation(comp_id, atom_name, s)) |ann| return ann;
    }

    // Check side-chain annotations first
    if (getSideChainSlice(comp_id)) |sc_entries| {
        if (lookupName(sc_entries, atom_name)) |ann| return ann;
    }

    // Fall back to backbone annotations
    return lookupName(&backbone_annotations, atom_name);
}

fn getOverrideAnnotation(comp_id: []const u8, atom_name: [4]u8, state: protonation.ResidueState) ?ChemAnnotation {
    if (std.mem.eql(u8, comp_id, "HIS") and state == .his) {
        if (nameEql(&atom_name, &n("ND1 "))) {
            return switch (state.his) {
                .auto => .{ .atom_type = .Nacc, .flags = ARDA },
                .hid => .{ .atom_type = .N, .flags = ARD },
                .hie => .{ .atom_type = .Nacc, .flags = ARA },
                .hip => .{ .atom_type = .N, .flags = ARP },
            };
        }
        if (nameEql(&atom_name, &n("NE2 "))) {
            return switch (state.his) {
                .auto => .{ .atom_type = .Nacc, .flags = ARDA },
                .hid => .{ .atom_type = .Nacc, .flags = ARA },
                .hie => .{ .atom_type = .N, .flags = ARD },
                .hip => .{ .atom_type = .N, .flags = ARP },
            };
        }
    }

    if (std.mem.eql(u8, comp_id, "ASP") and state == .asp) {
        if (nameEql(&atom_name, &n("OD1 "))) {
            return switch (state.asp) {
                .deprotonated => .{ .atom_type = .O, .flags = NEG_A },
                .atom1, .protonated => .{ .atom_type = .O, .flags = D },
                .atom2 => .{ .atom_type = .O, .flags = A },
            };
        }
        if (nameEql(&atom_name, &n("OD2 "))) {
            return switch (state.asp) {
                .deprotonated => .{ .atom_type = .O, .flags = NEG_A },
                .atom1, .protonated => .{ .atom_type = .O, .flags = A },
                .atom2 => .{ .atom_type = .O, .flags = D },
            };
        }
    }

    if (std.mem.eql(u8, comp_id, "GLU") and state == .glu) {
        if (nameEql(&atom_name, &n("OE1 "))) {
            return switch (state.glu) {
                .deprotonated => .{ .atom_type = .O, .flags = NEG_A },
                .atom1, .protonated => .{ .atom_type = .O, .flags = D },
                .atom2 => .{ .atom_type = .O, .flags = A },
            };
        }
        if (nameEql(&atom_name, &n("OE2 "))) {
            return switch (state.glu) {
                .deprotonated => .{ .atom_type = .O, .flags = NEG_A },
                .atom1, .protonated => .{ .atom_type = .O, .flags = A },
                .atom2 => .{ .atom_type = .O, .flags = D },
            };
        }
    }

    if (std.mem.eql(u8, comp_id, "LYS") and state == .lys and nameEql(&atom_name, &n(" NZ "))) {
        return switch (state.lys) {
            .charged => .{ .atom_type = .N, .flags = POS_D },
            .neutral => .{ .atom_type = .N, .flags = D },
        };
    }

    if (std.mem.eql(u8, comp_id, "CYS") and state == .cys and nameEql(&atom_name, &n(" SG "))) {
        return switch (state.cys) {
            .thiol => .{ .atom_type = .S, .flags = A },
            .thiolate => .{ .atom_type = .S, .flags = NEG_A },
        };
    }

    return null;
}

/// Returns terminal-specific annotation for an atom based on its terminal state.
/// These annotations provide additional flags (charge) to be OR-merged with
/// standard annotations, not used as replacements.
pub fn getTerminalAnnotation(atom_name: [4]u8, is_nterm: bool, is_cterm: bool) ?ChemAnnotation {
    if (is_nterm) {
        if (nameEql(&atom_name, &n(" N  "))) {
            return .{ .atom_type = .N, .flags = POS_D };
        }
    }
    if (is_cterm) {
        if (nameEql(&atom_name, &n("OXT ")) or nameEql(&atom_name, &n(" OXT"))) {
            return .{ .atom_type = .O, .flags = NEG_A };
        }
        if (nameEql(&atom_name, &n(" O  "))) {
            return .{ .atom_type = .O, .flags = NEG_A };
        }
    }
    return null;
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "backbone C is carbonyl type" {
    const ann = getAnnotation("ALA", n(" C  "));
    try std.testing.expect(ann != null);
    try std.testing.expectEqual(element.AtomType.C_eq_O, ann.?.atom_type);
}

test "backbone O is acceptor" {
    const ann = getAnnotation("ALA", n(" O  "));
    try std.testing.expect(ann != null);
    try std.testing.expect(ann.?.flags.acceptor == true);
}

test "backbone N is donor" {
    const ann = getAnnotation("ALA", n(" N  "));
    try std.testing.expect(ann != null);
    try std.testing.expect(ann.?.flags.donor == true);
}

test "ASP OD1 is negative acceptor" {
    const ann = getAnnotation("ASP", n("OD1 "));
    try std.testing.expect(ann != null);
    try std.testing.expect(ann.?.flags.negative == true);
    try std.testing.expect(ann.?.flags.acceptor == true);
}

test "LYS NZ is positive donor" {
    const ann = getAnnotation("LYS", n(" NZ "));
    try std.testing.expect(ann != null);
    try std.testing.expect(ann.?.flags.positive == true);
    try std.testing.expect(ann.?.flags.donor == true);
}

test "PHE CD1 is aromatic carbon" {
    const ann = getAnnotation("PHE", n("CD1 "));
    try std.testing.expect(ann != null);
    try std.testing.expectEqual(element.AtomType.Car, ann.?.atom_type);
    try std.testing.expect(ann.?.flags.aromatic == true);
}

test "HIS NE2 is aromatic donor/acceptor nitrogen" {
    const ann = getAnnotation("HIS", n("NE2 "));
    try std.testing.expect(ann != null);
    try std.testing.expectEqual(element.AtomType.Nacc, ann.?.atom_type);
    try std.testing.expect(ann.?.flags.aromatic == true);
    try std.testing.expect(ann.?.flags.acceptor == true);
    try std.testing.expect(ann.?.flags.donor == true);
}

test "HIS override maps ND1 and NE2 chemistry" {
    const hid_nd1 = getAnnotationWithOverride("HIS", n("ND1 "), .{ .his = .hid }).?;
    try std.testing.expectEqual(element.AtomType.N, hid_nd1.atom_type);
    try std.testing.expect(hid_nd1.flags.donor);
    try std.testing.expect(!hid_nd1.flags.acceptor);

    const hid_ne2 = getAnnotationWithOverride("HIS", n("NE2 "), .{ .his = .hid }).?;
    try std.testing.expectEqual(element.AtomType.Nacc, hid_ne2.atom_type);
    try std.testing.expect(hid_ne2.flags.acceptor);
    try std.testing.expect(!hid_ne2.flags.donor);
}

test "ASP protonated override removes negative flag from protonated oxygen" {
    const od2 = getAnnotationWithOverride("ASP", n("OD2 "), .{ .asp = .atom2 }).?;
    try std.testing.expect(od2.flags.donor);
    try std.testing.expect(!od2.flags.negative);

    const od1 = getAnnotationWithOverride("ASP", n("OD1 "), .{ .asp = .atom2 }).?;
    try std.testing.expect(od1.flags.acceptor);
    try std.testing.expect(!od1.flags.negative);
}

test "ALA CB has no annotation" {
    const ann = getAnnotation("ALA", n(" CB "));
    try std.testing.expect(ann == null);
}

test "unknown residue returns null" {
    const ann = getAnnotation("XYZ", n(" C  "));
    try std.testing.expect(ann == null);
}

test "N-terminal N gets positive flag" {
    const ann = getTerminalAnnotation(n(" N  "), true, false);
    try std.testing.expect(ann != null);
    try std.testing.expect(ann.?.flags.positive);
}

test "C-terminal O gets negative flag" {
    const ann = getTerminalAnnotation(n(" O  "), false, true);
    try std.testing.expect(ann != null);
    try std.testing.expect(ann.?.flags.negative);
}

test "OXT gets negative acceptor flags" {
    const ann = getTerminalAnnotation(n("OXT "), false, true);
    try std.testing.expect(ann != null);
    try std.testing.expectEqual(element.AtomType.O, ann.?.atom_type);
    try std.testing.expect(ann.?.flags.negative);
    try std.testing.expect(ann.?.flags.acceptor);
}

test "internal residue atoms get no terminal annotation" {
    const ann = getTerminalAnnotation(n(" N  "), false, false);
    try std.testing.expect(ann == null);
}

test "C-terminal N gets no annotation" {
    const ann = getTerminalAnnotation(n(" N  "), false, true);
    try std.testing.expect(ann == null);
}
