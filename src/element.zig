//! Element type system with VDW radii and atom flags.

const std = @import("std");

pub const AtomFlags = packed struct {
    donor: bool = false,
    acceptor: bool = false,
    aromatic: bool = false,
    positive: bool = false,
    negative: bool = false,
    metallic: bool = false,
    hb_only_dummy: bool = false,
    bonded_inter_residue: bool = false,
};

pub fn mergeFlags(a: AtomFlags, b: AtomFlags) AtomFlags {
    return @bitCast(@as(u8, @bitCast(a)) | @as(u8, @bitCast(b)));
}

pub const MoverHint = enum(u3) {
    none,
    rotate, // OH, SH
    rotate_nh3, // NH3+
    rotate_methyl, // CH3: placed staggered and not optimized (no mover)
    rotate_nh2, // NH2 (N-terminal PRO NH2+, neutral primary amine)
    flip_amide, // Asn/Gln
    flip_his, // His
    flip_carboxyl, // protonated Asp/Glu (proton on either oxygen)
};

pub const AtomTypeInfo = struct {
    explicit_radius: f32,
    implicit_radius: f32,
    covalent_radius: f32,
    flags: AtomFlags,
};

pub const AtomType = enum(u8) {
    // Hydrogen variants
    H, // non-polar H (1.22 Å)
    Har, // aromatic H (1.05 Å)
    Hpol, // polar H, donor (1.05 Å)
    Ha_p, // aromatic + polar H (1.05 Å)
    HOd, // H-bond only dummy (1.05 Å)

    // Carbon variants
    C, // sp3 carbon (1.70 Å)
    Car, // aromatic carbon (1.75 Å, acceptor)
    C_eq_O, // carbonyl carbon (1.65 Å)

    // Nitrogen variants
    N, // nitrogen (1.55 Å)
    Nacc, // nitrogen acceptor (1.55 Å)

    // Others
    O, // oxygen (1.40 Å, acceptor)
    P, // phosphorus (1.80 Å)
    S, // sulfur (1.80 Å, acceptor)
    Se, // selenium (1.90 Å)
    F, // fluorine (1.30 Å, acceptor)
    Cl, // chlorine (1.77 Å, acceptor)
    Br, // bromine (1.95 Å, acceptor)
    I, // iodine (2.10 Å, acceptor)

    // Metals
    Li, Na, Mg, K, Ca, Mn, Fe, Co, Ni, Cu, Zn, As, Rb, Sr, Mo, Ag, Cd, Sn, Cs, Ba, W, Pt, Au, Hg, Pb, U,

    unknown,

    const count = @typeInfo(AtomType).@"enum".fields.len;

    pub fn info(self: AtomType) AtomTypeInfo {
        return atom_type_table[@intFromEnum(self)];
    }
};

const NONE = AtomFlags{};
const D = AtomFlags{ .donor = true };
const A = AtomFlags{ .acceptor = true };
const AR = AtomFlags{ .aromatic = true };
const ARA = AtomFlags{ .aromatic = true, .acceptor = true };
const DAR = AtomFlags{ .donor = true, .aromatic = true };
const M = AtomFlags{ .metallic = true };
const HBD = AtomFlags{ .donor = true, .hb_only_dummy = true };

const atom_type_table: [AtomType.count]AtomTypeInfo = build: {
    var table: [AtomType.count]AtomTypeInfo = undefined;

    // Hydrogen variants
    table[@intFromEnum(AtomType.H)] = .{ .explicit_radius = 1.22, .implicit_radius = 0.00, .covalent_radius = 0.30, .flags = NONE };
    table[@intFromEnum(AtomType.Har)] = .{ .explicit_radius = 1.05, .implicit_radius = 0.00, .covalent_radius = 0.30, .flags = AR };
    table[@intFromEnum(AtomType.Hpol)] = .{ .explicit_radius = 1.05, .implicit_radius = 0.00, .covalent_radius = 0.30, .flags = D };
    table[@intFromEnum(AtomType.Ha_p)] = .{ .explicit_radius = 1.05, .implicit_radius = 0.00, .covalent_radius = 0.30, .flags = DAR };
    table[@intFromEnum(AtomType.HOd)] = .{ .explicit_radius = 1.05, .implicit_radius = 0.00, .covalent_radius = 0.30, .flags = HBD };

    // Carbon variants
    table[@intFromEnum(AtomType.C)] = .{ .explicit_radius = 1.70, .implicit_radius = 1.90, .covalent_radius = 0.77, .flags = NONE };
    table[@intFromEnum(AtomType.Car)] = .{ .explicit_radius = 1.75, .implicit_radius = 1.75, .covalent_radius = 0.77, .flags = ARA };
    table[@intFromEnum(AtomType.C_eq_O)] = .{ .explicit_radius = 1.65, .implicit_radius = 1.65, .covalent_radius = 0.80, .flags = NONE };

    // Nitrogen variants
    table[@intFromEnum(AtomType.N)] = .{ .explicit_radius = 1.55, .implicit_radius = 1.70, .covalent_radius = 0.70, .flags = NONE };
    table[@intFromEnum(AtomType.Nacc)] = .{ .explicit_radius = 1.55, .implicit_radius = 1.70, .covalent_radius = 0.70, .flags = A };

    // Others
    table[@intFromEnum(AtomType.O)] = .{ .explicit_radius = 1.40, .implicit_radius = 1.50, .covalent_radius = 0.66, .flags = A };
    table[@intFromEnum(AtomType.P)] = .{ .explicit_radius = 1.80, .implicit_radius = 1.80, .covalent_radius = 1.80, .flags = NONE };
    table[@intFromEnum(AtomType.S)] = .{ .explicit_radius = 1.80, .implicit_radius = 1.90, .covalent_radius = 1.04, .flags = A };
    table[@intFromEnum(AtomType.Se)] = .{ .explicit_radius = 1.90, .implicit_radius = 1.90, .covalent_radius = 1.17, .flags = A };
    table[@intFromEnum(AtomType.F)] = .{ .explicit_radius = 1.30, .implicit_radius = 1.30, .covalent_radius = 0.58, .flags = A };
    table[@intFromEnum(AtomType.Cl)] = .{ .explicit_radius = 1.77, .implicit_radius = 1.77, .covalent_radius = 0.99, .flags = A };
    table[@intFromEnum(AtomType.Br)] = .{ .explicit_radius = 1.95, .implicit_radius = 1.95, .covalent_radius = 1.14, .flags = A };
    table[@intFromEnum(AtomType.I)] = .{ .explicit_radius = 2.10, .implicit_radius = 2.10, .covalent_radius = 1.33, .flags = A };

    // Metals
    table[@intFromEnum(AtomType.Li)] = .{ .explicit_radius = 1.82, .implicit_radius = 1.82, .covalent_radius = 1.23, .flags = M };
    table[@intFromEnum(AtomType.Na)] = .{ .explicit_radius = 2.27, .implicit_radius = 2.27, .covalent_radius = 1.54, .flags = M };
    table[@intFromEnum(AtomType.Mg)] = .{ .explicit_radius = 1.73, .implicit_radius = 1.73, .covalent_radius = 1.36, .flags = M };
    table[@intFromEnum(AtomType.K)] = .{ .explicit_radius = 2.75, .implicit_radius = 2.75, .covalent_radius = 1.96, .flags = M };
    table[@intFromEnum(AtomType.Ca)] = .{ .explicit_radius = 2.31, .implicit_radius = 2.31, .covalent_radius = 1.74, .flags = M };
    table[@intFromEnum(AtomType.Mn)] = .{ .explicit_radius = 1.73, .implicit_radius = 1.73, .covalent_radius = 1.17, .flags = M };
    table[@intFromEnum(AtomType.Fe)] = .{ .explicit_radius = 1.65, .implicit_radius = 1.65, .covalent_radius = 1.17, .flags = M };
    table[@intFromEnum(AtomType.Co)] = .{ .explicit_radius = 1.67, .implicit_radius = 1.67, .covalent_radius = 1.16, .flags = M };
    table[@intFromEnum(AtomType.Ni)] = .{ .explicit_radius = 1.50, .implicit_radius = 1.50, .covalent_radius = 1.15, .flags = M };
    table[@intFromEnum(AtomType.Cu)] = .{ .explicit_radius = 1.52, .implicit_radius = 1.52, .covalent_radius = 1.17, .flags = M };
    table[@intFromEnum(AtomType.Zn)] = .{ .explicit_radius = 1.65, .implicit_radius = 1.65, .covalent_radius = 1.25, .flags = M };
    table[@intFromEnum(AtomType.As)] = .{ .explicit_radius = 1.85, .implicit_radius = 1.85, .covalent_radius = 1.21, .flags = NONE };
    table[@intFromEnum(AtomType.Rb)] = .{ .explicit_radius = 2.75, .implicit_radius = 2.75, .covalent_radius = 2.11, .flags = M };
    table[@intFromEnum(AtomType.Sr)] = .{ .explicit_radius = 2.49, .implicit_radius = 2.49, .covalent_radius = 1.92, .flags = M };
    table[@intFromEnum(AtomType.Mo)] = .{ .explicit_radius = 1.90, .implicit_radius = 1.90, .covalent_radius = 1.30, .flags = M };
    table[@intFromEnum(AtomType.Ag)] = .{ .explicit_radius = 1.72, .implicit_radius = 1.72, .covalent_radius = 1.34, .flags = M };
    table[@intFromEnum(AtomType.Cd)] = .{ .explicit_radius = 1.58, .implicit_radius = 1.58, .covalent_radius = 1.48, .flags = M };
    table[@intFromEnum(AtomType.Sn)] = .{ .explicit_radius = 2.17, .implicit_radius = 2.17, .covalent_radius = 1.40, .flags = M };
    table[@intFromEnum(AtomType.Cs)] = .{ .explicit_radius = 3.01, .implicit_radius = 3.01, .covalent_radius = 2.25, .flags = M };
    table[@intFromEnum(AtomType.Ba)] = .{ .explicit_radius = 2.68, .implicit_radius = 2.68, .covalent_radius = 1.98, .flags = M };
    table[@intFromEnum(AtomType.W)] = .{ .explicit_radius = 1.80, .implicit_radius = 1.80, .covalent_radius = 1.30, .flags = M };
    table[@intFromEnum(AtomType.Pt)] = .{ .explicit_radius = 1.75, .implicit_radius = 1.75, .covalent_radius = 1.28, .flags = M };
    table[@intFromEnum(AtomType.Au)] = .{ .explicit_radius = 1.66, .implicit_radius = 1.66, .covalent_radius = 1.34, .flags = M };
    table[@intFromEnum(AtomType.Hg)] = .{ .explicit_radius = 1.55, .implicit_radius = 1.55, .covalent_radius = 1.49, .flags = M };
    table[@intFromEnum(AtomType.Pb)] = .{ .explicit_radius = 2.02, .implicit_radius = 2.02, .covalent_radius = 1.47, .flags = M };
    table[@intFromEnum(AtomType.U)] = .{ .explicit_radius = 1.86, .implicit_radius = 1.86, .covalent_radius = 1.42, .flags = M };

    // Unknown
    table[@intFromEnum(AtomType.unknown)] = .{ .explicit_radius = 1.80, .implicit_radius = 1.80, .covalent_radius = 1.00, .flags = NONE };

    break :build table;
};

/// Normalize a symbol to 2-char: uppercase first letter, lowercase second letter (or space).
fn normalizeSymbol(symbol: []const u8) [2]u8 {
    var buf: [2]u8 = .{ ' ', ' ' };
    if (symbol.len == 0) return buf;
    buf[0] = std.ascii.toUpper(symbol[0]);
    if (symbol.len >= 2) {
        buf[1] = std.ascii.toLower(symbol[1]);
    }
    return buf;
}

const symbol_map = std.StaticStringMap(AtomType).initComptime(.{
    .{ "H ", .H },
    .{ "D ", .H },
    .{ "C ", .C },
    .{ "N ", .N },
    .{ "O ", .O },
    .{ "P ", .P },
    .{ "S ", .S },
    .{ "Se", .Se },
    .{ "F ", .F },
    .{ "Cl", .Cl },
    .{ "Br", .Br },
    .{ "I ", .I },
    .{ "Li", .Li },
    .{ "Na", .Na },
    .{ "Mg", .Mg },
    .{ "K ", .K },
    .{ "Ca", .Ca },
    .{ "Mn", .Mn },
    .{ "Fe", .Fe },
    .{ "Co", .Co },
    .{ "Ni", .Ni },
    .{ "Cu", .Cu },
    .{ "Zn", .Zn },
    .{ "As", .As },
    .{ "Rb", .Rb },
    .{ "Sr", .Sr },
    .{ "Mo", .Mo },
    .{ "Ag", .Ag },
    .{ "Cd", .Cd },
    .{ "Sn", .Sn },
    .{ "Cs", .Cs },
    .{ "Ba", .Ba },
    .{ "W ", .W },
    .{ "Pt", .Pt },
    .{ "Au", .Au },
    .{ "Hg", .Hg },
    .{ "Pb", .Pb },
    .{ "U ", .U },
});

pub fn elementFromSymbol(symbol: []const u8) AtomType {
    const normalized = normalizeSymbol(symbol);
    return symbol_map.get(&normalized) orelse .unknown;
}

test "hydrogen VDW radius" {
    const h_info = AtomType.H.info();
    try std.testing.expectApproxEqAbs(1.22, h_info.explicit_radius, 1e-6);
}

test "polar hydrogen is donor" {
    const hpol_info = AtomType.Hpol.info();
    try std.testing.expect(hpol_info.flags.donor == true);
    try std.testing.expectApproxEqAbs(1.05, hpol_info.explicit_radius, 1e-6);
}

test "oxygen is acceptor" {
    const o_info = AtomType.O.info();
    try std.testing.expect(o_info.flags.acceptor == true);
    try std.testing.expectApproxEqAbs(1.40, o_info.explicit_radius, 1e-6);
}

test "elementFromSymbol" {
    try std.testing.expectEqual(AtomType.C, elementFromSymbol("C"));
    try std.testing.expectEqual(AtomType.Fe, elementFromSymbol("FE"));
    try std.testing.expectEqual(AtomType.Fe, elementFromSymbol("Fe"));
    try std.testing.expectEqual(AtomType.unknown, elementFromSymbol("Xx"));
}

test "mergeFlags combines flags via OR" {
    const a = AtomFlags{ .donor = true };
    const b = AtomFlags{ .acceptor = true, .negative = true };
    const merged = mergeFlags(a, b);
    try std.testing.expect(merged.donor);
    try std.testing.expect(merged.acceptor);
    try std.testing.expect(merged.negative);
    try std.testing.expect(!merged.aromatic);
}

test "bonded_inter_residue flag" {
    var flags = AtomFlags{};
    try std.testing.expect(!flags.bonded_inter_residue);
    flags.bonded_inter_residue = true;
    try std.testing.expect(flags.bonded_inter_residue);

    // mergeFlags preserves bonded_inter_residue
    const a = AtomFlags{ .donor = true };
    const b = AtomFlags{ .bonded_inter_residue = true };
    const merged = mergeFlags(a, b);
    try std.testing.expect(merged.donor);
    try std.testing.expect(merged.bonded_inter_residue);
}
