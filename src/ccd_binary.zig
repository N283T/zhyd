//! Binary format for CCD (Chemical Component Dictionary).
//! Provides fast load/save of ComponentDict without CIF text parsing.
//!
//! Format:
//!   Header: MAGIC(4) + version(1) + reserved(3) + component_count(u32 LE)
//!   Per component:
//!     comp_id_len(u8) + comp_id(bytes)
//!     comp_type_len(u8) + comp_type(bytes)
//!     atom_count(u16 LE) + PackedAtom * atom_count
//!     bond_count(u16 LE) + PackedBond * bond_count
//!
//! PackedAtom/PackedBond are written as raw struct bytes assuming little-endian
//! layout. A comptime assertion rejects big-endian targets.

const std = @import("std");
const builtin = @import("builtin");

comptime {
    // Binary format uses native byte order for packed structs; only LE is supported.
    std.debug.assert(builtin.target.cpu.arch.endian() == .little);
}
const Allocator = std.mem.Allocator;
const ccd = @import("ccd.zig");

pub const ComponentDict = ccd.ComponentDict;
pub const Component = ccd.Component;
pub const CompAtom = ccd.CompAtom;
pub const CompBond = ccd.CompBond;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

pub const MAGIC: [4]u8 = .{ 'Z', 'R', 'D', 'C' };
pub const FORMAT_VERSION: u8 = 1;
pub const HEADER_SIZE: usize = 12; // 4 (magic) + 1 (version) + 3 (reserved) + 4 (count)

// ---------------------------------------------------------------------------
// isBinaryDict
// ---------------------------------------------------------------------------

/// Returns true if data begins with the binary dict magic bytes and is at
/// least HEADER_SIZE bytes long.
pub fn isBinaryDict(data: []const u8) bool {
    if (data.len < HEADER_SIZE) return false;
    return std.mem.eql(u8, data[0..4], &MAGIC);
}

// ---------------------------------------------------------------------------
// PackedAtom  (extern struct, 24 bytes)
// ---------------------------------------------------------------------------

pub const PackedAtom = extern struct {
    name: [4]u8,
    name_len: u8, // stores u4 value; mask with 0x0F on read
    element_symbol: [2]u8,
    charge: i8,
    /// bit0 = leaving, bit1 = aromatic
    flags: u8,
    _pad: [3]u8,
    ideal_x: f32,
    ideal_y: f32,
    ideal_z: f32,

    comptime {
        std.debug.assert(@sizeOf(PackedAtom) == 24);
    }

    pub fn fromCompAtom(a: CompAtom) PackedAtom {
        var flags: u8 = 0;
        if (a.leaving) flags |= 0x01;
        if (a.aromatic) flags |= 0x02;
        return .{
            .name = a.name,
            .name_len = @as(u8, a.name_len),
            .element_symbol = a.element_symbol,
            .charge = a.charge,
            .flags = flags,
            ._pad = .{ 0, 0, 0 },
            .ideal_x = a.ideal_x,
            .ideal_y = a.ideal_y,
            .ideal_z = a.ideal_z,
        };
    }

    pub fn toCompAtom(self: PackedAtom) CompAtom {
        return .{
            .name = self.name,
            .name_len = @truncate(self.name_len & 0x0F),
            .element_symbol = self.element_symbol,
            .charge = self.charge,
            .leaving = (self.flags & 0x01) != 0,
            .aromatic = (self.flags & 0x02) != 0,
            .ideal_x = self.ideal_x,
            .ideal_y = self.ideal_y,
            .ideal_z = self.ideal_z,
        };
    }
};

// ---------------------------------------------------------------------------
// PackedBond  (extern struct, 6 bytes)
// ---------------------------------------------------------------------------

pub const PackedBond = extern struct {
    atom_idx_1: u16,
    atom_idx_2: u16,
    order: u8, // BondOrder tag value
    /// bit0 = aromatic
    flags: u8,

    comptime {
        std.debug.assert(@sizeOf(PackedBond) == 6);
    }

    pub fn fromCompBond(b: CompBond) PackedBond {
        var flags: u8 = 0;
        if (b.aromatic) flags |= 0x01;
        return .{
            .atom_idx_1 = b.atom_idx_1,
            .atom_idx_2 = b.atom_idx_2,
            .order = @backingInt(b.order),
            .flags = flags,
        };
    }

    pub fn toCompBond(self: PackedBond) CompBond {
        return .{
            .atom_idx_1 = self.atom_idx_1,
            .atom_idx_2 = self.atom_idx_2,
            .order = if (self.order <= @backingInt(ccd.BondOrder.unknown))
                @fromBackingInt(@intCast(self.order))
            else
                .unknown,
            .aromatic = (self.flags & 0x01) != 0,
        };
    }
};

// ---------------------------------------------------------------------------
// writeDict
// ---------------------------------------------------------------------------

/// Serialize dict to writer in binary format.
pub fn writeDict(writer: anytype, dict: *const ComponentDict) !void {
    // Header
    try writer.writeAll(&MAGIC);
    try writer.writeByte(FORMAT_VERSION);
    try writer.writeAll(&[_]u8{ 0, 0, 0 }); // reserved
    const count: u32 = @intCast(dict.components.count());
    try writer.writeInt(u32, count, .little);

    // Components
    var iter = dict.components.iterator();
    while (iter.next()) |entry| {
        const comp = entry.value_ptr.*;

        // comp_id
        const id_len: u8 = @intCast(comp.comp_id.len);
        try writer.writeByte(id_len);
        try writer.writeAll(comp.comp_id);

        // comp_type
        const type_len: u8 = @intCast(comp.comp_type.len);
        try writer.writeByte(type_len);
        try writer.writeAll(comp.comp_type);

        // atoms
        const atom_count: u16 = @intCast(comp.atoms.len);
        try writer.writeInt(u16, atom_count, .little);
        for (comp.atoms) |atom| {
            const pa = PackedAtom.fromCompAtom(atom);
            try writer.writeAll(std.mem.asBytes(&pa));
        }

        // bonds
        const bond_count: u16 = @intCast(comp.bonds.len);
        try writer.writeInt(u16, bond_count, .little);
        for (comp.bonds) |bond| {
            const pb = PackedBond.fromCompBond(bond);
            try writer.writeAll(std.mem.asBytes(&pb));
        }
    }
}

// ---------------------------------------------------------------------------
// Tests: isBinaryDict
// ---------------------------------------------------------------------------

const testing = std.testing;

test "isBinaryDict: valid magic" {
    var header: [HEADER_SIZE]u8 = undefined;
    @memcpy(header[0..4], &MAGIC);
    header[4] = FORMAT_VERSION;
    @memset(header[5..], 0);
    try testing.expect(isBinaryDict(&header));
}

test "isBinaryDict: wrong magic" {
    var header: [HEADER_SIZE]u8 = undefined;
    @memcpy(header[0..4], "XXXX");
    @memset(header[4..], 0);
    try testing.expect(!isBinaryDict(&header));
}

test "isBinaryDict: empty data" {
    try testing.expect(!isBinaryDict(&[_]u8{}));
}

test "isBinaryDict: too short" {
    try testing.expect(!isBinaryDict(&[_]u8{ 'Z', 'R' }));
}

// ---------------------------------------------------------------------------
// readDict
// ---------------------------------------------------------------------------

pub const ReadError = error{
    InvalidMagic,
    UnsupportedVersion,
    UnexpectedEof,
    OutOfMemory,
    CountTooLarge,
};

/// Deserialize a ComponentDict from reader. Caller owns the returned dict.
pub fn readDict(allocator: Allocator, reader: anytype) ReadError!ComponentDict {
    // Read header
    var magic: [4]u8 = undefined;
    reader.readSliceAll(&magic) catch return error.UnexpectedEof;
    if (!std.mem.eql(u8, &magic, &MAGIC)) return error.InvalidMagic;

    const version = reader.takeByte() catch return error.UnexpectedEof;
    if (version != FORMAT_VERSION) return error.UnsupportedVersion;

    var reserved: [3]u8 = undefined;
    reader.readSliceAll(&reserved) catch return error.UnexpectedEof;

    const count = reader.takeInt(u32, .little) catch return error.UnexpectedEof;
    // Sanity check: real CCD has ~40,000 components; 500,000 is a generous upper bound.
    if (count > 500_000) return error.CountTooLarge;

    var dict = ComponentDict{
        .components = std.StringHashMap(Component).init(allocator),
        .allocator = allocator,
    };
    errdefer dict.deinit();

    try dict.components.ensureTotalCapacity(count);

    var i: u32 = 0;
    while (i < count) : (i += 1) {
        // comp_id
        const id_len = reader.takeByte() catch return error.UnexpectedEof;
        const comp_id = allocator.alloc(u8, id_len) catch return error.OutOfMemory;
        errdefer allocator.free(comp_id);
        reader.readSliceAll(comp_id) catch return error.UnexpectedEof;

        // comp_type
        const type_len = reader.takeByte() catch return error.UnexpectedEof;
        const comp_type = allocator.alloc(u8, type_len) catch return error.OutOfMemory;
        errdefer allocator.free(comp_type);
        reader.readSliceAll(comp_type) catch return error.UnexpectedEof;

        // atoms
        const atom_count = reader.takeInt(u16, .little) catch return error.UnexpectedEof;
        // Sanity check: no real CCD component has more than 1000 atoms.
        if (atom_count > 1000) return error.CountTooLarge;
        const atoms = allocator.alloc(CompAtom, atom_count) catch return error.OutOfMemory;
        errdefer allocator.free(atoms);
        for (atoms) |*atom| {
            var pa: PackedAtom = undefined;
            reader.readSliceAll(std.mem.asBytes(&pa)) catch return error.UnexpectedEof;
            atom.* = pa.toCompAtom();
        }

        // bonds
        const bond_count = reader.takeInt(u16, .little) catch return error.UnexpectedEof;
        // Sanity check: no real CCD component has more than 1000 bonds.
        if (bond_count > 1000) return error.CountTooLarge;
        const bonds = allocator.alloc(CompBond, bond_count) catch return error.OutOfMemory;
        errdefer allocator.free(bonds);
        for (bonds) |*bond| {
            var pb: PackedBond = undefined;
            reader.readSliceAll(std.mem.asBytes(&pb)) catch return error.UnexpectedEof;
            bond.* = pb.toCompBond();
        }

        // comp_id serves as both HashMap key and Component.comp_id (shared allocation).
        dict.components.putAssumeCapacity(comp_id, .{
            .comp_id = comp_id,
            .comp_type = comp_type,
            .atoms = atoms,
            .bonds = bonds,
        });
    }

    return dict;
}

// ---------------------------------------------------------------------------
// loadDict — auto-detect text vs binary
// ---------------------------------------------------------------------------

/// Load a ComponentDict from raw bytes. Detects binary vs CIF text automatically.
pub fn loadDict(allocator: Allocator, data: []const u8) !ComponentDict {
    if (isBinaryDict(data)) {
        var reader: std.Io.Reader = .fixed(data);
        return readDict(allocator, &reader);
    }
    return ccd.parseComponentDict(allocator, data);
}

// ---------------------------------------------------------------------------
// Tests: writeDict, readDict, loadDict
// ---------------------------------------------------------------------------

test "writeDict + readDict: round-trip" {
    const allocator = testing.allocator;

    // Build a dict manually
    var dict = ComponentDict{
        .components = std.StringHashMap(Component).init(allocator),
        .allocator = allocator,
    };
    defer dict.deinit();

    // Atoms
    const atoms = try allocator.dupe(CompAtom, &[_]CompAtom{
        .{
            .name = .{ 'C', '1', ' ', ' ' },
            .name_len = 2,
            .element_symbol = .{ 'C', ' ' },
            .charge = 0,
            .leaving = false,
            .aromatic = true,
            .ideal_x = 1.0,
            .ideal_y = 2.0,
            .ideal_z = 3.0,
        },
        .{
            .name = .{ 'O', '1', ' ', ' ' },
            .name_len = 2,
            .element_symbol = .{ 'O', ' ' },
            .charge = -1,
            .leaving = true,
            .aromatic = false,
            .ideal_x = 4.5,
            .ideal_y = 5.5,
            .ideal_z = 6.5,
        },
    });

    // Bonds
    const bonds = try allocator.dupe(CompBond, &[_]CompBond{
        .{
            .atom_idx_1 = 0,
            .atom_idx_2 = 1,
            .order = .double,
            .aromatic = true,
        },
    });

    // Single allocation shared as both HashMap key and Component.comp_id.
    const key = try allocator.dupe(u8, "TST");
    const comp_type = try allocator.dupe(u8, "NON-POLYMER");
    try dict.components.put(key, .{
        .comp_id = key,
        .comp_type = comp_type,
        .atoms = atoms,
        .bonds = bonds,
    });

    // Write to buffer
    var buf: std.Io.Writer.Allocating = .init(allocator);
    defer buf.deinit();
    try writeDict(&buf.writer, &dict);

    // Read back
    var reader: std.Io.Reader = .fixed(buf.writer.buffered());
    var dict2 = try readDict(allocator, &reader);
    defer dict2.deinit();

    const comp = dict2.get("TST") orelse return error.TestUnexpectedResult;

    try testing.expectEqualStrings("TST", comp.comp_id);
    try testing.expectEqualStrings("NON-POLYMER", comp.comp_type);
    try testing.expectEqual(@as(usize, 2), comp.atoms.len);
    try testing.expectEqual(@as(usize, 1), comp.bonds.len);

    // Atom 0
    const a0 = comp.atoms[0];
    try testing.expectEqualSlices(u8, "C1", a0.nameSlice());
    try testing.expectEqual(@as([2]u8, .{ 'C', ' ' }), a0.element_symbol);
    try testing.expectEqual(@as(i8, 0), a0.charge);
    try testing.expect(!a0.leaving);
    try testing.expect(a0.aromatic);
    try testing.expectApproxEqAbs(@as(f32, 1.0), a0.ideal_x, 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 2.0), a0.ideal_y, 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 3.0), a0.ideal_z, 1e-6);

    // Atom 1
    const a1 = comp.atoms[1];
    try testing.expectEqualSlices(u8, "O1", a1.nameSlice());
    try testing.expectEqual(@as([2]u8, .{ 'O', ' ' }), a1.element_symbol);
    try testing.expectEqual(@as(i8, -1), a1.charge);
    try testing.expect(a1.leaving);
    try testing.expect(!a1.aromatic);
    try testing.expectApproxEqAbs(@as(f32, 4.5), a1.ideal_x, 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 5.5), a1.ideal_y, 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 6.5), a1.ideal_z, 1e-6);

    // Bond 0
    const b0 = comp.bonds[0];
    try testing.expectEqual(@as(u16, 0), b0.atom_idx_1);
    try testing.expectEqual(@as(u16, 1), b0.atom_idx_2);
    try testing.expectEqual(ccd.BondOrder.double, b0.order);
    try testing.expect(b0.aromatic);
}

test "writeDict: empty dictionary" {
    const allocator = testing.allocator;

    var dict = ComponentDict{
        .components = std.StringHashMap(Component).init(allocator),
        .allocator = allocator,
    };
    defer dict.deinit();

    var buf: std.Io.Writer.Allocating = .init(allocator);
    defer buf.deinit();
    try writeDict(&buf.writer, &dict);

    try testing.expectEqual(HEADER_SIZE, buf.writer.buffered().len);
    // Verify component count field is 0
    const count = std.mem.readInt(u32, buf.writer.buffered()[8..12], .little);
    try testing.expectEqual(@as(u32, 0), count);
}

test "readDict: version mismatch" {
    var buf: [HEADER_SIZE]u8 = undefined;
    @memcpy(buf[0..4], &MAGIC);
    buf[4] = 99; // wrong version
    @memset(buf[5..8], 0);
    std.mem.writeInt(u32, buf[8..12], 0, .little);

    var reader: std.Io.Reader = .fixed(&buf);
    const result = readDict(testing.allocator, &reader);
    try testing.expectError(error.UnsupportedVersion, result);
}

test "readDict: invalid magic" {
    var buf: [HEADER_SIZE]u8 = undefined;
    @memcpy(buf[0..4], "XXXX");
    buf[4] = FORMAT_VERSION;
    @memset(buf[5..8], 0);
    std.mem.writeInt(u32, buf[8..12], 0, .little);

    var reader: std.Io.Reader = .fixed(&buf);
    const result = readDict(testing.allocator, &reader);
    try testing.expectError(error.InvalidMagic, result);
}

test "readDict: truncated data" {
    // Header claims 1 component but provides no component data
    var buf: [HEADER_SIZE]u8 = undefined;
    @memcpy(buf[0..4], &MAGIC);
    buf[4] = FORMAT_VERSION;
    @memset(buf[5..8], 0);
    std.mem.writeInt(u32, buf[8..12], 1, .little);

    var reader: std.Io.Reader = .fixed(&buf);
    const result = readDict(testing.allocator, &reader);
    try testing.expectError(error.UnexpectedEof, result);
}

test "loadDict: CIF text input" {
    const allocator = testing.allocator;

    const cif_text =
        \\data_TST
        \\_chem_comp.id TST
        \\_chem_comp.type "NON-POLYMER"
        \\loop_
        \\_chem_comp_atom.comp_id
        \\_chem_comp_atom.atom_id
        \\_chem_comp_atom.type_symbol
        \\_chem_comp_atom.charge
        \\_chem_comp_atom.pdbx_leaving_atom_flag
        \\_chem_comp_atom.model_Cartn_x_ideal
        \\_chem_comp_atom.model_Cartn_y_ideal
        \\_chem_comp_atom.model_Cartn_z_ideal
        \\TST C1 C 0 N 1.0 2.0 3.0
        \\
    ;

    var dict = try loadDict(allocator, cif_text);
    defer dict.deinit();

    const comp = dict.get("TST") orelse return error.TestUnexpectedResult;
    try testing.expectEqualStrings("TST", comp.comp_id);
    try testing.expectEqual(@as(usize, 1), comp.atoms.len);
}

test "loadDict: binary input round-trip" {
    const allocator = testing.allocator;

    // Build a minimal dict
    var dict = ComponentDict{
        .components = std.StringHashMap(Component).init(allocator),
        .allocator = allocator,
    };
    defer dict.deinit();

    const atoms = try allocator.dupe(CompAtom, &[_]CompAtom{
        .{
            .name = .{ 'N', ' ', ' ', ' ' },
            .name_len = 1,
            .element_symbol = .{ 'N', ' ' },
            .charge = 0,
            .leaving = false,
            .aromatic = false,
            .ideal_x = 0.0,
            .ideal_y = 0.0,
            .ideal_z = 0.0,
        },
    });
    const bonds = try allocator.dupe(CompBond, &[_]CompBond{});
    // Single allocation shared as both HashMap key and Component.comp_id.
    const key = try allocator.dupe(u8, "ATP");
    const comp_type = try allocator.dupe(u8, "ATP");
    try dict.components.put(key, .{
        .comp_id = key,
        .comp_type = comp_type,
        .atoms = atoms,
        .bonds = bonds,
    });

    // Write binary
    var buf: std.Io.Writer.Allocating = .init(allocator);
    defer buf.deinit();
    try writeDict(&buf.writer, &dict);

    // loadDict should detect binary and parse it
    var dict2 = try loadDict(allocator, buf.writer.buffered());
    defer dict2.deinit();

    const comp = dict2.get("ATP") orelse return error.TestUnexpectedResult;
    try testing.expectEqualStrings("ATP", comp.comp_id);
    try testing.expectEqual(@as(usize, 1), comp.atoms.len);
    try testing.expectEqual(@as(usize, 0), comp.bonds.len);
}
