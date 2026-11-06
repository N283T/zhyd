//! Internal hydrogens of the empirical pKa engine ("Internal hydrogens").
//!
//! These hydrogens exist only for the energy terms; they are never written
//! out and do not depend on zreduce's hydrogen placement.

const std = @import("std");
const math = @import("../math.zig");

pub const Vec = math.Vec3(f64);

/// Hydrogen on N with two heavy neighbours A and B, along -(â + b̂).
pub fn bisector(n: Vec, a: Vec, b: Vec, length: f64) Vec {
    const ua = a.sub(n).normalize();
    const ub = b.sub(n).normalize();
    const dir = ua.add(ub).negate().normalize();
    return n.add(dir.scale(length));
}

/// Two hydrogens on N bonded to one trigonal atom X, each 120° from the N→X
/// direction, in the plane through N and X that contains X's other
/// neighbours. When X has two other neighbours Y₁ and Y₂ that are not
/// coplanar with N and X, the in-plane direction is the unit component of
/// Y₂→Y₁ perpendicular to N→X (proposed amendment A4). The first hydrogen
/// is on Y₁'s side (cis to Y₁), the second on the opposite side.
pub fn planarPair(n: Vec, x: Vec, y1: Vec, y2: ?Vec, length: f64) [2]Vec {
    const u = x.sub(n).normalize();
    const v = perpendicular(if (y2) |second| y1.sub(second) else y1.sub(x), u);
    const along = u.scale(-0.5);
    const across = v.scale(@sqrt(3.0) / 2.0);
    return .{
        n.add(along.add(across).scale(length)),
        n.add(along.sub(across).scale(length)),
    };
}

/// Unit component of `w` perpendicular to the unit vector `u`.
fn perpendicular(w: Vec, u: Vec) Vec {
    return w.sub(u.scale(w.dot(u))).normalize();
}

const testing = std.testing;

test "bisector hydrogen points away from both neighbours" {
    // N at origin, neighbours along +x and +y: H along -(1,1)/sqrt2.
    const h = bisector(.{ .x = 0, .y = 0, .z = 0 }, .{ .x = 1.5, .y = 0, .z = 0 }, .{ .x = 0, .y = 1.3, .z = 0 }, 1.01);
    const s = 1.01 / @sqrt(2.0);
    try testing.expectApproxEqAbs(-s, h.x, 1e-12);
    try testing.expectApproxEqAbs(-s, h.y, 1e-12);
    try testing.expectApproxEqAbs(@as(f64, 0), h.z, 1e-12);
}

test "planar pair hydrogens sit at 120 degrees in the plane" {
    // N at origin, X on +x, Y above X in +y: H at 1.01*(-0.5, ±sqrt3/2, 0).
    const hs = planarPair(.{ .x = 0, .y = 0, .z = 0 }, .{ .x = 1.33, .y = 0, .z = 0 }, .{ .x = 2.0, .y = 1.2, .z = 0 }, null, 1.01);
    try testing.expectApproxEqAbs(@as(f64, -0.505), hs[0].x, 1e-12);
    try testing.expectApproxEqAbs(1.01 * @sqrt(3.0) / 2.0, hs[0].y, 1e-12);
    try testing.expectApproxEqAbs(@as(f64, -0.505), hs[1].x, 1e-12);
    try testing.expectApproxEqAbs(-1.01 * @sqrt(3.0) / 2.0, hs[1].y, 1e-12);
    try testing.expectApproxEqAbs(@as(f64, 0), hs[1].z, 1e-12);
}

test "planar pair with two non-coplanar neighbours uses the Y2→Y1 direction" {
    // N at origin, X on +x. Y1 = (2, 1.2, 0), Y2 = (2, 0, 1.2):
    // Y1 − Y2 = (0, 1.2, −1.2), already perpendicular to x, so
    // v = (0, 1, −1)/√2 and H = (−0.5, ±(√3/2)/√2, ∓(√3/2)/√2).
    const hs = planarPair(.{ .x = 0, .y = 0, .z = 0 }, .{ .x = 1.33, .y = 0, .z = 0 }, .{ .x = 2.0, .y = 1.2, .z = 0 }, .{ .x = 2.0, .y = 0, .z = 1.2 }, 1.0);
    const c = @sqrt(3.0) / 2.0 / @sqrt(2.0);
    try testing.expectApproxEqAbs(@as(f64, -0.5), hs[0].x, 1e-12);
    try testing.expectApproxEqAbs(c, hs[0].y, 1e-12);
    try testing.expectApproxEqAbs(-c, hs[0].z, 1e-12);
    try testing.expectApproxEqAbs(-c, hs[1].y, 1e-12);
    try testing.expectApproxEqAbs(c, hs[1].z, 1e-12);
}
