//! Dot sphere generation for probe-based scoring.
//! Matches the original reduce algorithm (DotSph.cpp).
//! Dots are placed on a sphere surface using concentric rings
//! with alternating half-dot offset (staggered grid).

const std = @import("std");
const math_mod = @import("../math.zig");
const Vec3 = math_mod.Vec3;

pub const DotSphere = struct {
    points: []Vec3(f32),
    density: f32,
    radius: f32,
    allocator: std.mem.Allocator,

    /// Generate a dot sphere matching original reduce's concentric ring algorithm.
    /// density: dots per square Angstrom (default 16.0)
    /// radius: sphere radius in Angstroms
    pub fn generate(allocator: std.mem.Allocator, radius: f32, density: f32) !DotSphere {
        const pi = std.math.pi;

        // Number of latitude rings: n_rings = π * radius * sqrt(density)
        const sqrt_density = @sqrt(density);
        const n_rings_f = pi * radius * sqrt_density;
        const n_rings: u32 = @max(1, @as(u32, @intFromFloat(@round(n_rings_f))));

        // First pass: count total dots to allocate exactly
        var total: u32 = 0;
        const delta_theta = pi / @as(f32, @floatFromInt(n_rings));
        for (0..n_rings) |ring_idx| {
            const theta = -pi / 2.0 + (@as(f32, @floatFromInt(ring_idx)) + 0.5) * delta_theta;
            const ring_r = radius * @cos(theta);
            const n_dots: u32 = @max(1, @as(u32, @intFromFloat(@round(2.0 * pi * ring_r * sqrt_density))));
            total += n_dots;
        }

        var points = try allocator.alloc(Vec3(f32), total);
        errdefer allocator.free(points);

        // Second pass: fill points
        var idx: u32 = 0;
        for (0..n_rings) |ring_idx| {
            const theta = -pi / 2.0 + (@as(f32, @floatFromInt(ring_idx)) + 0.5) * delta_theta;
            const ring_r = radius * @cos(theta);
            const z = radius * @sin(theta);
            const n_dots: u32 = @max(1, @as(u32, @intFromFloat(@round(2.0 * pi * ring_r * sqrt_density))));
            const delta_phi = 2.0 * pi / @as(f32, @floatFromInt(n_dots));
            // Alternate rings are offset by half a dot spacing (staggered grid)
            const offset: f32 = if (ring_idx % 2 == 1) delta_phi / 2.0 else 0.0;

            for (0..n_dots) |dot_idx| {
                const phi = @as(f32, @floatFromInt(dot_idx)) * delta_phi + offset;
                points[idx] = .{
                    .x = ring_r * @cos(phi),
                    .y = ring_r * @sin(phi),
                    .z = z,
                };
                idx += 1;
            }
        }

        return DotSphere{
            .points = points,
            .density = density,
            .radius = radius,
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *DotSphere) void {
        self.allocator.free(self.points);
    }
};

// ── Tests ────────────────────────────────────────────────────────────────────

const testing = std.testing;

test "dot sphere count at density 16" {
    var sphere = try DotSphere.generate(testing.allocator, 1.40, 16.0);
    defer sphere.deinit();
    // Expected: ~4π * 16 * 1.4² ≈ 394 dots
    try testing.expect(sphere.points.len > 350);
    try testing.expect(sphere.points.len < 450);
}

test "dot sphere all points at correct radius" {
    var sphere = try DotSphere.generate(testing.allocator, 1.70, 16.0);
    defer sphere.deinit();
    for (sphere.points) |p| {
        const dist = @sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        try testing.expectApproxEqAbs(dist, 1.70, 0.01);
    }
}

test "dot sphere uniform distribution" {
    // Check that dots are roughly uniformly distributed across octants
    var sphere = try DotSphere.generate(testing.allocator, 1.0, 16.0);
    defer sphere.deinit();
    // Count dots in each octant
    var octants: [8]u32 = @splat(0);
    for (sphere.points) |p| {
        const idx: u3 = @intCast(
            (@as(u3, if (p.x > 0) 1 else 0)) |
                (@as(u3, if (p.y > 0) 1 else 0) << 1) |
                (@as(u3, if (p.z > 0) 1 else 0) << 2),
        );
        octants[idx] += 1;
    }
    // Each octant should have roughly 1/8 of the dots
    const total: f32 = @floatFromInt(sphere.points.len);
    for (octants) |count| {
        const frac = @as(f32, @floatFromInt(count)) / total;
        try testing.expect(frac > 0.05); // at least 5% in each octant
        try testing.expect(frac < 0.25); // at most 25% in each octant
    }
}
