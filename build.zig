const std = @import("std");

const version = "0.3.0";

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const with_propka = b.option(
        bool,
        "propka",
        "Enable the PROPKA 3.5.1 pKa engine via zpropka (LGPL-2.1; default off)",
    ) orelse false;

    const options = b.addOptions();
    options.addOption([]const u8, "version", version);
    options.addOption(bool, "propka", with_propka);
    const options_module = options.createModule();

    // Library module
    const mod = b.addModule("zhyd", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .imports = &.{.{ .name = "build_options", .module = options_module }},
    });
    if (with_propka) {
        if (b.lazyDependency("zpropka", .{ .target = target, .optimize = optimize })) |dep| {
            mod.addImport("zpropka", dep.module("zpropka"));
        }
    }

    // CLI executable
    const exe_module = b.createModule(.{
        .root_source_file = b.path("src/main.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "zhyd", .module = mod },
            .{ .name = "build_options", .module = options_module },
        },
    });

    const exe = b.addExecutable(.{
        .name = "zhyd",
        .root_module = exe_module,
    });
    b.installArtifact(exe);

    // Run step
    const run_step = b.step("run", "Run zhyd");
    const run_cmd = b.addRunArtifact(exe);
    run_step.dependOn(&run_cmd.step);
    run_cmd.step.dependOn(b.getInstallStep());
    if (b.args) |args| {
        run_cmd.addArgs(args);
    }

    // Test step
    const mod_tests = b.addTest(.{ .root_module = mod });
    const exe_tests = b.addTest(.{ .root_module = exe.root_module });
    const test_step = b.step("test", "Run tests");
    test_step.dependOn(&b.addRunArtifact(mod_tests).step);
    test_step.dependOn(&b.addRunArtifact(exe_tests).step);
}
