const std = @import("std");
const zreduce = @import("zreduce");
const build_options = @import("build_options");

const Allocator = std.mem.Allocator;

fn defaultIo() std.Io {
    return std.Io.Threaded.global_single_threaded.io();
}

const RunConfig = struct {
    input_path: []const u8,
    output_path: ?[]const u8 = null,
    dict_path: ?[]const u8 = null,
    sdf_path: ?[]const u8 = null,
    json_path: ?[]const u8 = null,
    state_json_path: ?[]const u8 = null,
    no_opt: bool = false,
    no_flip: bool = false,
    validate: bool = false,
    water: zreduce.place.WaterConfig = .{},
    bond_policy: zreduce.place.BondPolicy = .{},
    nterm_mode: zreduce.place.NtermMode = .auto,
    protonation_path: ?[]const u8 = null,
    fix_path: ?[]const u8 = null,
    dump_movers_path: ?[]const u8 = null,
    strip_h: bool = false,
    model_filter: zreduce.run.ModelFilter = .all,
    show_timing: bool = false,
};

const PkaConfig = struct {
    input_path: []const u8,
    output_path: ?[]const u8 = null,
    model_filter: zreduce.run.ModelFilter = .all,
    model_only: bool = false,
    parameters_path: ?[]const u8 = null,
    engine: ?PkaEngine = null,
    propka_keep_waters: bool = false,
    empirical_mode: ?zreduce.pka_empirical.Mode = null,
    empirical_fixes: ?zreduce.pka_empirical.Fixes = null,
    show_timing: bool = false,
    // Directory input only:
    n_threads: ?u32 = null,
    jsonl_path: ?[]const u8 = null,
    quiet: bool = false,
};

const PkaEngine = zreduce.pka_run.Engine;

/// Fields shared between RunConfig and BatchConfig, parsed by parseCommonOption.
const CommonOptions = struct {
    no_opt: bool = false,
    no_flip: bool = false,
    water: zreduce.place.WaterConfig = .{},
    bond_policy: zreduce.place.BondPolicy = .{},
    nterm_mode: zreduce.place.NtermMode = .auto,
    protonation_path: ?[]const u8 = null,
    fix_path: ?[]const u8 = null,
    sdf_path: ?[]const u8 = null,
    strip_h: bool = false,
    model_filter: zreduce.run.ModelFilter = .all,
};

fn parseBondModeValue(s: []const u8) ?zreduce.place.BondLengthMode {
    if (std.ascii.eqlIgnoreCase(s, "neutron")) return .neutron;
    if (std.ascii.eqlIgnoreCase(s, "xray")) return .xray;
    return null;
}

fn parseOutputIsotopeValue(s: []const u8) ?zreduce.place.OutputIsotope {
    if (std.ascii.eqlIgnoreCase(s, "hydrogen") or std.ascii.eqlIgnoreCase(s, "h")) return .hydrogen;
    if (std.ascii.eqlIgnoreCase(s, "deuterium") or std.ascii.eqlIgnoreCase(s, "d")) return .deuterium;
    return null;
}

/// Attempts to parse one common option from args[i.*].
/// Returns true if the argument was consumed (recognized as a common option).
/// Advances i.* by one extra if the option takes a value argument.
fn parseCommonOption(args: []const []const u8, i: *usize, common: *CommonOptions) bool {
    const arg = args[i.*];
    if (std.mem.eql(u8, arg, "--no-opt")) {
        common.no_opt = true;
        return true;
    } else if (std.mem.eql(u8, arg, "--no-flip")) {
        common.no_flip = true;
        return true;
    } else if (std.mem.eql(u8, arg, "--water")) {
        common.water.enabled = true;
        return true;
    } else if (std.mem.eql(u8, arg, "--no-water")) {
        common.water.enabled = false;
        return true;
    } else if (std.mem.eql(u8, arg, "--water-phantom")) {
        common.water.enabled = true;
        common.water.phantom = true;
        return true;
    } else if (std.mem.eql(u8, arg, "--water-occ-cutoff")) {
        i.* += 1;
        if (i.* >= args.len) {
            std.debug.print("Error: --water-occ-cutoff requires a numeric argument\n", .{});
            std.process.exit(1);
        }
        const val = std.fmt.parseFloat(f32, args[i.*]) catch {
            std.debug.print("Error: invalid water occupancy cutoff '{s}'\n", .{args[i.*]});
            std.process.exit(1);
        };
        if (!std.math.isFinite(val) or val < 0.0 or val > 1.0) {
            std.debug.print("Error: --water-occ-cutoff must be between 0.0 and 1.0\n", .{});
            std.process.exit(1);
        }
        common.water.occupancy_cutoff = val;
        common.water.enabled = true;
        return true;
    } else if (std.mem.eql(u8, arg, "--water-b-cutoff")) {
        i.* += 1;
        if (i.* >= args.len) {
            std.debug.print("Error: --water-b-cutoff requires a numeric argument\n", .{});
            std.process.exit(1);
        }
        const val = std.fmt.parseFloat(f32, args[i.*]) catch {
            std.debug.print("Error: invalid water B-factor cutoff '{s}'\n", .{args[i.*]});
            std.process.exit(1);
        };
        if (!std.math.isFinite(val) or val < 0.0) {
            std.debug.print("Error: --water-b-cutoff must be a non-negative number\n", .{});
            std.process.exit(1);
        }
        common.water.b_factor_cutoff = val;
        common.water.enabled = true;
        return true;
    } else if (std.mem.eql(u8, arg, "--bond-mode")) {
        i.* += 1;
        if (i.* >= args.len) {
            std.debug.print("Error: --bond-mode requires 'neutron' or 'xray'\n", .{});
            std.process.exit(1);
        }
        const parsed = parseBondModeValue(args[i.*]) orelse {
            std.debug.print("Error: invalid --bond-mode '{s}' (expected 'neutron' or 'xray')\n", .{args[i.*]});
            std.process.exit(1);
        };
        common.bond_policy.mode = parsed;
        return true;
    } else if (std.mem.eql(u8, arg, "--nterm")) {
        i.* += 1;
        if (i.* >= args.len) {
            std.debug.print("Error: --nterm requires 'auto', 'aggressive', or 'neutral'\n", .{});
            std.process.exit(1);
        }
        const parsed = zreduce.place.NtermMode.fromString(args[i.*]) orelse {
            std.debug.print("Error: invalid --nterm '{s}' (expected auto|aggressive|neutral)\n", .{args[i.*]});
            std.process.exit(1);
        };
        common.nterm_mode = parsed;
        return true;
    } else if (std.mem.eql(u8, arg, "--isotope")) {
        i.* += 1;
        if (i.* >= args.len) {
            std.debug.print("Error: --isotope requires 'hydrogen'/'h' or 'deuterium'/'d'\n", .{});
            std.process.exit(1);
        }
        const parsed = parseOutputIsotopeValue(args[i.*]) orelse {
            std.debug.print("Error: invalid --isotope '{s}' (expected hydrogen|h|deuterium|d)\n", .{args[i.*]});
            std.process.exit(1);
        };
        common.bond_policy.output_isotope = parsed;
        return true;
    } else if (std.mem.eql(u8, arg, "--deuterium")) {
        common.bond_policy.output_isotope = .deuterium;
        return true;
    } else if (std.mem.eql(u8, arg, "--protonation")) {
        i.* += 1;
        if (i.* >= args.len) {
            std.debug.print("Error: --protonation requires a path argument\n", .{});
            std.process.exit(1);
        }
        common.protonation_path = args[i.*];
        return true;
    } else if (std.mem.eql(u8, arg, "--fix")) {
        i.* += 1;
        if (i.* >= args.len) {
            std.debug.print("Error: --fix requires a path argument\n", .{});
            std.process.exit(1);
        }
        common.fix_path = args[i.*];
        return true;
    } else if (std.mem.eql(u8, arg, "--sdf") or std.mem.eql(u8, arg, "-s")) {
        i.* += 1;
        if (i.* >= args.len) {
            std.debug.print("Error: {s} requires a path argument\n", .{arg});
            std.process.exit(1);
        }
        common.sdf_path = args[i.*];
        return true;
    } else if (std.mem.eql(u8, arg, "--strip-h")) {
        common.strip_h = true;
        return true;
    } else if (std.mem.eql(u8, arg, "--model")) {
        i.* += 1;
        if (i.* >= args.len) {
            std.debug.print("Error: --model requires 'all' or a model number\n", .{});
            std.process.exit(1);
        }
        const val = args[i.*];
        if (std.ascii.eqlIgnoreCase(val, "all")) {
            common.model_filter = .all;
        } else {
            const num = std.fmt.parseInt(u32, val, 10) catch {
                std.debug.print("Error: invalid --model value '{s}' (expected 'all' or a number)\n", .{val});
                std.process.exit(1);
            };
            common.model_filter = .{ .specific = num };
        }
        return true;
    }
    return false;
}

fn parseRunArgs(args: []const []const u8) ?RunConfig {
    var config = RunConfig{ .input_path = undefined };
    var common = CommonOptions{};
    var input_set = false;
    var i: usize = 0;

    while (i < args.len) : (i += 1) {
        const arg = args[i];

        if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {
            printRunUsage();
            return null;
        } else if (std.mem.eql(u8, arg, "-o") or std.mem.eql(u8, arg, "--output")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: {s} requires a path argument\n", .{arg});
                std.process.exit(1);
            }
            config.output_path = args[i];
        } else if (std.mem.eql(u8, arg, "-d") or std.mem.eql(u8, arg, "--dict")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: {s} requires a path argument\n", .{arg});
                std.process.exit(1);
            }
            config.dict_path = args[i];
        } else if (std.mem.eql(u8, arg, "--json")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: --json requires a path argument\n", .{});
                std.process.exit(1);
            }
            config.json_path = args[i];
        } else if (std.mem.eql(u8, arg, "--state-json")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: --state-json requires a path argument\n", .{});
                std.process.exit(1);
            }
            config.state_json_path = args[i];
        } else if (std.mem.eql(u8, arg, "--dump-movers")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: --dump-movers requires a path argument\n", .{});
                std.process.exit(1);
            }
            config.dump_movers_path = args[i];
        } else if (std.mem.eql(u8, arg, "--validate")) {
            config.validate = true;
        } else if (std.mem.eql(u8, arg, "--timing")) {
            config.show_timing = true;
        } else if (parseCommonOption(args, &i, &common)) {
            // consumed by parseCommonOption
        } else if (arg.len > 0 and arg[0] == '-') {
            std.debug.print("Error: unknown option '{s}'\n", .{arg});
            std.process.exit(1);
        } else {
            if (input_set) {
                std.debug.print("Error: unexpected positional argument '{s}'\n", .{arg});
                std.process.exit(1);
            }
            config.input_path = arg;
            input_set = true;
        }
    }

    if (!input_set) {
        std.debug.print("Error: missing input file\n", .{});
        printRunUsage();
        std.process.exit(1);
    }

    config.no_opt = common.no_opt;
    config.no_flip = common.no_flip;
    config.water = common.water;
    config.bond_policy = common.bond_policy;
    config.nterm_mode = common.nterm_mode;
    config.protonation_path = common.protonation_path;
    config.fix_path = common.fix_path;
    config.sdf_path = common.sdf_path;
    config.strip_h = common.strip_h;
    config.model_filter = common.model_filter;

    return config;
}

fn parsePkaArgs(args: []const []const u8) ?PkaConfig {
    var config = PkaConfig{ .input_path = undefined };
    var input_set = false;
    var i: usize = 0;
    while (i < args.len) : (i += 1) {
        const arg = args[i];
        if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {
            printPkaUsage();
            return null;
        } else if (std.mem.eql(u8, arg, "-o") or std.mem.eql(u8, arg, "--output")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: {s} requires a path argument\n", .{arg});
                std.process.exit(1);
            }
            config.output_path = args[i];
        } else if (std.mem.eql(u8, arg, "--model-only")) {
            config.model_only = true;
        } else if (std.mem.eql(u8, arg, "--timing")) {
            config.show_timing = true;
        } else if (std.mem.eql(u8, arg, "--quiet")) {
            config.quiet = true;
        } else if (std.mem.eql(u8, arg, "--jsonl")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: --jsonl requires a path\n", .{});
                std.process.exit(1);
            }
            config.jsonl_path = args[i];
        } else if (std.mem.eql(u8, arg, "-j") or std.mem.eql(u8, arg, "--threads")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: {s} requires a number\n", .{arg});
                std.process.exit(1);
            }
            config.n_threads = std.fmt.parseInt(u32, args[i], 10) catch {
                std.debug.print("Error: invalid thread count '{s}'\n", .{args[i]});
                std.process.exit(1);
            };
        } else if (std.mem.eql(u8, arg, "--propka-keep-waters")) {
            config.propka_keep_waters = true;
        } else if (std.mem.eql(u8, arg, "--engine")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: --engine requires 'builtin', 'empirical' or 'propka'\n", .{});
                std.process.exit(1);
            }
            config.engine = std.meta.stringToEnum(PkaEngine, args[i]) orelse {
                std.debug.print("Error: unknown --engine '{s}' (use 'builtin', 'empirical' or 'propka')\n", .{args[i]});
                std.process.exit(1);
            };
        } else if (std.mem.eql(u8, arg, "--empirical-mode")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: --empirical-mode requires 'parity' or 'zhyd'\n", .{});
                std.process.exit(1);
            }
            config.empirical_mode = std.meta.stringToEnum(zreduce.pka_empirical.Mode, args[i]) orelse {
                std.debug.print("Error: unknown --empirical-mode '{s}' (use 'parity' or 'zhyd')\n", .{args[i]});
                std.process.exit(1);
            };
        } else if (std.mem.eql(u8, arg, "--empirical-fixes")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: --empirical-fixes requires a list such as 'Q2,Q7' or 'none'\n", .{});
                std.process.exit(1);
            }
            config.empirical_fixes = zreduce.pka_empirical.Fixes.parse(args[i]) catch {
                std.debug.print("Error: invalid --empirical-fixes '{s}' (comma-separated Q1..Q8, or 'none')\n", .{args[i]});
                std.process.exit(1);
            };
        } else if (std.mem.eql(u8, arg, "--parameters")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: --parameters requires a path argument\n", .{});
                std.process.exit(1);
            }
            config.parameters_path = args[i];
        } else if (std.mem.eql(u8, arg, "--model")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: --model requires 'all' or a model number\n", .{});
                std.process.exit(1);
            }
            if (std.ascii.eqlIgnoreCase(args[i], "all")) {
                config.model_filter = .all;
            } else {
                const num = std.fmt.parseInt(u32, args[i], 10) catch {
                    std.debug.print("Error: invalid --model value '{s}'\n", .{args[i]});
                    std.process.exit(1);
                };
                config.model_filter = .{ .specific = num };
            }
        } else if (arg.len > 0 and arg[0] == '-') {
            std.debug.print("Error: unknown option '{s}'\n", .{arg});
            std.process.exit(1);
        } else {
            if (input_set) {
                std.debug.print("Error: unexpected positional argument '{s}'\n", .{arg});
                std.process.exit(1);
            }
            config.input_path = arg;
            input_set = true;
        }
    }
    if (!input_set) {
        std.debug.print("Error: missing input file\n", .{});
        printPkaUsage();
        std.process.exit(1);
    }
    return config;
}

fn printUsage(program_name: []const u8) void {
    std.debug.print(
        \\zreduce {s} - Hydrogen placement for mmCIF/PDB structures
        \\
        \\USAGE:
        \\    {s} <command> [OPTIONS] <args>
        \\
        \\COMMANDS:
        \\    run           Process a single mmCIF or PDB file
        \\    batch         Process all mmCIF/PDB files in a directory
        \\    pka           Experimental pKa prediction report for a file or a directory
        \\                  (does not alter coordinates)
        \\    compile-dict  Pre-compile CCD dictionary to binary format
        \\
        \\GLOBAL OPTIONS:
        \\    -h, --help       Show this help message
        \\    -V, --version    Show version
        \\
        \\Use '{s} <command> --help' for more information.
        \\
    , .{ build_options.version, program_name, program_name });
}

fn printPkaUsage() void {
    std.debug.print(
        \\USAGE:
        \\    zreduce pka [OPTIONS] <input.cif|input.pdb>
        \\    zreduce pka [OPTIONS] <input_dir> [-o <output_dir>]
        \\
        \\EXPERIMENTAL:
        \\    Predict pKa values and write a versioned JSON report. This command
        \\    does not alter coordinates or claim validated charge assignment.
        \\
        \\OPTIONS:
        \\    -h, --help         Show this help message
        \\    -o, --output PATH  JSON output (default: stdout); for a directory input, the
        \\                       output directory (default: <input_dir>_pka/)
        \\    --model VALUE      Model selection: 'all' (default) or a number
        \\    --model-only       Emit intrinsic model-pKa baseline without structure terms
        \\    --parameters PATH  JSON parameter set (calibration/research; needs its own "version").
        \\                       builtin: any built-in parameter; empirical (zhyd mode only):
        \\                       "model_pka" and "desolvation_scale" per site (e.g. "CYS")
        \\    --engine NAME      'propka' (PROPKA 3.5.1 via zpropka; default when built with
        \\                       -Dpropka=true), 'empirical' (experimental MIT engine
        \\                       based on the PROPKA 3 model; default otherwise) or
        \\                       'builtin' (earlier experimental zhyd model; implied by
        \\                       --model-only)
        \\    --empirical-mode MODE  Empirical engine: 'zhyd' (default; adopted corrections
        \\                       and recalibrated parameters) or 'parity' (PROPKA 3
        \\                       model, quirks included)
        \\    --empirical-fixes LIST  Override the empirical corrections: comma-separated
        \\                       Q1..Q8 (e.g. 'Q2,Q7') or 'none' (ablation)
        \\    --propka-keep-waters  Pass crystal waters to PROPKA (exact PROPKA-on-file
        \\                       reproduction; default drops them like PDB2PQR)
        \\    --timing           Print a per-stage timing breakdown to stderr
        \\
        \\DIRECTORY INPUT (batch):
        \\    Every .cif/.pdb/.ent file (optionally .gz) in <input_dir> gets
        \\    <output_dir>/<stem>.json, identical to a single-file run on
        \\    <input_dir>/<file>. A failing file does not stop the batch; the exit
        \\    status is non-zero if any file failed.
        \\    -j, --threads N    Thread count (default: auto-detect). --engine propka
        \\                       always runs on 1 thread
        \\    --jsonl PATH       Per-file JSONL log (input, output, status, n_sites,
        \\                       seconds, error)
        \\    --quiet            Suppress progress output
        \\
    , .{});
}

fn printRunUsage() void {
    std.debug.print(
        \\USAGE:
        \\    zreduce run [OPTIONS] <input.cif|input.pdb>
        \\
        \\OPTIONS:
        \\    -h, --help         Show this help message
        \\    -d, --dict PATH    CCD dictionary
        \\    -s, --sdf PATH     SDF/MOL file with ligand topology (for non-CCD compounds)
        \\    -o, --output PATH  Output file (default: stdout)
        \\    --json PATH        Write JSON log to file
        \\    --state-json PATH  Write versioned applied protonation/formal-charge report
        \\    --protonation PATH  Residue protonation override file
        \\    --fix PATH         Force mover states from control file
        \\    --dump-movers PATH  Write available mover IDs/states to file
        \\    --no-opt           Skip optimization
        \\    --no-flip          Disable Asn/Gln/His flips
        \\    --validate         Print validation diagnostics
        \\    --timing           Print a per-stage timing breakdown to stderr
        \\    --water            Add water hydrogens (default: off)
        \\    --water-phantom    Allow zero-occupancy phantom water H when orientation is underdetermined
        \\    --water-occ-cutoff N  Skip waters with occupancy below N (default: 0.66)
        \\    --water-b-cutoff N    Skip waters with B-factor above N (default: 40.0)
        \\    --bond-mode MODE   Bond-length mode: neutron|xray (default: neutron)
        \\    --isotope NAME     Output isotope for added H: hydrogen|h|deuterium|d (default: hydrogen)
        \\    --deuterium        Shortcut for --isotope deuterium
        \\    --nterm MODE       N-terminal protonation: auto|aggressive|neutral (default: auto)
        \\                         auto       NH3+/NH2+ on real chain-first residues only (ChimeraX-compatible)
        \\                         aggressive also NH3+/NH2+ on chain-break residues (reduce2 first_in_chain)
        \\                         neutral    NH2 (no + flag) on non-PRO real N-termini; PRO keeps NH2+
        \\                                    chain-break residues keep the single break-amide H in auto/neutral
        \\    --strip-h          Remove existing H atoms before placement
        \\    --model VALUE      Model selection: 'all' (default) or a model number
        \\
    , .{});
}

pub fn main(init: std.process.Init) !void {
    const allocator = init.gpa;

    const args = init.minimal.args.toSlice(init.arena.allocator()) catch {
        std.debug.print("Fatal: cannot read process arguments\n", .{});
        std.process.exit(2);
    };

    if (args.len < 2) {
        printUsage(args[0]);
        std.process.exit(1);
    }

    const subcmd = args[1];

    if (std.mem.eql(u8, subcmd, "-h") or std.mem.eql(u8, subcmd, "--help")) {
        printUsage(args[0]);
        return;
    }
    if (std.mem.eql(u8, subcmd, "-V") or std.mem.eql(u8, subcmd, "--version")) {
        std.debug.print("zreduce {s}\n", .{build_options.version});
        return;
    }

    if (std.mem.eql(u8, subcmd, "run")) {
        runSubcommand(allocator, args[2..]);
    } else if (std.mem.eql(u8, subcmd, "batch")) {
        batchSubcommand(allocator, args[2..]);
    } else if (std.mem.eql(u8, subcmd, "pka")) {
        pkaSubcommand(allocator, args[2..]);
    } else if (std.mem.eql(u8, subcmd, "compile-dict")) {
        compileDictSubcommand(allocator, args[2..]);
    } else {
        std.debug.print("Error: unknown subcommand '{s}'\n", .{subcmd});
        printUsage(args[0]);
        std.process.exit(1);
    }
}

fn parseBatchArgs(args: []const []const u8) ?zreduce.batch.BatchConfig {
    var config = zreduce.batch.BatchConfig{ .input_dir = undefined };
    var common = CommonOptions{};
    var input_set = false;
    var i: usize = 0;

    while (i < args.len) : (i += 1) {
        const arg = args[i];
        if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {
            printBatchUsage();
            return null;
        } else if (std.mem.eql(u8, arg, "-o") or std.mem.eql(u8, arg, "--output")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: {s} requires a path\n", .{arg});
                std.process.exit(1);
            }
            config.output_dir = args[i];
        } else if (std.mem.eql(u8, arg, "-d") or std.mem.eql(u8, arg, "--dict")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: {s} requires a path\n", .{arg});
                std.process.exit(1);
            }
            config.dict_path = args[i];
        } else if (std.mem.eql(u8, arg, "--jsonl")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: --jsonl requires a path\n", .{});
                std.process.exit(1);
            }
            config.jsonl_path = args[i];
        } else if (std.mem.eql(u8, arg, "-j") or std.mem.eql(u8, arg, "--threads")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: {s} requires a number\n", .{arg});
                std.process.exit(1);
            }
            config.n_threads = std.fmt.parseInt(u32, args[i], 10) catch {
                std.debug.print("Error: invalid thread count '{s}'\n", .{args[i]});
                std.process.exit(1);
            };
        } else if (std.mem.eql(u8, arg, "--quiet")) {
            config.quiet = true;
        } else if (std.mem.eql(u8, arg, "--timing")) {
            config.show_timing = true;
        } else if (std.mem.eql(u8, arg, "--gz")) {
            config.gzip_output = true;
        } else if (parseCommonOption(args, &i, &common)) {
            // consumed by parseCommonOption
        } else if (arg.len > 0 and arg[0] == '-') {
            std.debug.print("Error: unknown option '{s}'\n", .{arg});
            std.process.exit(1);
        } else {
            if (input_set) {
                std.debug.print("Error: unexpected argument '{s}'\n", .{arg});
                std.process.exit(1);
            }
            config.input_dir = arg;
            input_set = true;
        }
    }

    if (!input_set) {
        std.debug.print("Error: missing input directory\n", .{});
        printBatchUsage();
        std.process.exit(1);
    }

    config.no_opt = common.no_opt;
    config.no_flip = common.no_flip;
    config.water = common.water;
    config.bond_policy = common.bond_policy;
    config.nterm_mode = common.nterm_mode;
    config.protonation_path = common.protonation_path;
    config.fix_path = common.fix_path;
    config.sdf_path = common.sdf_path;
    config.strip_h = common.strip_h;
    config.model_filter = common.model_filter;
    config.json_version = build_options.version;
    return config;
}

fn printBatchUsage() void {
    std.debug.print(
        \\USAGE:
        \\    zreduce batch [OPTIONS] <input_dir>
        \\
        \\OPTIONS:
        \\    -h, --help         Show this help message
        \\    -d, --dict PATH    CCD dictionary (loaded once)
        \\    -s, --sdf PATH     SDF/MOL file with ligand topology (for non-CCD compounds)
        \\    -o, --output PATH  Output directory (default: <input>_reduced/)
        \\    -j, --threads N    Thread count (default: auto-detect)
        \\    --jsonl PATH       Aggregated JSONL log file
        \\    --protonation PATH  Residue protonation override file
        \\    --fix PATH         Force mover states from control file
        \\    --no-opt           Skip optimization
        \\    --no-flip          Disable flips
        \\    --quiet            Suppress progress output
        \\    --timing           Print per-stage timing (summed over files), wall time
        \\                       and throughput to stderr
        \\    --water            Add water hydrogens (default: off)
        \\    --water-phantom    Allow zero-occupancy phantom water H when orientation is underdetermined
        \\    --water-occ-cutoff N  Skip waters with occupancy below N (default: 0.66)
        \\    --water-b-cutoff N    Skip waters with B-factor above N (default: 40.0)
        \\    --bond-mode MODE   Bond-length mode: neutron|xray (default: neutron)
        \\    --isotope NAME     Output isotope for added H: hydrogen|h|deuterium|d (default: hydrogen)
        \\    --deuterium        Shortcut for --isotope deuterium
        \\    --nterm MODE       N-terminal protonation: auto|aggressive|neutral (default: auto)
        \\    --strip-h          Remove existing H atoms before placement
        \\    --model VALUE      Model selection: 'all' (default) or a model number
        \\    --gz               Write gzip-compressed output (.cif.gz)
        \\
    , .{});
}

const CompileDictConfig = struct {
    input_path: []const u8,
    output_path: ?[]const u8 = null,
};

fn parseCompileDictArgs(args: []const []const u8) ?CompileDictConfig {
    var config = CompileDictConfig{ .input_path = undefined };
    var input_set = false;
    var i: usize = 0;

    while (i < args.len) : (i += 1) {
        const arg = args[i];
        if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {
            printCompileDictUsage();
            return null;
        } else if (std.mem.eql(u8, arg, "-o") or std.mem.eql(u8, arg, "--output")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: {s} requires a path argument\n", .{arg});
                std.process.exit(1);
            }
            config.output_path = args[i];
        } else if (arg.len > 0 and arg[0] == '-') {
            std.debug.print("Error: unknown option '{s}'\n", .{arg});
            std.process.exit(1);
        } else {
            if (input_set) {
                std.debug.print("Error: unexpected argument '{s}'\n", .{arg});
                std.process.exit(1);
            }
            config.input_path = arg;
            input_set = true;
        }
    }

    if (!input_set) {
        std.debug.print("Error: missing input CIF path\n", .{});
        printCompileDictUsage();
        std.process.exit(1);
    }

    if (config.output_path == null) {
        std.debug.print("Error: -o/--output is required\n", .{});
        std.process.exit(1);
    }

    return config;
}

fn printCompileDictUsage() void {
    std.debug.print(
        \\USAGE:
        \\    zreduce compile-dict [OPTIONS] <input.cif>
        \\
        \\OPTIONS:
        \\    -h, --help         Show this help message
        \\    -o, --output PATH  Output binary dictionary file (required)
        \\
    , .{});
}

fn compileDictSubcommand(allocator: Allocator, args: []const []const u8) void {
    const io = defaultIo();
    const config = parseCompileDictArgs(args) orelse return;

    // Read input
    const source = zreduce.run.readFile(allocator, config.input_path) catch |err| {
        std.debug.print("Error: cannot read '{s}': {s}\n", .{ config.input_path, @errorName(err) });
        std.process.exit(1);
    };
    defer allocator.free(source);

    // Reject if already binary
    if (zreduce.ccd_binary.isBinaryDict(source)) {
        std.debug.print("Error: input is already a compiled dictionary\n", .{});
        std.process.exit(1);
    }

    // Parse CIF
    var dict = zreduce.ccd.parseComponentDict(allocator, source) catch |err| {
        std.debug.print("Error: failed to parse CCD dictionary: {s}\n", .{@errorName(err)});
        std.process.exit(1);
    };
    defer dict.deinit();

    // Write binary
    const output_path = config.output_path.?;
    const file = std.Io.Dir.cwd().createFile(io, output_path, .{}) catch |err| {
        std.debug.print("Error: cannot create '{s}': {s}\n", .{ output_path, @errorName(err) });
        std.process.exit(1);
    };
    defer file.close(io);

    var write_buf: [65536]u8 = undefined;
    var file_writer = file.writer(io, &write_buf);
    zreduce.ccd_binary.writeDict(&file_writer.interface, &dict) catch |err| {
        std.debug.print("Error: failed to write binary dictionary: {s}\n", .{@errorName(err)});
        std.Io.Dir.cwd().deleteFile(io, output_path) catch {};
        std.process.exit(1);
    };
    file_writer.interface.flush() catch |err| {
        std.debug.print("Error: failed to flush output: {s}\n", .{@errorName(err)});
        std.Io.Dir.cwd().deleteFile(io, output_path) catch {};
        std.process.exit(1);
    };

    std.debug.print("Compiled {d} components to '{s}'\n", .{ dict.components.count(), output_path });
}

fn batchSubcommand(allocator: Allocator, args: []const []const u8) void {
    const config = parseBatchArgs(args) orelse return;
    zreduce.batch.run(allocator, config) catch |err| {
        switch (err) {
            error.SomeFilesFailed, error.NoFilesFound => {},
            else => std.debug.print("Error: batch processing failed: {s}\n", .{@errorName(err)}),
        }
        std.process.exit(1);
    };
}

fn runSubcommand(allocator: Allocator, args: []const []const u8) void {
    const config = parseRunArgs(args) orelse return;
    var wall = zreduce.timing.Stopwatch.start();

    // Load CCD dictionary (once, before processFile)
    var ccd_dict: ?zreduce.ccd.ComponentDict = null;
    if (config.dict_path) |dict_path| {
        const dict_source = zreduce.run.readFile(allocator, dict_path) catch |err| {
            std.debug.print("Error: cannot read dictionary '{s}': {s}\n", .{ dict_path, @errorName(err) });
            std.process.exit(1);
        };
        defer allocator.free(dict_source);
        ccd_dict = zreduce.ccd_binary.loadDict(allocator, dict_source) catch |err| {
            std.debug.print("Error: failed to load CCD dictionary: {s}\n", .{@errorName(err)});
            std.process.exit(1);
        };
    }
    defer if (ccd_dict) |*d| d.deinit();

    // Load SDF topology dictionary (optional)
    var sdf_dict: ?zreduce.ccd.ComponentDict = null;
    if (config.sdf_path) |sdf_path| {
        const sdf_source = zreduce.run.readFile(allocator, sdf_path) catch |err| {
            std.debug.print("Error: cannot read SDF file '{s}': {s}\n", .{ sdf_path, @errorName(err) });
            std.process.exit(1);
        };
        defer allocator.free(sdf_source);
        sdf_dict = zreduce.sdf.parseSdf(allocator, sdf_source) catch |err| {
            std.debug.print("Error: failed to parse SDF file: {s}\n", .{@errorName(err)});
            std.process.exit(1);
        };
    }
    defer if (sdf_dict) |*d| d.deinit();
    const load_ns = wall.read();

    const proc_config = zreduce.run.ProcessConfig{
        .input_path = config.input_path,
        .output_path = config.output_path,
        .dict = if (ccd_dict) |*d| d else null,
        .sdf_dict = if (sdf_dict) |*d| d else null,
        .json_path = config.json_path,
        .state_json_path = config.state_json_path,
        .json_version = build_options.version,
        .no_opt = config.no_opt,
        .no_flip = config.no_flip,
        .validate_flag = config.validate,
        .water = config.water,
        .bond_policy = config.bond_policy,
        .nterm_mode = config.nterm_mode,
        .protonation_path = config.protonation_path,
        .fix_path = config.fix_path,
        .dump_movers_path = config.dump_movers_path,
        .strip_h = config.strip_h,
        .format = zreduce.run.detectFormat(config.input_path),
        .model_filter = config.model_filter,
    };

    const result = zreduce.run.processFile(allocator, proc_config) catch |err| {
        std.debug.print("Error: processing failed: {s}\n", .{@errorName(err)});
        std.process.exit(1);
    };

    if (result.n_movers > 0) {
        std.debug.print("  Movers: {d} ({d} singletons, {d} brute-force, {d} greedy)\n", .{
            result.n_movers,
            result.n_singletons,
            result.n_brute_force,
            result.n_vertex_cut,
        });
    }

    std.debug.print("zreduce: placed {d} H atoms on {d} residues ({d} skipped)\n", .{
        result.n_placed,
        result.n_residues,
        result.totalSkipped(),
    });
    if (result.n_distance_derived > 0) {
        std.debug.print("  note: {d} residues used distance-based bond inference (no dictionary entry)\n", .{result.n_distance_derived});
    }
    if (result.n_skipped_missing_ref > 0) {
        std.debug.print("  warning: {d} H skipped due to missing reference atoms (potential plan bug)\n", .{result.n_skipped_missing_ref});
    }
    if (config.show_timing) {
        var stages = result.timing;
        stages.addNs(.load, load_ns);
        zreduce.timing.printStageProfile("RUN", &zreduce.batch.run_stages, stages, wall.read(), null);
    }
}

/// Prediction options resolved from the command line, plus the parsed
/// parameter files they borrow from.
const ResolvedPka = struct {
    options: zreduce.pka_run.Options,
    parsed_params: ?std.json.Parsed(zreduce.pka.Parameters) = null,
    parsed_overrides: ?std.json.Parsed(zreduce.pka_empirical.params.Overrides) = null,

    fn deinit(self: *ResolvedPka) void {
        if (self.parsed_params) |*p| p.deinit();
        if (self.parsed_overrides) |*p| p.deinit();
    }
};

/// Validate engine options, load the parameter file (once) and print the
/// engine banner. Exits on invalid input.
fn resolvePkaOptions(allocator: Allocator, config: PkaConfig) ResolvedPka {
    // --model-only is a builtin-engine feature, so it selects that engine.
    const engine = config.engine orelse if (config.model_only) PkaEngine.builtin else PkaEngine.default();
    if (engine == .propka and !zreduce.pka_propka.available) {
        std.debug.print("Error: this zreduce was built without PROPKA; rebuild with -Dpropka=true\n", .{});
        std.process.exit(1);
    }
    const empirical_mode = config.empirical_mode orelse .zhyd;
    if (config.model_only and engine != .builtin) {
        std.debug.print("Error: --model-only applies to --engine builtin only\n", .{});
        std.process.exit(1);
    }
    if (config.parameters_path != null and !(engine == .builtin or (engine == .empirical and empirical_mode == .zhyd))) {
        std.debug.print("Error: --parameters applies to --engine builtin, or --engine empirical with --empirical-mode zhyd\n", .{});
        std.process.exit(1);
    }
    if (engine != .empirical and (config.empirical_mode != null or config.empirical_fixes != null)) {
        std.debug.print("Error: --empirical-mode and --empirical-fixes apply to --engine empirical only\n", .{});
        std.process.exit(1);
    }

    var resolved = ResolvedPka{ .options = .{ .engine = engine } };
    if (config.parameters_path) |path| {
        const json = zreduce.run.readFile(allocator, path) catch |err| {
            std.debug.print("Error: cannot read parameters '{s}': {s}\n", .{ path, @errorName(err) });
            std.process.exit(1);
        };
        defer allocator.free(json);
        if (engine == .empirical) {
            resolved.parsed_overrides = zreduce.pka_empirical.params.parseOverrides(allocator, json) catch |err| {
                std.debug.print("Error: invalid parameters '{s}': {s}\n", .{ path, @errorName(err) });
                std.process.exit(1);
            };
        } else {
            resolved.parsed_params = zreduce.pka.parseParameters(allocator, json) catch |err| {
                std.debug.print("Error: invalid parameters '{s}': {s}\n", .{ path, @errorName(err) });
                std.process.exit(1);
            };
        }
    }
    const parameters: zreduce.pka.Parameters = if (resolved.parsed_params) |p| p.value else .{};
    const empirical_options = zreduce.pka_empirical.Options{
        .mode = empirical_mode,
        .fixes = config.empirical_fixes,
        .overrides = if (resolved.parsed_overrides) |p| p.value else null,
    };
    switch (engine) {
        .builtin => std.debug.print("warning: pKa model {s} is experimental and not validated for production charge assignment\n", .{parameters.version}),
        .propka => std.debug.print("note: pKa engine {s}; states are reported, not applied to coordinates\n", .{zreduce.pka_propka.engine_version}),
        .empirical => std.debug.print("warning: pKa engine {s} ({s}) is experimental; states are reported, not applied to coordinates\n", .{
            empirical_options.engineName(),
            empirical_options.parameters().version,
        }),
    }
    resolved.options = .{
        .engine = engine,
        .model_filter = config.model_filter,
        .prediction_config = .{
            .parameters = parameters,
            .include_structure_terms = !config.model_only,
            .include_coupling = !config.model_only,
        },
        .propka_options = .{ .keep_waters = config.propka_keep_waters },
        .empirical_options = empirical_options,
    };
    return resolved;
}

fn isDirectory(path: []const u8) bool {
    const io = defaultIo();
    var dir = std.Io.Dir.cwd().openDir(io, path, .{}) catch return false;
    dir.close(io);
    return true;
}

fn pkaSubcommand(allocator: Allocator, args: []const []const u8) void {
    const config = parsePkaArgs(args) orelse return;
    if (isDirectory(config.input_path)) return pkaBatchSubcommand(allocator, config);
    if (config.n_threads != null or config.jsonl_path != null or config.quiet) {
        std.debug.print("Error: -j/--threads, --jsonl and --quiet apply to a directory input only\n", .{});
        std.process.exit(1);
    }

    const total = zreduce.timing.Stopwatch.start();
    var wall = total;
    var times = zreduce.timing.StageTimes{};
    const source = zreduce.pka_run.readInput(allocator, config.input_path) catch |err| {
        std.debug.print("Error: cannot read input '{s}': {s}\n", .{ config.input_path, @errorName(err) });
        std.process.exit(1);
    };
    defer allocator.free(source);
    wall.lapInto(&times, .read);

    var resolved = resolvePkaOptions(allocator, config);
    defer resolved.deinit();
    wall.lapInto(&times, .load);

    var failure = zreduce.pka_run.Failure{};
    _ = zreduce.pka_run.processSource(allocator, source, config.input_path, config.output_path, &resolved.options, &failure, &times) catch |err| {
        switch (failure.stage) {
            .parse => if (zreduce.run.detectFormat(config.input_path) == .mmcif)
                std.debug.print("Error: failed to parse mmCIF: {s}\n", .{@errorName(err)})
            else
                std.debug.print("Error: failed to parse PDB: {s}\n", .{@errorName(err)}),
            .no_models => std.debug.print("Error: no matching models\n", .{}),
            .predict => std.debug.print("Error: pKa prediction failed: {s}\n", .{@errorName(err)}),
            .write => std.debug.print("Error: failed to write pKa report: {s}\n", .{@errorName(err)}),
            .read => std.debug.print("Error: cannot read input '{s}': {s}\n", .{ config.input_path, @errorName(err) }),
        }
        std.process.exit(1);
    };
    if (config.show_timing) {
        zreduce.timing.printStageProfile("PKA", &zreduce.pka_run.stages, times, total.read(), null);
    }
}

fn pkaBatchSubcommand(allocator: Allocator, config: PkaConfig) void {
    var load = zreduce.timing.Stopwatch.start();
    var resolved = resolvePkaOptions(allocator, config);
    defer resolved.deinit();
    var n_threads = config.n_threads orelse 0;
    if (resolved.options.engine == .propka and n_threads != 1) {
        // zpropka's thread safety is not established; run PROPKA serially.
        std.debug.print("note: --engine propka runs on 1 thread\n", .{});
        n_threads = 1;
    }
    zreduce.pka_batch.run(allocator, .{
        .input_dir = config.input_path,
        .output_dir = config.output_path,
        .jsonl_path = config.jsonl_path,
        .n_threads = n_threads,
        .quiet = config.quiet,
        .show_timing = config.show_timing,
        .options = resolved.options,
        .load_ns = load.read(),
    }) catch |err| {
        switch (err) {
            error.SomeFilesFailed, error.NoFilesFound => {},
            else => std.debug.print("Error: pKa batch failed: {s}\n", .{@errorName(err)}),
        }
        std.process.exit(1);
    };
}
