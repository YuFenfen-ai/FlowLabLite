// FlowLabLite — CLI Format Dispatch
//
// Parses --format, --solver, --list-formats from @env.args() and dispatches
// to the appropriate Tier-1 format function (or the NetCDF-4 stub).
//
// Usage (via moon run):
//   moon run cmd/main --target wasm -- --list-formats
//   moon run cmd/main --target wasm -- --format csv
//   moon run cmd/main --target wasm -- --format vtk [--solver chorin|simple|pcg|mac]
//   moon run cmd/main --target wasm -- --format tecplot [--solver all|chorin|...]
//   moon run cmd/main --target wasm -- --format netcdf   (stub: prints error)
//
// Output goes to stdout.  Redirect with shell:
//   moon run cmd/main --target wasm -- --format vtk --solver pcg > cavity_pcg.vtk

// ─────────────────────────────────────────────────────────────────────────────
// Argument helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Return true if any format-related flag is present in argv.
/// Scans from index 1 (skipping argv[0] = program name).
fn cli_wants_format(args : Array[String]) -> Bool {
  let mut i = 1
  while i < args.length() {
    if args[i] == "--format" || args[i] == "--list-formats" {
      return true
    }
    i = i + 1
  }
  false
}

/// Return the value that follows `flag` in argv, or `default` if absent.
/// E.g. cli_get_string(args, "--format", "json") → "vtk" when args has "--format" "vtk".
fn cli_get_string(
  args    : Array[String],
  flag    : String,
  default : String,
) -> String {
  let mut i = 1
  while i < args.length() - 1 {
    if args[i] == flag {
      return args[i + 1]
    }
    i = i + 1
  }
  default
}

// ─────────────────────────────────────────────────────────────────────────────
// Format availability listing
// ─────────────────────────────────────────────────────────────────────────────

fn print_list_formats() -> Unit {
  println("Available formats in this build:")
  println("  json      JSON (default, with local_viewer.html support)     \u2713")
  println("  csv       RFC 4180 CSV (.csv)                                \u2713")
  println("  vtk       VTK Legacy ASCII Format (.vtk, ParaView / VisIt)   \u2713")
  println("  vti       VTK XML ImageData (.vti, ParaView preferred)       \u2713")
  println("  tecplot   Tecplot ASCII Data Format (.dat)                   \u2713")
  if netcdf_available {
    println("  netcdf    CF-1.11 NetCDF-4 (.nc, ncview / Panoply / CDO)    \u2713")
  } else {
    println("  netcdf    CF-1.11 NetCDF-4 (.nc, ncview / Panoply / CDO)    \u2717  (libnetcdf not linked)")
  }
  println("")
  println("To redirect output to a file:")
  println("  moon run cmd/main --target wasm -- --format vtk > cavity.vtk")
  println("  moon run cmd/main --target wasm -- --format vtk --solver pcg > pcg.vtk")
}

// ─────────────────────────────────────────────────────────────────────────────
// Data collection from global solver state
// ─────────────────────────────────────────────────────────────────────────────

/// Build a cell-centred SolverExport from the MAC staggered global state.
/// Interpolates staggered u/v to cell centres using get_u_mac_at / get_v_mac_at.
fn collect_mac_export() -> SolverExport {
  let u_cc = create_zeros_2d(mac_nc, mac_nc)
  let v_cc = create_zeros_2d(mac_nc, mac_nc)
  let p_cc = create_zeros_2d(mac_nc, mac_nc)
  for i = 0; i < mac_nc; i = i + 1 {
    for j = 0; j < mac_nc; j = j + 1 {
      u_cc[i][j] = get_u_mac_at(i, j)
      v_cc[i][j] = get_v_mac_at(i, j)
      p_cc[i][j] = g_p_mac[i][j]
    }
  }
  {
    name:  "mac",
    label: "MAC Staggered-Grid PCG",
    nc:    mac_nc,
    h:     dx,
    x0:    dx / 2.0,
    y0:    dy / 2.0,
    steps: local_mac_nt,
    u:     u_cc,
    v:     v_cc,
    p:     p_cc,
  }
}

/// Collect post-run global state into SolverExport structs.
///
/// solver_filter: "all" | "chorin" | "simple" | "pcg" | "mac"
/// Only solvers whose run_* flag is true AND whose name matches the filter
/// are included.
fn collect_solver_exports(solver_filter : String) -> Array[SolverExport] {
  let result : Array[SolverExport] = []
  let want_all = solver_filter == "all"
  if (want_all || solver_filter == "chorin") && run_chorin {
    result.push({
      name:  "chorin",
      label: "Chorin (Projection + Gauss-Seidel)",
      nc:    nx,
      h:     dx,
      x0:    0.0,
      y0:    0.0,
      steps: local_nt,
      u:     g_u,
      v:     g_v,
      p:     g_p,
    })
  }
  if (want_all || solver_filter == "simple") && run_simple {
    result.push({
      name:  "simple",
      label: "SIMPLE (Pressure-Correction)",
      nc:    nx,
      h:     dx,
      x0:    0.0,
      y0:    0.0,
      steps: local_simple_n,
      u:     g_u_s,
      v:     g_v_s,
      p:     g_p_s,
    })
  }
  if (want_all || solver_filter == "pcg") && run_pcg {
    result.push({
      name:  "pcg",
      label: "Chorin-PCG (Collocated)",
      nc:    nx,
      h:     dx,
      x0:    0.0,
      y0:    0.0,
      steps: local_pcg_nt,
      u:     g_u_pcg,
      v:     g_v_pcg,
      p:     g_p_pcg,
    })
  }
  if (want_all || solver_filter == "mac") && run_mac {
    result.push(collect_mac_export())
  }
  result
}

// ─────────────────────────────────────────────────────────────────────────────
// Main dispatch
// ─────────────────────────────────────────────────────────────────────────────

/// Dispatch CLI format output.  Called after all solvers have been run.
///
/// Reads --format, --solver from args.
/// All output goes to stdout (redirect with > for file output).
fn dispatch_format_cli(args : Array[String]) -> Unit {
  // --list-formats: no solver run needed
  if args.contains("--list-formats") {
    print_list_formats()
    return
  }

  let fmt    = cli_get_string(args, "--format",  "json")
  let solver = cli_get_string(args, "--solver",  "all")

  if fmt == "json" {
    println("[io] Hint: omit --format (or all flags) for the full JSON output with timing.")
    println("[io] Explicit --format json is not supported; re-run without --format.")
    return
  }

  // For VTK / VTI (one file per solver): default solver is "chorin"
  if fmt == "vtk" || fmt == "vti" {
    let single_solver = if solver == "all" { "chorin" } else { solver }
    let exports = collect_solver_exports(single_solver)
    if exports.is_empty() {
      println("[io] Solver '\{single_solver}' is not available (check --solver or run_* flags).")
      println("[io] Available: chorin | simple | pcg | mac")
      return
    }
    if fmt == "vtk" {
      println(format_vtk(exports[0]))
    } else {
      println(format_vti(exports[0]))
    }
    return
  }

  // For CSV / Tecplot: collect all (or filtered) solvers
  let solvers = collect_solver_exports(solver)
  if solvers.is_empty() {
    println("[io] No solver data to export (check --solver or run_* flags).")
    return
  }

  if fmt == "csv" {
    println(format_csv(solvers))
  } else if fmt == "tecplot" {
    println(format_tecplot(solvers))
  } else if fmt == "netcdf" {
    match format_netcdf(solvers[0]) {
      Ok(s)    => println(s)
      Err(msg) => println("[io] \{msg}")
    }
  } else {
    println("[io] Unknown format: '\{fmt}'.")
    println("[io] Run with --list-formats to see supported formats.")
  }
}
