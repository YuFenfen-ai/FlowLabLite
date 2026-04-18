// FlowLabLite — I/O Format Writers  (Tier-1, zero external dependencies)
//
// Supported formats:
//   CSV     – RFC 4180 Comma-Separated Values
//   VTK     – VTK Legacy ASCII Format v2.0 (STRUCTURED_GRID)
//   VTI     – VTK XML ImageData Format v0.1 (PointData, ascii)
//   Tecplot – Tecplot ASCII Data Format (.dat, structured ZONE)
//
// All functions are pure: they take a SolverExport value and return a String.
// No global state is read; no side-effects beyond memory allocation.
//
// Tier-2 format (NetCDF-4) lives in io_netcdf_stub.mbt.

// ─────────────────────────────────────────────────────────────────────────────
// Data model
// ─────────────────────────────────────────────────────────────────────────────

/// Exported field data for one solver, sampled on a regular nc×nc grid.
///
/// Grid coordinates:
///   x[col] = x0 + col * h    (col = 0 .. nc-1)
///   y[row] = y0 + row * h    (row = 0 .. nc-1)
///
/// Node-centred solvers (Chorin / SIMPLE / PCG):  x0 = y0 = 0.0
/// Cell-centred MAC staggered solver:              x0 = y0 = h/2
struct SolverExport {
  name  : String                   // "chorin" | "simple" | "pcg" | "mac"
  label : String                   // human-readable description
  nc    : Int                      // grid points per dimension
  h     : Double                   // uniform grid spacing
  x0    : Double                   // x-coordinate of the first grid point
  y0    : Double                   // y-coordinate of the first grid point
  steps : Int                      // time steps (or iterations) completed
  u     : Array[Array[Double]]     // u[row][col]  (y × x indexing)
  v     : Array[Array[Double]]
  p     : Array[Array[Double]]
}

// ─────────────────────────────────────────────────────────────────────────────
// CSV  (RFC 4180)
// ─────────────────────────────────────────────────────────────────────────────

/// Produce a RFC 4180 CSV document containing fields for all solvers.
///
/// Columns: solver, i (row/y-index), j (col/x-index), x, y, u, v, p, vel_mag
///
/// All solvers are concatenated in the same file; the "solver" column
/// identifies each row's origin.
fn format_csv(solvers : Array[SolverExport]) -> String {
  let re = 2.0 / nu
  let buf = StringBuilder::new()
  buf.write_string("# FlowLabLite export | Re=\{re} | solvers: \{solvers.length()}\n")
  buf.write_string("solver,i,j,x,y,u,v,p,vel_mag\n")
  for si = 0; si < solvers.length(); si = si + 1 {
    let s = solvers[si]
    for row = 0; row < s.nc; row = row + 1 {
      for col = 0; col < s.nc; col = col + 1 {
        let x   = s.x0 + col.to_double() * s.h
        let y   = s.y0 + row.to_double() * s.h
        let u_  = s.u[row][col]
        let v_  = s.v[row][col]
        let p_  = s.p[row][col]
        let mag = (u_ * u_ + v_ * v_).sqrt()
        buf.write_string("\{s.name},\{row},\{col},\{x},\{y},\{u_},\{v_},\{p_},\{mag}\n")
      }
    }
  }
  buf.to_string()
}

// ─────────────────────────────────────────────────────────────────────────────
// VTK Legacy ASCII  (VTK File Formats v4.2, STRUCTURED_GRID)
// ─────────────────────────────────────────────────────────────────────────────

/// Produce a VTK Legacy ASCII file for one solver.
///
/// Dataset type: STRUCTURED_GRID with DIMENSIONS nc nc 1.
/// Point ordering: i (x/col) varies fastest, then j (y/row) — standard VTK convention.
/// Fields written:
///   SCALARS pressure  double 1
///   VECTORS velocity  double      (u, v, 0)
///   SCALARS vel_magnitude  double 1
///
/// Compatible with ParaView ≥ 5.0 and VisIt.
fn format_vtk(s : SolverExport) -> String {
  let re   = 2.0 / nu
  let npts = s.nc * s.nc
  let buf  = StringBuilder::new()

  // ── Header ────────────────────────────────────────────────────────────────
  buf.write_string("# vtk DataFile Version 2.0\n")
  buf.write_string("FlowLabLite \{s.label} Re=\{re} step=\{s.steps}\n")
  buf.write_string("ASCII\n")
  buf.write_string("DATASET STRUCTURED_GRID\n")
  buf.write_string("DIMENSIONS \{s.nc} \{s.nc} 1\n")
  buf.write_string("POINTS \{npts} double\n")

  // ── Point coordinates (i=col fastest, j=row outer) ───────────────────────
  for row = 0; row < s.nc; row = row + 1 {
    for col = 0; col < s.nc; col = col + 1 {
      let x = s.x0 + col.to_double() * s.h
      let y = s.y0 + row.to_double() * s.h
      buf.write_string("\{x} \{y} 0.0\n")
    }
  }

  // ── Point data ────────────────────────────────────────────────────────────
  buf.write_string("\nPOINT_DATA \{npts}\n")

  // Pressure scalar field
  buf.write_string("SCALARS pressure double 1\n")
  buf.write_string("LOOKUP_TABLE default\n")
  for row = 0; row < s.nc; row = row + 1 {
    for col = 0; col < s.nc; col = col + 1 {
      buf.write_string("\{s.p[row][col]}\n")
    }
  }

  // Velocity vector field
  buf.write_string("\nVECTORS velocity double\n")
  for row = 0; row < s.nc; row = row + 1 {
    for col = 0; col < s.nc; col = col + 1 {
      let u_ = s.u[row][col]
      let v_ = s.v[row][col]
      buf.write_string("\{u_} \{v_} 0.0\n")
    }
  }

  // Velocity magnitude scalar field
  buf.write_string("\nSCALARS vel_magnitude double 1\n")
  buf.write_string("LOOKUP_TABLE default\n")
  for row = 0; row < s.nc; row = row + 1 {
    for col = 0; col < s.nc; col = col + 1 {
      let u_  = s.u[row][col]
      let v_  = s.v[row][col]
      let mag = (u_ * u_ + v_ * v_).sqrt()
      buf.write_string("\{mag}\n")
    }
  }

  buf.to_string()
}

// ─────────────────────────────────────────────────────────────────────────────
// Tecplot ASCII Data Format  (Tecplot Data Format Guide, structured ZONE)
// ─────────────────────────────────────────────────────────────────────────────

/// Produce a Tecplot ASCII .dat file containing all solvers as separate ZONEs.
///
/// Variables: X Y U V P VelMag
/// Each solver appears as one ZONE block with I=nc, J=nc, DATAPACKING=POINT.
/// SOLUTIONTIME = steps × dt (physical time; for SIMPLE, iterations × dt).
///
/// Compatible with Tecplot 360 EX and ParaView (Tecplot reader plug-in).
fn format_tecplot(solvers : Array[SolverExport]) -> String {
  let re  = 2.0 / nu
  let buf = StringBuilder::new()

  buf.write_string("TITLE = \"FlowLabLite Re=\{re} all solvers\"\n")
  buf.write_string("VARIABLES = \"X\", \"Y\", \"U\", \"V\", \"P\", \"VelMag\"\n")

  for si = 0; si < solvers.length(); si = si + 1 {
    let s    = solvers[si]
    let time = s.steps.to_double() * dt
    buf.write_string("ZONE T=\"\{s.name}\", I=\{s.nc}, J=\{s.nc}, K=1, ")
    buf.write_string("DATAPACKING=POINT, SOLUTIONTIME=\{time}\n")
    for row = 0; row < s.nc; row = row + 1 {
      for col = 0; col < s.nc; col = col + 1 {
        let x   = s.x0 + col.to_double() * s.h
        let y   = s.y0 + row.to_double() * s.h
        let u_  = s.u[row][col]
        let v_  = s.v[row][col]
        let p_  = s.p[row][col]
        let mag = (u_ * u_ + v_ * v_).sqrt()
        buf.write_string("\{x} \{y} \{u_} \{v_} \{p_} \{mag}\n")
      }
    }
  }

  buf.to_string()
}

// ─────────────────────────────────────────────────────────────────────────────
// Centerline velocity CSV
// ─────────────────────────────────────────────────────────────────────────────

/// Produce a CSV of centerline velocity profiles for one solver.
///
/// Extracts:
///   u_vcl — u along the vertical centerline   (col = nc/2, all rows)
///   v_hcl — v along the horizontal centerline (row = nc/2, all cols)
///
/// Output columns: component, index, coord, value
///   component: "u_vcl" or "v_hcl"
///   index    : grid index (0..nc-1)
///   coord    : physical coordinate (y for u_vcl, x for v_hcl)
///   value    : velocity value at that point
///
/// Total rows = 2 × nc (nc for u_vcl + nc for v_hcl).
fn format_centerline_csv(s : SolverExport) -> String {
  let re  = 2.0 / nu
  let buf = StringBuilder::new()
  buf.write_string("# FlowLabLite centerline velocity | solver=\{s.name} | Re=\{re} | step=\{s.steps}\n")
  buf.write_string("component,index,coord,value\n")

  // u along vertical centerline (x = x0 + (nc/2)*h, y varies)
  let jmid = s.nc / 2
  for row = 0; row < s.nc; row = row + 1 {
    let y = s.y0 + row.to_double() * s.h
    let u_val = s.u[row][jmid]
    buf.write_string("u_vcl,\{row},\{y},\{u_val}\n")
  }

  // v along horizontal centerline (y = y0 + (nc/2)*h, x varies)
  let imid = s.nc / 2
  for col = 0; col < s.nc; col = col + 1 {
    let x = s.x0 + col.to_double() * s.h
    let v_val = s.v[imid][col]
    buf.write_string("v_hcl,\{col},\{x},\{v_val}\n")
  }

  buf.to_string()
}

// ─────────────────────────────────────────────────────────────────────────────
// Monitor history CSV  (MonitorRecord — unified residual + physics monitor)
// ─────────────────────────────────────────────────────────────────────────────

/// One monitoring snapshot at a given step.
///
/// Fields:
///   step      — time step number (or SIMPLE iteration)
///   time      — physical time = step × dt (pseudo-time for SIMPLE)
///   div_norm  — divergence norm (mass conservation indicator)
///   r_u       — relative u velocity change between snapshots (-1.0 = not available)
///   r_v       — relative v velocity change between snapshots (-1.0 = not available)
///   ke        — total kinetic energy (0.5 · Σ(u²+v²)·dx·dy)
///   cfl_max   — maximum CFL number (-1.0 = not applicable for steady-state)
///   pcg_iters — pressure solver iteration count (50 for GS-fixed; actual for PCG)
struct MonitorRecord {
  step      : Int
  time      : Double
  div_norm  : Double
  r_u       : Double
  r_v       : Double
  ke        : Double
  cfl_max   : Double
  pcg_iters : Int
}

/// Produce a CSV of monitoring history records.
///
/// Header comment includes solver name and key parameters.
/// Columns: step, time, div_norm, r_u, r_v, ke, cfl_max, pcg_iters
/// Values of -1.0 are rendered as "-" to indicate "not available".
fn format_monitor_csv(
  records     : Array[MonitorRecord],
  solver_name : String,
  re          : Double,
) -> String {
  let buf = StringBuilder::new()
  buf.write_string("# FlowLabLite monitor | solver=\{solver_name} | Re=\{re} | dt=\{dt}\n")
  buf.write_string("step,time,div_norm,r_u,r_v,ke,cfl_max,pcg_iters\n")
  for k = 0; k < records.length(); k = k + 1 {
    let r = records[k]
    // Render -1.0 sentinel as "-"
    let r_u_str   = if r.r_u    < 0.0 { "-" } else { "\{r.r_u}" }
    let r_v_str   = if r.r_v    < 0.0 { "-" } else { "\{r.r_v}" }
    let cfl_str   = if r.cfl_max < 0.0 { "-" } else { "\{r.cfl_max}" }
    buf.write_string("\{r.step},\{r.time},\{r.div_norm},\{r_u_str},\{r_v_str},\{r.ke},\{cfl_str},\{r.pcg_iters}\n")
  }
  buf.to_string()
}

// ─────────────────────────────────────────────────────────────────────────────
// VTK XML ImageData (.vti)  (VTK XML Formats v0.1, ImageData)
// ─────────────────────────────────────────────────────────────────────────────

/// Produce a VTK XML ImageData (.vti) file for one solver.
///
/// Format: VTK XML 4.2 ImageData, ascii encoding, PointData.
/// Compatible with ParaView ≥ 3.0 (preferred format over VTK Legacy).
///
/// Grid mapping:
///   nc=41 (node-centred) → WholeExtent="0 40 0 40 0 0", Origin=(x0, y0, 0)
///   nc=40 (MAC cell-centred) → WholeExtent="0 39 0 39 0 0", Origin=(x0, y0, 0)
///
/// Three DataArrays (PointData):
///   1. pressure       — Float64, 1 component
///   2. velocity       — Float64, 3 components (u, v, 0)
///   3. vel_magnitude  — Float64, 1 component
///
/// Data ordering: x (col) varies fastest, then y (row) — VTK standard.
fn format_vti(s : SolverExport) -> String {
  let re     = 2.0 / nu
  let nc1    = s.nc - 1                 // extent upper bound per axis
  let buf    = StringBuilder::new()

  // ── XML header ────────────────────────────────────────────────────────────
  buf.write_string("<?xml version=\"1.0\"?>\n")
  buf.write_string("<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n")
  buf.write_string("  <!-- FlowLabLite \{s.label} Re=\{re} step=\{s.steps} -->\n")
  buf.write_string("  <ImageData WholeExtent=\"0 \{nc1} 0 \{nc1} 0 0\"")
  buf.write_string(" Origin=\"\{s.x0} \{s.y0} 0\"")
  buf.write_string(" Spacing=\"\{s.h} \{s.h} 1\">\n")
  buf.write_string("    <Piece Extent=\"0 \{nc1} 0 \{nc1} 0 0\">\n")
  buf.write_string("      <PointData>\n")

  // ── DataArray 1: pressure ─────────────────────────────────────────────────
  buf.write_string("        <DataArray type=\"Float64\" Name=\"pressure\" format=\"ascii\">\n")
  for row = 0; row < s.nc; row = row + 1 {
    for col = 0; col < s.nc; col = col + 1 {
      buf.write_string("          \{s.p[row][col]}\n")
    }
  }
  buf.write_string("        </DataArray>\n")

  // ── DataArray 2: velocity (3-component, w=0) ──────────────────────────────
  buf.write_string("        <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n")
  for row = 0; row < s.nc; row = row + 1 {
    for col = 0; col < s.nc; col = col + 1 {
      let u_ = s.u[row][col]
      let v_ = s.v[row][col]
      buf.write_string("          \{u_} \{v_} 0.0\n")
    }
  }
  buf.write_string("        </DataArray>\n")

  // ── DataArray 3: vel_magnitude ────────────────────────────────────────────
  buf.write_string("        <DataArray type=\"Float64\" Name=\"vel_magnitude\" format=\"ascii\">\n")
  for row = 0; row < s.nc; row = row + 1 {
    for col = 0; col < s.nc; col = col + 1 {
      let u_  = s.u[row][col]
      let v_  = s.v[row][col]
      let mag = (u_ * u_ + v_ * v_).sqrt()
      buf.write_string("          \{mag}\n")
    }
  }
  buf.write_string("        </DataArray>\n")

  // ── Closing tags ──────────────────────────────────────────────────────────
  buf.write_string("      </PointData>\n")
  buf.write_string("    </Piece>\n")
  buf.write_string("  </ImageData>\n")
  buf.write_string("</VTKFile>\n")

  buf.to_string()
}
