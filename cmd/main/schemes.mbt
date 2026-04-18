// FlowLabLite — Convection Schemes
//
// Provides alternative convection discretisations for the momentum equations.
// Currently implemented:
//   - First-order upwind  (baseline, used by default in main.mbt)
//   - QUICK               (Quadratic Upstream Interpolation for Convective
//                          Kinematics, Leonard 1979)
//
// Reference: Leonard, B.P. (1979). "A stable and accurate convective
//   modelling procedure based on quadratic upstream interpolation."
//   Comput. Methods Appl. Mech. Eng. 19, 59–98.
//
// Interface contract
// ------------------
// Both functions return a 2D array of size [n_rows][n_cols] containing the
// convective contribution  u * (∂φ/∂x) * dt  or  v * (∂φ/∂y) * dt  at every
// node (i, j).  The caller subtracts this array from the velocity field, exactly
// as the existing conv_u_x / conv_u_y terms are used in main.mbt.
//
// Interior nodes (i = 1..nr-2, j = 1..nc-2) use QUICK interpolation.
// Boundary rows/columns are clamped to 0.0 (BCs are applied elsewhere).
// Nodes adjacent to boundaries that lack the full QUICK stencil fall back to
// first-order upwind to avoid out-of-bounds access.

// ─────────────────────────────────────────────────────────────────────────────
// QUICK helper — interpolate φ at east face (j + ½) using upwind quadratic
// ─────────────────────────────────────────────────────────────────────────────

// East face (j + ½) for a positive-u flow (upstream = j-1, j, j+1):
//   φ_e = (6 φ_j + 3 φ_{j+1} - φ_{j-1}) / 8
// East face (j + ½) for a negative-u flow (upstream = j+2, j+1, j):
//   φ_e = (6 φ_{j+1} + 3 φ_j - φ_{j+2}) / 8
fn quick_east(
  phi  : Array[Array[Double]],
  i    : Int,
  j    : Int,
  nc   : Int,
  u_ij : Double,
) -> Double {
  if u_ij >= 0.0 {
    // Positive flow: upstream stencil needs j-1
    if j == 0 {
      return phi[i][j]   // boundary: degenerate to node value
    }
    (6.0 * phi[i][j] + 3.0 * phi[i][j + 1] - phi[i][j - 1]) / 8.0
  } else {
    // Negative flow: upstream stencil needs j+2
    if j + 2 >= nc {
      return phi[i][j + 1]   // boundary: degenerate to next node value
    }
    (6.0 * phi[i][j + 1] + 3.0 * phi[i][j] - phi[i][j + 2]) / 8.0
  }
}

// West face (j − ½):
//   Positive u: upstream stencil is j-2, j-1, j  → (6 φ_{j-1} + 3 φ_j - φ_{j-2}) / 8
//   Negative u: upstream stencil is j+1, j, j-1  → (6 φ_j + 3 φ_{j-1} - φ_{j+1}) / 8
fn quick_west(
  phi  : Array[Array[Double]],
  i    : Int,
  j    : Int,
  nc   : Int,
  u_ij : Double,
) -> Double {
  if u_ij >= 0.0 {
    if j - 2 < 0 {
      return phi[i][j - 1]   // boundary: degenerate
    }
    (6.0 * phi[i][j - 1] + 3.0 * phi[i][j] - phi[i][j - 2]) / 8.0
  } else {
    if j + 1 >= nc {
      return phi[i][j]   // boundary: degenerate
    }
    (6.0 * phi[i][j] + 3.0 * phi[i][j - 1] - phi[i][j + 1]) / 8.0
  }
}

// North face (i + ½) for y-direction:
//   Positive v: upstream stencil i-1, i, i+1 → (6 φ_i + 3 φ_{i+1} - φ_{i-1}) / 8
//   Negative v: upstream stencil i+2, i+1, i → (6 φ_{i+1} + 3 φ_i - φ_{i+2}) / 8
fn quick_north(
  phi  : Array[Array[Double]],
  i    : Int,
  j    : Int,
  nr   : Int,
  v_ij : Double,
) -> Double {
  if v_ij >= 0.0 {
    if i == 0 {
      return phi[i][j]
    }
    (6.0 * phi[i][j] + 3.0 * phi[i + 1][j] - phi[i - 1][j]) / 8.0
  } else {
    if i + 2 >= nr {
      return phi[i + 1][j]
    }
    (6.0 * phi[i + 1][j] + 3.0 * phi[i][j] - phi[i + 2][j]) / 8.0
  }
}

// South face (i − ½):
fn quick_south(
  phi  : Array[Array[Double]],
  i    : Int,
  j    : Int,
  nr   : Int,
  v_ij : Double,
) -> Double {
  if v_ij >= 0.0 {
    if i - 2 < 0 {
      return phi[i - 1][j]
    }
    (6.0 * phi[i - 1][j] + 3.0 * phi[i][j] - phi[i - 2][j]) / 8.0
  } else {
    if i + 1 >= nr {
      return phi[i][j]
    }
    (6.0 * phi[i][j] + 3.0 * phi[i - 1][j] - phi[i + 1][j]) / 8.0
  }
}

// ─────────────────────────────────────────────────────────────────────────────
// Public API
// ─────────────────────────────────────────────────────────────────────────────

/// QUICK convective term in the x-direction: u * ∂φ/∂x * dt.
///
/// Returns a [nr][nc] array.  Interior nodes (1..nr-2, 1..nc-2) use the QUICK
/// quadratic upstream interpolation; nodes that lack a full 3-point stencil
/// (next to walls) fall back to first-order upwind automatically via the
/// quick_east/west helpers.  Boundary rows/columns return 0.0.
///
/// Parameters
///   phi   – field being advected (e.g. u-velocity)
///   u_vel – advecting velocity in x-direction (same layout as phi)
///   nr    – number of rows
///   nc    – number of columns
///   dx    – x-cell size
///   dt    – time-step
fn quick_convection_x(
  phi   : Array[Array[Double]],
  u_vel : Array[Array[Double]],
  nr    : Int,
  nc    : Int,
  dx    : Double,
  dt    : Double,
) -> Array[Array[Double]] {
  let result = create_zeros_2d(nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      let u_ij  = u_vel[i][j]
      let phi_e = quick_east(phi, i, j, nc, u_ij)
      let phi_w = quick_west(phi, i, j, nc, u_ij)
      result[i][j] = u_ij * (phi_e - phi_w) / dx * dt
    }
  }
  result
}

/// QUICK convective term in the y-direction: v * ∂φ/∂y * dt.
///
/// Same conventions as quick_convection_x but operates on rows (i-direction).
fn quick_convection_y(
  phi   : Array[Array[Double]],
  v_vel : Array[Array[Double]],
  nr    : Int,
  nc    : Int,
  dy    : Double,
  dt    : Double,
) -> Array[Array[Double]] {
  let result = create_zeros_2d(nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      let v_ij  = v_vel[i][j]
      let phi_n = quick_north(phi, i, j, nr, v_ij)
      let phi_s = quick_south(phi, i, j, nr, v_ij)
      result[i][j] = v_ij * (phi_n - phi_s) / dy * dt
    }
  }
  result
}

/// First-order upwind convective term in x: u * ∂φ/∂x * dt  (backward difference).
/// Provided for comparison with QUICK; used when --scheme upwind (default).
fn upwind_convection_x(
  phi   : Array[Array[Double]],
  u_vel : Array[Array[Double]],
  nr    : Int,
  nc    : Int,
  dx    : Double,
  dt    : Double,
) -> Array[Array[Double]] {
  let result = create_zeros_2d(nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      result[i][j] = u_vel[i][j] * dt / dx * (phi[i][j] - phi[i][j - 1])
    }
  }
  result
}

/// First-order upwind convective term in y: v * ∂φ/∂y * dt  (backward difference).
fn upwind_convection_y(
  phi   : Array[Array[Double]],
  v_vel : Array[Array[Double]],
  nr    : Int,
  nc    : Int,
  dy    : Double,
  dt    : Double,
) -> Array[Array[Double]] {
  let result = create_zeros_2d(nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      result[i][j] = v_vel[i][j] * dt / dy * (phi[i][j] - phi[i - 1][j])
    }
  }
  result
}
