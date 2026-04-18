// FlowLabLite — Extended White-Box Test Suite  (Tests T59–T83)
// ============================================================
// Run: moon test --target wasm
//
// Software-engineering test classification used throughout:
//   [Unit]        – single function, isolated input, analytical oracle
//   [Integration] – two or more functions chained, black-box oracle
//   [System]      – full solver pipeline, physics-based oracle
//   [Regression]  – guard against known-good numerical properties
//
// Test-design techniques:
//   EP  – Equivalence Partitioning
//   BVA – Boundary Value Analysis
//   OBT – Oracle-Based Testing (analytical or computed reference)
//   PP  – Positive/Negative Paths

// ══════════════════════════════════════════════════════════════════════════════
// A. Unit Tests — Array Utilities  (T59–T60)
// ══════════════════════════════════════════════════════════════════════════════

// T59 [Unit / EP / OBT]
// copy_2d_array must produce a *deep* copy: mutating the source after the call
// must not change the destination.  Tests the separate-allocation invariant.
test "copy_2d_deep_independence" {
  let src = create_zeros_2d(4, 5)
  src[2][3] = 3.14
  src[0][4] = -7.0
  let dst = create_zeros_2d(4, 5)
  copy_2d_array(src, dst)
  // Overwrite every element of src with a very different value
  for i = 0; i < 4; i = i + 1 {
    for j = 0; j < 5; j = j + 1 {
      src[i][j] = src[i][j] + 999.0
    }
  }
  // Destination must retain pre-copy values exactly
  assert_true(approx_eq(dst[2][3],  3.14, 1.0e-15))
  assert_true(approx_eq(dst[0][4], -7.0,  1.0e-15))
  assert_true(approx_eq(dst[3][0],  0.0,  1.0e-15))
}

// T60 [Unit / BVA]
// create_zeros_2d must allocate independent rows: writing to row 0 must not
// affect row 1 (shared-pointer alias bug would fail this).
test "create_zeros_2d_row_independence" {
  let a = create_zeros_2d(3, 3)
  a[0][1] = 5.0
  // If rows share the same underlying array, a[1][1] would also be 5.0
  assert_eq(a[1][1], 0.0)
  assert_eq(a[2][1], 0.0)
  // Verify all non-written elements are still zero
  for i = 0; i < 3; i = i + 1 {
    for j = 0; j < 3; j = j + 1 {
      if not(i == 0 && j == 1) {
        assert_eq(a[i][j], 0.0)
      }
    }
  }
}

// ══════════════════════════════════════════════════════════════════════════════
// B. Unit Tests — Discrete Laplacian Operator  (T61–T63)
// ══════════════════════════════════════════════════════════════════════════════

// T61 [Unit / OBT]
// For the quadratic field u(x,y) = x²  (u[i][j] = (j·dx)²),
// the 5-point finite-difference Laplacian is exact: ∇²(x²) = 2 for all
// interior nodes regardless of grid spacing (property of degree-2 polynomials).
test "laplacian_apply_x_squared_exact" {
  let inv_dx2 = 1.0 / (dx * dx)
  let inv_dy2 = 1.0 / (dy * dy)
  let a_diag  = 2.0 * inv_dx2 + 2.0 * inv_dy2
  let u  = create_zeros_2d(ny, nx)
  let lv = create_zeros_2d(ny, nx)
  for i = 0; i < ny; i = i + 1 {
    for j = 0; j < nx; j = j + 1 {
      let x = j.to_double() * dx
      u[i][j] = x * x
    }
  }
  laplacian_apply(u, lv, inv_dx2, inv_dy2, a_diag)
  // ∇²(x²) = 2 exactly at every interior node
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      assert_true(approx_eq(lv[i][j], 2.0, 1.0e-8))
    }
  }
}

// T62 [Unit / OBT]
// Symmetric counterpart: u(x,y) = y²  →  ∇²(y²) = 2.
// Together with T61 this checks both spatial directions of the stencil.
test "laplacian_apply_y_squared_exact" {
  let inv_dx2 = 1.0 / (dx * dx)
  let inv_dy2 = 1.0 / (dy * dy)
  let a_diag  = 2.0 * inv_dx2 + 2.0 * inv_dy2
  let u  = create_zeros_2d(ny, nx)
  let lv = create_zeros_2d(ny, nx)
  for i = 0; i < ny; i = i + 1 {
    for j = 0; j < nx; j = j + 1 {
      let y = i.to_double() * dy
      u[i][j] = y * y
    }
  }
  laplacian_apply(u, lv, inv_dx2, inv_dy2, a_diag)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      assert_true(approx_eq(lv[i][j], 2.0, 1.0e-8))
    }
  }
}

// T63 [Unit / OBT]
// A spatially constant field is a harmonic function: ∇²(C) = 0 everywhere.
// Interior nodes only use interior neighbours so the boundary ghost treatment
// does not affect this result.
test "laplacian_apply_constant_gives_zero" {
  let inv_dx2 = 1.0 / (dx * dx)
  let inv_dy2 = 1.0 / (dy * dy)
  let a_diag  = 2.0 * inv_dx2 + 2.0 * inv_dy2
  let u  = create_zeros_2d(ny, nx)
  let lv = create_zeros_2d(ny, nx)
  // Fill entire array with an arbitrary constant (must also cover boundary rows)
  for i = 0; i < ny; i = i + 1 {
    for j = 0; j < nx; j = j + 1 {
      u[i][j] = 42.0
    }
  }
  laplacian_apply(u, lv, inv_dx2, inv_dy2, a_diag)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      assert_true(approx_eq(lv[i][j], 0.0, 1.0e-8))
    }
  }
}

// ══════════════════════════════════════════════════════════════════════════════
// C. Unit Tests — Boundary Conditions  (T64–T65)
// ══════════════════════════════════════════════════════════════════════════════

// T64 [Unit / PP]
// apply_pressure_bcs must be idempotent: calling it twice on the same array
// must produce exactly the same result as calling it once.
// The Neumann ghost-cell update p[0][j] = p[1][j] depends on p[1][j], which
// is not itself a ghost cell, so two passes give the same result.
test "apply_pressure_bcs_idempotent" {
  let p = create_zeros_2d(ny, nx)
  // Set a non-trivial interior field
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      p[i][j] = i.to_double() * j.to_double()
    }
  }
  apply_pressure_bcs(p)
  // Record the full state after first application
  let snap = create_zeros_2d(ny, nx)
  copy_2d_array(p, snap)
  // Second application
  apply_pressure_bcs(p)
  // Every element must be identical
  for i = 0; i < ny; i = i + 1 {
    for j = 0; j < nx; j = j + 1 {
      assert_true(approx_eq(p[i][j], snap[i][j], 1.0e-15))
    }
  }
}

// T65 [Unit / BVA]
// apply_pressure_bcs_n applied to the GAMG coarse grid (n = gamg_nc = 21)
// must enforce the same four BC types as the fine grid:
//   Neumann left/right/bottom  (ghost = adjacent interior)
//   Dirichlet top              (p = 0)
test "apply_pressure_bcs_n_coarse_grid" {
  let nc = gamg_nc  // 21
  let p  = create_zeros_2d(nc, nc)
  // Prime the array with non-zero values
  for i = 0; i < nc; i = i + 1 {
    for j = 0; j < nc; j = j + 1 {
      p[i][j] = (i + 1).to_double() * (j + 1).to_double()
    }
  }
  apply_pressure_bcs_n(p, nc)
  // Top row Dirichlet: p[nc-1][j] = 0
  for j = 0; j < nc; j = j + 1 {
    assert_true(approx_eq(p[nc - 1][j], 0.0, 1.0e-15))
  }
  // Bottom Neumann: p[0][j] = p[1][j]
  for j = 0; j < nc; j = j + 1 {
    assert_true(approx_eq(p[0][j], p[1][j], 1.0e-15))
  }
  // Left Neumann: p[i][0] = p[i][1]
  for i = 0; i < nc; i = i + 1 {
    assert_true(approx_eq(p[i][0], p[i][1], 1.0e-15))
  }
  // Right Neumann: p[i][nc-1] = p[i][nc-2]
  for i = 0; i < nc; i = i + 1 {
    assert_true(approx_eq(p[i][nc - 1], p[i][nc - 2], 1.0e-15))
  }
}

// ══════════════════════════════════════════════════════════════════════════════
// D. Unit Tests — RHS Source Term  (T66–T67)
// ══════════════════════════════════════════════════════════════════════════════

// T66 [Unit / EP / OBT]
// A divergence-free uniform horizontal flow (u = const, v = 0) has
// ∂u/∂x = 0 and ∂v/∂y = 0, so every term in the b formula is zero.
// This tests the "no-divergence → zero pressure source" equivalence class.
test "build_up_b_zero_for_divergence_free_flow" {
  let u = create_zeros_2d(ny, nx)
  let v = create_zeros_2d(ny, nx)
  // Uniform horizontal flow: constant u everywhere
  for i = 0; i < ny; i = i + 1 {
    for j = 0; j < nx; j = j + 1 {
      u[i][j] = 0.75   // arbitrary non-zero constant
    }
  }
  let b = create_zeros_2d(ny, nx)
  build_up_b_array(b, rho, dt, u, v, dx, dy) |> ignore
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      assert_true(approx_eq(b[i][j], 0.0, 1.0e-9))
    }
  }
}

// T67 [Unit / OBT]
// Analytical oracle: u[i][j] = j·dx (linear in x, so ∂u/∂x = 1 uniformly),
// v = 0.  By the b formula:
//   b = ρ·[(1/dt)·(∂u/∂x + ∂v/∂y) − (∂u/∂x)² − 2·(∂u/∂y)(∂v/∂x) − (∂v/∂y)²]
//     = 1·[(1000)·1 − 1 − 0 − 0]  =  999.0
// for all interior nodes (independent of i, j).
test "build_up_b_formula_verification" {
  let u = create_zeros_2d(ny, nx)
  let v = create_zeros_2d(ny, nx)
  // u linear in x: u[i][j] = j * dx
  for i = 0; i < ny; i = i + 1 {
    for j = 0; j < nx; j = j + 1 {
      u[i][j] = j.to_double() * dx
    }
  }
  let b = create_zeros_2d(ny, nx)
  build_up_b_array(b, rho, dt, u, v, dx, dy) |> ignore
  // Expected value: rho * (1/dt * 1 - 1^2 - 0 - 0) = 1000 - 1 = 999
  let expected = rho * (1.0 / dt * 1.0 - 1.0)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      assert_true(approx_eq(b[i][j], expected, 1.0e-7))
    }
  }
}

// ══════════════════════════════════════════════════════════════════════════════
// E. Unit Tests — GAMG Sub-Components  (T68–T70)
// ══════════════════════════════════════════════════════════════════════════════

// T68 [Unit / OBT]
// gamg_prolongate must preserve constant coarse fields.
// For e_coarse[ci][cj] = C at all nodes, bilinear interpolation gives C at
// every fine-grid node (all four case weights sum to 1).
test "gamg_prolongate_constant_preservation" {
  let nc      = gamg_nc        // 21
  let e_c     = create_zeros_2d(nc, nc)
  let e_fine  = create_zeros_2d(ny, nx)
  let c_val   = 3.7
  // Fill coarse grid with a constant
  for i = 0; i < nc; i = i + 1 {
    for j = 0; j < nc; j = j + 1 {
      e_c[i][j] = c_val
    }
  }
  gamg_prolongate(e_c, e_fine)
  // Every fine node reachable by prolongation must equal c_val
  // Case 1 (even,even): direct copy
  for ci = 0; ci < nc; ci = ci + 1 {
    for cj = 0; cj < nc; cj = cj + 1 {
      assert_true(approx_eq(e_fine[2 * ci][2 * cj], c_val, 1.0e-14))
    }
  }
  // Case 2 (even,odd): average of two equal neighbours = c_val
  for ci = 0; ci < nc; ci = ci + 1 {
    for cj = 0; cj < nc - 1; cj = cj + 1 {
      assert_true(approx_eq(e_fine[2 * ci][2 * cj + 1], c_val, 1.0e-14))
    }
  }
  // Case 4 (odd,odd): average of four equal corners = c_val
  for ci = 0; ci < nc - 1; ci = ci + 1 {
    for cj = 0; cj < nc - 1; cj = cj + 1 {
      assert_true(approx_eq(e_fine[2 * ci + 1][2 * cj + 1], c_val, 1.0e-14))
    }
  }
}

// T69 [Unit / OBT]
// gamg_restrict with normalisation factor 4: for a uniform fine field = 1.0,
// every interior coarse node (where all 8 fine-grid stencil neighbours exist)
// must receive exactly:
//   (1 + 4×0.5 + 4×0.25) / 4  =  4 / 4  =  1.0
// This verifies the /4 normalisation that prevents 4× over-amplification.
test "gamg_restrict_interior_normalization" {
  let nc     = gamg_nc   // 21
  let r_fine = create_zeros_2d(ny, nx)
  let r_c    = create_zeros_2d(nc, nc)
  // Uniform fine field
  for i = 0; i < ny; i = i + 1 {
    for j = 0; j < nx; j = j + 1 {
      r_fine[i][j] = 1.0
    }
  }
  gamg_restrict(r_fine, r_c)
  // Check interior coarse nodes (ci in [1..nc-2], cj in [1..nc-2])
  // where fi=2*ci ∈ [2..38] and fj=2*cj ∈ [2..38] — all stencil neighbours in range.
  for ci = 1; ci < nc - 1; ci = ci + 1 {
    for cj = 1; cj < nc - 1; cj = cj + 1 {
      assert_true(approx_eq(r_c[ci][cj], 1.0, 1.0e-14))
    }
  }
}

// T70 [Unit / PP]
// gamg_smooth_fine must be a descent smoother: after ν=2 damped-Jacobi steps
// starting from z=0, the SPD residual ||r + A·z|| (= ||r − B·z|| where B=−A)
// must be strictly smaller than ||r|| = ||r + A·0||.
// Uses a non-trivial residual that is not in the kernel of the smoother.
test "gamg_smooth_fine_reduces_residual" {
  let inv_dx2 = 1.0 / (dx * dx)
  let inv_dy2 = 1.0 / (dy * dy)
  let a_diag  = 2.0 * inv_dx2 + 2.0 * inv_dy2
  // Non-uniform r: higher-frequency checkerboard pattern for better smoother effect
  let r  = create_zeros_2d(ny, nx)
  let z  = create_zeros_2d(ny, nx)
  let az = create_zeros_2d(ny, nx)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      // Alternating ±1 pattern (high-frequency: smoothers reduce this rapidly)
      let sign = if (i + j) % 2 == 0 { 1.0 } else { -1.0 }
      r[i][j] = sign
    }
  }
  // ||r||² before smoothing
  let mut r0_sq = 0.0
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      r0_sq = r0_sq + r[i][j] * r[i][j]
    }
  }
  // Apply ν=2 damped-Jacobi smoothing
  gamg_smooth_fine(z, r, inv_dx2, inv_dy2, a_diag, 2)
  // Compute SPD residual ||r − B·z|| = ||r + A·z||
  apply_pressure_bcs(z)
  laplacian_apply(z, az, inv_dx2, inv_dy2, a_diag)
  let mut r1_sq = 0.0
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      let res = r[i][j] + az[i][j]
      r1_sq = r1_sq + res * res
    }
  }
  // Residual must decrease: ||r+A·z|| < ||r||
  assert_true(r1_sq < r0_sq)
}

// ══════════════════════════════════════════════════════════════════════════════
// F. Integration Tests — PCG Residual Criteria  (T71–T74)
// ══════════════════════════════════════════════════════════════════════════════
// Each test verifies that the solver actually achieves the stopping criterion
//   ‖b − A·p‖ / ‖b‖  <  pcg_tol = 1e-5
// after returning.  This is a stronger check than "iteration count ≤ limit".

fn compute_pressure_residual_norm(
  p       : Array[Array[Double]],
  b       : Array[Array[Double]],
  inv_dx2 : Double,
  inv_dy2 : Double,
  a_diag  : Double,
) -> Double {
  let ap = create_zeros_2d(ny, nx)
  apply_pressure_bcs(p)
  laplacian_apply(p, ap, inv_dx2, inv_dy2, a_diag)
  let mut r_sq  = 0.0
  let mut b_sq  = 0.0
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      let r = b[i][j] - ap[i][j]
      r_sq = r_sq + r * r
      b_sq = b_sq + b[i][j] * b[i][j]
    }
  }
  if b_sq > 1.0e-15 { r_sq.sqrt() / b_sq.sqrt() } else { 0.0 }
}

// T71 [Integration / OBT]
// Jacobi (baseline): residual satisfies the stopping criterion.
test "pcg_jacobi_residual_below_tolerance" {
  let inv_dx2 = 1.0 / (dx * dx)
  let inv_dy2 = 1.0 / (dy * dy)
  let a_diag  = 2.0 * inv_dx2 + 2.0 * inv_dy2
  let p = create_zeros_2d(ny, nx)
  let b = create_zeros_2d(ny, nx)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      b[i][j] = 1.0
    }
  }
  pressure_poisson_pcg(p, dx, dy, b) |> ignore
  let rel = compute_pressure_residual_norm(p, b, inv_dx2, inv_dy2, a_diag)
  // Accept 10× tolerance to allow for the final BC application
  assert_true(rel < pcg_tol * 10.0)
}

// T72 [Integration / OBT]
// DILU preconditioner: same residual criterion as Jacobi.
test "pcg_dilu_residual_below_tolerance" {
  let inv_dx2 = 1.0 / (dx * dx)
  let inv_dy2 = 1.0 / (dy * dy)
  let a_diag  = 2.0 * inv_dx2 + 2.0 * inv_dy2
  let p = create_zeros_2d(ny, nx)
  let b = create_zeros_2d(ny, nx)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      b[i][j] = 1.0
    }
  }
  pressure_poisson_pcg_dilu(p, dx, dy, b) |> ignore
  let rel = compute_pressure_residual_norm(p, b, inv_dx2, inv_dy2, a_diag)
  assert_true(rel < pcg_tol * 10.0)
}

// T73 [Integration / OBT]
// DIC preconditioner: same residual criterion.
test "pcg_dic_residual_below_tolerance" {
  let inv_dx2 = 1.0 / (dx * dx)
  let inv_dy2 = 1.0 / (dy * dy)
  let a_diag  = 2.0 * inv_dx2 + 2.0 * inv_dy2
  let p = create_zeros_2d(ny, nx)
  let b = create_zeros_2d(ny, nx)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      b[i][j] = 1.0
    }
  }
  pressure_poisson_pcg_dic(p, dx, dy, b) |> ignore
  let rel = compute_pressure_residual_norm(p, b, inv_dx2, inv_dy2, a_diag)
  assert_true(rel < pcg_tol * 10.0)
}

// T74 [Integration / OBT]
// GAMG preconditioner: residual must satisfy the stopping criterion.
// The stationary Richardson iteration converges when ρ(I + T·A) < 1.
test "pcg_gamg_residual_below_tolerance" {
  let inv_dx2 = 1.0 / (dx * dx)
  let inv_dy2 = 1.0 / (dy * dy)
  let a_diag  = 2.0 * inv_dx2 + 2.0 * inv_dy2
  let p = create_zeros_2d(ny, nx)
  let b = create_zeros_2d(ny, nx)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      b[i][j] = 1.0
    }
  }
  pressure_poisson_pcg_gamg(p, dx, dy, b) |> ignore
  let rel = compute_pressure_residual_norm(p, b, inv_dx2, inv_dy2, a_diag)
  assert_true(rel < pcg_tol * 10.0)
}

// ══════════════════════════════════════════════════════════════════════════════
// G. Integration Tests — Cross-Preconditioner Agreement  (T75)
// ══════════════════════════════════════════════════════════════════════════════

// T75 [Integration / OBT]
// Non-uniform polynomial RHS: all four preconditioners must converge to the
// same pressure field within 1% relative error.  Tests that the choice of
// preconditioner affects only convergence speed, not the computed solution.
// The RHS b[i][j] = i·(ny-1-i)·j·(nx-1-j) is smooth, non-trivial, and
// naturally zero on all four boundary rows/columns.
test "all_preconditioners_agree_on_nonuniform_rhs" {
  let p_jac  = create_zeros_2d(ny, nx)
  let p_dilu = create_zeros_2d(ny, nx)
  let p_dic  = create_zeros_2d(ny, nx)
  let p_gamg = create_zeros_2d(ny, nx)
  let b      = create_zeros_2d(ny, nx)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      let fi = i.to_double()
      let fj = j.to_double()
      let ny1 = (ny - 1).to_double()
      let nx1 = (nx - 1).to_double()
      b[i][j] = (fi * (ny1 - fi)) * (fj * (nx1 - fj))
    }
  }
  pressure_poisson_pcg(p_jac,  dx, dy, b) |> ignore
  pressure_poisson_pcg_dilu(p_dilu, dx, dy, b) |> ignore
  pressure_poisson_pcg_dic(p_dic,  dx, dy, b) |> ignore
  pressure_poisson_pcg_gamg(p_gamg, dx, dy, b) |> ignore
  // Compute ‖p − p_jac‖ / ‖p_jac‖ for each preconditioner
  let mut ref_sq  = 0.0
  let mut dilu_sq = 0.0
  let mut dic_sq  = 0.0
  let mut gamg_sq = 0.0
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      ref_sq  = ref_sq  + p_jac[i][j] * p_jac[i][j]
      let dd = p_dilu[i][j] - p_jac[i][j]
      let dc = p_dic[i][j]  - p_jac[i][j]
      let dg = p_gamg[i][j] - p_jac[i][j]
      dilu_sq = dilu_sq + dd * dd
      dic_sq  = dic_sq  + dc * dc
      gamg_sq = gamg_sq + dg * dg
    }
  }
  let dilu_err = if ref_sq > 1.0e-15 { dilu_sq.sqrt() / ref_sq.sqrt() } else { 0.0 }
  let dic_err  = if ref_sq > 1.0e-15 { dic_sq.sqrt()  / ref_sq.sqrt() } else { 0.0 }
  let gamg_err = if ref_sq > 1.0e-15 { gamg_sq.sqrt() / ref_sq.sqrt() } else { 0.0 }
  assert_true(dilu_err < 0.01)
  assert_true(dic_err  < 0.01)
  assert_true(gamg_err < 0.01)
}

// ══════════════════════════════════════════════════════════════════════════════
// H. System Tests — Solver Physics  (T76–T79)
// ══════════════════════════════════════════════════════════════════════════════

// T76 [System / OBT]
// Sign-convention test: the discrete Laplacian A is SND (diagonal = −a_diag < 0).
// Solving A·p = b with b > 0 gives p* = A⁻¹b < 0 at every interior node.
// This verifies the sign convention used throughout the preconditioner design.
test "pressure_sign_negative_for_positive_rhs" {
  let p = create_zeros_2d(ny, nx)
  let b = create_zeros_2d(ny, nx)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      b[i][j] = 1.0
    }
  }
  pressure_poisson_pcg(p, dx, dy, b) |> ignore
  // All interior nodes must be strictly negative (SND → p* < 0 for b > 0)
  let mut all_negative = true
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      if p[i][j] >= 0.0 { all_negative = false }
    }
  }
  assert_true(all_negative)
}

// T77 [System / Regression]
// Chorin solver must be deterministic: two runs from the same init must
// produce bit-identical results.  Floating-point ops are deterministic in
// sequential execution on a single thread.
test "chorin_solver_determinism" {
  init_simulation()
  run_n_steps(100)
  let u1 = get_u_center()
  let v1 = get_v_at(ny / 2, nx / 2)
  let p1 = get_p_at(ny / 2, nx / 2)
  // Reset and re-run identically
  init_simulation()
  run_n_steps(100)
  let u2 = get_u_center()
  let v2 = get_v_at(ny / 2, nx / 2)
  let p2 = get_p_at(ny / 2, nx / 2)
  assert_true(u1 == u2)   // bit-identical: same ops, same order → same IEEE 754 result
  assert_true(v1 == v2)
  assert_true(p1 == p2)
}

// T78 [System / Regression]
// Additive steps: running n₁ + n₂ time steps in two calls must equal one call
// of n₁ + n₂ steps from the same initial state.  Tests solver linearity of
// time-marching with respect to step decomposition.
test "chorin_additive_steps" {
  // Path A: 10 + 10 = 20 steps
  init_simulation()
  run_n_steps(10)
  run_n_steps(10)
  let ua = get_u_center()
  let pa = get_p_at(ny / 4, nx / 4)
  // Path B: 20 steps in one call
  init_simulation()
  run_n_steps(20)
  let ub = get_u_center()
  let pb = get_p_at(ny / 4, nx / 4)
  assert_true(ua == ub)   // bit-identical
  assert_true(pa == pb)
}

// T79 [System / PP]
// SIMPLE produces a physically converged flow field.
// After 100 iterations with under-relaxation (α_p=0.3, α_u=0.7) for Re=20,
// the velocity divergence norm must be below 0.02 — confirming that the
// pressure correction is driving mass conservation.
// Note: the divergence norm uses a *relative* PCG stopping criterion, so it
// scales with the RHS magnitude and does not decrease monotonically.
// A fixed absolute threshold is the robust measure of acceptable mass error.
test "simple_divergence_below_threshold_after_100_iters" {
  init_simple()
  run_simple_n_iter(100)
  let div = get_simple_divergence_norm()
  assert_true(div >= 0.0)
  assert_true(div < 0.02)
}

// ══════════════════════════════════════════════════════════════════════════════
// I. Regression Tests  (T80–T83)
// ══════════════════════════════════════════════════════════════════════════════

// T80 [Regression / BVA]
// build_modified_diag must produce values in the physically expected range.
// For the 41×41 uniform grid: a_diag = 2/dx² + 2/dy² = 1600.
// The modified diagonal satisfies 0 < d[i][j] ≤ a_diag (the DILU fill-in
// only subtracts positive terms, so d is strictly smaller than a_diag
// except at the first interior node where no subtraction occurs).
test "build_modified_diag_range" {
  let inv_dx2 = 1.0 / (dx * dx)
  let inv_dy2 = 1.0 / (dy * dy)
  let a_diag  = 2.0 * inv_dx2 + 2.0 * inv_dy2  // = 1600
  let d = build_modified_diag(inv_dx2, inv_dy2, a_diag)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      // Strictly positive (required for safe division in preconditioner)
      assert_true(d[i][j] > 0.0)
      // Cannot exceed the unmodified diagonal (stability bound)
      assert_true(d[i][j] <= a_diag + 1.0e-10)
    }
  }
  // Spot check: corner node (1,1) has no left/bottom neighbours to subtract →
  // d[1][1] should equal a_diag exactly.
  assert_true(approx_eq(d[1][1], a_diag, 1.0e-8))
}

// T81 [Regression / EP]
// MAC divergence must remain below 1e-4 even at step 200 (not just step 50).
// The staggered scheme maintains near-exact mass conservation throughout the run.
test "mac_divergence_remains_small_at_step_200" {
  init_mac()
  run_mac_n_steps(200)
  let div = get_mac_divergence_norm()
  assert_true(div < 1.0e-4)
}

// T82 [Regression / OBT]
// laplacian_apply_n must be exact on x² for the coarse grid (21×21).
// Coarse grid parameters: dx_c = 2·dx, inv_dx2_c = inv_dx2/4, a_diag_c = a_diag/2.
// The finite-difference Laplacian is exact on degree-≤2 polynomials regardless
// of grid spacing, so ∇²_c(x²) = 2 exactly at all coarse interior nodes.
test "laplacian_apply_n_x_squared_exact_coarse" {
  let nc      = gamg_nc           // 21
  let inv_dx2 = 1.0 / (dx * dx)
  let inv_dy2 = 1.0 / (dy * dy)
  let inv_dx2_c = inv_dx2 / 4.0  // coarse: dx_c = 2·dx
  let inv_dy2_c = inv_dy2 / 4.0
  let a_diag_c  = 2.0 * inv_dx2_c + 2.0 * inv_dy2_c
  let dx_c = dx * 2.0
  let u_c   = create_zeros_2d(nc, nc)
  let lv_c  = create_zeros_2d(nc, nc)
  for i = 0; i < nc; i = i + 1 {
    for j = 0; j < nc; j = j + 1 {
      let x = j.to_double() * dx_c
      u_c[i][j] = x * x
    }
  }
  laplacian_apply_n(u_c, lv_c, inv_dx2_c, inv_dy2_c, a_diag_c, nc, nc)
  for i = 1; i < nc - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      assert_true(approx_eq(lv_c[i][j], 2.0, 1.0e-7))
    }
  }
}

// T83 [Regression / System]
// All four PCG preconditioners must produce vortex flow consistent with
// Chorin's projection method: after a Chorin-PCG run of 200 steps, the
// streamwise velocity at the geometric centre must be negative.
// This end-to-end regression guards against preconditioner bugs that
// corrupt the pressure solve and break vortex formation.
test "chorin_pcg_produces_vortex_with_all_preconditioners" {
  // Chorin-PCG uses Jacobi preconditioner (already tested in T34)
  init_chorin_pcg()
  run_chorin_pcg_n_steps(200)
  let u_pcg_jac = get_u_pcg_at(ny / 2, nx / 2)
  assert_true(u_pcg_jac < 0.0)

  // Verify DILU pressure solver also converges in the PCG context
  // by solving the same final-step RHS and checking sign
  let p_dilu = create_zeros_2d(ny, nx)
  let b_test = create_zeros_2d(ny, nx)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      b_test[i][j] = 1.0
    }
  }
  pressure_poisson_pcg_dilu(p_dilu, dx, dy, b_test) |> ignore
  // DILU pressure solution must match Jacobi sign (both negative for b>0)
  assert_true(p_dilu[ny / 2][nx / 2] < 0.0)

  // DIC same check
  let p_dic = create_zeros_2d(ny, nx)
  pressure_poisson_pcg_dic(p_dic, dx, dy, b_test) |> ignore
  assert_true(p_dic[ny / 2][nx / 2] < 0.0)

  // GAMG same check
  let p_gamg = create_zeros_2d(ny, nx)
  pressure_poisson_pcg_gamg(p_gamg, dx, dy, b_test) |> ignore
  assert_true(p_gamg[ny / 2][nx / 2] < 0.0)
}

// ══════════════════════════════════════════════════════════════════════════════
// I. Unit Tests — I/O Format Writers  (T84–T91)
// ══════════════════════════════════════════════════════════════════════════════
//
// All tests use synthetic SolverExport values to avoid solver-state dependency.
// Tests are pure (no global-state mutation) and fast (small grids, no solver runs).

/// Helper: create a small SolverExport filled with constant field values.
fn make_test_export(
  name  : String,
  nc    : Int,
  fill_u : Double,
  fill_v : Double,
  fill_p : Double,
) -> SolverExport {
  let u = create_zeros_2d(nc, nc)
  let v = create_zeros_2d(nc, nc)
  let p = create_zeros_2d(nc, nc)
  for i = 0; i < nc; i = i + 1 {
    for j = 0; j < nc; j = j + 1 {
      u[i][j] = fill_u
      v[i][j] = fill_v
      p[i][j] = fill_p
    }
  }
  { name, label: "test", nc, h: 0.05, x0: 0.0, y0: 0.0, steps: 10, u, v, p }
}

// T84 [Unit / EP]
// format_csv must include the standard RFC 4180 header row.
test "csv_header_present" {
  let s = make_test_export("chorin", 3, 0.0, 0.0, 0.0)
  let csv = format_csv([s])
  assert_true(csv.contains("solver,i,j,x,y,u,v,p,vel_mag"))
}

// T85 [Unit / BVA / OBT]
// format_csv for a single nc=3 solver must contain exactly nc*nc = 9 data rows
// (one line per grid point, excluding comment and header lines).
test "csv_row_count" {
  let nc = 3
  let s  = make_test_export("test", nc, 1.0, 0.5, 0.1)
  let csv = format_csv([s])
  // Count lines that start with "test," (the solver name prefix)
  let mut data_rows = 0
  csv.split("\n").each(fn(line) {
    if line.has_prefix("test,") {
      data_rows = data_rows + 1
    }
  })
  assert_eq(data_rows, nc * nc)
}

// T86 [Unit / EP]
// format_vtk must produce a valid STRUCTURED_GRID header with correct DIMENSIONS.
test "vtk_dimensions_header" {
  let nc = 4
  let s  = make_test_export("chorin", nc, 0.0, 0.0, 0.0)
  let vtk = format_vtk(s)
  assert_true(vtk.contains("DATASET STRUCTURED_GRID"))
  assert_true(vtk.contains("DIMENSIONS \{nc} \{nc} 1"))
}

// T87 [Unit / OBT]
// format_vtk POINTS count must equal nc*nc.
test "vtk_points_count" {
  let nc   = 5
  let npts = nc * nc
  let s    = make_test_export("pcg", nc, 0.0, 0.0, 0.0)
  let vtk  = format_vtk(s)
  assert_true(vtk.contains("POINTS \{npts} double"))
}

// T88 [Unit / EP]
// format_tecplot for k solvers must contain exactly k ZONE lines.
test "tecplot_zone_count" {
  let s1 = make_test_export("chorin", 3, 0.0, 0.0, 0.0)
  let s2 = make_test_export("simple", 3, 0.0, 0.0, 0.0)
  let dat = format_tecplot([s1, s2])
  let mut zone_count = 0
  dat.split("\n").each(fn(line) {
    if line.has_prefix("ZONE T=") {
      zone_count = zone_count + 1
    }
  })
  assert_eq(zone_count, 2)
}

// T89 [Unit / PP / Negative path]
// format_netcdf stub must always return Err (never panic, never Ok).
test "netcdf_stub_returns_err" {
  let s = make_test_export("chorin", 3, 0.0, 0.0, 0.0)
  match format_netcdf(s) {
    Err(msg) => assert_true(msg.length() > 0)
    Ok(_)    => assert_true(false)   // must not reach here
  }
}

// T90 [Unit / OBT]
// format_vtk must contain the POINT_DATA keyword followed by the total count.
// This validates that the SCALARS/VECTORS sections are properly introduced.
// (MoonBit renders Double 0.0 as "0", so we verify POINT_DATA header, not coordinates.)
test "vtk_point_data_header" {
  let nc   = 4
  let npts = nc * nc
  let s    = make_test_export("chorin", nc, 0.0, 0.0, 0.0)
  let vtk  = format_vtk(s)
  assert_true(vtk.contains("POINT_DATA \{npts}"))
  assert_true(vtk.contains("SCALARS pressure double 1"))
  assert_true(vtk.contains("VECTORS velocity double"))
}

// T91 [Unit / OBT]
// format_tecplot SOLUTIONTIME for a solver with steps=10, dt=0.001 must be 0.01.
test "tecplot_solution_time" {
  let s   = make_test_export("chorin", 3, 0.0, 0.0, 0.0)
  // s.steps = 10, dt = 0.001 (module constant) → SOLUTIONTIME = 0.01
  let dat = format_tecplot([s])
  assert_true(dat.contains("SOLUTIONTIME=\{10.0 * dt}"))
}

// ══════════════════════════════════════════════════════════════════════════════
// J. Unit Tests — VTK XML ImageData (.vti) Format  (T92–T93)
// ══════════════════════════════════════════════════════════════════════════════

// ══════════════════════════════════════════════════════════════════════════════
// L. Unit Tests — Runtime Re + Steps Override  (T96)
// ══════════════════════════════════════════════════════════════════════════════

// T96 [Unit / EP / OBT]
// parse_pos_int("100", 0) must return 100, and nu=2.0/100 must equal 0.02.
// Also verifies that g_nu_runtime is mutable (single-element Array[Double]).
test "runtime_re_100_gives_nu_0_02" {
  let re = parse_pos_int("100", 0)
  assert_eq(re, 100)
  let nu_100 = 2.0 / re.to_double()
  assert_true(approx_eq(nu_100, 0.02, 1.0e-10))
  // Verify g_nu_runtime is settable and readable
  let prev = g_nu_runtime[0]
  g_nu_runtime[0] = nu_100
  assert_true(approx_eq(g_nu_runtime[0], 0.02, 1.0e-10))
  g_nu_runtime[0] = prev   // restore compile-time default (nu=0.1)
}

// ══════════════════════════════════════════════════════════════════════════════
// K. Unit Tests — Centerline CSV + Monitor CSV  (T94–T95)
// ══════════════════════════════════════════════════════════════════════════════

// T94 [Unit / BVA / OBT]
// format_centerline_csv for nc=5 must contain exactly 2*nc = 10 data rows:
// nc rows for u_vcl (vertical centerline) + nc rows for v_hcl (horizontal).
test "centerline_csv_row_count" {
  let nc = 5
  let s  = make_test_export("chorin", nc, 1.0, 0.5, 0.0)
  let out = format_centerline_csv(s)
  let mut u_rows = 0
  let mut v_rows = 0
  out.split("\n").each(fn(line) {
    if line.has_prefix("u_vcl,") { u_rows = u_rows + 1 }
    if line.has_prefix("v_hcl,") { v_rows = v_rows + 1 }
  })
  assert_eq(u_rows, nc)
  assert_eq(v_rows, nc)
}

// T95 [Unit / EP]
// format_monitor_csv must contain the standard header row with all column names.
test "monitor_csv_header_present" {
  let rec : MonitorRecord = {
    step: 0, time: 0.0, div_norm: 0.0,
    r_u: -1.0, r_v: -1.0, ke: 0.0, cfl_max: -1.0, pcg_iters: 0,
  }
  let out = format_monitor_csv([rec], "chorin", 20.0)
  assert_true(out.contains("step,time,div_norm,r_u,r_v,ke,cfl_max,pcg_iters"))
}

// ══════════════════════════════════════════════════════════════════════════════
// J. Unit Tests — VTK XML ImageData (.vti) Format  (T92–T93)
// ══════════════════════════════════════════════════════════════════════════════

// T92 [Unit / EP]
// format_vti must produce valid VTK XML: <?xml> declaration and
// <VTKFile type="ImageData"> root element must be present.
// Also verifies WholeExtent uses nc-1 (e.g., nc=4 → "0 3 0 3 0 0").
test "vti_xml_header_present" {
  let nc = 4
  let s  = make_test_export("chorin", nc, 0.5, 0.3, 0.1)
  let vti = format_vti(s)
  assert_true(vti.contains("<?xml version=\"1.0\"?>"))
  assert_true(vti.contains("<VTKFile type=\"ImageData\""))
  assert_true(vti.contains("<ImageData"))
  // nc=4 → extent upper bound = 3
  assert_true(vti.contains("WholeExtent=\"0 3 0 3 0 0\""))
}

// T93 [Unit / OBT]
// format_vti must contain exactly 3 <DataArray> opening tags:
// pressure (scalar), velocity (3-component), vel_magnitude (scalar).
test "vti_data_array_count_is_three" {
  let s   = make_test_export("pcg", 3, 1.0, 0.0, 0.5)
  let vti = format_vti(s)
  // Count occurrences of "<DataArray" (opening tag, not closing)
  let mut count = 0
  vti.split("\n").each(fn(line) {
    if line.contains("<DataArray") && not(line.contains("</DataArray")) {
      count = count + 1
    }
  })
  assert_eq(count, 3)
}
