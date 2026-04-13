// White-box tests for the cavity flow solver.
// Run: moon test --target wasm

// ── Helper: approximate floating-point equality ──
fn approx_eq(a : Double, b : Double, tol : Double) -> Bool {
  (a - b).abs() < tol
}

// ── Test 1: create_zeros_2d returns all-zero array of correct size ──
test "create_zeros_2d" {
  let arr = create_zeros_2d(5, 7)
  assert_eq(arr.length(), 5)
  assert_eq(arr[0].length(), 7)
  for i = 0; i < 5; i = i + 1 {
    for j = 0; j < 7; j = j + 1 {
      assert_eq(arr[i][j], 0.0)
    }
  }
}

// ── Test 2: copy_2d_array copies values correctly ──
test "copy_2d_array" {
  let src = create_zeros_2d(3, 3)
  src[0][0] = 1.0
  src[1][2] = 5.5
  src[2][1] = -3.0
  let dst = create_zeros_2d(3, 3)
  copy_2d_array(src, dst)
  assert_eq(dst[0][0], 1.0)
  assert_eq(dst[1][2], 5.5)
  assert_eq(dst[2][1], -3.0)
  // Verify it's a deep copy: modifying src should not affect dst
  src[0][0] = 99.0
  assert_eq(dst[0][0], 1.0)
}

// ── Test 3: generate_mesh_grid produces correct coordinate ranges ──
test "generate_mesh_grid" {
  let (x, y) = generate_mesh_grid()
  assert_eq(x.length(), ny)
  assert_eq(x[0].length(), nx)
  // Bottom-left corner should be (0, 0)
  assert_true(approx_eq(x[0][0], 0.0, 1.0e-10))
  assert_true(approx_eq(y[0][0], 0.0, 1.0e-10))
  // Top-right corner should be (2, 2)
  assert_true(approx_eq(x[ny - 1][nx - 1], 2.0, 1.0e-10))
  assert_true(approx_eq(y[ny - 1][nx - 1], 2.0, 1.0e-10))
  // X should increase along columns, Y along rows
  assert_true(x[0][1] > x[0][0])
  assert_true(y[1][0] > y[0][0])
}

// ── Test 4: init_simulation resets global state ──
test "init_simulation_resets_state" {
  // Dirty the global arrays
  g_u[ny / 2][nx / 2] = 42.0
  g_v[0][0] = -1.0
  g_p[ny - 1][nx - 1] = 99.0
  g_steps[0] = 999
  // Reset
  init_simulation()
  assert_eq(g_steps[0], 0)
  assert_eq(g_u[ny / 2][nx / 2], 0.0)
  assert_eq(g_v[0][0], 0.0)
  assert_eq(g_p[ny - 1][nx - 1], 0.0)
}

// ── Test 5: boundary conditions after solver run ──
test "boundary_conditions_after_run" {
  init_simulation()
  run_n_steps(10)
  // Top lid: u = 1
  for j = 0; j < nx; j = j + 1 {
    assert_true(approx_eq(get_u_at(ny - 1, j), 1.0, 1.0e-10))
  }
  // Bottom wall: u = 0, v = 0
  for j = 0; j < nx; j = j + 1 {
    assert_true(approx_eq(get_u_at(0, j), 0.0, 1.0e-10))
    assert_true(approx_eq(get_v_at(0, j), 0.0, 1.0e-10))
  }
  // Left wall: u = 0, v = 0 (exclude top lid row i=ny-1 where u=1)
  for i = 0; i < ny - 1; i = i + 1 {
    assert_true(approx_eq(get_u_at(i, 0), 0.0, 1.0e-10))
    assert_true(approx_eq(get_v_at(i, 0), 0.0, 1.0e-10))
  }
  // Right wall: u = 0, v = 0 (exclude top lid row)
  for i = 0; i < ny - 1; i = i + 1 {
    assert_true(approx_eq(get_u_at(i, nx - 1), 0.0, 1.0e-10))
    assert_true(approx_eq(get_v_at(i, nx - 1), 0.0, 1.0e-10))
  }
  // Top/bottom wall: v = 0
  for j = 0; j < nx; j = j + 1 {
    assert_true(approx_eq(get_v_at(ny - 1, j), 0.0, 1.0e-10))
  }
}

// ── Test 6: step counter increments correctly ──
test "step_counter" {
  init_simulation()
  assert_eq(get_step_count(), 0)
  run_n_steps(5)
  assert_eq(get_step_count(), 5)
  run_n_steps(3)
  assert_eq(get_step_count(), 8)
}

// ── Test 7: get_*_at returns 0 for out-of-range indices ──
test "out_of_range_returns_zero" {
  assert_eq(get_u_at(-1, 0), 0.0)
  assert_eq(get_u_at(0, -1), 0.0)
  assert_eq(get_u_at(ny, 0), 0.0)
  assert_eq(get_u_at(0, nx), 0.0)
  assert_eq(get_v_at(-1, -1), 0.0)
  assert_eq(get_p_at(999, 999), 0.0)
  assert_eq(get_velocity_magnitude_at(-1, 0), 0.0)
}

// ── Test 8: constant accessors return correct values ──
test "constant_accessors" {
  assert_eq(get_nx(), 41)
  assert_eq(get_ny(), 41)
  assert_eq(get_nt(), 500)
  assert_eq(get_nit(), 50)
  assert_true(approx_eq(get_rho(), 1.0, 1.0e-10))
  assert_true(approx_eq(get_nu(), 0.1, 1.0e-10))
  assert_true(approx_eq(get_dt(), 0.001, 1.0e-10))
  assert_true(approx_eq(get_dx(), 0.05, 1.0e-10))
  assert_true(approx_eq(get_dy(), 0.05, 1.0e-10))
  assert_true(approx_eq(get_re(), 20.0, 1.0e-10))
}

// ── Test 9: velocity magnitude consistency ──
test "velocity_magnitude_consistency" {
  init_simulation()
  run_n_steps(20)
  let mid_i = ny / 2
  let mid_j = nx / 2
  let u = get_u_at(mid_i, mid_j)
  let v = get_v_at(mid_i, mid_j)
  let expected_mag = (u * u + v * v).sqrt()
  let actual_mag = get_velocity_magnitude_at(mid_i, mid_j)
  assert_true(approx_eq(actual_mag, expected_mag, 1.0e-12))
}

// ── Test 10: max velocity magnitude is non-negative and >= max_u ──
test "max_velocity_magnitude_bounds" {
  init_simulation()
  run_n_steps(50)
  let max_mag = get_max_velocity_magnitude()
  let max_u = get_max_u()
  assert_true(max_mag >= 0.0)
  assert_true(max_mag >= max_u)
}

// ── Test 11: pressure Poisson solver produces bounded pressure ──
test "pressure_bounded" {
  init_simulation()
  run_n_steps(50)
  let max_p = get_max_p()
  let min_p = get_min_p()
  // Pressure should be finite and within reasonable physical bounds
  assert_true(max_p < 1000.0)
  assert_true(min_p > -1000.0)
  assert_true(max_p >= min_p)
}

// ── Test 12: center point getters match get_*_at(ny/2, nx/2) ──
test "center_getters_consistent" {
  init_simulation()
  run_n_steps(30)
  let mid_i = ny / 2
  let mid_j = nx / 2
  assert_true(approx_eq(get_u_center(), get_u_at(mid_i, mid_j), 1.0e-15))
  assert_true(approx_eq(get_v_center(), get_v_at(mid_i, mid_j), 1.0e-15))
  assert_true(approx_eq(get_p_center(), get_p_at(mid_i, mid_j), 1.0e-15))
}

// ── Test 13: divergence norm is non-negative after simulation ──
test "divergence_norm_nonnegative" {
  init_simulation()
  run_n_steps(50)
  let div = get_divergence_norm()
  assert_true(div >= 0.0)
}

// ── Test 14: build_up_b_array produces non-zero source term with non-trivial input ──
test "build_up_b_nonzero" {
  let u = create_zeros_2d(ny, nx)
  let v = create_zeros_2d(ny, nx)
  // Set top lid velocity
  for j = 0; j < nx; j = j + 1 {
    u[ny - 1][j] = 1.0
  }
  // Add some interior velocity to produce non-trivial b
  u[ny / 2][nx / 2] = 0.5
  let b = create_zeros_2d(ny, nx)
  build_up_b_array(b, rho, dt, u, v, dx, dy) |> ignore
  // b should have at least some non-zero entries in the interior
  let mut has_nonzero = false
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      if b[i][j].abs() > 1.0e-15 { has_nonzero = true }
    }
  }
  assert_true(has_nonzero)
}

// ── Test 15: run_all_steps produces lid-driven vortex ──
test "full_simulation_produces_vortex" {
  init_simulation()
  run_all_steps()
  assert_eq(get_step_count(), 500)
  // After 500 steps, the interior should have non-zero velocities
  let u_center = get_u_center()
  let v_center = get_v_center()
  // Center u should be negative (backflow) for lid-driven cavity
  assert_true(u_center < 0.0)
  // v_center should be small but non-zero
  assert_true(v_center.abs() > 0.0)
  // Max velocity should be >= 1.0 (lid velocity)
  assert_true(get_max_u() >= 1.0 - 1.0e-10)
}

// ── Test 16: pressure boundary conditions satisfied ──
test "pressure_boundary_dp_zero" {
  init_simulation()
  run_n_steps(50)
  // dp/dy = 0 at y = 0: p[0,:] = p[1,:] (approximately)
  for j = 0; j < nx; j = j + 1 {
    assert_true(approx_eq(get_p_at(0, j), get_p_at(1, j), 1.0e-10))
  }
  // p = 0 at y = 2 (top): p[ny-1,:] = 0
  for j = 0; j < nx; j = j + 1 {
    assert_true(approx_eq(get_p_at(ny - 1, j), 0.0, 1.0e-10))
  }
}
