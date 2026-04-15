# FlowLabLite

FlowLabLite is a lightweight, AI-native CFD solver for 2D lid-driven cavity flow,
written in [MoonBit](https://www.moonbitlang.com/).
It compiles to WebAssembly (wasm-gc) and ships an interactive browser-based
visualization page (`cmd/main/main.html`) with no runtime dependencies.

Two independent solvers are provided:

| Solver | Algorithm | Use case |
|---|---|---|
| **Chorin** | Projection method, explicit time-marching | Transient simulation, 500 time steps |
| **SIMPLE** | Semi-Implicit Method for Pressure-Linked Equations | Steady-state iteration with under-relaxation |

---

## Code Structure

```
FlowLabLite/
├── cmd/main/
│   ├── main.mbt          # Core CFD solver + exported WASM API (pub fn …)
│   ├── main_bench.mbt    # Micro-benchmark (fib baseline)
│   ├── main_wbtest.mbt   # White-box unit tests (26 tests)
│   ├── main.html         # Browser visualization (loads WASM, draws canvases)
│   ├── local_viewer.html # Browser viewer for locally computed JSON results
│   ├── moon.pkg.json     # Package config: imports, link/exports for both wasm targets
│   ├── moon.pkg          # Legacy pkg marker (is-main: true) – required by moon tool
│   └── now_clock.c       # Optional native timing shim (not needed for WASM builds)
├── lib/
│   └── moon.pkg.json     # Library package placeholder
├── moon.mod.json         # Module manifest (package name, version)
├── build_wasm.sh         # Build script that re-links WASM with -exported_functions
├── run_local.sh          # Run solver locally + extract JSON results for viewer
└── README.md
```

### Key source file: `cmd/main/main.mbt`

| Section | Description |
|---|---|
| Constants (`nx`, `ny`, `nt`, …) | Grid dimensions 41×41, 500 time steps, Re = 20 |
| Global arrays (`g_u`, `g_v`, `g_p`) | Chorin simulation state, allocated once |
| `init_simulation()` | Zeroes Chorin state; JS calls this before running |
| `run_n_steps(n)` | Advances Chorin solver by *n* time steps; used for progressive rendering |
| `run_all_steps()` | Runs all `nt` steps in one call |
| `get_u_at(i,j)` / `get_v_at` / `get_p_at` | Point access to Chorin velocity/pressure fields |
| `get_max_u()` … `get_min_p()` | Chorin field statistics |
| `get_divergence_norm()` | Mean \|div u\| over interior cells (convergence indicator) |
| `cavity_flow_array()` | Navier-Stokes solver core (Chorin projection, finite-difference) |
| `pressure_poisson_array()` | Iterative pressure solve (Gauss-Seidel, `nit` iterations) |
| `output_json()` | Outputs full grid data as JSON to stdout (for local_viewer.html) |
| **SIMPLE constants** | `simple_alpha_p = 0.3`, `simple_alpha_u = 0.7` (under-relaxation) |
| **SIMPLE state** (`g_u_s`, `g_v_s`, `g_p_s`) | Independent 41×41 arrays for SIMPLE solver |
| `init_simple()` | Zeroes SIMPLE state; must be called before `run_simple_n_iter` |
| `run_simple_n_iter(n)` | Runs *n* SIMPLE iterations; caches convergence residual |
| `get_simple_step_count()` | Number of SIMPLE iterations completed |
| `get_simple_residual()` | Latest continuity residual ‖div u‖ (convergence indicator) |
| `get_u_simple_at(i,j)` / `get_v_simple_at` / `get_p_simple_at` | Point access to SIMPLE fields |
| `get_max_u_simple()` | Maximum \|u\| in SIMPLE field |
| `get_simple_divergence_norm()` | Mean \|div u\| over SIMPLE interior cells |
| `simple_one_iter()` | One SIMPLE sweep: momentum predictor → pressure correction → update |

---

## Prerequisites

| Tool | Minimum version | Install |
|---|---|---|
| `moon` (MoonBit build tool) | 0.1.20260309 | https://www.moonbitlang.com/download/ |
| `moonc` (MoonBit compiler) | 0.8.3 | bundled with moon |
| `node` (for export verification only) | 18+ | https://nodejs.org/ |
| Chrome / Edge / Brave | 115+ | required to run wasm-gc with stringref |
| Any HTTP server | — | see "Serving the page" below |

---

## Quick Start — From Source to Browser

This section shows the complete workflow from a clean checkout to viewing results
in the browser.

### 1. Install dependencies

```bash
# MoonBit toolchain (moon + moonc)
# Download and install from https://www.moonbitlang.com/download/
# After install, verify:
moon version
moonc -v

# Node.js (only needed for WASM export verification in build_wasm.sh)
node --version   # 18+ required
```

### 2. Clone and enter the project

```bash
git clone https://github.com/YuFenfen-ai/FlowLabLite.git
cd FlowLabLite
```

### 3. Run the tests

```bash
moon test --target wasm    # runs all 26 tests (Chorin + SIMPLE)
```

Expected output:
```
Total tests: 26, passed: 26, failed: 0.
```

### 4. Build the WASM binary

```bash
bash build_wasm.sh release    # optimised build  → _build/wasm-gc/release/…/main.wasm
# or
bash build_wasm.sh            # debug build       → _build/wasm-gc/debug/…/main.wasm
```

The script prints the export list (38 functions + `_start` + `memory` = 40 symbols).

### 5. Serve over HTTP

The browser uses `fetch()` to load the WASM, so the page must be served over HTTP
(not opened as `file://`). Serve from the **project root**:

```bash
# Python (built-in, no install needed)
python -m http.server 8080

# Node.js
npx serve .

# VS Code Live Server extension — right-click project root → "Open with Live Server"
```

### 6. Open in Chrome 115+

```
http://localhost:8080/cmd/main/main.html
```

### 7. Run the simulation

1. Click **Load WASM** — status shows "WebAssembly module loaded successfully"
2. Click **Run Full Simulation** (or press `F`) — canvases update progressively
3. To run SIMPLE in the browser, open the browser DevTools console and type:

```javascript
// SIMPLE solver API (all functions are on the wasm instance)
wasm.init_simple();
// Run 200 SIMPLE iterations in chunks of 50 for progressive updates
for (let i = 0; i < 4; i++) {
  wasm.run_simple_n_iter(50);
  console.log('iter', wasm.get_simple_step_count(),
              'residual', wasm.get_simple_residual());
}
```

---

## How to Build

### Step 1 – Compile MoonBit source to `.core` files

```bash
moon build --target wasm-gc          # debug build (includes source-map)
moon build --target wasm-gc --release  # release build (smaller, faster)
```

> **Known issue in moon 0.1.20260309:** `moon build` alone does not export
> the API functions — the compiled WASM only has `_start`.
> Use `build_wasm.sh` (step 2) to fix this.

### Step 2 – Re-link with exported functions (required)

```bash
bash build_wasm.sh          # debug WASM
bash build_wasm.sh release  # release WASM
```

**Output:**
```
_build/wasm-gc/release/build/cmd/main/main.wasm   ← used by main.html
_build/wasm-gc/debug/build/cmd/main/main.wasm     ← fallback
```

### Step 3 – Run the tests

```bash
moon test                     # all targets
moon test --target wasm-gc    # wasm-gc only
moon test --target wasm       # wasm (Wasmtime) only
```

### Step 4 – Run benchmarks (optional)

```bash
moon bench
```

---

## How to Use – Local Run + Browser Viewer

Run the solver locally and view results in `local_viewer.html`.

### Step 1 – Run the simulation locally

```bash
bash run_local.sh                  # → results.json
bash run_local.sh my_results.json  # → custom filename
```

### Step 2 – Open the local viewer

Open `cmd/main/local_viewer.html` directly in a browser (no HTTP server needed —
it uses `FileReader`, not `fetch()`).

### Step 3 – Load results

- **Drag & drop** `results.json` onto the page, or click the drop zone to browse.
- Or click **Paste JSON from clipboard**.

---

## How to Use – WASM in Browser (Interactive)

### Serving the page

```bash
python -m http.server 8080
```

Then open:

```
http://localhost:8080/cmd/main/main.html
```

### Controls

| Action | Button / Key |
|---|---|
| Load WebAssembly module | **Load WASM** button |
| Run full simulation (500 steps, progressive) | **Run Full Simulation** / `F` |
| Run quick test (50 steps) | **Quick Test** / `Q` |
| Clear canvases | **Clear** / `C` |
| Download grid data as JSON | **Download Data** / `D` |

### Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| "WASM file not found" | Build not run or wrong serve root | Run `bash build_wasm.sh release` then serve from project root |
| Canvases show mock rainbow patterns | WASM loaded but `init_simulation` not found | Rebuild with `build_wasm.sh` |
| Blank page / console error about `0x77` | Browser too old | Use Chrome 115+ or Edge 115+ |
| "Failed to load WebAssembly module" | CORS / file:// protocol | Must serve over HTTP, not file:// |

---

## SIMPLE Algorithm

### Background

The **SIMPLE** (Semi-Implicit Method for Pressure-Linked Equations, Patankar & Spalding, 1972)
algorithm solves the steady incompressible Navier-Stokes equations iteratively.
Unlike the Chorin projection method (which time-marches to equilibrium), SIMPLE
uses a pressure-correction loop with under-relaxation to converge directly to
the steady state.

### Algorithm (per iteration)

```
1. Momentum predictor
   Solve u*, v* explicitly using current pressure p*:
     u* = u - (convection) - (1/ρ)·∂p*/∂x + ν∇²u    (with α_u = 0.7 under-relaxation)
     v* = v - (convection) - (1/ρ)·∂p*/∂y + ν∇²v

2. Pressure-correction source
   b[i,j] = ρ/dt · (∂u*/∂x + ∂v*/∂y)   (same form as Chorin pressure RHS)

3. Pressure-correction Poisson solve
   ∇²p' = b    (Gauss-Seidel, 50 inner iterations)

4. Field update
   p  ← p* + α_p · p'          (α_p = 0.3)
   u  ← u* − (dt/ρ)·∂p'/∂x
   v  ← v* − (dt/ρ)·∂p'/∂y

5. Boundary conditions
   Lid: u[top,:] = 1,  walls: u=v=0
   Pressure: Neumann on walls, p=0 at top

6. Convergence check: residual = mean|∂u/∂x + ∂v/∂y|
```

### Parameters

| Parameter | Value | Purpose |
|---|---|---|
| `simple_alpha_p` | 0.3 | Pressure under-relaxation (stabilises p-correction) |
| `simple_alpha_u` | 0.7 | Velocity under-relaxation (prevents oscillation) |
| Inner iterations | 50 (nit) | Gauss-Seidel sweeps per pressure solve |

### WASM API (11 new exports)

| Function | Description |
|---|---|
| `init_simple()` | Reset SIMPLE state (u_s, v_s, p_s = 0, counter = 0) |
| `run_simple_n_iter(n)` | Advance by n SIMPLE iterations |
| `get_simple_step_count()` | Number of iterations completed |
| `get_simple_residual()` | Last cached ‖div u‖ residual |
| `get_u_simple_at(i,j)` | u-velocity at grid point (i,j) |
| `get_v_simple_at(i,j)` | v-velocity at grid point (i,j) |
| `get_p_simple_at(i,j)` | pressure at grid point (i,j) |
| `get_max_u_simple()` | Maximum \|u\| in SIMPLE field |
| `get_simple_divergence_norm()` | Mean \|∂u/∂x + ∂v/∂y\| over interior |

---

## Test Report

### Running tests

```bash
moon test --target wasm      # all 26 tests on wasm target
moon test                    # all targets
```

### Test summary (26 tests, all passing)

#### Chorin solver tests (1–16)

| # | Test name | Category | Validates |
|---|---|---|---|
| 1 | `create_zeros_2d` | Array utility | Correct dimensions and all-zero init |
| 2 | `copy_2d_array` | Array utility | Deep copy, source independence |
| 3 | `generate_mesh_grid` | Mesh | Coordinate range [0,2]×[0,2], monotonicity |
| 4 | `init_simulation_resets_state` | State mgmt | Global arrays zeroed, step counter reset |
| 5 | `boundary_conditions_after_run` | Solver | Lid u=1, no-slip walls after 10 steps |
| 6 | `step_counter` | State mgmt | Increments correctly across multiple run_n_steps |
| 7 | `out_of_range_returns_zero` | Safety | get_*_at returns 0 for out-of-bounds indices |
| 8 | `constant_accessors` | API | nx=41, ny=41, Re=20, rho=1, nu=0.1, dt=0.001 |
| 9 | `velocity_magnitude_consistency` | Solver | get_velocity_magnitude_at = sqrt(u²+v²) |
| 10 | `max_velocity_magnitude_bounds` | Solver | max_mag >= 0, max_mag >= max_u |
| 11 | `pressure_bounded` | Solver | Pressure finite within ±1000 after 50 steps |
| 12 | `center_getters_consistent` | API | get_u/v/p_center = get_*_at(ny/2, nx/2) |
| 13 | `divergence_norm_nonnegative` | Solver | Mean \|div u\| >= 0 |
| 14 | `build_up_b_nonzero` | Solver | Pressure source term non-trivial with lid velocity |
| 15 | `full_simulation_produces_vortex` | Physics | After 500 steps: negative u at centre, max_u >= 1 |
| 16 | `pressure_boundary_dp_zero` | Physics | dp/dy=0 at y=0, p=0 at y=2 |

#### SIMPLE solver tests (17–26)

| # | Test name | Category | Validates |
|---|---|---|---|
| 17 | `simple_init_resets_state` | State mgmt | All SIMPLE fields zeroed, counter reset |
| 18 | `simple_step_counter` | State mgmt | Counter increments correctly across batches |
| 19 | `simple_boundary_conditions` | Solver | Lid u=1, no-slip walls u=v=0, v=0 at top after 20 iters |
| 20 | `simple_out_of_range_returns_zero` | Safety | get_*_simple_at returns 0 for out-of-bounds |
| 21 | `simple_residual_nonnegative` | Solver | get_simple_residual() >= 0, divergence_norm >= 0 |
| 22 | `simple_produces_flow` | Physics | Non-zero interior u, v after 50 iterations |
| 23 | `simple_max_u_bounds` | Solver | 0 <= max_u_simple <= 2 (under-relaxation keeps field bounded) |
| 24 | `simple_pressure_boundary` | Physics | p=0 at top, p[0,:]=p[1,:] at bottom (Neumann) |
| 25 | `simple_vs_chorin_qualitative` | Physics | Both solvers show negative u at centre (backflow) after 200 iters |
| 26 | `simple_independent_of_chorin` | Isolation | SIMPLE iterations do not modify Chorin global state |

### Physical validation

**Chorin solver (Re = 20, 500 time steps):**

- Centre u-velocity ≈ −0.06 (negative = backflow, confirms clockwise primary vortex)
- Max u-velocity ≥ 1.0 (lid velocity maintained)
- Pressure Neumann BC at walls and Dirichlet p=0 at top satisfied
- Mean divergence norm decreases to O(10⁻²)

**SIMPLE solver (Re = 20, 200–500 iterations):**

- Interior velocities are non-zero after 50 iterations (lid-driven flow established)
- Both u and v are non-trivial in interior cells
- Negative u at centre after 200 iterations confirms vortex formation
- Pressure boundary conditions (p=0 at top, Neumann at walls) satisfied at each iteration
- Convergence residual ‖div u‖ is non-negative and decreases with iterations
- Under-relaxation keeps max \|u\| bounded within [0, 2] (lid = 1, alpha_u = 0.7)

---

## Simulation Parameters

| Parameter | Value | Description |
|---|---|---|
| `nx`, `ny` | 41 × 41 | Grid points in x and y |
| `nt` | 500 | Total Chorin time steps |
| `nit` | 50 | Pressure-Poisson inner iterations per step/SIMPLE iter |
| `dx`, `dy` | 0.05 | Grid spacing (domain = 2 × 2) |
| `dt` | 0.001 | Time step (Chorin) / pseudo-time step (SIMPLE) |
| `rho` | 1.0 | Fluid density |
| `nu` | 0.1 | Kinematic viscosity |
| Re | 20 | Reynolds number (u_lid × L / nu = 1 × 2 / 0.1) |
| `simple_alpha_p` | 0.3 | SIMPLE pressure under-relaxation factor |
| `simple_alpha_u` | 0.7 | SIMPLE velocity under-relaxation factor |

---

## License

Apache-2.0 — see [LICENSE](LICENSE).

Developed by Fenfen Yu (余芬芬) in collaboration with
Beihang University (北京航空航天大学) and
Ezhou Hi-Modeling Technology Co., Ltd. (鄂州海慕科技有限公司).
