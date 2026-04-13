# FlowLabLite

FlowLabLite is a lightweight, AI-native CFD solver for 2D lid-driven cavity flow,
written in [MoonBit](https://www.moonbitlang.com/).
It compiles to WebAssembly (wasm-gc) and ships an interactive browser-based
visualization page (`cmd/main/main.html`) with no runtime dependencies.

---

## Code Structure

```
FlowLabLite/
├── cmd/main/
│   ├── main.mbt          # Core CFD solver + exported WASM API (pub fn …)
│   ├── main_bench.mbt    # Micro-benchmark (fib baseline)
│   ├── main_wbtest.mbt   # White-box unit tests
│   ├── main.html         # Browser visualization (loads WASM, draws canvases)
│   ├── local_viewer.html # Browser viewer for locally computed JSON results
│   ├── moon.pkg.json     # Package config: imports, link/exports for both wasm targets
│   ├── moon.pkg          # Legacy pkg marker (is-main: true) – required by moon tool
│   └── now_clock.c       # Optional native timing shim (not needed for WASM builds)
├── lib/
│   └── moon.pkg.json     # Library package placeholder
├── moon.mod.json         # Module manifest (package name, version)
├── build_wasm.sh         # Build script that re-links WASM with -exported_functions
│                         # (workaround for moon 0.1.20260309 not passing exports to moonc)
├── run_local.sh          # Run solver locally + extract JSON results for viewer
└── README.md
```

### Key source file: `cmd/main/main.mbt`

| Section | Description |
|---|---|
| Constants (`nx`, `ny`, `nt`, …) | Grid dimensions 41×41, 500 time steps, Re = 20 |
| Global arrays (`g_u`, `g_v`, `g_p`, `g_b`) | Simulation state, allocated once at module init |
| `init_simulation()` | Zeroes all global arrays; JS calls this before running |
| `run_n_steps(n)` | Advances the solver by *n* time steps; used for progressive rendering |
| `run_all_steps()` | Runs all `nt` steps in one call |
| `get_u_at(i, j)` / `get_v_at` / `get_p_at` | Point access to velocity/pressure fields |
| `get_max_u()` … `get_min_p()` | Field statistics |
| `get_divergence_norm()` | Mean |div u| over interior cells (convergence indicator) |
| `cavity_flow_array()` | Navier-Stokes solver core (Chorin projection, finite-difference) |
| `pressure_poisson_array()` | Iterative pressure solve (Gauss-Seidel, `nit` iterations) |
| `output_json()` | Outputs full grid data as JSON to stdout (for local_viewer.html) |

### Browser visualization: `cmd/main/main.html`

A single-file, zero-dependency HTML/JS application.

- Fetches the WASM binary (tries release build, falls back to debug)
- Uses a `Proxy`-based import object to tolerate any MoonBit runtime imports
- Calls `init_simulation()` then `run_n_steps(50)` in a loop, yielding to the
  browser between chunks for progressive canvas updates
- Renders velocity magnitude (jet colormap) and pressure field (blue→orange→red)
  with bilinear interpolation to fill every canvas pixel from the 41×41 grid

---

## Prerequisites

| Tool | Minimum version | Install |
|---|---|---|
| `moon` (MoonBit build tool) | 0.1.20260309 | https://www.moonbitlang.com/download/ |
| `moonc` (MoonBit compiler) | 0.8.3 | bundled with moon |
| `node` (for export verification only) | 18+ | https://nodejs.org/ |
| Chrome / Edge / Brave | 115+ | required to run wasm-gc with stringref |
| Any HTTP server | — | see "Serving the page" below |

> **Note:** Node.js cannot load this WASM directly because it uses wasm-gc stringref
> (type code `0x77`), which requires a browser with full WasmGC support (Chrome 115+).
> Node.js is only used for verifying the export section.

---

## How to Build

### Step 1 – Compile MoonBit source to `.core` files

```bash
moon build --target wasm-gc          # debug build (includes source-map)
moon build --target wasm-gc --release  # release build (smaller, faster)
```

This generates `.core` intermediate files under `_build/wasm-gc/`.

> **Known issue in moon 0.1.20260309:** `moon build` alone does not export
> the 27 API functions — the compiled WASM only has `_start`.
> Use `build_wasm.sh` (step 2) to fix this.

### Step 2 – Re-link with exported functions (required)

```bash
bash build_wasm.sh          # debug WASM (with source-map, -O0)
bash build_wasm.sh release  # release WASM (default moonc optimisation)
```

This script calls `moon build` first, then re-runs `moonc link-core` with
`-exported_functions init_simulation,run_all_steps,...` so all 27 API
functions appear in the WASM export section.

**Output:**
```
_build/wasm-gc/release/build/cmd/main/main.wasm   ← used by main.html
_build/wasm-gc/debug/build/cmd/main/main.wasm     ← fallback
```

The script prints a confirmation line showing all 29 WASM exports
(27 API functions + `_start` + `memory`).

### Step 3 – Run the tests (optional)

```bash
moon test                     # all targets
moon test --target wasm-gc    # wasm-gc only
```

### Step 4 – Run benchmarks (optional)

```bash
moon bench
```

---

## How to Use – Local Run + Browser Viewer

This mode runs the solver **locally** (not in the browser) and produces a JSON
file with velocity/pressure results, which you then view in `local_viewer.html`.

### Step 1 – Run the simulation locally

```bash
# One-step script (builds, runs, extracts JSON):
bash run_local.sh                  # → results.json
bash run_local.sh my_results.json  # → custom filename

# Or manually:
moon run cmd/main --target wasm 2>&1 \
  | sed -n '/===JSON_DATA_START===/,/===JSON_DATA_END===/p' \
  | grep -v '===JSON_DATA' \
  > results.json
```

### Step 2 – Open the local viewer

Open `cmd/main/local_viewer.html` **directly in a browser** (no HTTP server needed —
it uses `FileReader`, not `fetch()`).

> Alternatively, serve with any HTTP server and navigate to the page.

### Step 3 – Load results

- **Drag & drop** `results.json` onto the page, or click the drop zone to browse.
- Or click **Paste JSON from clipboard**.
- The velocity magnitude field and pressure field render immediately.

### Expected output

| Field | Description |
|---|---|
| Left canvas | Velocity magnitude with jet colormap + vector arrows |
| Right canvas | Pressure field (blue→cyan→orange→red) |
| Info panel | Grid size, Re, time steps, max u/v/p, center-point data, divergence norm |

---

## How to Use – WASM in Browser (Interactive)

### Serving the page

The HTML file uses `fetch()` with relative paths to load the WASM binary,
so it **must be served over HTTP** (not opened as `file://`).

Serve the project root (not the `cmd/main/` subdirectory):

```bash
# Python (built-in)
python -m http.server 8080

# Node.js (npx)
npx serve .

# VS Code Live Server
# Right-click index or use "Open with Live Server" from the project root
```

Then open in Chrome 115+:

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

### Expected workflow

1. Open the page → click **Load WASM** → status shows "WebAssembly module loaded successfully"
2. Click **Run Full Simulation** → canvases update progressively as steps complete
3. Left canvas: velocity magnitude (blue = slow, red = fast)
4. Right canvas: pressure field (blue = low, red = high)
5. Parameter panel below shows Re, grid size, divergence norm, max velocity

### Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| "WASM file not found" | Build not run or wrong serve root | Run `bash build_wasm.sh release` then serve from project root |
| Canvases show mock rainbow patterns | WASM loaded but `init_simulation` not found | Rebuild with `build_wasm.sh` (old build missing exports) |
| Blank page / console error about `0x77` | Browser too old | Use Chrome 115+ or Edge 115+ |
| "Failed to load WebAssembly module" | CORS / file:// protocol | Must serve over HTTP, not file:// |

---

## Test Report

### Running tests

```bash
moon test --target wasm      # run all 16 tests on wasm target
moon test                    # all targets
```

### Test summary (16 tests, all passing)

| # | Test name | Category | Validates |
|---|---|---|---|
| 1 | `create_zeros_2d` | Array utility | Correct dimensions and all-zero init |
| 2 | `copy_2d_array` | Array utility | Deep copy, source independence |
| 3 | `generate_mesh_grid` | Mesh | Coordinate range [0, 2] × [0, 2], monotonicity |
| 4 | `init_simulation_resets_state` | State mgmt | Global arrays zeroed, step counter reset |
| 5 | `boundary_conditions_after_run` | Solver | Lid u=1, no-slip walls u=v=0 after 10 steps |
| 6 | `step_counter` | State mgmt | Increments correctly across multiple run_n_steps |
| 7 | `out_of_range_returns_zero` | Safety | get_*_at returns 0 for negative/overflow indices |
| 8 | `constant_accessors` | API | nx=41, ny=41, Re=20, rho=1, nu=0.1, dt=0.001 |
| 9 | `velocity_magnitude_consistency` | Solver | get_velocity_magnitude_at = sqrt(u²+v²) |
| 10 | `max_velocity_magnitude_bounds` | Solver | max_mag >= 0, max_mag >= max_u |
| 11 | `pressure_bounded` | Solver | Pressure finite and within ±1000 after 50 steps |
| 12 | `center_getters_consistent` | API | get_u/v/p_center = get_*_at(ny/2, nx/2) |
| 13 | `divergence_norm_nonnegative` | Solver | Mean |div u| >= 0 |
| 14 | `build_up_b_nonzero` | Solver | Pressure source term non-trivial with lid velocity |
| 15 | `full_simulation_produces_vortex` | Physics | After 500 steps: negative u at center (backflow), max_u >= 1 |
| 16 | `pressure_boundary_dp_zero` | Physics | dp/dy=0 at y=0, p=0 at y=2 satisfied |

### Physical validation

The lid-driven cavity at Re = 20 produces a primary clockwise vortex centered
slightly above and to the right of the geometric center. After 500 time steps:

- Center u-velocity is negative (≈ −0.06), confirming backflow
- Center v-velocity is near zero (order 10⁻⁵)
- Pressure boundary conditions (Neumann on walls, Dirichlet p=0 at top) are satisfied
- Mean divergence norm decreases to O(10⁻²), indicating reasonable mass conservation
  for the 50-iteration Gauss-Seidel pressure solver

---

## Simulation Parameters

| Parameter | Value | Description |
|---|---|---|
| `nx`, `ny` | 41 × 41 | Grid points in x and y |
| `nt` | 500 | Total time steps |
| `nit` | 50 | Pressure-Poisson iterations per step |
| `dx`, `dy` | 0.05 | Grid spacing (domain = 2 × 2) |
| `dt` | 0.001 | Time step |
| `rho` | 1.0 | Fluid density |
| `nu` | 0.1 | Kinematic viscosity |
| Re | 20 | Reynolds number (u_lid × L / nu = 1 × 2 / 0.1) |

---

## License

Apache-2.0 — see [LICENSE](LICENSE).

Developed by Fenfen Yu (余芬芬) in collaboration with
Beihang University (北京航空航天大学) and
Ezhou Hi-Modeling Technology Co., Ltd. (鄂州海慕科技有限公司).
