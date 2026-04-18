# FlowLabLite — todo3 实施计划

> 前置状态：todo2 全部完成（Steps A-H, M），108 个测试全部通过（2026-04-18）。
> 理论公式见 `docs/theory.md`，架构决策见 `decides.md`。
> 每步实施后必须运行 `moon test --target wasm` 确认 108 个测试继续通过。
> 新增测试从 T109 开始编号。

---

## 步骤总览

| 步骤 | 内容 | 难度 | 依赖 | 优先级 | 状态 |
|---|---|---|---|---|---|
| **I** | TVD Van Leer / Superbee 对流格式 | 高 | Step E | P1 | 待实施 |
| **J** | RK3 / 隐式 Euler 时间推进 | 高 | Step I | P1 | 待实施 |
| **2a** | 等值线数值标注（HTML Canvas）| 中 | 无 | P1 | 待实施 |
| **2b** | 图例单位显示 | 低 | 无 | P1 | 待实施 |
| **3** | main.html ↔ local_viewer.html 互链 | 低 | 无 | P1 | 待实施 |
| **4** | 颜色主题（Dark/Light/Transparent/Print）| 低 | 无 | P1 | 待实施 |
| **1** | FVM 基础求解器（5th SolverExport）| 高 | 无 | P2 | 待实施 |
| **K** | 运行时 nx/ny（动态网格尺寸）| 极高 | 无 | P2 | 待实施 |
| **6** | 被动标量输运 | 高 | 速度场 | P2 | 待实施 |
| **7** | 能量方程 + Boussinesq 浮力 | 高 | Step 6 | P2 | 待实施 |
| **8** | BC 库扩展（Neumann/对称/周期）| 高 | Step 1 | P2 | 待实施 |
| **5** | 3D 顶盖驱动流（独立包）| 极高 | Step K | P3 | 待实施 |
| **9** | 拉伸网格（tanh 压缩）| 中高 | 无 | P3 | 待实施 |
| **10** | 曲线网格基础（ξ-η 变换）| 极高 | Step 9 | P3 | 待实施 |

---

## Step I — TVD Van Leer / Superbee 对流格式

**理论：** `docs/theory.md` §2

**实现位置：** `cmd/main/schemes.mbt`（追加）

### 新增函数

```moonbit
// 比例因子（返回 Option 以处理分母为零）
fn tvd_ratio(phi_up2: Double, phi_up1: Double, phi_down1: Double) -> Double

// 限制函数
fn limiter_van_leer(r: Double) -> Double  // (r + |r|) / (1 + |r|)
fn limiter_superbee(r: Double) -> Double  // max(0, min(2r,1), min(r,2))

// TVD 面通量（东/西/北/南）
fn tvd_east(phi: Array[Array[Double]], u_vel: Array[Array[Double]],
            i: Int, j: Int, limiter: (Double) -> Double) -> Double
fn tvd_west(phi, u_vel, i, j, limiter) -> Double
fn tvd_north(phi, v_vel, i, j, limiter) -> Double
fn tvd_south(phi, v_vel, i, j, limiter) -> Double

// TVD 对流项（全场）
fn tvd_convection_x(phi, u_vel, nr, nc, dx, dt, limiter) -> Array[Array[Double]]
fn tvd_convection_y(phi, v_vel, nr, nc, dy, dt, limiter) -> Array[Array[Double]]
```

### CLI 集成

```
--scheme van_leer  → 使用 Van Leer 限制函数
--scheme superbee  → 使用 Superbee 限制函数
--scheme quick     → 已实现（Step E）
--scheme upwind    → 默认（已实现）
```

### 测试（T109–T113）

| 编号 | 测试名 | 验证内容 |
|---|---|---|
| T109 | `van_leer_r_zero` | r=0 → ψ=0（一阶迎风极限）|
| T110 | `van_leer_r_one` | r=1 → ψ=1（二阶中心极限）|
| T111 | `superbee_r_gt2` | r≥2 → ψ=2（TVD 上界）|
| T112 | `tvd_uniform_field_zero` | 常数场 → TVD 对流项 = 0 |
| T113 | `tvd_better_than_upwind_sine` | TVD 误差 < 迎风误差（sine profile）|

**验证：** 与 Ghia Re=100 对比（TVD 误差应 ≤ 迎风误差）；与 QUICK 对比（TVD 无振荡）。

**提交：** `feat(schemes): add TVD Van Leer/Superbee convection limiters (Step I)`

---

## Step J — RK3 / 隐式 Euler 时间推进

**理论：** `docs/theory.md` §3

**实现位置：** `cmd/main/main.mbt`（追加新求解器函数）

### 设计决策（D-J1）

**RK3 实现策略：** 基于现有 Chorin 求解器，抽取右端项计算为独立函数：

```moonbit
// 单步右端项（对流 + 粘性，不含压力梯度）
fn rhs_u(u, v, nu_, dt, dx, dy, nr, nc) -> Array[Array[Double]]
fn rhs_v(u, v, nu_, dt, dx, dy, nr, nc) -> Array[Array[Double]]

// RK3 子步（包含压力投影）
fn chorin_rk3_step(u, v, p, nu_, rho, dt, dx, dy, nr, nc, nit_) -> Unit
```

**隐式 Euler：** 每时间步用 Gauss-Seidel 迭代求解线性化动量方程：

```moonbit
fn chorin_implicit_euler_step(u, v, p, nu_, rho, dt, dx, dy, nr, nc, nit_) -> Unit
```

### 全局状态（新增）

```moonbit
let g_u_rk3 : Array[Array[Double]] = ...   // RK3 求解器状态
let g_v_rk3 : Array[Array[Double]] = ...
let g_p_rk3 : Array[Array[Double]] = ...
```

### CLI 集成

```
--time-scheme rk3        → 使用 RK3 时间推进
--time-scheme implicit   → 使用隐式 Euler
--time-scheme euler      → 默认（显式 Euler，已实现）
```

### 测试（T114–T117）

| 编号 | 测试名 | 验证内容 |
|---|---|---|
| T114 | `rk3_state_independent` | RK3 状态不影响 Chorin/SIMPLE/PCG/MAC |
| T115 | `rk3_conservation_better` | RK3 散度范数 < 显式 Euler（相同步数）|
| T116 | `implicit_stable_large_dt` | dt=0.01 隐式 Euler 稳定；显式 Euler 发散 |
| T117 | `rk3_ghia_match` | Re=100, 2000步, RK3 Ghia 误差 < 迎风结果 |

**提交：** `feat(solver): add RK3 and implicit Euler time integration (Step J)`

---

## Step 2a — 等值线数值标注（HTML Canvas）

**实现位置：** `cmd/main/main.html`

### 功能描述

在压力等值线图中，在等值线上每隔 N 个点标注数值（Canvas `fillText`）：

```javascript
function drawContourLabels(ctx, contourData, colorScale, labelInterval) {
  // 在等值线路径上等间距放置数值标签
  // 标签背景：半透明白色矩形（避免与等值线颜色混淆）
  // 字体：monospace 11px
}
```

### 等值线算法（Marching Squares 简化版）

从当前网格数据计算等值线路径，在 Canvas 上绘制并标注。
可选：使用 wasm 导出的 `get_p_at(i, j)` 读取压力值，在 JS 中计算等值线。

---

## Step 2b — 图例单位显示

**实现位置：** `cmd/main/main.html`、`cmd/main/local_viewer.html`

```javascript
const fieldUnits = {
  pressure: '[Pa]',
  u: '[m/s]',
  v: '[m/s]',
  magnitude: '[m/s]',
  temperature: '[°C]',
  scalar: '[-]'
}
// 在颜色条标题旁追加单位字符串
```

---

## Step 3 — main.html ↔ local_viewer.html 互链

**实现位置：** `cmd/main/main.html`、`cmd/main/local_viewer.html`

```html
<!-- main.html header 添加 -->
<a href="local_viewer.html" class="nav-link">📂 Viewer</a>

<!-- local_viewer.html header 添加 -->
<a href="main.html" class="nav-link">🔴 Live Solver</a>
```

CSS 样式：`.nav-link { color: #aaf; text-decoration: none; margin: 0 8px; }`

---

## Step 4 — 颜色主题

**实现位置：** `cmd/main/main.html`、`cmd/main/local_viewer.html`

### 四种主题

```css
:root { --c-bg: #1a1a1a; --c-surface: #2a2a2a; --c-text: #eee; }  /* Dark（默认）*/
.theme-light     { --c-bg: #fff; --c-surface: #f5f5f5; --c-text: #222; }
.theme-transparent { --c-bg: transparent; --c-surface: rgba(255,255,255,0.05); }
.theme-print     { --c-bg: #fff; --c-surface: #fff; --c-text: #000; }  /* 打印友好 */
```

切换按钮：`<button onclick="cycleTheme()">🎨 Theme</button>`

---

## Step 1 — FVM 基础求解器（5th SolverExport）

**理论：** `docs/theory.md` §4

**设计决策（D-1）：**
- 实现为独立文件 `cmd/main/solver_fvm.mbt`
- 作为第 5 个 SolverExport 类型（`--solver fvm`）
- 不修改现有 4 个求解器的任何代码
- 使用 SIMPLE 算法在有限体积法框架下实现

**实现位置：** `cmd/main/solver_fvm.mbt`（新建）

### 全局状态（新增）

```moonbit
let g_u_fvm : Array[Array[Double]] = ...
let g_v_fvm : Array[Array[Double]] = ...
let g_p_fvm : Array[Array[Double]] = ...
```

### 新增 WASM 导出（计划）

`init_fvm`, `run_fvm_n_iter`, `get_fvm_step_count`,
`get_u_fvm_at`, `get_v_fvm_at`, `get_p_fvm_at`, `get_fvm_divergence_norm`

### 测试（T118–T122）

| 编号 | 测试名 | 验证内容 |
|---|---|---|
| T118 | `fvm_state_independent` | FVM 不影响其他 4 个求解器 |
| T119 | `fvm_mass_conservation` | 100 步后散度范数 < 0.01 |
| T120 | `fvm_bc_correct` | 顶盖 u=1，壁面 u=v=0 |
| T121 | `fvm_vs_chorin_pressure_sign` | 压力分布符号与 Chorin 一致 |
| T122 | `fvm_ghia_l2_error` | Re=100, FVM Ghia L2 误差 < 0.1 |

**提交：** `feat(solver): add FVM-SIMPLE solver as 5th SolverExport (Step 1)`

---

## Step K — 运行时 nx/ny（动态网格尺寸）

**理论：** 动态内存分配，`Array::make(n, ...)` 运行时尺寸

**难度：极高**（全部 2D 数组须动态分配，所有 `nx`/`ny` 编译期常量须改为函数参数传递）

**设计决策（D-K）：**
- 将 `nx`, `ny` 从模块顶层常量改为全局状态数组：
  ```moonbit
  let g_nx : Array[Int] = [41]
  let g_ny : Array[Int] = [41]
  ```
- 所有求解器函数接受 `nr: Int, nc: Int` 参数（已在部分函数中实现）
- 全局状态数组（`g_u`, `g_v` 等）改为二级包装：`Array[Array[Array[Double]]]`（外层 size-1 数组包含实际 nr×nc 数组）
- 或改用懒初始化模式（首次调用 `init_*` 时分配）

**CLI 集成：**
```
--grid 81x81  → 运行时设置 g_nx[0]=81, g_ny[0]=81
```

**测试（T123–T126）：**
- T123: `grid_81x81_valid` — 81×81 初始化后 g_nx[0]=81
- T124: `grid_switch_independent` — 41→81 切换后两套结果独立
- T125: `grid_41x41_matches_baseline` — 41×41 数值与原 Chorin 结果一致
- T126: `grid_161x161_stable` — 161×161 运行 100 步无 NaN

**提交：** `feat(core): add runtime nx/ny dynamic grid size (Step K)`

---

## Step 6 — 被动标量输运

**理论：** `docs/theory.md` §6

**实现位置：** `cmd/main/solver_scalar.mbt`（新建）

### 编译期开关

```moonbit
let run_scalar : Array[Bool] = [false]  // CLI --scalar 激活
```

### 全局状态（新增）

```moonbit
let g_phi : Array[Array[Double]] = ...   // 标量场
```

### CLI 集成

```
--scalar --alpha 0.01    → 激活标量输运，扩散系数 α=0.01
--scalar-bc left=1,right=0,top=0,bottom=0
```

### 测试（T127–T130）

| 编号 | 测试名 | 验证内容 |
|---|---|---|
| T127 | `scalar_state_independent` | 标量不影响速度场 |
| T128 | `scalar_mass_conservation_pure_diffusion` | 纯扩散（u=v=0）质量守恒 |
| T129 | `scalar_pe_limit_convection` | Pe→∞：标量沿流线守恒 |
| T130 | `scalar_steady_state_harmonic` | Pe→0 稳态 = 调和函数（∇²φ≈0）|

**提交：** `feat(physics): add passive scalar transport (Step 6)`

---

## Step 7 — 能量方程 + Boussinesq 近似

**理论：** `docs/theory.md` §7

**实现位置：** `cmd/main/solver_energy.mbt`（新建）

### 编译期开关

```moonbit
let run_energy : Array[Bool] = [false]
```

### 全局状态（新增）

```moonbit
let g_theta : Array[Array[Double]] = ...   // 无量纲温度场
```

### CLI 集成

```
--energy --ra 1e5 --pr 0.71
--energy-bc left=1,right=0,top=zero-grad,bottom=zero-grad
```

### 测试（T131–T135）

| 编号 | 测试名 | 验证内容 |
|---|---|---|
| T131 | `energy_state_independent` | 温度场不影响速度场（当 Ra=0）|
| T132 | `energy_bc_hot_cold` | 左壁 θ=1，右壁 θ=0 正确施加 |
| T133 | `boussinesq_buoyancy_direction` | Ra>0 时 v 分量浮力源项 > 0（热壁附近）|
| T134 | `nusselt_number_ra1e3` | Ra=1000, Pr=0.71 → Nu ∈ [1.0, 1.3] |
| T135 | `de_vahl_davis_ra1e4` | Ra=1e4, 41×41 → |Nu-2.243|/2.243 < 10% |

**提交：** `feat(physics): add energy equation with Boussinesq buoyancy (Step 7)`

---

## Step 8 — BC 库扩展

**理论：** `docs/theory.md` §8

**实现位置：** `cmd/main/bc_library.mbt`（新建）

### 新增 BC 类型

```moonbit
enum BCType {
  Dirichlet(Double)    // 给定值
  Neumann              // 零梯度
  Symmetric            // 法向分量为零
  Periodic             // 周期性
}

struct BoundaryConditions {
  left: BCType; right: BCType; top: BCType; bottom: BCType
}

fn apply_bc(phi: Array[Array[Double]], bc: BoundaryConditions, nr: Int, nc: Int) -> Unit
```

### 测试（T136–T139）

| 编号 | 测试名 | 验证内容 |
|---|---|---|
| T136 | `bc_neumann_zero_gradient` | Neumann: `phi[i][nc-1] == phi[i][nc-2]` |
| T137 | `bc_symmetric_tangential` | 对称: 法向为零，切向梯度为零 |
| T138 | `bc_periodic_wraparound` | 周期: `phi[i][0] == phi[i][nc-1]` |
| T139 | `bc_dirichlet_exact` | Dirichlet: 边界值精确等于设定值 |

**提交：** `feat(bc): add Neumann/symmetric/periodic BC library (Step 8)`

---

## Step 5 — 3D 顶盖驱动流

**理论：** `docs/theory.md` §5

**设计决策（D-5）：**
- 独立包 `cmd/main3d/`（零回归风险）
- 独立 `moon.pkg.json`、`build_wasm3d.sh`
- 网格：`nz = 41`（可运行时配置，依赖 Step K）

**依赖：** Step K（运行时网格尺寸）

**参考：** Albensoeder & Kuhlmann (2002) Re=100 基准

**测试（T140–T144）：**
- T140: `3d_state_init` — 初始化后 nz=41
- T141: `3d_bc_6faces` — 6 个面 BC 正确施加
- T142: `3d_divergence_initial` — 初始散度 = 0
- T143: `3d_vortex_center` — 500 步后 u/v/w 涡心符号正确
- T144: `3d_velocity_profile` — Re=100, 中面速度与 Albensoeder 定性一致

**提交：** `feat(3d): add 3D lid-driven cavity solver in cmd/main3d/ (Step 5)`

---

## Step 9 — 拉伸网格

**理论：** `docs/theory.md` §9

**实现位置：** `cmd/main/grid_stretched.mbt`（新建）

```moonbit
fn stretched_grid_tanh(n: Int, L: Double, beta_s: Double) -> Array[Double]
// → 返回 n 个节点坐标

fn local_dx(x_nodes: Array[Double]) -> Array[Double]
// → 返回 n-1 个局部间距

fn laplacian_nonuniform(phi, x_nodes, y_nodes, i, j) -> Double
// → 非均匀网格拉普拉斯算子（二阶精度）
```

**测试（T145–T147）：**
- T145: `stretched_grid_endpoints` — x[0]=0, x[n-1]=L 精确
- T146: `stretched_grid_monotone` — x[j+1] > x[j]（单调）
- T147: `laplacian_nonuniform_quadratic_exact` — $\phi=x^2$ 时拉普拉斯精确等于 2

**提交：** `feat(grid): add tanh stretched grid with non-uniform Laplacian (Step 9)`

---

## Step 10 — 曲线网格基础

**理论：** `docs/theory.md` §10

**依赖：** Step 9（拉伸网格）+ Step 1（FVM）

**实现位置：** `cmd/main/grid_curvilinear.mbt`（新建）

```moonbit
struct JacobianCell {
  x_xi: Double; x_eta: Double
  y_xi: Double; y_eta: Double
  det_J: Double  // |J|
  g11: Double; g12: Double; g22: Double  // 度量张量
}

fn compute_jacobian(x_phys, y_phys, xi, eta) -> JacobianCell
fn curvilinear_laplacian(phi, jacobians, i, j) -> Double
fn curvilinear_convection(phi, u_contra, v_contra, jacobians, i, j) -> Double
```

**测试（T148–T151）：**
- T148: `jacobian_uniform_is_identity` — 均匀网格 → J=I, det=1
- T149: `curvilinear_laplacian_harmonic` — 调和函数 → 拉普拉斯 ≈ 0
- T150: `curvilinear_vs_uniform_square` — 均匀曲线网格 ≈ 标准 FDM
- T151: `curvilinear_conservation` — 散度定理（面积分 = 体积分）

**提交：** `feat(grid): add curvilinear grid coordinate transformation (Step 10)`

---

## 回归保护协议

每步实施后执行：
```bash
moon test --target wasm               # 必须：原 108 个测试全部通过
bash run_local.sh --solver chorin --format csv > /tmp/chorin_new.csv
diff docs/chorin_baseline.csv /tmp/chorin_new.csv    # 数值零差异
```

新增测试总计划：
- Step I: T109–T113（5 个）
- Step J: T114–T117（4 个）
- Step 1: T118–T122（5 个）
- Step K: T123–T126（4 个）
- Step 6: T127–T130（4 个）
- Step 7: T131–T135（5 个）
- Step 8: T136–T139（4 个）
- Step 5: T140–T144（5 个）
- Step 9: T145–T147（3 个）
- Step 10: T148–T151（4 个）
- **计划总计：151 个测试**

---

## 关键设计决策汇总

| 编号 | 决策 | 选择 | 理由 |
|---|---|---|---|
| D-1 | FVM 实现方式 | 独立文件 `solver_fvm.mbt` | 零回归风险 |
| D-5 | 3D 策略 | 独立包 `cmd/main3d/` | 零回归风险 |
| D-J1 | RK3 集成方式 | 抽取 RHS 函数 + 新全局状态 | 不修改现有 Chorin 路径 |
| D-K | 动态网格策略 | 全局 `Array[Int]` + 懒初始化 | 最小化改动量 |
| D-6 | 多物理场开关 | 编译期 `Array[Bool]` + CLI 激活 | 零性能损失（未激活时）|

---

*作者：Fenfen Yu（余芬芬），AI 协作：Claude Sonnet 4.6*
*日期：2026-04-18*
