# FlowLabLite — Claude Code 项目简报

## 项目概述

2D 顶盖驱动方腔流（lid-driven cavity flow）CFD 求解器，MoonBit 实现，编译为 WebAssembly (wasm-gc)，通过浏览器可视化结果。无运行时依赖。

- 作者：Fenfen Yu（余芬芬）
- 合作单位：北京航空航天大学、鄂州海慕科技有限公司
- 许可证：Apache-2.0

---

## 两个独立求解器

| 求解器 | 算法 | 全局状态 | 用途 |
|---|---|---|---|
| **Chorin** | 投影法，显式时间推进 | `g_u`, `g_v`, `g_p` | 瞬态模拟，500 时间步 |
| **SIMPLE** | 压力修正迭代，欠松弛 | `g_u_s`, `g_v_s`, `g_p_s` | 稳态迭代求解 |

两套状态数组完全独立，运行 SIMPLE 不会影响 Chorin 状态，反之亦然。

---

## 关键参数

```
网格：nx = ny = 41（41×41 节点）
时间步：dt = 0.001，总步数 nt = 500
压力 Poisson 内迭代：nit = 50（Gauss-Seidel）
流体：rho = 1.0，nu = 0.1，Re = 20
域：2×2，dx = dy = 0.05

SIMPLE 欠松弛：alpha_p = 0.3（压力），alpha_u = 0.7（速度）
```

---

## 重要实现约定

### 边界条件顺序（必须遵守）
侧壁 BC **先设**，顶盖 BC **最后设**，否则角点 `u[ny-1][0]` 和 `u[ny-1][nx-1]` 会被侧壁循环覆盖为 0：

```moonbit
// 正确顺序
for i = 0; i < ny; i = i + 1 { u[i][0] = 0.0; u[i][nx-1] = 0.0 }  // 侧壁
for j = 0; j < nx; j = j + 1 { u[ny-1][j] = 1.0 }                  // 顶盖（最后）
```

### MoonBit 语法限制
不支持行首 `+` 续行，多行表达式必须拆成中间变量：

```moonbit
// 错误（编译报错）
let visc = nu * dt / dx2 * (...)
          + nu * dt / dy2 * (...)   // ← 不合法

// 正确
let visc_x = nu * dt / dx2 * (...)
let visc_y = nu * dt / dy2 * (...)
let visc   = visc_x + visc_y
```

---

## 构建与测试命令

```bash
# 测试（必须加 --target wasm，有 26 个测试）
moon test --target wasm

# 构建 WASM（必须用脚本，moon build 单独不导出函数）
bash build_wasm.sh release   # 生产版本
bash build_wasm.sh           # 调试版本

# 本地运行求解器 + 生成 results.json
bash run_local.sh

# 启动本地 HTTP 服务（浏览器查看）
python -m http.server 8080
# 然后打开 http://localhost:8080/cmd/main/main.html
```

---

## 已知问题

- **moon 0.1.20260309 不传 `-exported_functions` 给链接器**：`moon build` 生成的 WASM 只有 `_start`，必须用 `build_wasm.sh` 重新链接，才能导出全部 API 函数。

- **浏览器要求 Chrome 115+**：WASM 使用 wasm-gc stringref（type code `0x77`），旧版浏览器不支持。Node.js 也不能直接加载，仅用于验证导出段。

---

## WASM 导出函数（共 38 个）

**Chorin 求解器（27 个）：**
`init_simulation`, `run_all_steps`, `run_n_steps`,
`get_nx`, `get_ny`, `get_nt`, `get_nit`, `get_re`, `get_dx`, `get_dy`, `get_dt`, `get_rho`, `get_nu`,
`get_step_count`, `get_u_at`, `get_v_at`, `get_p_at`,
`get_velocity_magnitude_at`, `get_max_velocity_magnitude`,
`get_u_center`, `get_v_center`, `get_p_center`,
`get_divergence_norm`, `get_max_u`, `get_max_v`, `get_max_p`, `get_min_p`

**SIMPLE 求解器（11 个）：**
`init_simple`, `run_simple_n_iter`,
`get_simple_step_count`, `get_simple_residual`,
`get_u_simple_at`, `get_v_simple_at`, `get_p_simple_at`,
`get_max_u_simple`, `get_simple_divergence_norm`

---

## 文件结构要点

```
cmd/main/main.mbt         # 全部求解器代码（Chorin + SIMPLE）
cmd/main/main_wbtest.mbt  # 26 个白盒测试
cmd/main/moon.pkg.json    # WASM 导出列表（两套求解器函数都在这里）
cmd/main/main.html        # 浏览器可视化页面
cmd/main/local_viewer.html# 本地 JSON 结果查看器
build_wasm.sh             # WASM 构建脚本（含导出函数列表）
run_local.sh              # 本地运行 + 提取 JSON
```

---

## 测试状态

26 个测试，全部通过（`moon test --target wasm`）：
- 测试 1–16：Chorin 求解器（数组工具、边界条件、物理验证）
- 测试 17–26：SIMPLE 求解器（状态重置、边界条件、收敛残差、涡旋形成、独立性）
