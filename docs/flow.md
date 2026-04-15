# FlowLabLite — Execution Flow Diagrams

## 1. `main()` Function

> **更新说明（2026-04-15）：** main() 已重构，使用 `@bench` 实现真实计时，  
> 同时运行 Chorin 和 SIMPLE 两个求解器，输出包含双结果的统一 JSON。  
> 已移除：`now_seconds()`（存根）、`generate_mesh_grid()`（仅本地使用）、  
> `write_velocity_to_file()`（无效占位实现）。

```mermaid
flowchart TD
    A([main 入口]) --> B[@bench.monotonic_clock_start\n记录程序启动时间戳]
    B --> C[打印求解器配置\n网格 41×41 / Re=20\nChorin local_nt 步 / SIMPLE local_simple_n 次]

    C --> D{run_chorin = true?}
    D -- 是 --> E["timed('Chorin init')\n→ init_simulation()"]
    D -- 否 --> G
    E --> F["timed('Chorin N steps')\n→ run_n_steps(local_nt)"]
    F --> G{run_simple = true?}

    G -- 是 --> H["timed('SIMPLE init')\n→ init_simple()"]
    G -- 否 --> J
    H --> I["timed('SIMPLE N iters')\n→ run_simple_n_iter(local_simple_n)"]
    I --> J[@bench.monotonic_clock_end\n计算 total_ms]

    J --> K{timing_enabled?}
    K -- 是 --> L[打印 Timing Summary\n各阶段 ms + 总计]
    K -- 否 --> M
    L --> M[output_json\ng_u/v/p + g_u_s/v_s/p_s\n+ timing → JSON_DATA 标记块]
    M --> N([结束])
```

### 本地运行配置常量（`main.mbt` 顶部）

| 常量 | 类型 | 默认值 | 说明 |
|---|---|---|---|
| `timing_enabled` | Bool | true | 计时总开关 |
| `time_chorin_phase` | Bool | true | Chorin 分项计时 |
| `time_simple_phase` | Bool | true | SIMPLE 分项计时 |
| `run_chorin` | Bool | true | 是否运行 Chorin |
| `local_nt` | Int | nt=500 | Chorin 步数 |
| `run_simple` | Bool | true | 是否运行 SIMPLE |
| `local_simple_n` | Int | 100 | SIMPLE 迭代次数 |

---

## 2. `cavity_flow_array()` — Chorin 投影法

```mermaid
flowchart TD
    A([cavity_flow_array\nnt, u, v, dt, dx, dy, p, rho, nu]) --> B[分配工作数组\nun, vn, b]
    B --> C{n = 0; n < nt?}
    C -- 否，循环结束 --> Z([返回 u, v, p])
    C -- 是 --> D[copy_2d_array\nu → un,  v → vn\n保存旧时间层]

    D --> E[build_up_b_array\n计算压力 Poisson 右端项 b\nb = ρ·1/dt·∇·u − 对流项²]

    E --> F[pressure_poisson_array\nGauss-Seidel 迭代 nit=50 次\n求解 ∇²p = b]

    F --> G[更新 u 动量方程\n内部节点 i=1..ny-2, j=1..nx-2\nu = un − 对流 − ∂p/∂x + 粘性]

    G --> H[更新 v 动量方程\n内部节点 i=1..ny-2, j=1..nx-2\nv = vn − 对流 − ∂p/∂y + 粘性]

    H --> I[施加速度边界条件\n底壁/左壁/右壁: u=v=0\n顶盖: u=1, v=0]

    I --> J{n+1 mod 10 = 0?}
    J -- 是 --> K[打印进度\nTime step n+1/nt]
    J -- 否 --> L[n = n+1]
    K --> L
    L --> C
```

> **压力 Poisson 右端项 b（`build_up_b_array`）**
>
> $$b_{i,j} = \rho \left[\frac{1}{\Delta t}\left(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\right) - \left(\frac{\partial u}{\partial x}\right)^2 - 2\frac{\partial u}{\partial y}\frac{\partial v}{\partial x} - \left(\frac{\partial v}{\partial y}\right)^2\right]$$

---

## 3. SIMPLE 算法 — `simple_one_iter()` / `run_simple_n_iter()`

```mermaid
flowchart TD
    START([run_simple_n_iter\nn 次迭代]) --> INIT[分配工作数组\nu_star, v_star, b, p_corr]

    INIT --> LOOP{k = 0; k < n?}
    LOOP -- 否 --> POST[缓存收敛残差\ng_simple_residual =\nget_simple_divergence_norm]
    POST --> END([返回])

    LOOP -- 是 --> A

    subgraph ONE ["simple_one_iter — 单次 SIMPLE 扫描"]
        A[/"(a) 动量预测步\n内部节点 i=1..ny-2, j=1..nx-2"/]
        A --> A1["u_new = u − 对流_u − ∂p*/∂x + 粘性_u\nu* = α_u · u_new + (1−α_u) · u\n（α_u = 0.7 速度欠松弛）"]
        A1 --> A2["v_new = v − 对流_v − ∂p*/∂y + 粘性_v\nv* = α_u · v_new + (1−α_u) · v"]
        A2 --> B["施加 u*, v* 边界条件\n侧壁先设 → 顶盖最后设\nu*[top,:] = 1,  其余 = 0"]
        B --> C[/"(b) 构造压力修正源项\nb = build_up_b_array(u*, v*)"/]
        C --> D[/"(c) 求解压力修正 Poisson\np' = pressure_poisson_array(b)\nGauss-Seidel 50 次"/]
        D --> E[/"(d) 场量更新（欠松弛）\n内部节点"/]
        E --> E1["p  ← p* + α_p · p'\n（α_p = 0.3 压力欠松弛）"]
        E1 --> E2["u  ← u* − dt/(2ρ·dx) · (p'[j+1]−p'[j-1])\nv  ← v* − dt/(2ρ·dy) · (p'[i+1]−p'[i-1])"]
        E2 --> F["施加 u, v, p 边界条件\n速度：侧壁→顶盖\n压力：Neumann 底/侧壁，p=0 顶部"]
    end

    F --> G[g_simple_iters += 1]
    G --> LOOP
```

### 欠松弛因子

| 参数 | 值 | 作用 |
|---|---|---|
| `simple_alpha_u` | 0.7 | 速度欠松弛，抑制动量方程振荡 |
| `simple_alpha_p` | 0.3 | 压力欠松弛，稳定压力修正收敛 |

### 收敛判据

$$\text{residual} = \frac{1}{(N_x-2)(N_y-2)} \sum_{i,j} \left|\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y}\right|$$

每批 `run_simple_n_iter(n)` 结束后自动计算并缓存到 `g_simple_residual`，
可通过 `get_simple_residual()` 读取。
