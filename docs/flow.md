# FlowLabLite — Execution Flow Diagrams

## 1. `main()` Function

```mermaid
flowchart TD
    A([main 入口]) --> B[记录程序启动时间\nnow_seconds]
    B --> C[打印求解器信息\n网格 41×41 / 500步 / Re=20]
    C --> D[create_zeros_2d\n初始化 u, v, p 数组]
    D --> E[generate_mesh_grid\n生成坐标网格 x, y]
    E --> F[记录仿真开始时间\nsim_start]
    F --> G[cavity_flow_array\nnt_test=100 步]
    G --> H[记录仿真结束时间\nsim_end]
    H --> I{have_timer?}
    I -- 是 --> J[打印仿真耗时]
    I -- 否 --> K[打印 timing unavailable]
    J --> L[打印中心点结果\nu / v / p / 速度幅值]
    K --> L
    L --> M[打印顶盖速度剖面\nj=0,5,10,...,40]
    M --> N[统计全场最大值\nmax_u / max_v / max_p / min_p]
    N --> O[打印采样点速度幅值\n4个样本点]
    O --> P[write_velocity_to_file\nvelocity_field.txt]
    P --> Q[output_json\n生成 JSON 供 local_viewer.html]
    Q --> R([结束])
```

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
