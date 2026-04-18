# FlowLabLite — I/O 模块测试报告

**日期：** 2026-04-18  
**作者：** Fenfen Yu (余芬芬)  
**测试套件：** `moon test --target wasm` — **91 tests, passed: 91, failed: 0**

---

## 1. 范围

本报告覆盖 I/O 模块新增测试 T84–T91，及对现有 83 个测试（T1–T83）的回归验证。

---

## 2. 新增测试（T84–T91）

所有测试位于 `cmd/main/main_ext_wbtest.mbt`，分类遵循项目既有约定：

| 测试编号 | 测试名称 | 分类 | 被测函数 | 验证内容 | 结果 |
|---|---|---|---|---|---|
| T84 | `csv_header_present` | Unit / EP | `format_csv` | RFC 4180 标准列头行存在 | **PASS** |
| T85 | `csv_row_count` | Unit / BVA / OBT | `format_csv` | nc=3 时输出恰好 9 行数据 | **PASS** |
| T86 | `vtk_dimensions_header` | Unit / EP | `format_vtk` | DATASET STRUCTURED_GRID + DIMENSIONS 标头 | **PASS** |
| T87 | `vtk_points_count` | Unit / OBT | `format_vtk` | POINTS 计数 = nc² | **PASS** |
| T88 | `tecplot_zone_count` | Unit / EP | `format_tecplot` | k 个求解器 → k 个 ZONE 块 | **PASS** |
| T89 | `netcdf_stub_returns_err` | Unit / PP (负路径) | `format_netcdf` | Stub 必须返回 Err，不崩溃 | **PASS** |
| T90 | `vtk_point_data_header` | Unit / OBT | `format_vtk` | POINT_DATA + SCALARS + VECTORS 标头 | **PASS** |
| T91 | `tecplot_solution_time` | Unit / OBT | `format_tecplot` | SOLUTIONTIME = steps × dt | **PASS** |

### 测试辅助函数

`make_test_export(name, nc, fill_u, fill_v, fill_p)` — 创建常值填充的小网格 SolverExport，与全局求解器状态完全独立，测试运行速度快（<1 ms）。

---

## 3. 回归验证（T1–T83）

原有 83 个测试全部通过，无回归。

| 测试分组 | 范围 | 数量 | 状态 |
|---|---|---|---|
| Chorin 求解器 | T1–T16 | 16 | 全部 PASS |
| SIMPLE 求解器 | T17–T26 | 10 | 全部 PASS |
| Chorin-PCG 求解器 | T27–T34 | 8 | 全部 PASS |
| MAC 交错网格 | T35–T42 | 8 | 全部 PASS |
| DILU 预条件器 | T43–T48 | 6 | 全部 PASS |
| DIC 预条件器 | T49–T54 | 6 | 全部 PASS |
| GAMG 预条件器 | T55–T58 | 4 | 全部 PASS |
| 数组工具扩展 | T59–T60 | 2 | 全部 PASS |
| Laplacian 算子 | T61–T63 | 3 | 全部 PASS |
| 边界条件 | T64–T65 | 2 | 全部 PASS |
| RHS 源项 | T66–T67 | 2 | 全部 PASS |
| GAMG 子组件 | T68–T70 | 3 | 全部 PASS |
| PCG 残差集成 | T71–T74 | 4 | 全部 PASS |
| 跨预条件器一致性 | T75 | 1 | PASS |
| 系统物理测试 | T76–T79 | 4 | 全部 PASS |
| 回归测试 | T80–T83 | 4 | 全部 PASS |

---

## 4. I/O 功能集成验证（手工）

| 场景 | 命令 | 验证内容 | 结果 |
|---|---|---|---|
| --list-formats | `moon run cmd/main --target wasm -- --list-formats` | 5 种格式列表，NetCDF ✗ 正确标注 | **PASS** |
| VTK 输出 | `bash run_local.sh --format vtk --solver chorin` | 生成 6739 行 `.vtk`，DIMENSIONS 41 41 1 | **PASS** |
| Tecplot 输出 | `bash run_local.sh --format tecplot` | 4 个 ZONE（Chorin/SIMPLE/PCG/MAC），MAC I=40 J=40 | **PASS** |
| CSV 输出 | `bash run_local.sh --format csv` | 4 × 1681 + 2 头行 = 6726 行，RFC 4180 | **PASS** |
| NetCDF Stub | `moon run cmd/main --target wasm -- --format netcdf` | 打印 [io] 错误提示，退出码 0 | **PASS** |
| 无参数默认 JSON | `moon run cmd/main --target wasm` | 与修改前行为完全一致 | **PASS** |
| 数值完整性 | run_local.sh JSON 对比 | center_u / max_v / div_norm 未变化 | **PASS** |

---

## 5. 已知设计决策

### MoonBit Double 渲染
MoonBit 对 `0.0` 渲染为 `"0"`（整数样式），对 `0.05` 渲染为 `"0.05"`。VTK/Tecplot/CSV 解析器均能接受无小数点的整数形式，不影响 ParaView 读取。

### NetCDF-4 降级机制
`io_netcdf_stub.mbt` 永远编译通过（无外部依赖）。`netcdf_available = false` 常量触发 `--list-formats` 中的 ✗ 标注。升级路径见 `io.md §4.2`。

### 安静模式（quiet mode）
`g_quiet_mode[0] = true` 在 `--format` 模式下抑制三处进度 println：
- `cavity_flow_array`: "Time step N/500 completed"
- `cavity_flow_pcg`: "PCG step N/...: pressure iters=..."
- `cavity_flow_mac`: "MAC step N/...: pressure iters=..."

确保 stdout 输出为纯格式数据，可直接 `> file.vtk` 重定向。

---

## 6. 文件变更清单

| 文件 | 类型 | 说明 |
|---|---|---|
| `cmd/main/io_formats.mbt` | 新建 | SolverExport 结构体 + CSV/VTK/Tecplot 生成函数 |
| `cmd/main/io_netcdf_stub.mbt` | 新建 | NetCDF-4 stub（永远编译，返回 Err） |
| `cmd/main/io_cli.mbt` | 新建 | CLI 参数解析 + 格式分发 |
| `cmd/main/main.mbt` | 修改 | g_quiet_mode 全局；3 处进度保护；fn main 前置格式分发 |
| `cmd/main/moon.pkg.json` | 修改 | 添加 `moonbitlang/core/env` 导入 |
| `cmd/main/main_ext_wbtest.mbt` | 修改 | 新增 T84–T91（8 个测试） |
| `cmd/main/main.html` | 修改 | 工具栏添加 ⬇ VTK / ⬇ Tecplot 按钮；downloadVTK/downloadTecplot JS 方法 |
| `run_local.sh` | 修改 | 支持 --format / --solver 透传；非 JSON 格式直接保存文件 |
| `io.md` | 新建 | I/O 模块规划文档 |
| `docs/test_report_io_20260418.md` | 新建 | 本测试报告 |

---

**总计：91 tests, passed: 91, failed: 0**
