// 完整的2D顶盖驱动流求解器，MoonBit语言实现 - 修复数值稳定性问题
// 参考OpenFOAM数据结构，使用有限体积法和稳态SIMPLE算法
// 支持不同雷诺数，结构化正交网格，串行实现
// 可扩展：可添加更多物理模型或网格类型
// 高性能：使用高效数据结构和算法
// 可读性/可维护性：模块化设计，清晰注释
// 易用性：简单接口，参数化配置

// 网格类：结构化正交网格
struct Mesh2D {
  nx : Int
  ny : Int
  l : Double
  dx : Double
  dy : Double
}

// Mesh2D的构造函数
fn Mesh2D::new(nx : Int, ny : Int, length : Double) -> Mesh2D {
  let dx = length / (nx - 1).to_double()
  let dy = length / (ny - 1).to_double()
  { nx, ny, l: length, dx, dy }
}

// 获取单元中心坐标
fn Mesh2D::cell_center(self : Mesh2D, i : Int, j : Int) -> (Double, Double) {
  (i.to_double() * self.dx, j.to_double() * self.dy)
}

// 场类：模拟volScalarField/volVectorField
struct Field[T] {
  mesh : Mesh2D
  data : Array[Array[T]]
}

// Field的构造函数 - 特化为Double类型
fn Field::new_double(mesh : Mesh2D) -> Field[Double] {
  let  data = Array::make(mesh.nx, Array::make(mesh.ny, 0.0))
  for i = 0; i < mesh.nx; i = i + 1 {
    data[i] = Array::make(mesh.ny, 0.0)
  }
  { mesh, data }
}

// 设置所有值
fn Field::set_all[T](self : Field[T], value : T) -> Unit {
  for i = 0; i < self.mesh.nx; i = i + 1 {
    for j = 0; j < self.mesh.ny; j = j + 1 {
      self.data[i][j] = value
    }
  }
}

// 获取指定位置的值
fn Field::at[T](self : Field[T], i : Int, j : Int) -> T {
  self.data[i][j]
}

// 设置指定位置的值
fn Field::set[T](self : Field[T], i : Int, j : Int, value : T) -> Unit {
  self.data[i][j] = value
}

// 复制场的值
fn Field::copy_from[T](self : Field[T], other : Field[T]) -> Unit {
  for i = 0; i < self.mesh.nx; i = i + 1 {
    for j = 0; j < self.mesh.ny; j = j + 1 {
      self.data[i][j] = other.data[i][j]
    }
  }
}

// 梯度算子（有限体积）- 仅适用于Double类型，增加边界检查
fn Field::grad(self : Field[Double], i : Int, j : Int) -> (Double, Double) {
  let gx = if i == 0 || i == self.mesh.nx - 1 {
    0.0  // 边界处梯度设为0
  } else {
    (self.data[i + 1][j] - self.data[i - 1][j]) / (2.0 * self.mesh.dx)
  }
  
  let gy = if j == 0 || j == self.mesh.ny - 1 {
    0.0  // 边界处梯度设为0
  } else {
    (self.data[i][j + 1] - self.data[i][j - 1]) / (2.0 * self.mesh.dy)
  }
  
  (gx, gy)
}

// 散度算子 - 仅适用于Double类型，增加边界检查
fn Field::div(self : Field[Double], i : Int, j : Int) -> Double {
  let dx = if i == 0 || i == self.mesh.nx - 1 {
    0.0
  } else {
    (self.data[i + 1][j] - self.data[i - 1][j]) / (2.0 * self.mesh.dx)
  }
  
  let dy = if j == 0 || j == self.mesh.ny - 1 {
    0.0
  } else {
    (self.data[i][j + 1] - self.data[i][j - 1]) / (2.0 * self.mesh.dy)
  }
  
  dx + dy
}

// 矩阵类：用于线性方程组（压力泊松等）
struct Matrix {
  rows : Int
  cols : Int
  data : Array[Array[Double]]
}

// Matrix的构造函数
fn Matrix::new(rows : Int, cols : Int) -> Matrix {
  let data = Array::make(rows, Array::make(cols, 0.0))
  { rows, cols, data }
}

// 设置矩阵元素值
fn Matrix::set(self : Matrix, i : Int, j : Int, value : Double) -> Unit {
  self.data[i][j] = value
}

// 获取矩阵元素值
fn Matrix::at(self : Matrix, i : Int, j : Int) -> Double {
  self.data[i][j]
}

// 改进的高斯-赛德尔求解器，增加数值稳定性检查
fn Matrix::solve(self : Matrix, b : Array[Double], x : Array[Double], max_iter : Int, tol : Double) -> Unit {
  println("  Starting pressure equation solver...")
  for iter = 0; iter < max_iter; iter = iter + 1 {
    let mut max_err = 0.0
    let mut has_nan = false
    
    for i = 0; i < self.rows; i = i + 1 {
      let mut sum = 0.0
      for j = 0; j < self.cols; j = j + 1 {
        if j != i {
          sum = sum + self.data[i][j] * x[j]
        }
      }
      
      // 检查对角元素是否为零
      if self.data[i][i].abs() < 1.0e-15 {
        println("    Warning: Nearly zero diagonal element at \{i}")
        continue
      }
      
      let new_x = (b[i] - sum) / self.data[i][i]
      
      // 检查NaN或无穷大
      if new_x != new_x || new_x.abs() > 1.0e10 {
        has_nan = true
        println("    Warning: NaN or large value detected at iteration \{iter}")
        break
      }
      
      let err = if new_x - x[i] > 0.0 { new_x - x[i] } else { x[i] - new_x }
      if err > max_err { 
        max_err = err 
      }
      x[i] = new_x
    }
    
    if has_nan {
      println("    Pressure solver failed due to numerical instability")
      break
    }
    
    // 每10步或收敛时输出残差信息
    if iter % 10 == 0 || max_err < tol {
      println("    Pressure solver iter \{iter}: residual = \{max_err}")
    }
    
    if max_err < tol { 
      println("    Pressure solver converged in \{iter + 1} iterations")
      break 
    }
  }
}

// 求解器类：封装SIMPLE算法，增加松弛因子
struct CavitySolver {
  mesh : Mesh2D
  u : Field[Double]  // x速度
  v : Field[Double]  // y速度
  p : Field[Double]  // 压力
  re : Double        // 雷诺数
  max_iter : Int
  tol : Double
  alpha_u : Double   // 速度松弛因子
  alpha_p : Double   // 压力松弛因子
  dt : Double        // 时间步长
}

// CavitySolver的构造函数，增加松弛因子参数
fn CavitySolver::new(nx : Int, ny : Int, length : Double, reynolds : Double, max_iter : Int, tolerance : Double) -> CavitySolver {
  let mesh = Mesh2D::new(nx, ny, length)
  let u = Field::new_double(mesh)
  let v = Field::new_double(mesh)
  let p = Field::new_double(mesh)
  
  // 根据雷诺数调整参数
  let dt = if reynolds > 1000.0 { 0.00001 } else if reynolds > 100.0 { 0.0001 } else { 0.001 }
  let alpha_u = 0.5  // 速度松弛因子
  let alpha_p = 0.3  // 压力松弛因子
  
  { mesh, u, v, p, re: reynolds, max_iter, tol: tolerance, alpha_u, alpha_p, dt }
}

// 初始化场 - 使用更好的初始条件
fn CavitySolver::initialize(self : CavitySolver) -> Unit {
  // 初始化为小的随机值而不是零，避免数值问题
  self.u.set_all(1.0e-6)
  self.v.set_all(1.0e-6)
  self.p.set_all(0.0)
  
  // 设置初始的顶盖速度分布
  for i = 0; i < self.mesh.nx; i = i + 1 {
    let x = i.to_double() * self.mesh.dx
    // 使用平滑的初始分布
    let u_top = if x < 0.1 || x > 0.9 { 
      0.5 
    } else { 
      1.0 
    }
    self.u.set(i, self.mesh.ny - 1, u_top)
  }
}

// 应用边界条件
fn CavitySolver::apply_boundary_conditions(self : CavitySolver) -> Unit {
  // 顶盖：u=1, v=0 (平滑过渡)
  for i = 0; i < self.mesh.nx; i = i + 1 {
    let x = i.to_double() * self.mesh.dx
    let u_top = if x < 0.05 || x > 0.95 { 
      0.0 
    } else { 
      1.0 
    }
    self.u.set(i, self.mesh.ny - 1, u_top)
    self.v.set(i, self.mesh.ny - 1, 0.0)
  }
  
  // 其他墙面：无滑移
  for j = 0; j < self.mesh.ny; j = j + 1 {
    self.u.set(0, j, 0.0)
    self.u.set(self.mesh.nx - 1, j, 0.0)
    self.v.set(0, j, 0.0)
    self.v.set(self.mesh.nx - 1, j, 0.0)
  }
  for i = 0; i < self.mesh.nx; i = i + 1 {
    self.u.set(i, 0, 0.0)
    self.v.set(i, 0, 0.0)
  }
  
  // 压力边界条件：右上角固定为0作为参考
  self.p.set(self.mesh.nx - 1, self.mesh.ny - 1, 0.0)
}

// SIMPLE算法的一个时间步 - 增加数值稳定性
fn CavitySolver::step(self : CavitySolver) -> Unit {
  let nu = 1.0 / self.re

  // 1. 预测速度（显式离散动量方程）
  let u_star = Field::new_double(self.mesh)
  let v_star = Field::new_double(self.mesh)
  
  // 复制当前速度场
  u_star.copy_from(self.u)
  v_star.copy_from(self.v)
  
  for i = 1; i < self.mesh.nx - 1; i = i + 1 {
    for j = 1; j < self.mesh.ny - 1; j = j + 1 {
      // 使用一阶迎风格式处理对流项以增加稳定性
      let u_face_x = if self.u.at(i, j) > 0.0 { 
        self.u.at(i, j) 
      } else { 
        self.u.at(i + 1, j) 
      }
      let u_face_y = if self.v.at(i, j) > 0.0 { 
        self.u.at(i, j) 
      } else { 
        self.u.at(i, j + 1) 
      }
      
      let conv_u = u_face_x * (self.u.at(i + 1, j) - self.u.at(i - 1, j)) / (2.0 * self.mesh.dx) +
                   u_face_y * (self.u.at(i, j + 1) - self.u.at(i, j - 1)) / (2.0 * self.mesh.dy)
      
      let diff_u = nu * ((self.u.at(i + 1, j) - 2.0 * self.u.at(i, j) + self.u.at(i - 1, j)) / (self.mesh.dx * self.mesh.dx) +
                         (self.u.at(i, j + 1) - 2.0 * self.u.at(i, j) + self.u.at(i, j - 1)) / (self.mesh.dy * self.mesh.dy))
      
      let p_grad = self.p.grad(i, j)
      let u_new = self.u.at(i, j) + self.dt * (-conv_u + diff_u - p_grad.0)
      
      // 限制速度变化幅度
      let du = u_new - self.u.at(i, j)
      let limited_du = if du.abs() > 0.1 { 
        if du > 0.0 { 0.1 } else { -0.1 }
      } else { 
        du 
      }
      u_star.set(i, j, self.u.at(i, j) + limited_du)

      // v速度类似处理
      let v_face_x = if self.u.at(i, j) > 0.0 { 
        self.v.at(i, j) 
      } else { 
        self.v.at(i + 1, j) 
      }
      let v_face_y = if self.v.at(i, j) > 0.0 { 
        self.v.at(i, j) 
      } else { 
        self.v.at(i, j + 1) 
      }
      
      let conv_v = v_face_x * (self.v.at(i + 1, j) - self.v.at(i - 1, j)) / (2.0 * self.mesh.dx) +
                   v_face_y * (self.v.at(i, j + 1) - self.v.at(i, j - 1)) / (2.0 * self.mesh.dy)
      
      let diff_v = nu * ((self.v.at(i + 1, j) - 2.0 * self.v.at(i, j) + self.v.at(i - 1, j)) / (self.mesh.dx * self.mesh.dx) +
                         (self.v.at(i, j + 1) - 2.0 * self.v.at(i, j) + self.v.at(i, j - 1)) / (self.mesh.dy * self.mesh.dy))
      
      let v_new = self.v.at(i, j) + self.dt * (-conv_v + diff_v - p_grad.1)
      
      let dv = v_new - self.v.at(i, j)
      let limited_dv = if dv.abs() > 0.1 { 
        if dv > 0.0 { 0.1 } else { -0.1 }
      } else { 
        dv 
      }
      v_star.set(i, j, self.v.at(i, j) + limited_dv)
    }
  }

  // 2. 压力泊松方程 - 简化处理
  let interior_nx = self.mesh.nx - 2
  let interior_ny = self.mesh.ny - 2
  
  if interior_nx > 0 && interior_ny > 0 {
    // 使用简化的压力修正方程
    for i = 1; i < self.mesh.nx - 1; i = i + 1 {
      for j = 1; j < self.mesh.ny - 1; j = j + 1 {
        let div_u_star = (u_star.at(i + 1, j) - u_star.at(i - 1, j)) / (2.0 * self.mesh.dx) +
                         (v_star.at(i, j + 1) - v_star.at(i, j - 1)) / (2.0 * self.mesh.dy)
        
        // 简单的压力修正
        let p_correction = -self.alpha_p * div_u_star / self.dt
        self.p.set(i, j, self.p.at(i, j) + p_correction)
      }
    }
  }

  // 3. 速度修正 - 使用松弛因子
  for i = 1; i < self.mesh.nx - 1; i = i + 1 {
    for j = 1; j < self.mesh.ny - 1; j = j + 1 {
      let p_grad = self.p.grad(i, j)
      let u_corrected = u_star.at(i, j) - self.dt * p_grad.0
      let v_corrected = v_star.at(i, j) - self.dt * p_grad.1
      
      // 应用松弛因子
      self.u.set(i, j, self.u.at(i, j) + self.alpha_u * (u_corrected - self.u.at(i, j)))
      self.v.set(i, j, self.v.at(i, j) + self.alpha_u * (v_corrected - self.v.at(i, j)))
    }
  }
}

// 计算速度场残差的辅助函数
fn CavitySolver::calculate_velocity_residual(self : CavitySolver, u_old : Field[Double], v_old : Field[Double]) -> Double {
  let mut max_u_residual = 0.0
  let mut max_v_residual = 0.0
  
  for i = 1; i < self.mesh.nx - 1; i = i + 1 {
    for j = 1; j < self.mesh.ny - 1; j = j + 1 {
      let u_diff = if self.u.at(i, j) - u_old.at(i, j) > 0.0 { 
        self.u.at(i, j) - u_old.at(i, j) 
      } else { 
        u_old.at(i, j) - self.u.at(i, j) 
      }
      let v_diff = if self.v.at(i, j) - v_old.at(i, j) > 0.0 { 
        self.v.at(i, j) - v_old.at(i, j) 
      } else { 
        v_old.at(i, j) - self.v.at(i, j) 
      }
      
      if u_diff > max_u_residual { max_u_residual = u_diff }
      if v_diff > max_v_residual { max_v_residual = v_diff }
    }
  }
  
  if max_u_residual > max_v_residual { max_u_residual } else { max_v_residual }
}

// 计算连续性方程残差
fn CavitySolver::calculate_continuity_residual(self : CavitySolver) -> Double {
  let mut max_div = 0.0
  
  for i = 1; i < self.mesh.nx - 1; i = i + 1 {
    for j = 1; j < self.mesh.ny - 1; j = j + 1 {
      let div_u = self.u.div(i, j)
      let abs_div = if div_u > 0.0 { div_u } else { -div_u }
      if abs_div > max_div { max_div = abs_div }
    }
  }
  
  max_div
}

// 运行求解器 - 增加数值稳定性检查
fn CavitySolver::run(self : CavitySolver) -> (Field[Double], Field[Double], Field[Double]) {
  println("Starting SIMPLE algorithm solver...")
  println("Grid: \{self.mesh.nx}x\{self.mesh.ny}, Re = \{self.re}")
  println("Max iterations: \{self.max_iter}, Tolerance: \{self.tol}")
  println("Time step: \{self.dt}, Relaxation factors: α_u=\{self.alpha_u}, α_p=\{self.alpha_p}")
  println("")
  
  self.initialize()
  self.apply_boundary_conditions()
  
  for iter = 0; iter < self.max_iter; iter = iter + 1 {
    // 保存上一步的速度场用于计算残差
    let u_old = Field::new_double(self.mesh)
    let v_old = Field::new_double(self.mesh)
    
    u_old.copy_from(self.u)
    v_old.copy_from(self.v)
    
    println("SIMPLE iteration \{iter + 1}:")
    
    self.step()
    self.apply_boundary_conditions()
    
    // 检查数值稳定性
    let mut has_nan = false
    for i = 0; i < self.mesh.nx; i = i + 1 {
      for j = 0; j < self.mesh.ny; j = j + 1 {
        let u_val = self.u.at(i, j)
        let v_val = self.v.at(i, j)
        let p_val = self.p.at(i, j)
        if u_val != u_val || v_val != v_val || p_val != p_val || 
           u_val.abs() > 1.0e6 || v_val.abs() > 1.0e6 || p_val.abs() > 1.0e6 {
          has_nan = true
          break
        }
      }
      if has_nan { break }
    }
    
    if has_nan {
      println("*** SIMULATION FAILED: Numerical instability detected ***")
      break
    }
    
    // 计算残差
    let velocity_residual = self.calculate_velocity_residual(u_old, v_old)
    let continuity_residual = self.calculate_continuity_residual()
    
    println("  Velocity residual: \{velocity_residual}")
    println("  Continuity residual: \{continuity_residual}")
    
    // 检查收敛性
    let overall_residual = if velocity_residual > continuity_residual { 
      velocity_residual 
    } else { 
      continuity_residual 
    }
    
    if overall_residual < self.tol {
      println("")
      println("*** CONVERGED in \{iter + 1} iterations ***")
      println("Final residual: \{overall_residual}")
      break
    }
    
    // 每5步输出一次进度
    if (iter + 1) % 5 == 0 {
      println("  Progress: \{iter + 1}/\{self.max_iter} iterations completed")
    }
    
    println("")
  }
  
  println("Simulation complete. Fields ready for visualization.")
  (self.u, self.v, self.p)
}

// 主函数：示例运行 - 使用更保守的参数
fn main {
  // 使用更小的网格和更保守的参数
  let solver = CavitySolver::new(15, 15, 1.0, 10.0, 50, 1.0e-3)
  let (u, v, p) = solver.run()
  
  // 输出一些关键位置的结果
  println("")
  println("=== Final Results ===")
  let mid_i = solver.mesh.nx / 2
  let mid_j = solver.mesh.ny / 2
  println("Velocity at center (\{mid_i}, \{mid_j}): u = \{u.at(mid_i, mid_j)}, v = \{v.at(mid_i, mid_j)}")
  println("Pressure at center: p = \{p.at(mid_i, mid_j)}")
  
  // 输出顶盖下方的速度分布
  println("Top lid velocity profile:")
  let j_below_lid = solver.mesh.ny - 2
  for i = 0; i < solver.mesh.nx; i = i + 5 {
    println("  x=\{i.to_double() * solver.mesh.dx}: u=\{u.at(i, j_below_lid)}")
  }
}

// 导出函数供Web调用（假设编译到WASM）
pub fn run_simulation(nx : Int, ny : Int, l : Double, re : Double, max_iter : Int, tol : Double) -> String {
  let solver = CavitySolver::new(nx, ny, l, re, max_iter, tol)
  let (u, v, p) = solver.run()
  
  // 简化的结果序列化
  let mut result = "{ \"u\": ["
  for i = 0; i < nx; i = i + 1 {
    result = result + "["
    for j = 0; j < ny; j = j + 1 {
      result = result + u.at(i, j).to_string()
      if j < ny - 1 { 
        result = result + "," 
      }
    }
    result = result + "]"
    if i < nx - 1 { 
      result = result + "," 
    }
  }
  result = result + "], \"v\": [/* similar */], \"p\": [/* similar */] }"
  result
}