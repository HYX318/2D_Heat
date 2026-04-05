# CartesianTopology 实现摘要

## 概述
成功实现了 `CartesianTopology` 类及其单元测试，用于管理 MPI Cartesian 拓扑以进行域分解。

## 创建的文件

### 1. 头文件
**路径**: `/Users/galoishuang/Development/2D_Heat/src/mpi/cartesian_topology.hpp`
**大小**: 8.0 KB

**主要特性**:
- Direction 枚举 (X, Y) 和 Shift 枚举 (Forward, Backward)
- NeighborInfo 结构体，存储四个邻居的 ranks
- CartesianTopology 类，提供完整的 Cartesian 拓扑管理功能

**关键方法**:
- 构造函数：自动维度计算和用户指定维度
- 访问器：communicator(), rank(), size(), dims(), coords(), neighbors()
- 维度访问器：dim_x(), dim_y(), coord_x(), coord_y()
- 边界检查：is_on_boundary(), has_neighbor()
- 移动语义：支持移动构造和移动赋值

### 2. 实现文件
**路径**: `/Users/galoishuang/Development/2D_Heat/src/mpi/cartesian_topology.cpp`
**大小**: 9.6 KB

**实现细节**:
- 使用 `MPI_Cart_create` 创建 Cartesian communicator
- 使用 `MPI_Cart_coords` 获取进程坐标
- 使用 `MPI_Cart_shift` 计算邻居 ranks
- 自动计算最优维度（最接近正方形）
- 正确处理边界进程（MPI_PROC_NULL）
- 异常资源管理（RAII 模式）

### 3. 测试文件
**路径**: `/Users/galoishuang/Development/2D_Heat/tests/unit/test_cartesian_topology.cpp`
**大小**: 15 KB

## 测试覆盖

### 测试用例 (11 个测试)

1. **InitializationTest** - 测试拓扑创建
   - 验证可以成功创建 CartesianTopology
   - 验证 communicator 有效
   - 验证 rank 和 size 有效

2. **DimensionCalculationTest** - 测试自动维度计算
   - 验证最优维度计算正确
   - 验证维度乘积等于进程数
   - 验证维度接近正方形

3. **SpecifiedDimensionsTest** - 测试用户指定维度
   - 验证接受用户指定的维度
   - 验证维度匹配输入
   - 验证无效维度抛出异常

4. **CoordinateMappingTest** - 测试坐标映射
   - 验证每个进程获得唯一坐标
   - 验证坐标在有效范围内
   - 验证坐标访问器正确

5. **NeighborComputationTest** - 测试邻居计算
   - 验证每个进程获得正确的邻居 ranks
   - 验证边界进程使用 MPI_PROC_NULL
   - 验证进程不是自己的邻居

6. **BoundaryDetectionTest** - 测试边界检测
   - 验证 is_on_boundary() 正确识别边界
   - 验证 has_neighbor() 是 is_on_boundary 的反函数
   - 验证边界进程的邻居为 MPI_PROC_NULL

7. **EdgeCasesTest** - 测试边缘情况
   - 测试 1xN 拓扑
   - 测试 Nx1 拓扑
   - 验证所有进程在正确边界上

8. **SquareTopologyTest** - 测试正方形拓扑
   - 测试 2x2 和 3x3 拓扑
   - 验证角进程有 2 个邻居
   - 验证边界进程有 3 个邻居
   - 验证内部进程有 4 个邻居

9. **NoCopyTest** - 测试禁用复制语义
   - 验证复制构造函数被删除
   - 验证复制赋值被删除

10. **MoveTest** - 测试移动语义
    - 验证移动构造正确转移所有权
    - 验证移动赋值正确转移所有权
    - 验证移动后对象处于有效但无效状态

11. **DestructorTest** - 测试析构函数
    - 验证析构函数正确释放 communicator
    - 验证可以创建和销毁多个拓扑对象

12. **MultipleTopologiesTest** - 测试多个拓扑
    - 验证可以顺序创建多个拓扑
    - 验证每个拓扑独立

## 测试运行配置

### CMake 配置
已在 CMakeLists.txt 中添加以下测试配置：

- **CartesianTopology_2x2**: 4 个进程 (2x2 拓扑)
- **CartesianTopology_3x3**: 9 个进程 (3x3 拓扑)
- **CartesianTopology_2x3**: 6 个进程 (2x3 拓扑)

### 测试示例

#### 2x2 拓扑 (4 进程)
```
Rank 0: [0, 0] Dimensions: [2, 2]  Neighbors: S=-1 N=1 W=-1 E=2  Boundaries: W=Y E=N S=Y N=N
Rank 1: [0, 1] Dimensions: [2, 2]  Neighbors: S=0 N=-1 W=-1 E=3  Boundaries: W=Y E=N S=N N=Y
Rank 2: [1, 0] Dimensions: [2, 2]  Neighbors: S=-1 N=3 W=0 E=-1  Boundaries: W=N E=Y S=Y N=N
Rank 3: [1, 1] Dimensions: [2, 2]  Neighbors: S=2 N=-1 W=1 E=-1  Boundaries: W=N E=Y S=N N=Y
```

#### 3x3 拓扑 (9 进程)
```
Rank 0: [0, 0] Dimensions: [3, 3]  Neighbors: S=-1 N=1 W=-1 E=3  Boundaries: W=Y E=N S=Y N=N
Rank 4: [1, 1] Dimensions: [3, 3]  Neighbors: S=3 N=5 W=1 E=7  Boundaries: W=N E=N S=N N=N
```

## 主要特性

### 1. 自动维度计算
- 自动计算最接近正方形的维度
- 例如：4 进程 -> [2, 2], 6 进程 -> [2, 3], 9 进程 -> [3, 3]

### 2. 坐标管理
- 自动映射 rank 到坐标
- 提供便捷的访问器：coord_x(), coord_y()

### 3. 邻居识别
- 自动计算四个方向的邻居 ranks
- 边界进程使用 MPI_PROC_NULL
- 支持边界检测和邻居查询

### 4. 异常安全
- RAII 模式自动管理资源
- 移动语义支持资源转移
- 正确处理 MPI 错误

### 5. 内存管理
- 正确释放 Cartesian communicator
- 避免内存泄漏
- 支持多个拓扑对象

## 使用示例

```cpp
#include "mpi/cartesian_topology.hpp"

// 自动计算维度
CartesianTopology topology(MPI_COMM_WORLD);

// 指定维度
CartesianTopology topology(MPI_COMM_WORLD, {2, 3});

// 访问信息
int rank = topology.rank();
int size = topology.size();
int dim_x = topology.dim_x();
int dim_y = topology.dim_y();
int coord_x = topology.coord_x();
int coord_y = topology.coord_y();

// 邻居信息
const auto& neighbors = topology.neighbors();
int north = neighbors.north;
int south = neighbors.south;
int east = neighbors.east;
int west = neighbors.west;

// 边界检测
bool is_west_boundary = topology.is_on_boundary(Direction::X, Shift::Backward);
bool has_east_neighbor = topology.has_neighbor(Direction::X, Shift::Forward);
```

## 依赖

- MPI: 用于 Cartesian 拓扑功能
- MPIContext: 用于 MPI 初始化（可选）
- C++17: 用于标准库功能

## 构建和测试

### 使用 CMake
```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTS=ON -DENABLE_MPI=ON ..
make
./bin/test_cartesian_topology  # 单进程测试
mpirun -np 4 ./bin/test_cartesian_topology  # 多进程测试
```

### 使用 Make (简单测试)
```bash
mpic++ -Wall -g test_cartesian_topology_simple.cpp src/mpi/cartesian_topology.cpp -o test_cartesian_topology_simple -I src
export TMPDIR=/tmp
mpirun -np 4 ./test_cartesian_topology_simple
```

## 验证结果

### 手动测试结果
- ✅ 4 进程 (2x2): 所有测试通过
- ✅ 9 进程 (3x3): 所有测试通过
- ✅ 6 进程 (2x3): 所有测试通过

### 验证的功能
- ✅ 正确创建 Cartesian 拓扑
- ✅ 自动计算最优维度
- ✅ 正确映射 rank 到坐标
- ✅ 正确计算邻居 ranks
- ✅ 正确识别边界进程
- ✅ 正确处理 MPI_PROC_NULL
- ✅ 移动语义工作正常
- ✅ 异常安全
- ✅ 无内存泄漏

## 总结

CartesianTopology 类已成功实现并通过测试。该类提供了完整的 MPI Cartesian 拓扑管理功能，包括：
- 自动和手动维度配置
- 坐标映射和邻居识别
- 边界检测
- 异常安全的资源管理
- 全面的单元测试覆盖

该实现准备好用于 2D 热方程求解器的域分解功能。
