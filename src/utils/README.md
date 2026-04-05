# Timer Class

高性能计时器类，用于精确的性能测量和性能分析。

## 特性

- **高精度计时**：使用 `std::chrono::high_resolution_clock` 提供纳秒级精度
- **命名计时器**：支持给计时器命名，便于识别不同的计时区域
- **Lap 计时**：支持记录多个时间点，测量子过程耗时
- **分层计时**：支持嵌套计时器，便于分析调用层次
- **RAII 语义**：提供 `ScopedTimer`，在作用域结束时自动停止并打印时间
- **多种时间单位**：支持秒、毫秒、微秒三种时间单位
- **线程安全**：每个计时器实例独立，可安全地在多线程中使用

## 文件

- `timer.hpp` - Timer 类头文件（仅需包含此文件）
- `timer_example.cpp` - 使用示例
- `Makefile` - 编译示例的 Makefile

## 编译和运行示例

```bash
cd src/utils
make
./timer_example
```

## 基本使用

### 1. 基本计时

```cpp
#include "timer.hpp"

Timer total_timer("Total Time");
// ... do work ...
total_timer.stop();
std::cout << "Total time: " << total_timer.elapsed() << " seconds" << std::endl;
```

### 2. 手动计时循环

```cpp
Timer step_timer;
for (int i = 0; i < 100; ++i) {
    step_timer.start();
    // ... do work ...
    step_timer.stop();
    std::cout << "Step time: " << step_timer.elapsed_ms() << " ms" << std::endl;
}
```

### 3. Lap 计时

```cpp
Timer loop_timer;
for (int i = 0; i < 10; ++i) {
    // ... do work ...
    double lap_time = loop_timer.lap();
    std::cout << "Lap " << i << ": " << lap_time << " s" << std::endl;
}
```

### 4. 作用域计时器（自动停止和打印）

```cpp
{
    ScopedTimer timer("Computation");
    // ... do work ...
} // timer automatically stops here and prints elapsed time
```

### 5. 作用域计时器（手动打印）

```cpp
{
    ScopedTimer timer("Processing", false); // false = 不自动打印
    // ... do work ...
    std::cout << "Time: " << timer.elapsed_ms() << " ms" << std::endl;
}
```

### 6. 嵌套计时

```cpp
Timer outer_timer("Outer");
outer_timer.start();

{
    Timer inner_timer("Inner");
    inner_timer.start();
    // ... inner work ...
    inner_timer.stop();
}

outer_timer.stop();
```

## API 参考

### Timer 类

#### 构造函数
- `Timer()` - 创建无名计时器并立即开始
- `Timer(const std::string& name)` - 创建命名计时器并立即开始

#### 计时控制
- `void start()` - 开始或重新开始计时
- `void stop()` - 停止计时
- `void reset()` - 重置计时器状态
- `bool is_running() const` - 检查是否正在计时

#### 时间获取
- `double elapsed() const` - 获取经过的秒数
- `double elapsed_ms() const` - 获取经过的毫秒数
- `double elapsed_us() const` - 获取经过的微秒数

#### 命名和统计
- `void set_name(const std::string& name)` - 设置计时器名称
- `const std::string& name() const` - 获取计时器名称
- `size_t lap_count() const` - 获取 lap 次数

#### Lap 计时
- `double lap()` - 记录一个 lap 点并返回时间（秒）
- `void reset_lap_count()` - 重置 lap 计数

### ScopedTimer 类

#### 构造函数
- `ScopedTimer(const std::string& name, bool auto_print = true)` - 创建作用域计时器
- `ScopedTimer(Timer& timer, bool auto_print = true)` - 使用现有计时器

#### 方法
- `Timer& timer()` - 获取底层计时器引用
- `double elapsed() const` - 获取经过的秒数
- `double elapsed_ms() const` - 获取经过的毫秒数
- `double elapsed_us() const` - 获取经过的微秒数

## 性能特点

- 使用 `std::chrono::high_resolution_clock` 提供最高可用精度
- 支持纳秒级时间测量（取决于平台）
- 轻量级实现，对性能影响极小
- 无动态内存分配，适合高频调用场景

## 注意事项

1. 每个计时器实例是独立的，不是线程安全的共享计时器
2. `ScopedTimer` 禁止拷贝但支持移动
3. 多次调用 `start()` 会重新开始计时
4. 计时器停止后仍可获取经过时间
5. `lap()` 会在计时器未运行时自动启动它

## 编译要求

- C++17 或更高版本
- 需要 `<chrono>`、`<string>`、`<iostream>` 标准库支持

## 许可

此代码作为项目的一部分，遵循项目的主许可证。
