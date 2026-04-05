# GitHub Actions CI/CD Workflows

本项目配置了完整的 GitHub Actions CI/CD 流水线，用于自动化构建、测试和代码质量检查。

## 工作流文件

### 1. `ci.yml` - 主 CI 流水线

**触发条件：**
- Push 到 `master` 或 `main` 分支
- 创建或更新 Pull Request
- 手动触发（通过 GitHub Actions 页面）

**测试矩阵：**

| OS | Compiler | MPI | Build Type |
|----|----------|-----|-------------|
| ubuntu-latest | gcc-11 | openmpi | Release, Debug |
| ubuntu-latest | gcc-11 | mpich | Release, Debug |
| ubuntu-latest | clang-14 | openmpi | Release, Debug |
| ubuntu-latest | clang-14 | mpich | Release, Debug |
| macos-latest | apple-clang | openmpi | Release, Debug |

**工作步骤：**

1. **检出代码** - 使用 `actions/checkout@v3`
2. **环境设置** - 安装 MPI 和编译器
   - macOS: `brew install openmpi`
   - Ubuntu: 安装 OpenMPI、MPICH、G++、Clang、Google Test
3. **配置 CMake** - 支持自定义编译器和构建类型
4. **构建** - 使用 CMake 构建所有目标
5. **测试** - 运行单元测试和 MPI 测试（`mpirun -np 4`）
6. **代码覆盖率** - 对 GCC Debug 构建生成覆盖率报告
7. **静态分析** - 使用 clang-tidy 进行代码分析
8. **上传构建产物** - 保存 Release 版本的二进制文件

**特殊处理：**

- macOS: 设置 `TMPDIR=/tmp` 避免共享内存问题
- 超时: 30 分钟
- 缓存: 使用 `actions/cache@v3` 加速构建

### 2. `benchmark.yml` - 性能基准流水线

**触发条件：**
- 每月 1 日 00:00 UTC 自动运行
- 手动触发（支持自定义参数）

**参数（手动触发时可选）：**
- `grid_size`: 网格大小（默认：100x100）
- `time_steps`: 时间步数（默认：1000）
- `num_processes`: MPI 进程数（默认：4）

**工作步骤：**

1. **构建 Release 版本**
2. **运行基准测试** - 测试不同进程数（1, 2, 4, 8）的性能
3. **分析性能** - 提取并总结执行时间
4. **上传结果** - 保存基准测试结果为 artifact（保留 30 天）
5. **PR 评论** - 如果由 PR 触发，自动在 PR 中添加结果评论

### 3. `codeql.yml` - CodeQL 安全扫描

**触发条件：**
- Push 到 `master` 或 `main` 分支
- 创建或更新 Pull Request
- 每周日 00:00 UTC 自动运行
- 手动触发

**扫描类型：**
- security-extended: 扩展的安全查询
- security-and-quality: 安全和质量查询

**扫描语言：**
- C++

## 如何手动触发 CI

### 通过 GitHub Web 界面

1. 进入 GitHub 仓库页面
2. 点击 "Actions" 标签
3. 选择要运行的工作流
4. 点击 "Run workflow" 按钮
5. 选择分支并填写参数（如果需要）
6. 点击 "Run workflow" 确认

### 通过 GitHub CLI (gh)

```bash
# 触发 CI 工作流
gh workflow run ci.yml

# 触发性能基准测试（带参数）
gh workflow run benchmark.yml \
  -f grid_size="200x200" \
  -f time_steps="5000" \
  -f num_processes="8"

# 触发 CodeQL 扫描
gh workflow run codeql.yml
```

## 构建产物

### CI 工作流产物

每个 Release 构建会生成以下产物：
- 二进制文件：`Heat`
- 输出文件：`*.txt`
- 命名格式：`Heat-{os}-{compiler}-{mpi}`
- 保留时间：7 天

### 基准测试产物

每次运行会生成：
- `benchmark_1_processes.txt` - 单进程结果
- `benchmark_2_processes.txt` - 双进程结果
- `benchmark_4_processes.txt` - 四进程结果
- `benchmark_8_processes.txt` - 八进程结果
- `summary.md` - 性能总结报告
- 保留时间：30 天

## 代码覆盖率

CI 工作流会为 GCC Debug 构成生成代码覆盖率报告，并自动上传到 Codecov（如果配置）。

要配置 Codecov，需要：

1. 在 [codecov.io](https://codecov.io) 注册账户
2. 添加仓库并获取 token
3. 在 GitHub 仓库设置中添加 Secret：
   - Name: `CODECOV_TOKEN`
   - Value: 你的 Codecov token

## 故障排除

### macOS 构建失败

如果遇到共享内存错误，确保工作流中设置了 `TMPDIR=/tmp`。这已在 `ci.yml` 中配置。

### MPI 测试失败

确保 MPI 正确安装并配置：
- 检查 `mpicc` 和 `mpirun` 可用
- 验证 MPI 库路径正确
- 检查防火墙设置（本地测试）

### 超时问题

如果构建经常超时：
- 检查是否有无限循环
- 减少 `-j` 并行编译参数
- 优化测试代码执行时间

## 本地测试 CI 配置

使用 [act](https://github.com/nektos/act) 在本地测试 GitHub Actions：

```bash
# 安装 act
brew install act  # macOS
# 或
curl https://raw.githubusercontent.com/nektos/act/master/install.sh | sudo bash

# 运行所有工作流
act -l  # 列出工作流
act push  # 模拟 push 事件
act workflow_dispatch -W .github/workflows/ci.yml  # 手动触发 CI
```

## 自定义配置

### 添加新的编译器

在 `ci.yml` 的 `matrix.compiler` 中添加：

```yaml
compiler: [gcc-11, clang-14, gcc-12, clang-15]
```

### 添加新的 MPI 实现

在 `ci.yml` 的 `matrix.mpi` 中添加：

```yaml
mpi: [openmpi, mpich, intel-mpi]
```

### 修改缓存策略

调整缓存键以匹配你的需求：

```yaml
key: ${{ runner.os }}-${{ matrix.compiler }}-${{ hashFiles('**/*.cpp') }}
```

## 监控和维护

### 查看构建历史

GitHub Actions 页面显示所有工作流运行历史，包括：
- 构建状态（成功/失败/取消）
- 运行时间
- 触发者信息
- 日志输出

### 设置构建通知

在 GitHub 仓库设置中配置：
- Email 通知
- Slack/Discord 集成（通过第三方 Actions）
- Webhook 集成

### 清理旧构建产物

定期清理旧的构建产物以节省存储空间。GitHub 自动删除超过保留时间的产物。

## 最佳实践

1. **频繁提交** - 小步提交让 CI 更容易定位问题
2. **编写测试** - 每个功能都应该有对应的测试
3. **关注警告** - 修复编译警告，提高代码质量
4. **保持依赖更新** - 定期更新 Actions 和依赖版本
5. **监控性能** - 定期检查基准测试结果
6. **安全第一** - 定期运行 CodeQL 扫描并修复问题

## 联系和支持

如有问题或建议：
- 提交 Issue 到仓库
- 查看 [GitHub Actions 文档](https://docs.github.com/en/actions)
- 参考 [CMake 文档](https://cmake.org/documentation/)
