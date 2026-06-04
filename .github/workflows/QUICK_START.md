# GitHub Actions CI/CD - 快速参考

## 文件结构

```
.github/
└── workflows/
    ├── ci.yml           # 主 CI 流水线
    ├── benchmark.yml     # 性能基准测试
    ├── codeql.yml        # 安全扫描
    ├── README.md         # 详细文档
    └── QUICK_START.md    # 本文件
```

## 支持的平台和编译器

### 测试矩阵

| 平台 | 编译器 | MPI 实现 | 构建类型 |
|------|--------|----------|----------|
| Ubuntu Latest | GCC 11 | OpenMPI | Release, Debug |
| Ubuntu Latest | GCC 11 | MPICH | Release, Debug |
| Ubuntu Latest | Clang 14 | OpenMPI | Release, Debug |
| Ubuntu Latest | Clang 14 | MPICH | Release, Debug |
| macOS Latest | Apple Clang | OpenMPI | Release, Debug |

**总计：8 种组合 × 2 种构建类型 = 16 个作业**

## 工作流触发条件

### CI 流水线 (`ci.yml`)
- ✅ Push 到 `master`/`main` 分支
- ✅ Pull Requests
- ✅ 手动触发

### 性能基准 (`benchmark.yml`)
- 📅 每月 1 日 00:00 UTC
- ✅ 手动触发

### 安全扫描 (`codeql.yml`)
- ✅ Push 到 `master`/`main` 分支
- ✅ Pull Requests
- 📅 每周日 00:00 UTC
- ✅ 手动触发

## 手动触发 CI

### Web 界面

1. 进入仓库 → Actions 标签
2. 选择工作流 → 点击 "Run workflow"
3. 选择分支 → 填写参数 → 运行

### GitHub CLI

```bash
# 触发 CI
gh workflow run ci.yml

# 触发基准测试（带参数）
gh workflow run benchmark.yml \
  -f grid_size="200x200" \
  -f time_steps="5000" \
  -f num_processes="8"

# 触发安全扫描
gh workflow run codeql.yml
```

### GitHub API

```bash
curl -X POST \
  -H "Authorization: token YOUR_TOKEN" \
  -H "Accept: application/vnd.github.v3+json" \
  https://api.github.com/repos/OWNER/REPO/actions/workflows/ci.yml/dispatches \
  -d '{"ref":"master"}'
```

## 快速命令

```bash
# 查看所有工作流
gh workflow list

# 查看工作流运行历史
gh run list

# 查看特定工作流的运行
gh run list --workflow=ci.yml

# 查看运行详情
gh run view RUN_ID

# 下载运行日志
gh run view RUN_ID --log

# 取消运行
gh run cancel RUN_ID

# 重新运行失败的工作流
gh run rerun RUN_ID
```

## 本地测试

使用 `act` 在本地测试 GitHub Actions：

```bash
# 安装 act
brew install act  # macOS

# 列出工作流
act -l

# 运行 CI
act push -W .github/workflows/ci.yml

# 运行特定作业
act push -j build-and-test -W .github/workflows/ci.yml
```

## 关键特性

### ✅ 已配置的功能

- [x] 多平台测试 (Ubuntu, macOS)
- [x] 多编译器支持 (GCC, Clang, Apple Clang)
- [x] 多 MPI 实现 (OpenMPI, MPICH)
- [x] 构建缓存
- [x] 代码覆盖率 (gcov/lcov)
- [x] 静态分析 (clang-tidy)
- [x] 安全扫描 (CodeQL)
- [x] 性能基准测试
- [x] 构建产物上传
- [x] 超时保护 (30 分钟)
- [x] macOS 特殊处理 (TMPDIR)

### 📦 构建产物

- **CI 产物**：`Heat-{os}-{compiler}-{mpi}` (保留 7 天)
- **基准结果**：`benchmark-results-{os}-{run-number}` (保留 30 天)

### ⏱️ 超时设置

- CI: 30 分钟
- Benchmark: 60 分钟
- CodeQL: 30 分钟

## 故障排除

### macOS 构建失败

确保设置了 `TMPDIR=/tmp`（已在 ci.yml 中配置）

### MPI 测试失败

检查 MPI 安装和路径：
```bash
mpirun --version
mpicc --showme:version
```

### 超时问题

- 检查测试代码是否有无限循环
- 优化编译参数
- 减少 `-j` 并行数

## 更多信息

- 📖 详细文档：`.github/workflows/README.md`
- 🔄 GitHub Actions 文档：https://docs.github.com/en/actions
- 🔒 CodeQL 文档：https://docs.github.com/en/code-security
