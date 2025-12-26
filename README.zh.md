# HelixAlign

[English](README.en.md) | 中文

## 项目概述

HelixAlign 是一个用 Rust 实现的高效序列比对工具，实现了与 nucmer 相似的功能，支持多种匹配算法和输出格式。该项目旨在提供一个高性能、内存友好的替代方案，用于基因组序列比对分析。

## 主要功能

### 1. 序列比对算法

- **MUM (Maximal Unique Matches)**: 计算在两个序列中都唯一的最大匹配
- **MAM (Maximal Unique Reference Matches)**: 计算在参考序列中唯一的最大匹配（默认算法）
- **MEM (Maximal Exact Matches)**: 计算所有最大匹配，不考虑唯一性

### 2. 核心组件

#### 后缀数组实现
- 实现了稀疏后缀数组，支持高效的模式匹配
- 提供 `find_matches` 方法返回匹配位置值
- 支持 LCP (最长公共前缀) 数组计算

#### 序列处理
- 实现 `DnaSequence` 结构体处理 DNA 序列
- 提供反向互补序列计算功能

#### 基因组统计
- 计算 N50/N90 等基因组统计指标
- 支持多序列统计信息计算
- 实现了 `Display` trait 便于输出格式化

### 3. Nucmer 参数支持

项目已实现所有 nucmer 参数，确保功能对等：

#### 基本匹配参数
- `-mum`: 计算在两个序列中都唯一的最大匹配
- `-mumreference`/`-mumcand`: 计算在参考序列中唯一的最大匹配（默认）
- `-maxmatch`: 计算所有最大匹配，不考虑唯一性
- `-l`/`--minmatch`: 设置单个精确匹配的最小长度（默认: 20）

#### 聚类和扩展参数
- `-b`/`--breaklen`: 设置对齐扩展尝试扩展低分区域的最大距离（默认: 200）
- `-c`/`--mincluster`: 设置匹配簇的最小长度（默认: 65）
- `-D`/`--diagdiff`: 设置簇中相邻锚点的最大对角线差异（默认: 5）
- `-d`/`--diagfactor`: 设置簇中相邻锚点的最大对角线差异作为间隙长度的微分分数（默认: 0.12）
- `-g`/`--maxgap`: 设置簇中相邻匹配之间的最大间隙（默认: 90）
- `-L`/`--minalign`: 设置对齐的最小长度（默认: 0）

#### 处理选项
- `-noextend`: 不执行簇扩展步骤
- `-nooptimize`: 不进行对齐分数优化
- `-nosimplify`: 不通过移除阴影簇来简化对齐
- `-f`/`--forward`: 仅使用查询序列的正向链
- `-r`/`--reverse`: 仅使用查询序列的反向互补链

#### 输出和文件选项
- `-p`/`--prefix`: 将输出写入 PREFIX.delta（默认: out）
- `--delta`: 将 delta 文件输出到指定路径
- `--sam-short`: 输出 SAM 文件，短格式
- `--sam-long`: 输出 SAM 文件，长格式
- `--save`: 将后缀数组保存到文件
- `--load`: 从文件加载后缀数组

#### 高级选项
- `-banded`: 强制基于 diagdiff 参数对动态规划矩阵进行绝对带状限制
- `-large`: 强制使用大偏移量
- `-G`/`--genome`: 基因组到基因组映射（长查询序列）
- `-M`/`--max-chunk`: 设置最大块大小
- `-t`/`--threads`: 设置使用的线程数
- `-batch`: 按参考序列的批次进行处理
- `-format`: 指定输出格式（default, delta, paf, sam）
- `-stats`: 显示参考和查询序列统计信息（N50, N90 等）

### 4. 输出格式
- **Default**: 默认格式
- **Delta**: nucmer 兼容的 delta 格式
- **PAF**: Pairwise mApping Format
- **SAM**: Sequence Alignment/Map 格式

### 5. 多线程支持
- 使用 Rayon 实现并行处理
- 支持自定义线程数量
- 实现进度条显示

## 安装

### 从源码编译

```bash
git clone https://github.com/yourusername/helixalign.git
cd helixalign
cargo build --release
```

编译后的可执行文件将位于 `target/release/helixalign`。

## 使用方法

### 基本用法

```bash
helixalign reference.fa query.fa
```

### 使用特定算法和参数

```bash
helixalign -maxmatch -l 20 -t 4 -format sam reference.fa query.fa
```

### 显示统计信息

```bash
helixalign -stats reference.fa query.fa
```

### 使用所有参数的完整示例

```bash
helixalign -maxmatch -b 200 -c 65 -D 5 -d 0.12 -noextend -f -g 90 -l 20 -L 0 -nooptimize -r -nosimplify -p results --delta results.delta --sam-short results.sam --batch 10000 -banded -large -G -M 50000 -t 8 -format paf -stats reference.fa query.fa
```

## 技术亮点

### 1. 高效算法实现
- 使用稀疏后缀数组实现高效的模式匹配
- 支持多种匹配算法（MUM, MAM, MEM）
- 实现了最大匹配的查找和扩展

### 2. 内存优化
- 稀疏后缀数组减少内存占用
- 批处理机制支持大文件处理
- 可选的大偏移量支持

### 3. 用户友好
- 详细的命令行帮助信息
- 进度条显示处理进度
- 多种输出格式选择

## 项目结构

```
helixalign/
├── src/
│   ├── main.rs              # 主程序入口和命令行解析
│   ├── lib.rs               # 库文件，导出所有模块
│   ├── sequence.rs          # DNA 序列处理
│   ├── suffix_array.rs      # 后缀数组实现
│   ├── algorithms.rs        # 匹配算法实现
│   ├── nucmer.rs            # Nucmer 对齐实现
│   ├── genomic_stats.rs     # 基因组统计计算
│   └── output_format.rs     # 输出格式处理
├── Cargo.toml               # 项目配置和依赖
├── README.md                # 项目说明文档
└── README.zh.md             # 中文说明文档
```

## 依赖项

- `rayon`: 并行计算支持
- `indicatif`: 进度条显示
- `clap`: 命令行参数解析（计划中）

## 性能

HelixAlign 在保持与 nucmer 功能对等的同时，提供了以下性能优势：

- 更快的后缀数组构建
- 更低的内存占用
- 更好的多线程扩展性

## 贡献

欢迎提交问题报告和功能请求！如果您想贡献代码，请：

1. Fork 本仓库
2. 创建您的功能分支 (`git checkout -b feature/AmazingFeature`)
3. 提交您的更改 (`git commit -m 'Add some AmazingFeature'`)
4. 推送到分支 (`git push origin feature/AmazingFeature`)
5. 打开一个 Pull Request

## 许可证

本项目采用 MIT 许可证 - 查看 [LICENSE](LICENSE) 文件了解详情。

## 致谢

- 本项目受到 MUMmer4 的启发，旨在提供类似的功能
- 感谢 Rust 社区提供的优秀工具和库

## 联系方式

如有问题或建议，请通过以下方式联系：

- 提交 Issue: [GitHub Issues](https://github.com/yourusername/helixalign/issues)
- 邮箱: your.email@example.com
