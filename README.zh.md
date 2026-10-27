# HelixAlign 使用说明（正确性修复分支）

本分支实现 DNA 的 MUM、MAM、MEM **最大精确匹配种子**。它可以作为小规模数据的
额外独立 CPU 正确性对照。它不是 GPU 软件，也没有实现完整 NUCmer 的聚类和
有缺口延伸。原版问题、实际运行结果和最小复现见
[独立审核报告](audit/CORRECTNESS_AUDIT_20261027.zh.md)。

## 编译与基本用法

需要支持 Rust edition 2024 的工具链。manifest 的最低版本为 1.85，
本次实际使用 1.99.1，未额外验证每个更早版本。

~~~bash
cargo test --locked
cargo build --release --locked

# MUM：参考全部记录合计唯一，并在当前 query 记录内唯一。
target/release/helixalign -mum -l 20 -b --output-format matches reference.fa query.fa

# MAM：仅要求参考唯一；默认只处理 query 正向。
target/release/helixalign -mumreference -l 20 --output-format matches reference.fa query.fa

# MEM：保留所有左右均不可延伸的精确匹配，包括不同对角线的合法短匹配。
target/release/helixalign -maxmatch -l 20 -r --format paf reference.fa query.fa
~~~

-l/--minmatch 为正整数；-b/--both 为双向，-r/--reverse 为反向，
--forward 为正向。-t/--threads 设定 query 记录处理的线程数。
-f、--format、--output-format 为同一输出格式选项。
-stats/--stats 将统计输出到 stderr，不污染匹配文件。
-n、-F 被接受，ACGT-only 和参考 ID 输出始终启用。

## 记录、字符与坐标

参考唯一性跨全部 reference records 计数；query 唯一性仅在当前记录内判断，
重叠出现也计数。两个方向分别处理；回文匹配可在两条链上各输出一条。

只允许 A/C/G/T 参与匹配；N 和其他接受的 IUPAC 字符为不可跨越的屏障，
不作为通配符。FASTA 记录不直接连接，坐标不会跨 contig。
ID 使用 header 的第一个 token。空记录、缺失 header、非法字符、重复 ID 会报错；
多个 query 文件之间也要求 ID 唯一。

推荐输出 matches，其六列为：

~~~text
ref_id  query_id  strand  ref_start  query_start  length
~~~

实际使用 TAB 分隔。两个起点均为**原始正向记录上的零起始左端坐标**，
反向命中也保持这个坐标规则，并用 strand=- 标注。
例如 reference=AACCGGTA、query=CCTACCGGTTG、阈值8，
反向输出为 r q - 0 2 8。

PAF 使用真实 ID、每条记录的长度、零起始半开区间和正确方向；
MAPQ=255 表示未知。PAF 行表示精确种子，不代表完整 read alignment。
默认人类可读格式使用一起始左端位置，并明确打印方向。

## 明确不支持及兼容性变化

- 不输出原先错误的 SAM/delta；这些格式请求返回非零状态。
- 名为 nucmer/nucmer-rs 的 CLI 入口被拒绝，避免将种子当作完整比对。
- 不支持蛋白质 CLI、稀疏后缀采样或 NUCmer 聚类/延伸参数。
- SparseSuffixArray 只接受 k=1，名称保留不代表已经实现稀疏采样。
- print_matches_in_format 已替换为带记录上下文的 write_record_matches。
- NucmerAligner::align 返回 Vec<RecordMatch>，保留方向；
  align_parallel 返回 Result。历史聚类字段仅保留默认构造兼容性，
  改成其他值会明确报错；它们不执行相应算法。
- 正向使用 --forward，旧 -f 的输出格式意义保留。
- 错误输入或参数返回状态2；不会把未知格式静默改成默认格式。

## 已执行验证与未验证范围

原版8项测试通过，但432次正向差异比较失败298次。
修复版39项 Rust 测试全部通过；24对定向小样例和120对固定种子随机样例，
分别以3种算法、正向/反向/双向比较，共1296次，与独立穷举预期及官方
MUMmer4.0.1均一致，FP/FN均为0。

完整方法和证据见 [audit/REPRODUCE.md](audit/REPRODUCE.md)。
官方对照阈值至少2且不使用-c；该版本的长度1和反向-c边界问题另有记录，
长度1使用独立Rust穷举检查。不能把有限测试通过解释成所有参数完全等价。

仍使用完整后缀比较排序和朴素LCP构建，长重复可能产生很大的计算和内存代价。
尚未完成全基因组、资源耗尽或正式性能测试；也没有测量GPU加速比。
GPU性能主基线应使用固定版本、同样任务和输出规则的成熟官方CPU实现。
