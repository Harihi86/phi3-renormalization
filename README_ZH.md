# φ³ 重整化：实现、原理与精确验证

[English](README.md) · [原理与源码对应](docs/theory-zh.md) · [质量／动能项演示](docs/L2E3G4-zh.md)

本仓库对应爱丁堡大学 MPhys 理论物理项目 *Renormalization in phi-cubed theory*。
**0.2.0** 引入当前 Phi3Lab 计算核心：生成泛函计数、图同构、精确有理系数、BPH 递归减除、
反项作用量独立生成和逐格验证，同时保留已有数值核，补入数学研究源码与可编辑 LaTeX 推导。

## 安装并打开中文工作台

需要 Python **3.11+**。精确组合计算和一圈有限积分示例只需标准库。

```bash
python -m venv .venv
# 按操作系统的常规方式激活 .venv。
python -m pip install ".[dev]"
phi3-ct build --output runs/counterterms
phi3-ct verify --input runs/counterterms
phi3-ct serve --input runs/counterterms --language zh
```

打开打印的本机地址；英文界面使用 `--language en`。每次生成必须指定新目录。
这些命令不连接 SSH、不启动远程作业、不运行大型积分。

默认示例从零生成两条路径，输出两圈三点函数的 **7 行 × 11 列 = 77 格**完整对应表，包含结构零。
目录中有 `path-a.json`、`path-b.json`、`table.csv`、`verification.json`、完整 Wick 审查、
中英文 HTML、SVG 和可编辑 TikZ/LaTeX。[仓库内的紧凑示例](examples/counterterms/README.md)可直接查阅和重新验证。

![含二点子发散的图](examples/counterterms/figures/L2E3G4.svg)

## 当前实现与原理

| 环节 | 源码 |
|---|---|
| 图生成与同构 | `src/phi3lab/graphs.py`、`higher_graphs.py` |
| 生成泛函与类型化 Wick 展开 | `ct_wick.py`、`ct_typed.py` |
| 路径 A：实际递归 `C=−T R̄`、Henge 核验 | `ct_bph.py`、`ct_series_a.py` |
| 路径 B：独立作用量图、低圈来源与端口双射 | `ct_action.py`、`ct_series_b.py` |
| 结构、系数、符号及动量逐格核验 | `ct_verify.py`、`ct_series_verify.py` |
| 数值实现 | `automatic_bph.py` 及 `finite`、`primitive`、`nested_bubble`、`nested_triangle`、`overlap`、`two_loop` 等模块 |
| 研究源码与推导 | [326 个 Python 数学模块](research/README.md)、[168 个 C++／头文件](research/native-precision/README.md)、[53 份 LaTeX 源文件](docs/derivations/README.md) |

路径 B 不从 A 的原图收缩结果或 `insertion_columns` 汇总插入列。
回归测试禁止 B 调用这些生成器后，再独立生成并比较。两条路径共享图结构及同构工具，
图枚举和系数推导独立。

采用六维欧氏有质量理论，作用量含 `+gφ³/3!`，传播子为 `1/(q²+m²)`，`m²>0`，
圈测度为 `d⁶q/π³`。裸顶点为 `−g`，反项顶点为作用量系数的负值；
正组合权重 `1/S` 与振幅符号、耦合幂分开记录。

二点 Taylor 是外动量**向量总次数不超过 2**，分别对应质量 `P⁰` 和动能 `P²`，
不是对 `P²` 展开到二阶。三点为零阶。保留 tadpole，一点有限线性条件尚未指定。
[L2E3G4 逐项说明](docs/L2E3G4-zh.md)展示裸图 `−1/2`、质量插入 `+1/2`、动能插入 `+1/2`
及整体减除；该图子图外动量是宿主环变量，Taylor 只作用于子图因子。
[完整约定](docs/conventions.md)。

## 高圈支持与数值状态

```bash
phi3-ct family --loops 2 --legs 2 --output runs/L2E2
phi3-ct verify --input runs/L2E2
phi3-ct one-loop --output runs/one-loop
python -m pytest -q
```

完整族支持一至四圈、二点或三点；三／四圈由用户显式运行，成本较高。
历史完整族证据的规模如下；此次导入没有重跑四圈计算：

| 族 | 裸图行 | 对应列 | 精确单元格 |
|---|---:|---:|---:|
| L1E2 | 1 | 3 | 3 |
| L1E3 | 1 | 2 | 2 |
| L2E2 | 2 | 7 | 14 |
| L2E3 | 7 | 11 | 77 |
| L3E2 | 10 | 40 | 400 |
| L3E3 | 56 | 162 | 9,072 |
| L4E2 | 64 | 385 | 24,640 |
| L4E3 | 540 | 2,553 | 1,378,620 |

这些是固定外腿标签的截肢连通 1PI 图数，不能与全部生成泛函项数混用。
缺失的外部参考图仍待逐图比对，没有预设“10 张”或“11 张”。

[162 组积分研究总账](docs/evidence/precision-ledger.json)记录：**严格解析 24 组、经验双核验收 48 组、
未解决 90 组**。经验验收不是严格求积误差证明，精确组合相等也不等于积分精度通过。
重计算保持**暂停**；本次源码导入不增加验收数。[状态与边界](docs/status.md)。

可选数值依赖：`python -m pip install ".[numerics]"`。wheel 包含完整 55 模块计算核心；
研究源码和 TeX 随 Git 仓库及源码分发包提供。机器专用作业驱动、原生生成绑定和大型样本不在公开包中。

## 验证与来源

[实际验证记录](docs/validation.md)区分本次重新运行与历史证据。
[逐文件清单](docs/source-import.json)记录导入前后 SHA-256 和必要的界面文字修改。
原有 `phi3-demo` 命令继续可用。

```bash
python -m ruff check src tests
python -m ruff format --check src tests
python -m pytest -q
python -m build
```

研究核心保留原始格式，不为样式检查批量改写；其行为由精确回归测试覆盖。
[发布范围](docs/publication.md) · [后续工作](ROADMAP.md)。

原理参考 Kennedy 的 [BPH 证明](https://arxiv.org/abs/hep-th/9612113)、Kennedy、Binoth、Rippon 的
[重整化自动化](https://arxiv.org/abs/0712.1016)，以及 Srednicki 的[教材](https://www.physics.ucsb.edu/~mark/qft.html)。
本项目不宣称新的普遍全圈收敛定理。原创软件采用 MIT 许可，不打包外部教材、论文或依赖源码。
