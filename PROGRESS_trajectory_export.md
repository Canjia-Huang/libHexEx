# HexEx：六面体网格元素「轨迹」导出 —— 进度记录

> 用途：会话中断/重启后据此续接。
> **状态：已完成。** 四个测试网格（fiveFine / fiveadvanced / cylinderadvanced / sphere）全部通过结构不变量检查（0 failure）；轨迹以两个 geogram 网格输出，并在导出后回读逐项校验。

---

## 0. 结论速览（验收口径）

`demo/export_trajectory` 的检查分两级：

- **FAIL（必须为 0）**：结构不变量 —— 所有 piece 的体积之和 == 输入 tet 体积之和（体积守恒，相对误差 < 1e-9）；每个 refined tet 都有 piece 与输入 tet；所有 hex 元素索引在范围内；refined 网格无裂缝（`n_crack_edges == 0`，`n_vertex_weld_conflicts == 0`）；hex edge 折线上的 refined edge 顺序一致；写出 geogram 后回读，属性数量与数值完全一致。
- **warn（可由 extractor 的洞/折叠解释）**：没有轨迹的 hex 顶点/边/面/胞的数量。

四个网格的结果：

| | tets | pieces | groups / hex cells | refined v / e / f / tets | 体积相对误差 | 无轨迹的 hex 元素 | FAIL |
|---|---|---|---|---|---|---|---|
| **fiveFine** | 13492 | 89248 | 4418 / 4416 | 157021 / 407403 / 428881 / 836256 | 2.7e-15 | 0 v, 2 e, 1 f, 1 c | **0** |
| **fiveadvanced** | 13492 | 23846 | 96 / 96 | 34914 / 66348 / 79127 / 151796 | 4.4e-15 | **无** | **0** |
| **cylinderadvanced** | 18461 | 37843 | 140 / 140 | 55832 / 111604 / 131471 / 256970 | 6.5e-15 | **无** | **0** |
| **sphere** | 76 | 1475 | 224 / 216 | 2979 / 8307 / 8279 / 15392 | 3.8e-16 | 0 v, 31 e, 27 f, 8 c | **0** |

关键点：**体积守恒精确到 1e-15**，说明「tet ∩ 单位立方体」这一分解是真正的剖分。遗漏的少量元素集中在 sphere（极粗且严重折叠的参数化，76 个 tet 映射到 6×6×6 的网格，且 216 个格胞中约 59% 的候选立方体与 tet 的交为空）以及 fiveFine 的洞附近，均由诊断计数如实报告。

---

## 1. 任务目标（goal）

从 integer-grid parameterization 提取六面体网格的过程中，记录 hex vertices / hex edges / hex facets / hex cells 在输入四面体网格中的**轨迹**，最终以 geogram 网格输出（对标 libQEx fork 的 `Output (poly)`）：

1. 每个 hex 顶点 ← 哪些 refined 顶点（即哪些 tet / 参数位置）；
2. 每个 hex 边 ← 一条穿过 tet 网格的折线（polyline）；
3. 每个 hex 面 ← 一片 refined 三角面片（patch）；
4. 每个 hex 胞 ← 一组 refined 四面体（其并集就是该格胞在输入网格中的体积）。

---

## 2. 核心设计决定

### 2.1 refined 网格 = 「tet × 单位立方体」的剖分

对每个输入 tet，其在**自身参数空间**中的像是四面体 `P`。整数网格的单位立方体把参数空间剖分，因此

```
piece(tet, cube) = tet ∩ cube
```

就是轨迹的原子。参数空间里 cube 就是 `[I,I+1]×[J,J+1]×[K,K+1]`（`I,J,K` 为整数），**只用局部整数坐标**即可枚举——不需要任何全局坐标。

每个 piece 用 10 个半空间（4 个 tet 面 + 6 个 cube 面）裁剪得到凸多面体；体积用散度定理；边界多边形按**顶点集合**做规范三角化（从最小顶点扇形展开），再**从 piece 质心向每个三角形做锥**得到 refined 四面体。这样：

- refined facet（piece 边界三角形）在相邻 piece 之间共享（规范三角化保证一致）；
- refined tet 的体积和 == piece 体积（实测 1e-15 级）；
- 全部 piece 的体积和 == 输入 tet 体积（体积守恒）。

### 2.2 不用全局整数网格标架（重要）

最初（错误）方案：用 extractor 的 transition function 在 tet 邻接图上 BFS，给每个 tet 一个全局→局部的 `GridIsomorphism`，从而给每个 piece 一个全局格胞坐标 `(I,J,K)`，再用它对应 hex 元素。

**该方案在 `cylinderadvanced` 上失败**（140 个格胞只有 63 个能对上，2862 个 transition 冲突）。原因是**本质性的**：integer-grid map 在奇异边附近标架会旋转，transition 沿闭合回路复合后不是恒等（实测回路给出 90° 旋转 + 平移），因此**全局一致标架根本不存在**。这不是数值问题，也不是 transition 拟合错误（逐面验证：所有面的 transition 都把 3 个顶点精确映射过去，误差 < 1e-6）。

因此最终方案：

- 轨迹几何完全局部（参数空间整数平面），不依赖全局标架；
- **piece 分组**：同一个单位立方体的所有 piece 构成一个 group。piece 之间**只通过「落在 tet 面上」的 facet 连接**；对每条连接，用该面的 transition 把当前 piece 的局部 cube 映射到邻接 tet 的局部坐标，若与邻接 piece 的 cube 完全相等才连接。transition 逐面精确，所以即使全局有 holonomy，分组也正确。
- 沿 BFS 同时累计 `pieceToCanonical = pieceToCanonical[parent] * transition^{-1}`，把每个 piece 的 8 个角点映射到 group 的规范标架。
- hex 元素的识别：
  - hex 顶点：refined 网格中的格点顶点（3 个整数平面）按**3D 位置**匹配 extractor 的中间 hex 顶点（容差 1e-6），再取等价类 → 最终 hex 顶点；
  - group 的 8 个角点 = 各 piece 暴露的 cube 角点（3 整数平面）经规范标架换算后的 corner index → (final hex vertex)；
  - 8 个角都拿到后，用顶点集合在最终 hex 网格中查 cell；
  - hex face：cube 面（局部 axis/side）的 4 个角 → 4 个 hex 顶点 → 查 face；
  - hex edge：落在两条整数平面交线上的 refined edge，用 piece 局部的两条平面确定 cube 棱 → 2 个角 → 查 edge。
  - **tet 面与整数平面重合**的情况（域边界上的 hex 面）：facet 的支撑平面是 tet 面，但多边形所有顶点共享一个 cube 平面，因此同样识别为 hex face（否则 `fiveadvanced` 会丢 200 个边界 hex 面）。

### 2.3 输出形式

两个 geogram 网格：

- `<prefix>_refined.geogram`：vertices + edges + triangular facets + tetrahedra，全部带属性；
- `<prefix>_hex.geogram`：提取出的 hex 网格（hexes + quads + edges），带 `n_refined` 等属性。

refined 网格的 facet/edge 存储**只包含 piece 边界**（piece 之间和域边界），这正是 hex 面/边轨迹所在；piece 内部的三角形/边由四面体隐含。geogram 的 `.geogram` 格式同时保存顶点/边/面/胞及其属性（已实测 round-trip）。

---

## 3. 改动清单

| 文件 | 改动 |
|---|---|
| `src/HexEx/TrajectoryLayout.hh` | **新增**：`HexTrajectoryLayout`（refined mesh、pieces、groups、逐 hex 元素的对应关系与轨迹、诊断计数、`describe()`）；`TrajectoryMeshView`（与具体 mesh 类型解耦的连接关系） |
| `src/HexEx/TrajectoryLayout.cc` | **新增**：refined 网格构建（裁剪、规范三角化、质心锥化、按位置焊接）、piece 分组（transition 逐面传播）、hex 元素识别、折线排序、诊断 |
| `src/HexEx/HexExtractor.hh` | 新增 `getHexMeshWithTrajectory()`（模板，做 getHexMesh 的全部工作后构建 layout）与 `buildTrajectoryLayout()` 声明；`#include "TrajectoryLayout.hh"` |
| `CMakeLists.txt` | 源文件加入 `TrajectoryLayout.cc`；`HEXEX_WITH_GEOGRAM=AUTO/ON/OFF`、`GEOGRAM_ROOT`；找到 geogram 时加入 `interfaces/geogram` 与 `demo/export_trajectory` |
| `interfaces/geogram/hexex_geogram.hh/.cc` | **新增**：`GeogramBridge`（`toRefinedMesh` / `toHexMesh` / `save` / `describe`），全部信息写入 geogram attributes |
| `interfaces/geogram/CMakeLists.txt` | **新增**：`HexExGeogram` 库 |
| `cmake/FindGeogram.cmake` | **新增**（从 libQEx 移植）：支持已安装的 geogram 或源码树 + 构建目录 |
| `demo/export_trajectory/main.cc` + CMakeLists | **新增**：导出工具 + 两级验收 + geogram 回读校验 |
| `apps/hexex/main.cc` + CMakeLists | **新增**：`hexex` 命令行程序（CLI11），对齐 libQEx 的 `qex`：`--in/--out/--out-refined/--out-hex/--out-trajectory/--check/--refined-surface-only/--simple-merge-position/-v`；默认静默，只有请求 geogram 输出或 `--check` 时才构建轨迹 |
| `src/HexEx/Verbose.hh/.cc` | **新增**：`setVerbose()/verbose()`；`HEXEX_DEBUG_ONLY` 由编译期改为运行期开关（`apps/hexex` 默认关闭输出） |
| `third_party/CLI11` | **新增** git submodule（v2.7.2，与 libQEx 一致）；缺失时 CMake 用 FetchContent 拉同版本 |
| `tests/trajectory_test.cc` | **新增** gtest：立方体（含 transition）100% 覆盖；`testdata/sphere.hexex` 结构不变量 |
| `README.md` | 新增 “Trajectories of the extracted elements” 一节 |

---

## 4. 构建与运行

```bash
cmake -S . -B build            # 自动发现 ~/geogram
cmake --build build -j8

# 命令行程序（推荐）
./build/Build/bin/hexex --in tests/testdata/fiveadvanced.hexex \
    --out-trajectory /tmp/fiveadvanced --check
#  -> /tmp/fiveadvanced_refined.geogram
#  -> /tmp/fiveadvanced_hex.geogram
# 可只输出其一：--out-refined / --out-hex；--refined-surface-only 省略四面体；
# --out <file.ovm> 写 OpenVolumeMesh；不加 geogram 参数时只提取、不建轨迹。

# 完整校验 demo（写出 + 回读逐项校验）
./build/Build/bin/export_trajectory tests/testdata/fiveadvanced.hexex /tmp/fiveadvanced

# 单元测试
cd build/Build && ./unittests --gtest_filter='Trajectory.*'
# 或（ctest，工作目录 build/Build）
```

可选开关：`-DHEXEX_WITH_GEOGRAM=OFF` 关闭 geogram（核心库和原有测试不受影响）；`-DGEOGRAM_ROOT=/path/to/geogram` 指定源码树。

`hexex` 默认关闭库的进度输出（`HexEx::setVerbose(false)`），`-v` 打开。环境变量 **`HEXEX_TRAJECTORY_DIAGNOSTICS=1`** 打开轨迹诊断输出（拒绝的 piece 连接、缺角 group、非法 hex cell 角点等）。

---

## 5. 关键结论 / 踩过的坑

1. **不能用全局标架**：integer-grid map 有奇异边 → transition 沿回路不闭合（cylinderadvanced 2862 处）。凡是「给每个 tet 一个全局格子坐标」的做法都会在大约一半的网格上失效。正确做法是**逐面**使用 transition 做 piece 分组。
2. **transition 的方向**：`getTransitionFunction(hfh)` 把 `incident_cell(hfh)` 的参数映射到 `opposite(hfh)` 的参数。把当前 piece 的 cube 映射到邻接 tet 时用 **transition**（不是它的逆）；把规范标架从 parent 传到 child 时用 `canonicalToParent * transition^{-1}`。方向写反会表现为「相邻 piece 的 cube 完全对不上」（早期版本 group 碎裂成 322 个，用此诊断定位）。
3. **裁剪会留下未被任何面引用的顶点**（`poly.vertices` 里多于面实际使用的点）。必须先用面循环收集 `used` 顶点，否则：质心算在凸多面体外面 → 锥化体积放大数十倍；refined 网格里多出孤立顶点。修复后 `n_pieces_with_bad_tetrahedralization == 0`。
4. **体积要在同一个空间比较**：`poly.volume()` 是参数空间体积；乘上 `|det(paramTo3D 线性部分)|` 才是 3D 体积。输入 tet 体积也用 3D 坐标算。两者只在同一空间才守恒。
5. **tet 面可能与整数平面重合**（域边界）：这类 facet 的支撑平面是 tet 面，但它同时也是 hex 面。判定方法是「多边形所有顶点共享同一个 cube 平面（在该 piece 的局部标架中）」。不做这一步会丢掉大量边界 hex 面（fiveadvanced 少 200 个面）。
6. **facet 焊接用三角形关键字**：`facet_piece_b` 通过相同三角形顶点三元组识别邻居。曾被**共线退化三角形**错误地共享，导致假的 piece 连接（进而让 piece 分组错误拒绝连接）；由于 cube 面/ tet 面不重叠，面积非零的三角形不会跨两个 cube 共享。
7. **vertex 焊接容差**：refined 顶点按 3D 位置焊接（1e-7 网格单位）。同一几何点由不同 tet 的仿射映射算出，误差 ~1e-13；不同点若近于此容差视为数值退化。piece 内部两个不同角点被焊接在一起计为 `weld conflict` 并丢弃该 piece（体积损失 < 1e-14）。
8. **中间 hex cell 的顶点数可能 > 8**：extractor 的中间 cell 可能带额外顶点。必须像 `addCellsToHexMesh()` 那样，只取 `localCellUV ∈ {0,1}^3` 的 8 个角点，否则会误判为非法（sphere 早期 7 个 cell 因此丢失）。
9. **refined 边的「在网格线上」判定**：两端点共享的两个整数平面（同一 piece 的局部标架）→ 网格线；三个平面 → 格点。跨 piece 合并平面集合是错的（不同标架）。因此每 piece 单独存 `(vertex → 局部平面)`。
10. **折线要沿连接关系走**，不能按坐标排序：折叠参数化下同一条 hex 边可能被访问两次（坐标非单调、出现重复格点）。按邻接走链；分支时退化为排序并计入 `n_hex_edges_with_branched_polyline`。
11. **退化 tet**（参数体积为 0，实测 2/22/5/84 个）：裁剪会得到非零的伪体积，必须跳过（`tv > 1e-9`）。否则体积守恒破坏（曾导致 17% 误差）。
12. 几何是 3D 的还是参数的：refined 网格输出 3D 坐标（`getInverseParametrizationMatrix` 每 tet 一次）；参数空间只用于裁剪和体积。

---

## 6. 剩余可选项（非阻塞）

1. sphere 的 8 个无轨迹格胞：来自极粗/折叠参数化下 group 缺角（round 到 2/4/6 个角）。可在「7/8 角」时用 hex 网格反查补全（当前只实现到完整 8 角才查 cell），但 2/4/6 角的 group 无法唯一补全。
2. 若需要更小的输出：可为 `toRefinedMesh` 增加「不写 refined 四面体」的选项（只保留顶点/边/面），但那样 hex 胞轨迹就只剩属性、没有几何。
3. 可以把 refined 网格的 piece 多边形本身（而非其三角化）导出为 connector/polyhedron，但 geogram 只支持 tet/hex/prism/pyramid/connector（connector 固定 4 顶点 3 面），因此选择四面体化。
4. `edge_hex_edge_step` 对分支（折叠）折线只是列表顺序，不是几何顺序；已单独计数。
