# 刻意排除之物及其理由

生成于 2026-10-31 20:44:30 +0800

## 一 · 按名硬排除（隐患文件）

- **`step8c_compass_v3_ranking_SOLV.csv`** —— 名为 v3 而内含 28 行（V2 队列）。读者会据文件名误取为本稿队列。若要存档，须重命名为 …_main28_… 并在同目录加一份 README 说明其为被取代之队列。
- **`permutation_test_report_SOLV.txt`** —— 内含标为「PAPER-READY SENTENCE」之句，引 n = 28、r = 0.116。本稿之值为 ρ = +0.078、p = 0.676、n = 31；该报告会主动递给读者一句本稿已排除之话。
- **`phase3_cross_dG_newTier1_SOLV.csv`** —— **V2 跨臂产物，非本稿之物。** 六行，2026-06-23，分子为 `005` / `027` / `028`（`Cross_arm_decision_2026-10-29.md` §1.2 已按内容确证）。本稿跨臂为`004` / `012` / `057` / `074` / `279` / `001` 六个分子、12 个体系、36 条轨迹。二者同以 `phase3_cross_dG_` 起头，若并置于 `05_cross_ortholog/` 之下，读者会取它为本稿的一部分。本稿无任何论断系于它。

## 二 · 第三层：生成侧模型，不入本存档

理由：本稿不作生成阶段效力主张；将生成模型与之一同存档，会把主动学习的材料读入本稿，且这些构件是下一篇论文的证据基础，系于本稿 DOI 之下会混淆它们所支持的是哪一篇的主张。

- `model_best_V2.pth` —— 解码器权重 —— 属生成侧，为下一篇（PRISM-Gen 2.0）之证据基础
- `activity_predictor.pkl` —— pIC50 代理 —— 同上
- `activity_predictor_BBS_v*.pkl` —— BBS 代理各版本 —— 同上，且其版本比较属该篇
- `train_corrected.py` —— 训练脚本 —— 同上

**一处例外：** `fragments_V2.csv` 已存档。§2.1.1 以内容与哈希识别它（10,173 碎片 / 319,428 字节 / MD5 `bc8816cc`），而一个以哈希作出的识别若其所指文件未存档，便不是读者可用的识别。

## 二之二 · 单文件例外（覆盖了整类排除）

整类排除是钝器；下列各份落在被排除之树中而仍然入档，每一条附其理由。**一处不记录的破例，与一处沉默的缺口同类。**

- **`release/PRISM_Dyn/scripts/step7_pbsa_fix_v4.in`** —— 其参数与本稿所载逐项相符（startframe 5001 / endframe 10000 / interval 25 / radiopt 1），且与 `results/step7_gbsa_solv_v3/gbsa_v3_build.json` 所记的 `same_frames_as_pbsa: (5001,10000,25)` 与 `frames_per_replica_mmpbsa: 200` 独立相符。**此为「参数与所记运行相符的输入」，而非已确立为「该次运行所读的那个文件」**——后者须定位 v3 臂的 step7 驱动方可确立，README 须照此措辞陈述。

## 三 · 按目录整类排除

lang2 下存有多代被取代的产物树。以逐文件的方式排除挡不住它们，故按**路径段**整类排除。其中 `step7` 与 `step7_gbsa_solv` 两代产物属 **28 候选队列**，而该队列已由 May 于 2026-10-27 指示排除于本稿之外；若随存档发布，读者会得到与本稿 31 候选不相容的端点能量。

- **release 树为另一次发布的组装，与本存档的边界不同** —— 1 个文件，例：`release/FRATTVAE_V2/fragments_V2.csv`

### 三之一 · 整类排除之登记（本次未产生命中者）

下列目录已按路径段登记排除。本次组装中**没有任何模式指向它们**，故上表不会列出——但其排除理由仍须随存档可见，否则「已排除」一说便无从检验。

- **`15_全文修改前备份/`** —— 全文修改前备份
- **`_abandoned_20261015/`** —— 已废弃的跨臂首次尝试（§3.6.1 已载偏离记录）
- **`_archive_prebugfix_20260621/`** —— 修缺陷之前的产物
- **`_archive_release_old/`** —— 被取代的 release 旧档
- **`_s1_authoritative_inputs/`** —— ⚠ **其名曰 authoritative，而其内容不属本稿。** 2026-10-30 实读 70 个文件：它不是 MD 输入集（我先前的刻画有误，此处更正），而是**配体参数化（S1）阶段**的记录——逐分子的 `ligand_authoritative.mol` 与 `parameterization_manifest.json`、`failure_manifest.json`，另有 `_s1_bijection_audit.csv`、`_s1_provenance_frozen.json`、`_S1_FAILED.json`。内含 22 个 `pool_mol_*` 目录（日期 7 月 19–20 日），其中 002 / 018 / 023 / 025 / 026 / 028 六个**不在**本稿的 31 分子之列，而本稿的 003 / 017 / 022 / 029 / 030 / 031 / 034 / 057 / 061 / 073 / 074 / 120 / 146 / 243 / 279 十五个**不在其中**；交集仅 16 个。故它属被取代的队列。**正因其名具误导性，此项排除必须留痕**：若随档发布，读者会据目录名取它为本稿的权威输入。
- **`aud_cross_site_series/`** —— **30 条轨迹那一代的跨臂逐帧序列（2026-10-19 18:28）。** 本稿跨臂为 36 条（AMD-2，2026-10-31），其序列在 `aud_cross_site_occupancy_series/`（36 份）。二者只差目录名的一截，若并置入档，读者无从分辨哪一份对应稿件 §3.6。
- **`step7/`** —— 2026-06-01 的 step7 产物，属 28 候选队列（May 2026-10-27 指示排除本稿）
- **`step7_gbsa_solv/`** —— 2026-06-21 的 step7 产物，同属 28 候选队列
- **`test_single/`** —— 单体系试跑
- **`tmp/`** —— 临时产物

## 三之二 · 按内容去重（只存代表 + 哈希证明）

- `02_md_inputs` / `md_solv_v3/*/r*/03_min/em.mdp`：627 个文件仅 **9 种内容**，已只存代表；全部路径与其 sha256 见 `02_md_inputs/IDENTICAL.sha256`
- `02_md_inputs` / `cross_solv_v3/*/r*/03_min/em*.mdp`：304 个文件仅 **9 种内容**，已只存代表；全部路径与其 sha256 见 `02_md_inputs/IDENTICAL.sha256`
- `02_md_inputs` / `md_solv_v3/*/r*/06_md/ndx_groups.json`：145 个文件仅 **1 种内容**，已只存代表；全部路径与其 sha256 见 `02_md_inputs/IDENTICAL.sha256`

## 四 · 按扩展名排除

轨迹与二进制：.cpt, .edr, .gz, .h5, .log, .npz, .pkl, .pth, .tar, .tng, .trr, .xtc, .zip

