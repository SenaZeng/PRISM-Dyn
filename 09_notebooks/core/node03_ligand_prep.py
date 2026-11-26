# -*- coding: utf-8 -*-
"""
V3_clean/core/node03_ligand_prep.py
===================================
Node 03 — 配体制备（master SDF → docking-ready PDBQT）+ 正向原子映射。

消费 Node 02 的权威 3D 主分子（`identity_manifests.jsonl` + `sdf/<name>.sdf`），
为每个【已构建】分子产出一个对接就绪的 PDBQT，并建立一张【显式、带校验和】的
正向原子映射（pdbqt 原子 ↔ master 原子）。

【为什么本节点最危险】——V2 的身份级联缺陷正是在此发生：
  1) obabel `-p 7.4` 在转换时【重新质子化】，覆盖 Node 01 冻结的微观态；
  2) pdbqt 丢失键级（bond order），若下游从 pdbqt 反推化学 → 身份被改写；
  3) 依赖工具"保持原子顺序"而非显式映射 → 索引错配，坐标张冠李戴。

【本节点如何逐条封死】：
  1) 用 Meeko（读 SDF、保留显式 H 与形式电荷、【绝不】重新质子化）而非 obabel-带-p。
  2) 键级 / 立体 / 形式电荷【绝不】从 PDBQT 读取 —— 它们始终归 master SDF 所有，
     对接姿势只是 master 分子上的【新坐标】（"one molecule object" 原则的空间层落实）。
  3) 正向映射【不信任工具的原子顺序】，而是用【坐标恒等】重建：SDF→PDBQT 转换不移动
     原子，故 pdbqt 原子坐标与 master 原子坐标一一对应，按坐标匹配即得权威映射。

【不做化学判断】——质子化 / 电荷 / 立体已由 Node 01/02 冻结。Node 03 只加原子类型、
定义可旋转键、写坐标，并核验身份骨架未被改动。

QC 门（失败即 QC_FAIL，记录、不静默通过）：
  - 每个 PDBQT 原子都能匹配到一个 master 原子（无多余原子 → 防"加氢/重新质子化"）；
  - master 的【全部重原子 + 全部极性 H】都被覆盖（无缺失 → 防"脱氢/脱质子"）；
  - 映射在 pdbqt 侧双射（无一对多 → 防索引错配）；
  - （警告，不阻断）Gasteiger 电荷和 ≈ 冻结净电荷。
非极性 H（接在 C 上）被 Meeko 合并进母原子 —— 这是设计，不算缺失。

输入：`results/02_master_molecules/identity_manifests.jsonl`（Node 02 v0.2.0）
输出：每分子 `pdbqt/<name>.pdbqt` + `atom_maps/<name>.json`
      + `identity_manifests_node03.jsonl` + held / failures + manifest

DRY_RUN：默认 True（fail-safe），必须显式传入。write_outputs 带一致性守卫。
"""

import re
import json
import hashlib
import subprocess
from datetime import datetime, timezone
from pathlib import Path

NODE03_VERSION = "0.2.0"

# 坐标匹配容差（Å）。SDF→PDBQT 不移动原子，故匹配应精确；1e-2 给足浮点余量。
COORD_TOL = 1e-2

# v0.2.0：统一【刚性环】。Meeko 默认会把较大环拆开做柔性大环（macrocycle），
# 并插入 G# 胶水赝原子（与真实断点碳共坐标）。这会让被拆环的分子在对接中获得
# 额外构象自由度，而其余分子的环是刚性的 —— 于是 BBS 分布出现【表征不一致】，
# 恰好污染设定 MD 阈值的那个聚类。对接只是粗筛，环柔性在 MD 中充分采样，故这里
# 对全部分子统一施加刚性环，保证 BBS 可比。写入 Methods。
RIGID_MACROCYCLES = True

# Meeko 柔性大环的胶水赝原子类型：G0/G1/...（真实断点碳被重标为 CG0/CG1，仍是真原子）。
_GLUE_TYPE_RE = re.compile(r"^G\d+$")


def _is_glue(ad_type):
    """G# 赝原子（Meeko 柔性大环胶水）；与真实原子共坐标、无独立 master 对应，
    须排除出双射并单独记录。CG# 是真实碳，不算胶水。"""
    return bool(_GLUE_TYPE_RE.match(ad_type))
# Gasteiger 电荷和与冻结净电荷的允许偏差（仅警告，不阻断）。
GASTEIGER_SUM_TOL = 0.05


def utc_now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def file_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _str_sha256(s):
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


# ────────────────────────── 后端适配 ──────────────────────────
def meeko_available():
    try:
        import meeko  # noqa: F401
        return True
    except Exception:
        return False


def _meeko_prepare(mol, rigid_macrocycles=RIGID_MACROCYCLES):
    """用 Meeko 生成 PDBQT 字符串。跨 0.5–0.7 版本适配。
    返回 (pdbqt_str, ok, err)。Meeko 保留输入的显式 H 与形式电荷，合并非极性 H，
    分配 AutoDock 原子类型与扭转树 —— 【不】重新质子化。
    rigid_macrocycles=True 时不拆环、不产生 G# 胶水赝原子（全分子统一刚性环）。"""
    from meeko import MoleculePreparation
    try:
        prep = MoleculePreparation(rigid_macrocycles=rigid_macrocycles)
    except TypeError:
        # 旧版无该 kwarg：退回默认（会拆大环），由下游 glue 处理兜底
        prep = MoleculePreparation()
    setups = prep.prepare(mol)                       # 0.5+：返回 setup 列表
    setup = setups[0] if isinstance(setups, (list, tuple)) else setups
    # 0.5–0.7：PDBQTWriterLegacy.write_string -> (str, ok, err)
    try:
        from meeko import PDBQTWriterLegacy
        out = PDBQTWriterLegacy.write_string(setup)
        if isinstance(out, tuple):
            s = out[0]
            ok = bool(out[1]) if len(out) > 1 else True
            err = out[2] if len(out) > 2 else ""
            return s, ok, err
        return out, True, ""
    except ImportError:
        # 0.3–0.4 旧 API
        s = prep.write_pdbqt_string()
        return s, True, ""


def _obabel_prepare(sdf_path, out_pdbqt):
    """【次选，会告警】obabel 回退。绝不使用 -p / -h / --addhydrogens ——
    任何改变质子化状态的开关都被禁止。仅当 Meeko 不可用且用户显式选择时才走这里。
    obabel 的原子顺序不可信，但本节点用坐标匹配重建映射，故仍安全（只是扭转树质量差）。"""
    cmd = ["obabel", str(sdf_path), "-O", str(out_pdbqt), "--partialcharge", "gasteiger"]
    banned = {"-p", "-h", "--addhydrogens", "--addpolarh"}
    if banned & set(cmd):
        raise RuntimeError("obabel 命令包含被禁止的质子化开关 —— 会破坏冻结微观态")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0 or not Path(out_pdbqt).is_file():
        raise RuntimeError(f"obabel 失败: {r.stderr[:200]}")
    return Path(out_pdbqt).read_text()


# ────────────────────────── PDBQT 解析 + 坐标映射 ──────────────────────────
def parse_pdbqt_atoms(pdbqt_str):
    """解析 PDBQT 的 ATOM/HETATM 行，返回 [(x, y, z, ad_type, partial_charge)]。
    按固定列宽解析（PDBQT 规范）。"""
    atoms = []
    for ln in pdbqt_str.splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
            toks = ln.split()
            ad = toks[-1]                        # 完整 AD 类型（CG0 / G0 / HD …），非截断 2 字符
            try:
                q = float(toks[-2])              # 部分电荷（倒数第二个 token）
            except (ValueError, IndexError):
                q = 0.0
            atoms.append((x, y, z, ad, q))
    return atoms


def _polar_h_indices(mol):
    """极性 H = 接在非碳重原子上的 H（Meeko 会保留；非极性 H 接在 C 上，被合并）。"""
    out = set()
    for a in mol.GetAtoms():
        if a.GetAtomicNum() == 1:
            nbrs = a.GetNeighbors()
            if nbrs and nbrs[0].GetAtomicNum() != 6:
                out.add(a.GetIdx())
    return out


def build_forward_map(master_mol, pdbqt_atoms, tol=COORD_TOL):
    """按坐标恒等重建 pdbqt→master 原子映射（工具无关的权威映射）。
    G# 胶水赝原子（Meeko 柔性大环）与真实原子共坐标、无独立 master 对应，
    故排除出双射并单独返回。
    返回 (fwd, unmatched, glue)：
      fwd       : {pdbqt_index(0基, 仅真实原子): master_atom_idx(0基)}
      unmatched : [(pdbqt_index, ad_type, nearest_dist)]  —— 非空即真实原子有"多余"（疑似加氢）
      glue      : [pdbqt_index, ...]  —— G# 赝原子（记录，不参与双射）
    """
    conf = master_mol.GetConformer()
    mc = [(conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z)
          for i in range(master_mol.GetNumAtoms())]
    fwd, unmatched, glue = {}, [], []
    for pi, (x, y, z, ad, q) in enumerate(pdbqt_atoms):
        if _is_glue(ad):
            glue.append(pi)
            continue
        best, bd = None, 9e18
        for mi, (mx, my, mz) in enumerate(mc):
            d = (x - mx) ** 2 + (y - my) ** 2 + (z - mz) ** 2
            if d < bd:
                bd, best = d, mi
        if bd ** 0.5 <= tol:
            fwd[pi] = best
        else:
            unmatched.append((pi, ad, round(bd ** 0.5, 4)))
    return fwd, unmatched, glue


# ────────────────────────── 单分子制备 ──────────────────────────
def prepare_one(rec, Chem, out_dir, backend="meeko", dry_run=True,
                tol=COORD_TOL, rigid_macrocycles=RIGID_MACROCYCLES):
    """对一条 Node 02 记录制备 PDBQT + 原子映射，做 QC。
    返回 (status, manifest)：status ∈ {'built', 'held', 'qc_fail'}。"""
    name = rec["name"]
    sdf_path = rec.get("sdf_path")
    frozen_q = rec.get("frozen_net_charge", rec.get("formal_charge"))
    out_dir = Path(out_dir)
    pdbqt_path = out_dir / "pdbqt" / f"{name}.pdbqt"
    map_path = out_dir / "atom_maps" / f"{name}.json"

    base = {"name": name, "parent_id": rec.get("parent_id"),
            "node02_sdf_sha256": rec.get("sdf_sha256"),
            "frozen_net_charge": frozen_q, "backend": backend}

    # 读冻结 master SDF；核验它就是 Node 02 写出的那一个（sdf_sha256 一致）
    if not sdf_path or not Path(sdf_path).is_file():
        return "held", {**base, "reason": "master_sdf_missing"}
    if rec.get("sdf_sha256") and file_sha256(sdf_path) != rec["sdf_sha256"]:
        return "qc_fail", {**base, "reason": "node02_sdf_hash_drift",
                           "error": "master SDF 与 Node 02 记录的 sha256 不符"}
    mol = Chem.MolFromMolFile(str(sdf_path), sanitize=True, removeHs=False)
    if mol is None or mol.GetNumConformers() == 0:
        return "held", {**base, "reason": "master_sdf_unreadable"}

    # 生成 PDBQT
    try:
        if backend == "meeko":
            pdbqt_str, ok, err = _meeko_prepare(mol, rigid_macrocycles=rigid_macrocycles)
            if not ok or not pdbqt_str:
                return "held", {**base, "reason": "meeko_prepare_failed",
                                "error": str(err)[:200]}
        elif backend == "obabel":
            out_dir.joinpath("pdbqt").mkdir(parents=True, exist_ok=True)
            pdbqt_str = _obabel_prepare(sdf_path, pdbqt_path)
        else:
            return "held", {**base, "reason": f"unknown_backend:{backend}"}
    except Exception as e:
        return "held", {**base, "reason": "prepare_exception",
                        "error": f"{type(e).__name__}: {str(e)[:180]}"}

    # 解析 + 坐标映射
    pa = parse_pdbqt_atoms(pdbqt_str)
    if not pa:
        return "qc_fail", {**base, "reason": "empty_pdbqt", "error": "无 ATOM 行"}
    fwd, unmatched, glue = build_forward_map(mol, pa, tol)

    heavy = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1}
    polar_h = _polar_h_indices(mol)
    covered = set(fwd.values())

    # ---- QC 门 ----
    problems = []
    if unmatched:
        problems.append(f"unmatched_pdbqt_atoms={len(unmatched)} (疑似重新质子化/加氢)")
    missing_heavy = heavy - covered
    if missing_heavy:
        problems.append(f"missing_heavy_atoms={len(missing_heavy)}")
    missing_polar_h = polar_h - covered
    if missing_polar_h:
        problems.append(f"missing_polar_H={len(missing_polar_h)} (疑似脱质子)")
    if len(set(fwd.values())) != len(fwd):
        problems.append("non_bijective_map (一 master 原子被多个 pdbqt 原子占用)")

    if problems:
        return "qc_fail", {**base, "reason": "identity_qc_failed",
                           "error": "; ".join(problems),
                           "n_pdbqt_atoms": len(pa), "n_matched": len(fwd),
                           "n_unmatched": len(unmatched)}

    # Gasteiger 电荷和（仅警告；Vina 忽略部分电荷，权威电荷是冻结净电荷）
    gsum = round(sum(a[4] for a in pa), 3)
    charge_warn = None
    try:
        if frozen_q is not None and abs(gsum - float(frozen_q)) > GASTEIGER_SUM_TOL:
            charge_warn = f"gasteiger_sum={gsum} vs frozen_net_charge={frozen_q}"
    except (TypeError, ValueError):
        pass

    # 显式正向 + 反向映射（0 基索引；PDBQT ATOM serial = pdbqt_index + 1）
    fwd_map = {str(pi): mi for pi, mi in sorted(fwd.items())}
    rev_map = {str(mi): pi for pi, mi in fwd.items()}
    map_obj = {
        "name": name,
        "note": ("pdbqt_index(0-based, file order; PDBQT serial = index+1) -> "
                 "master_sdf_atom_idx(0-based). Built by coordinate identity, "
                 "NOT by trusting tool atom order."),
        "coord_tol_angstrom": tol,
        "n_pdbqt_atoms": len(pa),
        "n_master_atoms": mol.GetNumAtoms(),
        "n_master_heavy": len(heavy),
        "n_polar_h": len(polar_h),
        "n_nonpolar_h_merged": mol.GetNumAtoms() - len(heavy) - len(polar_h),
        "n_glue_pseudo_atoms": len(glue),
        "glue_pdbqt_indices": glue,
        "rigid_macrocycles": rigid_macrocycles,
        "pdbqt_to_master": fwd_map,
        "master_to_pdbqt": rev_map,
    }
    map_obj["map_checksum"] = _str_sha256(
        json.dumps(map_obj["pdbqt_to_master"], sort_keys=True))

    if not dry_run:
        pdbqt_path.parent.mkdir(parents=True, exist_ok=True)
        map_path.parent.mkdir(parents=True, exist_ok=True)
        pdbqt_path.write_text(pdbqt_str)
        map_path.write_text(json.dumps(map_obj, ensure_ascii=False, indent=2))

    man = {
        **base,
        "pdbqt_path": str(pdbqt_path),
        "pdbqt_sha256": _str_sha256(pdbqt_str),
        "atom_map_path": str(map_path),
        "map_checksum": map_obj["map_checksum"],
        "n_pdbqt_atoms": len(pa),
        "n_master_heavy": len(heavy),
        "n_polar_h": len(polar_h),
        "n_nonpolar_h_merged": map_obj["n_nonpolar_h_merged"],
        "n_glue_pseudo_atoms": len(glue),
        "rigid_macrocycles": rigid_macrocycles,
        "gasteiger_charge_sum": gsum,
        "charge_warning": charge_warn,
        # 溯源：把 Node 02 的身份钉进来，链路可追
        "canonical_isomeric_smiles": rec.get("canonical_isomeric_smiles"),
        "identity_hash": rec.get("identity_hash"),
        "written": (not dry_run),
    }
    return "built", man


# ────────────────────────── 批量编排 ──────────────────────────
def run_node03(node02_manifests, Chem, out_dir, backend="meeko",
               limit=None, dry_run=True, tol=COORD_TOL,
               rigid_macrocycles=RIGID_MACROCYCLES):
    """遍历 Node 02 的【已构建】记录，逐一制备 PDBQT + 映射。
    返回 (built_manifests, held, failures, stats)。"""
    out_dir = Path(out_dir)
    print(f"[node03] dry_run={dry_run} backend={backend} "
          f"rigid_macrocycles={rigid_macrocycles} -> "
          f"{'不写 PDBQT(预览)' if dry_run else '将写 PDBQT 到 ' + str(out_dir / 'pdbqt')}",
          flush=True)
    if backend == "meeko" and not meeko_available():
        raise ImportError(
            "Node 03 需要 Meeko（读 SDF、保留质子化、不重新加氢）。\n"
            "  安装：pip install meeko\n"
            "  若确需回退，显式传 backend='obabel'（次选，会告警；绝不用 -p/-h）。")

    built, held, failures = [], [], []
    n = 0
    for rec in node02_manifests:
        if limit is not None and n >= limit:
            break
        n += 1
        status, man = prepare_one(rec, Chem, out_dir, backend=backend,
                                  dry_run=dry_run, tol=tol,
                                  rigid_macrocycles=rigid_macrocycles)
        if status == "built":
            built.append(man)
        elif status == "held":
            held.append(man)
        else:
            failures.append(man)

    # 母体覆盖审计：某母体的全部子代都未制备成 PDBQT → 整体掉出 docking
    from collections import Counter, defaultdict
    built_names = {m["name"] for m in built}
    by_parent = defaultdict(list)
    for rec in node02_manifests[:n]:
        by_parent[rec.get("parent_id")].append(rec["name"])
    orphan_parents = [p for p, kids in by_parent.items()
                      if p is not None and not any(k in built_names for k in kids)]

    stats = {
        "n_processed": n,
        "n_built": len(built),
        "n_held": len(held),
        "n_qc_fail": len(failures),
        "backend": backend,
        "coord_tol_angstrom": tol,
        "rigid_macrocycles": rigid_macrocycles,
        "n_with_glue_pseudo_atoms": sum(1 for m in built if m.get("n_glue_pseudo_atoms")),
        "held_reasons": Counter(h.get("reason", "unknown") for h in held),
        "qc_fail_reasons": Counter(f.get("reason", "unknown") for f in failures),
        "n_charge_warnings": sum(1 for m in built if m.get("charge_warning")),
        "orphaned_parents": {"n_orphaned_parents": len(orphan_parents),
                             "orphan_parent_ids": orphan_parents[:50]},
        "failures": failures[:50],
        "held_examples": held[:20],
    }
    # Counter -> dict（JSON 友好）
    stats["held_reasons"] = dict(stats["held_reasons"])
    stats["qc_fail_reasons"] = dict(stats["qc_fail_reasons"])
    return built, held, failures, stats


def _read_node02_provenance(node02_dir):
    """钉住上游 Node 02 manifest（版本 + hash），链路完整可追。"""
    p = Path(node02_dir) / "Node02_manifest.json"
    if not p.is_file():
        return None
    try:
        m = json.loads(p.read_text(encoding="utf-8"))
        return {"node02_version": m.get("version"),
                "manifest_sha256": file_sha256(p),
                "node02_input_sha256": m.get("input_sha256")}
    except Exception:
        return None


def _tool_versions(backend):
    v = {}
    try:
        import rdkit
        v["rdkit"] = rdkit.__version__
    except Exception:
        v["rdkit"] = None
    if backend == "meeko":
        try:
            import meeko
            v["meeko"] = getattr(meeko, "__version__", "unknown")
        except Exception:
            v["meeko"] = None
    return v


def write_outputs(built, held, failures, stats, out_dir, node02_dir,
                  backend="meeko", dry_run=True):
    """写 identity_manifests_node03.jsonl + held.csv + qc_failures.json + manifest。"""
    out_dir = Path(out_dir)
    idman = out_dir / "identity_manifests_node03.jsonl"
    held_path = out_dir / "node03_held.csv"
    fail_path = out_dir / "node03_qc_failures.json"
    mani_path = out_dir / "Node03_manifest.json"

    manifest = {
        "module": "node03_ligand_prep",
        "version": NODE03_VERSION,
        "created_at_utc": utc_now(),
        "dry_run": dry_run,
        "backend": backend,
        "tool_versions": _tool_versions(backend),
        "upstream_node02": _read_node02_provenance(node02_dir),
        "invariant": {
            "internal_gate": (
                "master SDF 的【全部重原子 + 全部极性 H】按坐标恒等（tol=%g Å）双射到 "
                "PDBQT 原子；无多余 PDBQT 原子（防重新质子化），无缺失极性 H（防脱质子）。"
                "非极性 H 按设计被合并。" % COORD_TOL),
            "chemistry_ownership": (
                "键级 / 立体 / 形式电荷【不】从 PDBQT 读取 —— 始终归 master SDF 所有，"
                "对接姿势只是 master 分子上的新坐标（经 atom_maps/<name>.json 回映）。"
                "故 pdbqt 键级丢失（V2 已知缺陷）无法在此污染身份。"),
        },
        "policy": {
            "ligand_source": "Node 02 master SDF (frozen); read by sdf_sha256, not rebuilt",
            "protonation": "NEVER re-decided; Meeko preserves input H & formal charge",
            "banned": "obabel -p / -h / --addhydrogens (any protonation-altering switch)",
            "atom_map": "coordinate identity, tool-order-independent, checksummed per ligand",
            "rigid_macrocycles": (
                "True — 全分子统一刚性环，禁用 Meeko 柔性大环拆分，避免仅 28/1765 个"
                "分子获得额外构象自由度而污染 BBS 分布；环柔性在 MD 阶段充分采样"),
            "glue_pseudo_atoms": (
                "G# 赝原子（若启用柔性大环才会出现）与真实原子共坐标、不参与身份双射，"
                "单独记录于 atom_map；当前刚性环设定下应为 0"),
            "nonpolar_h": "merged into parent heavy atom (AutoDock convention)",
            "charge_note": "Vina ignores partial charges; frozen_net_charge is authoritative; "
                           "gasteiger_charge_sum recorded as a soft consistency check only",
        },
        "stats": stats,
    }

    if not dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)
        with open(idman, "w", encoding="utf-8") as f:
            for m in built:
                f.write(json.dumps(m, ensure_ascii=False) + "\n")
        import csv
        if held:
            keys = sorted({k for h in held for k in h.keys()})
            with open(held_path, "w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=keys)
                w.writeheader()
                for h in held:
                    w.writerow(h)
        fail_path.write_text(json.dumps(failures, indent=2, ensure_ascii=False),
                             encoding="utf-8")
        mani_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False),
                             encoding="utf-8")
    return manifest, {"identity_manifests_node03": idman, "held": held_path,
                      "failures": fail_path, "manifest": mani_path}
