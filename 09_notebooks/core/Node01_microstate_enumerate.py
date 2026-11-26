"""
core/Node01_microstate_enumerate.py
V3_clean · Node 01 — authoritative microstate assignment.

WHAT THIS DOES
--------------
Takes the Node 00 identity pool (all_candidates_v3_start.csv) and rules on ONE
authoritative protonation microstate per molecule at physiological pH, using
Dimorphite-DL. This is the deliberate, recorded charge decision the whole V3
rebuild depends on — the charge is NOT inherited from any stage-specific number
and NOT re-guessed later (and, pointedly, NOT set by OpenBabel).

CONFIRMED SCIENTIFIC POLICY (encoded below)
-------------------------------------------
1. Backend      : Dimorphite-DL (free, reproducible; not OpenBabel).
2. pH / count   : pH 7.4, keep the SINGLE dominant microstate. When >1 state is
                  plausible at 7.4 the choice is genuinely ambiguous, so we take
                  a deterministic placeholder and FLAG the molecule for review —
                  we never silently pick a contested state.
3. Protonation  : protonation ONLY. Tautomers are NOT auto-canonicalized and
   only          undefined stereocentres are NOT auto-assigned; such molecules
                  are FLAGGED for manual review, not silently resolved.

INVARIANT (fail-closed)
-----------------------
Protonation may add/remove H and change net charge, but must NOT change the
heavy-atom skeleton. Every microstate is checked against the input on:
  (a) heavy-atom formula  (always computable), and
  (b) HetAtomProtomer hash (RDKit rdMolHash, when available) —
and a microstate that fails is rejected, never carried as an "assigned" state.

BONUS — definitive QM verdict
-----------------------------
Node 00 could only tentatively flag the 24 explicitly-neutralized QM entries.
With the authoritative charge in hand, Node 01 emits `qm_final_status`:
QM is `valid` ONLY if the molecule has DFT AND is neutral here AND PySCF computed
it neutral; otherwise `recompute`. This also catches molecules PySCF treated as
neutral that are in fact charged at pH 7.4.

CONTRACT: Input -> Processing -> Output -> QC Gate -> Manifest.
Nothing is written unless dry_run=False and the QC gate passes.
Orchestrated from the notebook (import from core/; do not inline algorithms).
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from collections import Counter

NODE01_VERSION = "0.9.2"

DEFAULT_PH = 7.4
DEFAULT_PKA_PRECISION = 1.0
MAX_VARIANTS = 128

# columns expected from the Node 00 output that Node 01 reads for the QM verdict
QM_CHARGENOTE_COL = "PySCF_ChargeNote"
QM_CALCMETHOD_COL = "Calc_Method"


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def file_sha256(path):
    import hashlib
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# --------------------------------------------------------------------------- #
# Protonation backend adapter (Dimorphite-DL).
# If your installed Dimorphite-DL exposes a different API, edit ONLY this fn.
# --------------------------------------------------------------------------- #
def default_protonate(smiles, min_ph=DEFAULT_PH, max_ph=DEFAULT_PH,
                      precision=DEFAULT_PKA_PRECISION, max_variants=MAX_VARIANTS):
    """
    Return a list of protonated SMILES for `smiles` at [min_ph, max_ph].

    Dimorphite-DL's public API differs across versions, so this tries the known
    interfaces in order. If none match, the raised error lists the module's
    attributes so the exact adapter can be pinned in one place.
    """
    import importlib
    dm = importlib.import_module("dimorphite_dl")

    # API 1 — class:  DimorphiteDL(...).protonate(smi) -> list[str]
    if hasattr(dm, "DimorphiteDL"):
        eng = dm.DimorphiteDL(min_ph=min_ph, max_ph=max_ph, pka_precision=precision,
                              max_variants=max_variants, label_states=False)
        return [str(s) for s in eng.protonate(smiles)]

    # API 2 — function: protonate_smiles(smi, ...) -> list[str]
    if hasattr(dm, "protonate_smiles"):
        for kwargs in (dict(ph_min=min_ph, ph_max=max_ph, precision=precision),
                       dict(min_ph=min_ph, max_ph=max_ph, pka_precision=precision),
                       {}):
            try:
                return [str(s) for s in dm.protonate_smiles(smiles, **kwargs)]
            except TypeError:
                continue

    # API 3 — Durrant-lab classic: run_with_mol_list([mol], ...) -> list[RDKit Mol]
    if hasattr(dm, "run_with_mol_list"):
        from rdkit import Chem
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return []
        for kwargs in (dict(min_ph=min_ph, max_ph=max_ph, pka_precision=precision,
                            max_variants=max_variants, silent=True),
                       dict(min_ph=min_ph, max_ph=max_ph, pka_precision=precision)):
            try:
                outs = dm.run_with_mol_list([m], **kwargs)
                return [Chem.MolToSmiles(x) for x in outs if x is not None]
            except TypeError:
                continue

    raise ImportError(
        "Unrecognized dimorphite_dl API — none of DimorphiteDL / protonate_smiles / "
        "run_with_mol_list matched. Module exposes: "
        + ", ".join(a for a in dir(dm) if not a.startswith("_")))


def selftest_backend(protonate_fn=None, ph=DEFAULT_PH, precision=DEFAULT_PKA_PRECISION):
    """
    Probe the protonation backend on one known molecule BEFORE the batch.
    Acetic acid (CC(=O)O) should deprotonate to CC(=O)[O-] at pH 7.4.
    Returns a small dict; raises if the backend is unusable.
    """
    if protonate_fn is None:
        protonate_fn = default_protonate
    probe = "CC(=O)O"
    states = protonate_fn(probe, ph, ph, precision)  # raises if API mismatched
    states = [str(s) for s in states]
    ok = any("[O-]" in s for s in states)
    return {"probe": probe, "states": states, "deprotonated_as_expected": ok}


# --------------------------------------------------------------------------- #
# Structure helpers
# --------------------------------------------------------------------------- #
def heavy_formula(mol):
    """Heavy-atom composition string (H-independent). Invariant under protonation."""
    c = Counter(a.GetSymbol() for a in mol.GetAtoms() if a.GetSymbol() != "H")
    return "".join(f"{el}{c[el]}" for el in sorted(c))


def protomer_hash(mol, rdMolHash):
    """Protonation-invariant skeleton hash; None if unavailable."""
    if rdMolHash is None:
        return None
    try:
        return rdMolHash.MolHash(mol, rdMolHash.HashFunction.HetAtomProtomer)
    except Exception:
        return None


def count_undefined_stereo(mol, Chem):
    """Number of unassigned stereocentres (protonation policy = flag, not fix)."""
    try:
        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True,
                                            useLegacyImplementation=False)
        return sum(1 for _, lab in centers if lab in ("?", None))
    except Exception:
        try:
            centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
            return sum(1 for _, lab in centers if lab == "?")
        except Exception:
            return -1  # unknown


# --------------------------------------------------------------------------- #
# v0.9.0 — stereo perception correction
#
# DEFECT FIXED: stereo was assessed on the PROTONATED microstate while charge is
# assessed on the NEUTRALIZED scaffold (Rule 1). The inconsistency meant that the
# ammonium nitrogen created by our OWN protonation was counted as an "undefined
# stereocentre", inflating k, driving spurious 2^k expansion in finalize, and
# pushing molecules over the k>=4 hold threshold.
#
# Two atom classes are excluded from the undefined-stereocentre count:
#   (a) protonated amine N  ([NH+], [NH2+], [NH3+]) — configuration interconverts
#       by rapid proton exchange in water, so it is NOT configurationally stable
#       on the docking/MD timescale. NOTE: a quaternary ammonium with no N-H is
#       configurationally stable and is deliberately KEPT.
#   (b) bridgehead atoms of BRIDGED bicyclics (rings sharing >= 3 atoms, e.g.
#       bicyclo[1.1.1]pentane) — the bridges are equivalent, so there is no real
#       R/S, and an enumerated tag is geometrically unsatisfiable. FUSED ring
#       systems (rings sharing exactly 2 atoms, e.g. decalin) are genuine
#       stereocentres and are KEPT.
# --------------------------------------------------------------------------- #
def bridgehead_atoms(mol, Chem=None):
    """Atom indices that are bridgeheads of a BRIDGED bicyclic (not fused)."""
    out = set()
    if mol is None:
        return out
    try:
        ri = mol.GetRingInfo()
        rings = [set(r) for r in ri.AtomRings()]
        for i in range(len(rings)):
            for j in range(i + 1, len(rings)):
                shared = rings[i] & rings[j]
                if len(shared) >= 3:          # bridged (fused shares exactly 2)
                    for idx in shared:
                        a = mol.GetAtomWithIdx(idx)
                        ring_bonds = sum(1 for b in a.GetBonds()
                                         if ri.NumBondRings(b.GetIdx()) > 0)
                        if ring_bonds >= 3:
                            out.add(idx)
    except Exception:
        return set()
    return out


def labile_protonated_nitrogens(mol, Chem=None):
    """Atom indices of N+ bearing at least one H (configuration labile via
    proton exchange). Quaternary ammonium (no H) is excluded from this set."""
    out = set()
    if mol is None:
        return out
    try:
        for a in mol.GetAtoms():
            if (a.GetSymbol() == "N" and a.GetFormalCharge() == 1
                    and a.GetTotalNumHs() >= 1):
                out.add(a.GetIdx())
    except Exception:
        return set()
    return out


def count_undefined_stereo_corrected(mol, Chem):
    """
    Return (k_corrected, k_raw, n_excluded_N, n_excluded_bridgehead).

    k_raw       = what v0.8.0 counted (kept for audit/back-compat)
    k_corrected = genuine undefined stereocentres, after removing the two
                  artifact classes above. This is what gates enumeration.
    """
    if mol is None:
        return -1, -1, 0, 0
    try:
        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True,
                                            useLegacyImplementation=False)
    except Exception:
        try:
            centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        except Exception:
            return -1, -1, 0, 0

    unassigned = [idx for idx, lab in centers if lab in ("?", None)]
    k_raw = len(unassigned)
    bh = bridgehead_atoms(mol, Chem)
    npn = labile_protonated_nitrogens(mol, Chem)

    n_bh = sum(1 for i in unassigned if i in bh)
    n_np = sum(1 for i in unassigned if (i in npn and i not in bh))
    k_corr = k_raw - n_bh - n_np
    return max(k_corr, 0), k_raw, n_np, n_bh


def excluded_stereo_atoms(mol, Chem):
    """Union of atom indices whose stereochemistry must NOT be enumerated.
    Shared with Node01_finalize so gate and operation use the SAME atom set."""
    return bridgehead_atoms(mol, Chem) | labile_protonated_nitrogens(mol, Chem)


def _canon(mol, Chem):
    return Chem.MolToSmiles(mol, isomericSmiles=True)


# --------------------------------------------------------------------------- #
# Physiological (pH 7.4) formal-charge model  [SCIENTIFIC RULE TABLE — review]
# --------------------------------------------------------------------------- #
# Used ONLY to pick, among Dimorphite's enumerated states, the one whose NET
# charge matches the physiologically expected charge. It corrects the observed
# over-deprotonation of lactam/amide N-H, phenols, and pyridines by the
# lexicographic placeholder. Edit these SMARTS to change the policy — nothing
# is silent; the chosen rule is recorded in the manifest.
#
# Deprotonated at 7.4 (contribute -1 each): pKa well below 7.4
PHYSIO_ACID_SMARTS = {
    "carboxylic_acid":   "[CX3](=O)[OX2H1]",
    "sulfonic_acid":     "[SX4](=O)(=O)[OX2H1]",
    "phosphonic/phosphate": "[PX4](=O)([OX2H1])",
    "tetrazole":         "[$(c1[nH]nnn1),$(c1n[nH]nn1)]",
    "acyl_sulfonamide":  "[SX4](=O)(=O)[NX3H1][CX3]=O",
}
# Protonated at 7.4 (contribute +1 each): pKa well above 7.4
PHYSIO_BASE_SMARTS = {
    "guanidine":         "[NX3][CX3](=[NX2])[NX3]",
    "amidine":           "[NX3][CX3]=[NX2;!$([NX2]a)]",
    # aliphatic amine: match the N ATOM ONLY (sp3-C attachment as a recursive
    # constraint) so a ring N with two sp3 neighbours is counted ONCE, not twice.
    # v0.9.3 (2026-10-14): exclude nitrogens whose basicity is removed by an
    # adjacent electron-withdrawing substituent. A nitramine nitrogen delocalises
    # its lone pair into the nitro group and a cyanamide nitrogen into the nitrile,
    # so neither is protonated at pH 7.4 — cyanamide and N-cyanoguanidine are weak
    # ACIDS. Before this, both satisfied every clause of the aliphatic-amine
    # pattern and were assigned +1: 32 N-cyano and 2 N-nitro microstates in the
    # 1,790-row V3 table. N-nitroso is added pre-emptively (0 rows at present).
    "aliphatic_amine":   "[NX3;H0,H1,H2;!$(N=*);!$(N-a);!$(NC=O);!$(NS(=O)=O)"
                         ";!$(N[CX3]=[NX2]);!$(N[CX3](=[NX2])[NX3])"
                         ";!$(N[N+](=O)[O-]);!$(N[NX3](=O)=O)"
                         ";!$(N[CX2]#[NX1])"
                         ";!$(N[NX2]=[OX1])"
                         ";$([NX3][CX4])]",
}
# Everything else (phenol pKa~10, lactam/amide N-H pKa~9-11, pyridine pKa~5,
# primary amide) is treated as NEUTRAL at pH 7.4.


def neutralize_scaffold(mol, Chem):
    """恢复中性参考骨架：把简单质子得失的可电离基团（羧酸根/铵/酚氧等）恢复中性。
    v0.7.0 关键修复：expected_charge 必须从中性骨架算，否则输入若已含 [O-]/[NH+]，
    官能团 SMARTS（找 -OH/-NH）匹配不到，已电离基团被漏计（bug: 402/485/1218 得 0 而非 -1）。
    只中性化单质子得失，不动重原子骨架。"""
    patt = Chem.MolFromSmarts("[+1!H0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    if patt is None:
        return mol
    m = Chem.RWMol(mol)
    changed = False
    for match in m.GetSubstructMatches(patt):
        a = m.GetAtomWithIdx(match[0])
        ch = a.GetFormalCharge()
        try:
            a.SetFormalCharge(0)
            a.SetNumExplicitHs(a.GetTotalNumHs() - ch)
            changed = True
        except Exception:
            pass
    if not changed:
        return mol
    out = m.GetMol()
    try:
        Chem.SanitizeMol(out)
        return out
    except Exception:
        return mol


def construct_authoritative_microstate(mol, Chem):
    """v0.7.0：从中性骨架按化学规则【构造】pH7.4 权威微观态（规则=真值源）。
    不从 Dimorphite 候选里挑，而是直接应用规则：
      酸性基团（羧酸/磺酸/膦酸/酰基磺酰胺）→ 去质子 -1
      碱性基团（胍/脒/脂肪胺）→ 质子化 +1
      其它（酰胺/芳香N-H/酚/内酰胺/N-O）→ 保持中性
    返回 (canonical_smiles, net_charge)。sanitize 失败返回 (None, error)。
    注意：四氮唑不在此自动处理（其 SMARTS 需针对具体结构确认，留人工）。"""
    from rdkit.Chem import MolFromSmarts
    ref = neutralize_scaffold(mol, Chem)

    # 歧义检测（v0.7.0）：规则只对"确定"官能团有信心。若分子含规则覆盖【外】的
    # 可电离基团，构造层【拒绝】并退回人工——不强行给可能错的答案。
    # 覆盖外的模糊基团：四氮唑（1H/2H/C-取代难自动定）、烯醇、异羟肟酸等。
    AMBIGUOUS_SMARTS = {
        # v0.8.0: only a FREE-N-H tetrazole is tautomer-ambiguous (1H/2H). An
        # N-substituted tetrazole (no [nH]) is unambiguously neutral, so it must
        # NOT block construction — that over-triggering forced mol_1002/1225 to
        # 'conflict' when their true state is +1 (secondary amine).
        "tetrazole_NH": "[$(c1[nH]nnn1),$(c1n[nH]nn1),$(c1nn[nH]n1),$(c1nnn[nH]1),"
                        "$([nH]1nnnc1)]",
        "enol": "[CX3]([OX2H1,OX1-])=[CX3]",
        "hydroxamic": "[CX3](=O)[NX3][OX2H1,OX1-]",
    }
    for name, s in AMBIGUOUS_SMARTS.items():
        p = MolFromSmarts(s)
        if p is not None and ref.HasSubstructMatch(p):
            return None, f"ambiguous_group:{name}"

    rw = Chem.RWMol(ref)
    claimed = set()
    # 酸：去质子（该基团里带 H 的 O/N 设 -1）
    for smarts in PHYSIO_ACID_SMARTS.values():
        p = MolFromSmarts(smarts)
        if p is None:
            continue
        for match in ref.GetSubstructMatches(p):
            if set(match) & claimed:
                continue
            for idx in match:
                a = rw.GetAtomWithIdx(idx)
                if a.GetSymbol() in ("O", "N") and a.GetTotalNumHs() > 0:
                    a.SetFormalCharge(-1)
                    a.SetNumExplicitHs(max(0, a.GetTotalNumHs() - 1))
                    claimed.update(match)
                    break
    # 碱：质子化（N 设 +1）
    # v0.9.1: adjacent amine clusters take only ONE proton at pH 7.4 (see
    # adjacent_amine_correction). The constructed state must match the corrected
    # expected_charge, otherwise construction would be rejected as a mismatch.
    amine_cluster = _amine_cluster_map(ref, Chem)
    protonated_clusters = set()
    for smarts in PHYSIO_BASE_SMARTS.values():
        p = MolFromSmarts(smarts)
        if p is None:
            continue
        for match in ref.GetSubstructMatches(p):
            if set(match) & claimed:
                continue
            for idx in match:
                a = rw.GetAtomWithIdx(idx)
                if a.GetSymbol() != "N":
                    continue
                cid = amine_cluster.get(idx)
                if cid is not None and cid in protonated_clusters:
                    claimed.update(match)   # cluster already carries its proton
                    break
                a.SetFormalCharge(+1)
                a.SetNumExplicitHs(a.GetTotalNumHs() + 1)
                if cid is not None:
                    protonated_clusters.add(cid)
                claimed.update(match)
                break
    out = rw.GetMol()
    try:
        Chem.SanitizeMol(out)
        return Chem.MolToSmiles(out), Chem.GetFormalCharge(out)
    except Exception as e:
        return None, f"sanitize_failed:{e}"


def expected_charge_at_pH(mol, Chem):
    """
    Physiologically expected net formal charge at pH 7.4 from the rule table.

    Groups are counted with ATOM-OVERLAP de-duplication: a guanidine also matches
    the amidine SMARTS on shared atoms, so without this it would count twice.
    Bases are matched most-specific-first (guanidine, amidine, then amine); any
    later match that shares an atom with an already-counted group is skipped.
    """
    from rdkit.Chem import MolFromSmarts
    mol = neutralize_scaffold(mol, Chem)   # v0.7.0: 先中性化，修已电离基团漏计
    # v0.9.2 PERMANENT-CHARGE BASELINE.
    # Previously q started at 0 and only titratable groups were counted, so the
    # rule was structurally blind to charges that are NOT pH-dependent:
    #   quaternary ammonium  [N+](C)(C)(C)C   — no titratable proton, permanently +1
    #   hypervalent species  e.g. [Cl+2]
    # neutralize_scaffold only touches atoms bearing hydrogen (its SMARTS requires
    # !H0 on the cation branch), so such charges survive neutralisation and the
    # residual formal charge of the neutral scaffold IS the permanent charge.
    # Failure mode this fixes: pool_mol_707 / 1304 (quaternary N+, rule said 0 vs
    # true +1) and pool_mol_1589 ([Cl+2], rule said 0 vs true +2) -> false conflicts.
    q = Chem.GetFormalCharge(mol)
    claimed = set()

    def _count(smarts_dict, sign):
        nonlocal q
        for smarts in smarts_dict.values():
            patt = MolFromSmarts(smarts)
            if patt is None:
                continue
            for match in mol.GetSubstructMatches(patt):
                atoms = set(match)
                if atoms & claimed:      # overlaps an already-counted group
                    continue
                claimed.update(atoms)
                q += sign

    _count(PHYSIO_BASE_SMARTS, +1)       # order: guanidine -> amidine -> amine
    _count(PHYSIO_ACID_SMARTS, -1)
    q += adjacent_amine_correction(mol, Chem)   # v0.9.1
    return q


# --------------------------------------------------------------------------- #
# v0.9.1 — adjacent-amine dominant protonation state
#
# The per-group rule protonates each basic nitrogen independently. That is wrong
# for amines close enough to depress each other's successive pKa below 7.4:
# the SECOND protonation does not occur at physiological pH, so the cluster
# carries +1, not +1 per nitrogen.
#
# The correction is SPACING-DEPENDENT, because pKa2 rises with separation:
#   N-C-N     geminal (1,1)              pKa2 strongly depressed   -> cluster = +1
#   N-C-C-N   1,2-diamine / piperazine   pKa2 ~5.3 (pip) / ~6.9 (en) < 7.4 -> +1
#   N-C-C-C-N 1,3-diamine                pKa2 ~8.6 > 7.4  -> BOTH protonated,
#                                        NO correction (diprotonation is correct)
# So the cluster cutoff is a shortest path of <= 3 bonds between the nitrogens.
# --------------------------------------------------------------------------- #
ADJACENT_AMINE_MAX_BONDS = 3     # N-C-N (2) and N-C-C-N (3) cluster; N-C-C-C-N (4) does not


def _basic_amine_nitrogens(mol, Chem):
    """Indices of nitrogens the base rule would protonate (aliphatic amines only).
    Guanidine/amidine are excluded: their charge is delocalised, not an
    independent site, and they are not subject to this pairwise depression."""
    from rdkit.Chem import MolFromSmarts
    patt = MolFromSmarts(PHYSIO_BASE_SMARTS.get("aliphatic_amine", ""))
    if patt is None:
        return []
    idx = set()
    for match in mol.GetSubstructMatches(patt):
        for a in match:
            if mol.GetAtomWithIdx(a).GetSymbol() == "N":
                idx.add(a)
    return sorted(idx)


def _amine_clusters(mol, Chem, max_bonds=ADJACENT_AMINE_MAX_BONDS):
    """Group basic amine nitrogens that depress each other's pKa.
    Two nitrogens join the same cluster when their shortest path is
    <= max_bonds bonds. Returns a list of sets of atom indices."""
    nitro = _basic_amine_nitrogens(mol, Chem)
    if len(nitro) < 2:
        return [{n} for n in nitro]
    parent = {n: n for n in nitro}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for i, a in enumerate(nitro):
        for b in nitro[i + 1:]:
            try:
                path = Chem.GetShortestPath(mol, a, b)
            except Exception:
                continue
            if path and (len(path) - 1) <= max_bonds:   # len(path)-1 = bond count
                union(a, b)

    groups = {}
    for n in nitro:
        groups.setdefault(find(n), set()).add(n)
    return list(groups.values())


def _amine_cluster_map(mol, Chem, max_bonds=ADJACENT_AMINE_MAX_BONDS):
    """{atom_idx: cluster_id} for basic amine nitrogens (construction uses this
    to give each cluster exactly one proton)."""
    return {idx: cid
            for cid, cluster in enumerate(_amine_clusters(mol, Chem, max_bonds))
            for idx in cluster}


def adjacent_amine_correction(mol, Chem, max_bonds=ADJACENT_AMINE_MAX_BONDS):
    """Charge to SUBTRACT because adjacent amines are not both protonated at pH 7.4.
    A cluster of k nitrogens contributes +1 (not +k) -> correction -(k-1).
    Returns a value <= 0."""
    return -sum(len(c) - 1 for c in _amine_clusters(mol, Chem, max_bonds)
                if len(c) > 1)


def pick_by_expected_charge(valid_states, expected_q):
    """
    Among (canon_smiles, mol, charge) candidates, choose the one whose net charge
    equals expected_q. Returns (picked_tuple, resolved_by, candidate_charges).
    resolved_by: 'rule' (unique charge match) | 'rule_positional' (match but >1
    state at that charge) | 'conflict' (no state matches expected_q).
    """
    charges = sorted({c for _, _, c in valid_states})
    at_q = [t for t in valid_states if t[2] == expected_q]
    if len(at_q) == 1:
        return at_q[0], "rule", charges
    if len(at_q) > 1:
        at_q.sort(key=lambda t: t[0])
        return at_q[0], "rule_positional", charges
    # no candidate matches the expected charge -> genuine rule/tool conflict
    return None, "conflict", charges


DEFAULT_PER_MOL_TIMEOUT_S = 20   # a single molecule may not stall the whole batch


class _MicrostateTimeout(Exception):
    pass


def _alarm_handler(signum, frame):
    raise _MicrostateTimeout()


def _protonate_guarded(protonate_fn, smiles, min_ph, max_ph, precision, timeout_s):
    """Call the backend with a per-molecule wall-clock timeout (Unix main thread)."""
    import signal
    if timeout_s and timeout_s > 0 and hasattr(signal, "SIGALRM"):
        old = signal.signal(signal.SIGALRM, _alarm_handler)
        signal.setitimer(signal.ITIMER_REAL, float(timeout_s))
        try:
            return protonate_fn(smiles, min_ph, max_ph, precision)
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0.0)
            signal.signal(signal.SIGALRM, old)
    return protonate_fn(smiles, min_ph, max_ph, precision)


def assign_microstate(input_smiles, protonate_fn, Chem, rdMolHash,
                      ph=DEFAULT_PH, precision=DEFAULT_PKA_PRECISION,
                      timeout_s=DEFAULT_PER_MOL_TIMEOUT_S, dominant_rule=False):
    """Rule on one authoritative microstate for a single molecule."""
    rec = {
        "microstate_smiles": None, "net_charge": None, "protonation_status": None,
        "n_states_at_ph": 0, "n_states_rejected": 0, "skeleton_ok": None,
        "heavy_formula_in": None, "stereo_undefined_centers": None,
        "stereo_undefined_raw": None, "stereo_excluded_protonated_N": 0,
        "stereo_excluded_bridgehead": 0,
        "expected_charge": None, "candidate_charges": None, "resolved_by": None,
        "dimorphite_conflict": False,
        "possible_overprotonation": False, "review_tier": "B_manual",
        "tetrazole_found": False, "tetrazole_nh_present": None, "tetrazole_charge": 0,
        "needs_review": True, "review_reason": "",
    }
    in_mol = Chem.MolFromSmiles(str(input_smiles))
    if in_mol is None:
        rec.update(protonation_status="failed_parse", skeleton_ok=False,
                   review_reason="input SMILES did not parse")
        return rec

    in_charge = Chem.GetFormalCharge(in_mol)
    in_formula = heavy_formula(in_mol)
    in_prot = protomer_hash(in_mol, rdMolHash)
    rec["heavy_formula_in"] = in_formula

    # v0.8.0 auditable tetrazole check (records why the tetrazole is/ isn't acidic)
    t_found, t_nh, t_charge, t_type = tetrazole_proton_check(in_mol, Chem)
    rec["tetrazole_found"] = t_found
    rec["tetrazole_nh_present"] = t_nh
    rec["tetrazole_charge"] = t_charge

    try:
        states = _protonate_guarded(protonate_fn, str(input_smiles),
                                    ph, ph, precision, timeout_s)
    except _MicrostateTimeout:
        rec.update(microstate_smiles=_canon(in_mol, Chem), net_charge=in_charge,
                   protonation_status="failed_protonation", n_states_at_ph=0,
                   skeleton_ok=True,
                   stereo_undefined_centers=count_undefined_stereo(in_mol, Chem),
                   review_reason=f"protonation timed out after {timeout_s}s "
                                 f"(combinatorial blow-up on this scaffold)")
        return rec
    except Exception as e:
        # protonation failed -> keep input as neutral fallback, flag for review
        rec.update(microstate_smiles=_canon(in_mol, Chem), net_charge=in_charge,
                   protonation_status="failed_protonation", n_states_at_ph=0,
                   skeleton_ok=True, stereo_undefined_centers=count_undefined_stereo(in_mol, Chem),
                   review_reason=f"protonation backend error: {e}")
        return rec

    # keep only states whose heavy skeleton matches the input (reject mangling)
    valid = []
    rejected = 0
    for s in states:
        m = Chem.MolFromSmiles(str(s))
        if m is None:
            rejected += 1        # invalid protonation state (e.g. hypervalent N)
            continue
        formula_ok = (heavy_formula(m) == in_formula)
        prot_ok = True if (in_prot is None) else (protomer_hash(m, rdMolHash) == in_prot)
        if formula_ok and prot_ok:
            valid.append((_canon(m, Chem), m, Chem.GetFormalCharge(m)))
        else:
            rejected += 1        # skeleton would have been altered

    n = len(valid)
    rec["n_states_at_ph"] = n
    rec["n_states_rejected"] = rejected

    if n == 0:
        rec.update(microstate_smiles=_canon(in_mol, Chem), net_charge=in_charge,
                   protonation_status="failed_protonation", skeleton_ok=True,
                   stereo_undefined_centers=count_undefined_stereo(in_mol, Chem),
                   review_reason="no valid protonation state matched input skeleton")
        return rec

    valid.sort(key=lambda t: t[0])  # deterministic baseline order

    # v0.9.1 CRITICAL FIX — the rule must verify EVERY Dimorphite result.
    # Before v0.9.1, `n == 1` short-circuited straight to acceptance and the rule
    # was never consulted, so for those molecules Dimorphite silently became the
    # source of truth — the exact inversion this rebuild exists to prevent
    # (summary doc S2: "rule as source of truth, Dimorphite as cross-validation").
    # Failure mode it allowed: pool_mol_1096, three O-methyl hydroxamic acids
    # deprotonated to -3 while the rule says 0 (hydroxamic is inside the
    # construction guard precisely because the rule does not cover it, but
    # neither rule nor guard ran because n == 1).
    # Now: when dominant_rule is on, n == 1 goes through the same verification
    # path as n > 1. pick_by_expected_charge handles a single-candidate list.
    if n == 1 and not dominant_rule:
        micro_smi, micro_mol, charge = valid[0]
        status = "unchanged" if micro_smi == _canon(in_mol, Chem) else "assigned"
        reasons = []
    elif dominant_rule:
        # resolve among candidates by physiological expected charge
        exp_q = expected_charge_at_pH(in_mol, Chem)
        rec["expected_charge"] = exp_q
        picked, how, cand_charges = pick_by_expected_charge(valid, exp_q)
        rec["candidate_charges"] = ",".join(str(c) for c in cand_charges)
        if picked is not None:
            micro_smi, micro_mol, charge = picked
            # ---------------------------------------------------------------- #
            # v0.9.2 SITE-LEVEL VERIFICATION
            # pick_by_expected_charge verifies only the TOTAL charge. When several
            # Dimorphite candidates share the expected total, "rule_positional"
            # picks one deterministically by position — which can select a state
            # with the protons on chemically WRONG ATOMS while the total is right.
            # Observed in v0.9.1: all 19 held parents were rule_positional, e.g.
            # pool_mol_027 put a proton on a beta-lactam AMIDE nitrogen and kept an
            # adjacent diammonium, at the correct total of +2. Such site errors are
            # invisible in net charge and would propagate silently into MD.
            # Fix: when the choice is positional, the rule-constructed state is the
            # authoritative answer (the ambiguity guard already blocks construction
            # wherever the rule is not confident, so tetrazole/enol/hydroxamic
            # cases still fall back to the positional pick).
            # ---------------------------------------------------------------- #
            if how == "rule_positional":
                con_smi, con_chg = construct_authoritative_microstate(in_mol, Chem)
                if con_smi is not None and con_chg == exp_q:
                    con_mol = Chem.MolFromSmiles(con_smi)
                    if con_mol is not None and heavy_formula(con_mol) == in_formula:
                        micro_smi = _canon(con_mol, Chem)
                        micro_mol, charge = con_mol, con_chg
                        how = "rule_constructed_site"
            # v0.9.1: keep "unchanged" semantics — the rule confirming that the
            # input was already the correct state is not a reprotonation.
            status = ("unchanged" if micro_smi == _canon(in_mol, Chem)
                      else "assigned")
            rec["resolved_by"] = how
            reasons = ([] if how in ("rule", "rule_constructed_site") else
                       [f"charge {exp_q} fixed by rule; {n} states share it, "
                        f"position picked deterministically"])
        elif in_charge == exp_q:
            # Dimorphite offered no state at the expected charge, but the INPUT
            # already IS that state (e.g. it force-deprotonated a lactam N-H that
            # should stay neutral). Use the input as the authoritative microstate.
            micro_smi, micro_mol, charge = _canon(in_mol, Chem), in_mol, in_charge
            status = "assigned"
            rec["resolved_by"] = "rule_input"
            reasons = []
        else:
            # v0.7.0: Dimorphite 与输入都没给出 expected 电荷 ->
            # 规则作为真值源，直接【构造】权威微观态（不从工具错误候选里挑）。
            con_smi, con_chg = construct_authoritative_microstate(in_mol, Chem)
            if con_smi is not None and con_chg == exp_q:
                con_mol = Chem.MolFromSmiles(con_smi)
                # 构造态必须与输入重原子骨架一致（不改分子）
                if con_mol is not None and heavy_formula(con_mol) == in_formula:
                    micro_smi, micro_mol, charge = _canon(con_mol, Chem), con_mol, con_chg
                    status = "assigned"
                    rec["resolved_by"] = "rule_constructed"
                    rec["dimorphite_conflict"] = True   # 记录：Dimorphite 不一致
                    reasons = []
                else:
                    micro_smi, micro_mol, charge = valid[0]
                    status = "ambiguous"
                    rec["resolved_by"] = "conflict_construct_skeleton_mismatch"
                    reasons = [f"rule expects {exp_q}; constructed state failed "
                               f"skeleton check; needs manual review"]
            else:
                # 构造也无法达到 expected（真化学歧义：四氮唑/烯醇等）-> 人工
                micro_smi, micro_mol, charge = valid[0]
                status = "ambiguous"
                rec["resolved_by"] = "conflict"
                rec["dimorphite_conflict"] = True
                reasons = [f"rule expects net charge {exp_q}; candidates "
                           f"{{{rec['candidate_charges']}}}, input ({in_charge}), "
                           f"and rule-construction all fail to reach it -> "
                           f"genuine chemical ambiguity (manual review)"]
    else:
        micro_smi, micro_mol, charge = valid[0]
        rec["candidate_charges"] = ",".join(str(c) for c in sorted({c for _, _, c in valid}))
        status = "ambiguous"
        reasons = [f"{n} plausible states at pH {ph} (placeholder = lexicographically first)"]

    # v0.9.0: corrected perception — exclude protonation-created ammonium N and
    # bridged-bicyclic bridgeheads (see count_undefined_stereo_corrected).
    stereo_undef, stereo_raw, n_exN, n_exBH = count_undefined_stereo_corrected(
        micro_mol, Chem)
    rec["stereo_undefined_raw"] = stereo_raw
    rec["stereo_excluded_protonated_N"] = n_exN
    rec["stereo_excluded_bridgehead"] = n_exBH
    if stereo_undef and stereo_undef > 0:
        reasons.append(f"{stereo_undef} undefined stereocentre(s)")
    if (n_exN + n_exBH) > 0:
        reasons.append(f"{n_exN + n_exBH} phantom centre(s) excluded "
                       f"(protonated N={n_exN}, bridgehead={n_exBH})")

    # Rule 3 (v0.7.1): high-charge / over-protonation sanity check
    overprot, exclude, sanity_reasons = sanity_flags(micro_mol, charge, Chem)
    reasons.extend(sanity_reasons)

    # A/B/C triage. B is driven by GENUINE review reasons only — a mere
    # "position picked deterministically" note (rule_positional, correct charge)
    # is informational and stays in Class A.
    b_manual = bool(
        (stereo_undef and stereo_undef > 0)
        or (charge is not None and abs(charge) >= 3)
        or overprot
        or (status == "ambiguous")          # conflict / tetrazole / enol / genuine ambiguity
    )
    review_tier = ("C_exclude" if exclude
                   else "B_manual" if b_manual
                   else "A_resolved")
    needs_review = (review_tier != "A_resolved")

    rec.update(
        microstate_smiles=micro_smi, net_charge=charge, protonation_status=status,
        skeleton_ok=True, stereo_undefined_centers=stereo_undef,
        possible_overprotonation=overprot, review_tier=review_tier,
        needs_review=needs_review, review_reason="; ".join(reasons),
    )
    return rec


# --------------------------------------------------------------------------- #
# Rule 3 (v0.7.1) — high-charge / over-protonation sanity check
# The per-group rule protonates each amine independently; adjacent basic sites
# depress each other's successive pKa, so it can OVER-protonate (mirror of the
# earlier over-deprotonation). We do NOT model pKa depression — we flag the
# suspects and route them to review, keeping the pipeline interpretable.
# --------------------------------------------------------------------------- #
OVERPROT_SMARTS = {
    "hydrazinium":     "[NX4+]~[NX4+]",     # protonated N-N: non-physical at pH 7.4
    # v0.9.1 LABEL FIX: the previous names were off by one position.
    #   [NX4+]C[NX4+]  is GEMINAL (1,1) — both N on the same carbon
    #   [NX4+]CC[NX4+] is 1,2 (ethylenediamine / piperazine spacing)
    # The chemistry decisions were unaffected; only the labels were wrong.
    "1,1-geminal_diammonium": "[NX4+]C[NX4+]",
    "1,2-diammonium":         "[NX4+]CC[NX4+]",
}


def sanity_flags(micro_mol, net_charge, Chem):
    """Return (possible_overprotonation, exclude, reasons) for the chosen state."""
    from rdkit.Chem import MolFromSmarts
    reasons = []
    overprot = False
    exclude = False
    if micro_mol is not None:
        hz = MolFromSmarts(OVERPROT_SMARTS["hydrazinium"])
        if hz is not None and micro_mol.HasSubstructMatch(hz):
            overprot = True
            exclude = True          # hydrazinium is a clear anomaly -> Class C
            reasons.append("hydrazinium [N+]-[N+] (non-physical at pH 7.4)")
        for nm in ("1,1-geminal_diammonium", "1,2-diammonium"):
            p = MolFromSmarts(OVERPROT_SMARTS[nm])
            if p is not None and micro_mol.HasSubstructMatch(p):
                overprot = True
                reasons.append(f"adjacent {nm} (successive pKa depressed; "
                               f"likely over-protonated)")
    if net_charge is not None:
        if abs(net_charge) >= 5:
            exclude = True          # extreme charge -> Class C
            reasons.append(f"|net_charge|={abs(net_charge)} >= 5 (extreme)")
        elif abs(net_charge) >= 3:
            reasons.append(f"|net_charge|={abs(net_charge)} >= 3 (charge suspect)")
    return overprot, exclude, reasons


def tetrazole_proton_check(mol, Chem):
    """
    v0.8.0 — deterministic, auditable tetrazole charge decision.
    NOT 'saw a tetrazole -> -1'. The question is whether it carries a free N-H.
      free N-H tetrazole (unsubstituted / C5-linked) -> acidic, contributes -1
      N-substituted tetrazole (no N-H)               -> neutral, contributes  0
    Returns (found, nh_present, charge_contribution, ttype).
    """
    from rdkit.Chem import MolFromSmarts
    if mol is None:
        return False, None, 0, "none"
    ring_any = MolFromSmarts("[$(c1nnnn1),$(c1nnn[nH]1),$(c1nn[nH]n1),"
                             "$(c1n[nH]nn1),$(c1[nH]nnn1),$(n1nnnc1),$(n1nnnn1)]")
    if ring_any is None or not mol.HasSubstructMatch(ring_any):
        return False, None, 0, "none"
    nh = MolFromSmarts("[$(c1[nH]nnn1),$(c1n[nH]nn1),$(c1nn[nH]n1),"
                       "$(c1nnn[nH]1),$([nH]1nnnc1)]")
    nh_present = bool(nh is not None and mol.HasSubstructMatch(nh))
    if nh_present:
        return True, True, -1, "free_NH_acidic"
    return True, False, 0, "N_substituted_neutral"


def qm_final_verdict(net_charge, charge_note, calc_method):
    """Definitive QM validity given the authoritative charge."""
    has_qm = (calc_method is not None and str(calc_method) not in ("nan", ""))
    if not has_qm:
        return "none"
    neutral_here = (net_charge == 0)
    computed_neutral = (str(charge_note) == "neutral")
    return "valid" if (neutral_here and computed_neutral) else "recompute"


# --------------------------------------------------------------------------- #
# Batch driver
# --------------------------------------------------------------------------- #
def enumerate_microstates(df, pd, Chem, protonate_fn=None, rdMolHash=None,
                          ph=DEFAULT_PH, precision=DEFAULT_PKA_PRECISION,
                          smiles_col="canonical_smiles",
                          timeout_s=DEFAULT_PER_MOL_TIMEOUT_S, log_every=200,
                          quiet_rdkit=True, dominant_rule=False):
    """Assign one authoritative microstate per row; returns (out_df, stats)."""
    if protonate_fn is None:
        protonate_fn = default_protonate
    if smiles_col not in df.columns:
        raise KeyError(f"input must contain '{smiles_col}' (Node 00 output).")

    # Dimorphite emits chemically-invalid states on poly-aza scaffolds; RDKit
    # rejects them (correctly) but floods stderr. Mute during the batch — the
    # count of rejected states is captured per molecule instead.
    if quiet_rdkit:
        try:
            from rdkit import RDLogger
            RDLogger.DisableLog("rdApp.*")
        except Exception:
            pass

    import time
    t0 = time.time()
    n = len(df)
    recs = []
    for k, smi in enumerate(df[smiles_col].tolist(), 1):
        recs.append(assign_microstate(smi, protonate_fn, Chem, rdMolHash,
                                      ph, precision, timeout_s=timeout_s,
                                      dominant_rule=dominant_rule))
        if log_every and (k % log_every == 0 or k == n):
            done = [r["protonation_status"] for r in recs]
            n_ok = sum(1 for s in done if s in ("assigned", "unchanged", "ambiguous"))
            n_to = sum(1 for r in recs if "timed out" in (r["review_reason"] or ""))
            print(f"  [{k}/{n}] elapsed={time.time()-t0:5.0f}s  "
                  f"assigned={n_ok}  timeouts={n_to}", flush=True)

    add = pd.DataFrame(recs, index=df.index)
    out = pd.concat([df.reset_index(drop=True), add.reset_index(drop=True)], axis=1)

    # definitive QM verdict (uses authoritative net_charge + carried QM columns)
    note_col = QM_CHARGENOTE_COL if QM_CHARGENOTE_COL in out.columns else None
    calc_col = QM_CALCMETHOD_COL if QM_CALCMETHOD_COL in out.columns else None
    if calc_col:
        out["qm_final_status"] = [
            qm_final_verdict(nc,
                             out.at[i, note_col] if note_col else None,
                             out.at[i, calc_col])
            for i, nc in zip(out.index, out["net_charge"])
        ]

    def cnt(col):
        return out[col].value_counts(dropna=False).to_dict() if col in out.columns else {}

    ps = cnt("protonation_status")
    n_success = int(ps.get("assigned", 0) + ps.get("unchanged", 0) + ps.get("ambiguous", 0))
    n_failed = int(ps.get("failed_protonation", 0) + ps.get("failed_parse", 0))
    n_timeouts = int(out["review_reason"].fillna("").str.contains("timed out").sum())
    stats = {
        "n_total": len(out),
        "ph": ph, "pka_precision": precision,
        "protonation_status_counts": ps,
        "n_microstates_assigned": n_success,       # real protonation produced a state
        "n_failed": n_failed,                       # backend/parse failure (fallback used)
        "n_timeouts": n_timeouts,                   # subset of failed: combinatorial blow-up
        "failed_fraction": round(n_failed / max(len(out), 1), 4),
        "n_no_valid_protonation_state": int((out["n_states_at_ph"] == 0).sum()),
        "total_invalid_states_rejected": int(out.get("n_states_rejected",
                                                      pd.Series([0])).fillna(0).sum()),
        "n_needs_review": int(out["needs_review"].sum()),
        "review_tier_counts": (out["review_tier"].value_counts(dropna=False).to_dict()
                               if "review_tier" in out.columns else {}),
        "n_possible_overprotonation": (int(out["possible_overprotonation"].sum())
                                       if "possible_overprotonation" in out.columns else 0),
        "tetrazole_summary": ({
            "found": int(out["tetrazole_found"].sum()),
            "free_NH_acidic_-1": int(((out["tetrazole_found"]) & (out["tetrazole_charge"] == -1)).sum()),
            "N_substituted_neutral_0": int(((out["tetrazole_found"]) & (out["tetrazole_charge"] == 0)).sum()),
        } if "tetrazole_found" in out.columns else {}),
        "n_charge_nonzero": int((out["net_charge"].fillna(0) != 0).sum()),
        "net_charge_distribution": {str(k): int(v) for k, v in
                                    out["net_charge"].value_counts(dropna=False).items()},
        "n_stereo_undefined": int((out["stereo_undefined_centers"].fillna(0) > 0).sum()),
        "stereo_correction_v090": ({
            "n_raw_undefined": int((out["stereo_undefined_raw"].fillna(0) > 0).sum()),
            "n_corrected_undefined": int((out["stereo_undefined_centers"].fillna(0) > 0).sum()),
            "phantom_protonated_N_excluded": int(out["stereo_excluded_protonated_N"].fillna(0).sum()),
            "phantom_bridgehead_excluded": int(out["stereo_excluded_bridgehead"].fillna(0).sum()),
            "n_molecules_with_phantom_centres": int(
                ((out["stereo_excluded_protonated_N"].fillna(0)
                  + out["stereo_excluded_bridgehead"].fillna(0)) > 0).sum()),
        } if "stereo_undefined_raw" in out.columns else {}),
        "resolved_by_counts": (out["resolved_by"].value_counts(dropna=False).to_dict()
                               if "resolved_by" in out.columns else {}),
        "qm_final_status_counts": cnt("qm_final_status"),
    }
    return out, stats


def qc_gate(out_df, stats, max_failed_fraction=0.05):
    """Fail-closed checks before writing."""
    problems = []
    # 1) the node must actually have DONE something — reject a run where the
    #    backend produced no microstate for anything (e.g. an API mismatch that
    #    fell back to input for every molecule).
    if stats.get("n_microstates_assigned", 0) == 0:
        problems.append("no microstate was actually assigned — backend likely "
                        "failed for all rows (check review_reason / backend self-test)")
    # 2) even a partial systemic failure should not pass silently
    if stats.get("failed_fraction", 0.0) > max_failed_fraction:
        problems.append(f"failed_fraction {stats['failed_fraction']:.2%} exceeds "
                        f"{max_failed_fraction:.0%} — likely a backend/config problem")
    # 3) a skeleton-failing state must never wear an 'assigned'/'unchanged' badge
    bad = out_df[(out_df["protonation_status"].isin(["assigned", "unchanged"]))
                 & (out_df["skeleton_ok"] != True)]  # noqa: E712
    if len(bad):
        problems.append(f"{len(bad)} assigned rows with skeleton_ok != True")
    # 4) every row must carry a microstate SMILES (flagged rows keep a placeholder)
    if out_df["microstate_smiles"].isna().any():
        problems.append("rows with null microstate_smiles")
    # 5) no docking column may have leaked in from upstream
    DOCK = ["Broad_Spectrum_Score", "E_SARS_CoV_2", "E_SARS_CoV_1", "E_MERS_CoV",
            "Docking_Rank", "Is_Final_Top", "Selection_Mode"]
    leaked = [c for c in DOCK if c in out_df.columns]
    if leaked:
        problems.append(f"docking column leak: {leaked}")
    if stats["n_total"] == 0:
        problems.append("empty input")
    return {"passed": len(problems) == 0, "problems": problems}


def run(in_csv, out_csv, manifest_json, trace_json, pd, Chem,
        protonate_fn=None, rdMolHash=None, ph=DEFAULT_PH,
        precision=DEFAULT_PKA_PRECISION, dry_run=True,
        timeout_s=DEFAULT_PER_MOL_TIMEOUT_S, log_every=200, dominant_rule=False):
    """Full node: read -> enumerate -> QC -> (write if dry_run=False and QC passes)."""
    in_df = pd.read_csv(in_csv)
    out, stats = enumerate_microstates(
        in_df, pd, Chem, protonate_fn=protonate_fn, rdMolHash=rdMolHash,
        ph=ph, precision=precision, timeout_s=timeout_s, log_every=log_every,
        dominant_rule=dominant_rule)
    gate = qc_gate(out, stats)
    manifest = {
        "module": "Node01_microstate_enumerate",
        "version": NODE01_VERSION,
        "created_at_utc": utc_now(),
        "dry_run": dry_run,
        "policy": {
            "backend": "Dimorphite-DL",
            "pH": ph, "pka_precision": precision,
            "keep": "single dominant microstate (ambiguous flagged, not silently picked)",
            "dominant_rule": dominant_rule,
            "dominant_rule_acids": list(PHYSIO_ACID_SMARTS.keys()) if dominant_rule else [],
            "dominant_rule_bases": list(PHYSIO_BASE_SMARTS.keys()) if dominant_rule else [],
            "tautomers": "not auto-canonicalized (flagged for review)",
            "undefined_stereo": "not auto-assigned (flagged for review)",
        },
        "invariant": "heavy-atom formula + HetAtomProtomer hash preserved vs input",
        "input": str(in_csv),
        "qc_gate": gate,
        "stats": stats,
        "note": ("Authoritative charge is decided here, recorded, and never "
                 "re-guessed downstream. qm_final_status gives the definitive QM "
                 "validity now that the true charge is known."),
    }
    if not dry_run and gate["passed"]:
        Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(out_csv, index=False)
        Path(manifest_json).write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
        trace = {
            "created_at_utc": utc_now(),
            "input_csv": str(in_csv), "input_sha256": file_sha256(in_csv),
            "output_csv": str(out_csv), "output_sha256": file_sha256(out_csv),
            "module_version": NODE01_VERSION,
            "n_total": stats["n_total"],
            "protonation_status_counts": stats["protonation_status_counts"],
            "qm_final_status_counts": stats["qm_final_status_counts"],
        }
        Path(trace_json).write_text(
            json.dumps(trace, indent=2, ensure_ascii=False), encoding="utf-8")
    return out, manifest


# ============ v0.7.0: correction audit (checkagain24) ============
def build_correction_audit(df, pandas_mod):
    """从 Node 01 结果表提取 v0.7.0 校正审计：
    只列 resolved_by ∈ {rule_constructed} 或 dimorphite_conflict=True 的分子，
    记录 分子 / Dimorphite候选 / 规则最终电荷 / 原因。
    输出 Node01_microstate_correction_audit.csv 供人工签字与方法学引用。
    """
    pd = pandas_mod
    rows = []
    for _, r in df.iterrows():
        rb = str(r.get("resolved_by") or "")
        conflict = str(r.get("dimorphite_conflict")) in ("True", "true", "1")
        if rb == "rule_constructed" or conflict:
            rows.append({
                "name": r.get("name"),
                "dimorphite_candidates": r.get("candidate_charges"),
                "rule_expected_charge": r.get("expected_charge"),
                "final_net_charge": r.get("net_charge"),
                "resolved_by": rb,
                "status": r.get("protonation_status"),
                "needs_review": r.get("needs_review"),
                "review_reason": r.get("review_reason"),
            })
    return pd.DataFrame(rows)
