#!/bin/bash
# Node 05b · V3 生产驱动（协议逐字来自 V2 run_md2.sh；四处受控偏离见模块头注释）
set -euo pipefail
export PATH="/mnt/hpc/home/wanghongmei/anaconda3/envs/md_env/bin:$PATH"
set +u; source /mnt/hpc/home/wanghongmei/apps/gromacs-2025.4-cuda12.1/bin/GMXRC; set -u
cd /mnt/hpc/home/wanghongmei/lang2/cross_solv_v3/pool_mol_074__4YLU__primary/r2
LIG=MOL
echo "[MD2] $(pwd) replica=r2 seed=43 @ $(hostname) $(date -Is)"

# ---- 0) 由制备目录取输入（只读，不改制备产物）--------------------------------
mkdir -p 03_min 04_nvt 05_npt 06_md
cp -f /mnt/hpc/home/wanghongmei/lang2/cross_solv_v3/pool_mol_074__4YLU__primary/complex.gro 03_min/complex.gro
cp -f /mnt/hpc/home/wanghongmei/lang2/cross_solv_v3/pool_mol_074__4YLU__primary/topol.top   03_min/topol.top
cp -f /mnt/hpc/home/wanghongmei/lang2/cross_solv_v3/pool_mol_074__4YLU__primary/posre_prot.itp 03_min/
mkdir -p 03_min/ligand.acpype
cp -f /mnt/hpc/home/wanghongmei/lang2/cross_solv_v3/pool_mol_074__4YLU__primary/ligand.acpype/ligand_GMX.itp   03_min/ligand.acpype/
cp -f /mnt/hpc/home/wanghongmei/lang2/cross_solv_v3/pool_mol_074__4YLU__primary/ligand.acpype/posre_ligand.itp 03_min/ligand.acpype/
cd 03_min

# ---- 1) index：合并 Protein 与配体残基组；随后【按名】断言 -------------------
printf '"Protein" | "%s"\nq\n' "$LIG" | gmx make_ndx -f complex.gro -o index.ndx > make_ndx.log 2>&1
/mnt/hpc/home/wanghongmei/anaconda3/envs/md_env/bin/python -c "import sys; sys.path.insert(0,'/mnt/hpc/home/wanghongmei/lang2'); from core import node05_md_prep as N5; \
import json; g=N5.assert_ndx_groups('index.ndx', ['Protein','$LIG','Protein_$LIG','Water_and_ions','System','SOL']); \
open('ndx_groups.json','w').write(json.dumps(g, indent=1)); print('[NDX]', g)"

# ---- 2) EM（主路径）---------------------------------------------------------
cat > em.mdp <<'MDPEM'
integrator  = steep
emtol       = 500.0
emstep      = 0.01
nsteps      = 50000
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
MDPEM
gmx grompp -f em.mdp -c complex.gro -r complex.gro -p topol.top -o em.tpr -n index.ndx -maxwarn 5 > grompp_em.log 2>&1
gmx mdrun -deffnm em -ntmpi 1 -ntomp "${OMP_NUM_THREADS:-8}" >> mdrun_em.log 2>&1

# ---- 3) QC 闸（V2 阈值：fmax>500 或 mindist<0.12 nm → 修复分支）--------------
SOLG=$(/mnt/hpc/home/wanghongmei/anaconda3/envs/md_env/bin/python -c "import json;print(json.load(open('ndx_groups.json'))['SOL'])")
SYSG=$(/mnt/hpc/home/wanghongmei/anaconda3/envs/md_env/bin/python -c "import json;print(json.load(open('ndx_groups.json'))['System'])")
printf '%s\n%s\n' "$SYSG" "$SOLG" | gmx mindist -f em.gro -s em.tpr -n index.ndx -od mindist_em.xvg > mindist.log 2>&1 || true
/mnt/hpc/home/wanghongmei/anaconda3/envs/md_env/bin/python -c "
import re,json,math
fmax=None
for ln in open('em.log',errors='ignore'):
    m=re.search(r'Maximum force\s*=\s*([-\d.eE+]+)', ln)
    if m: fmax=float(m.group(1))
mind=None
try:
    rows=[l.split() for l in open('mindist_em.xvg') if l and l[0] not in '#@']
    if rows: mind=float(rows[-1][1])
except Exception: pass
legacy='repair' if ((mind is not None and mind<0.12) or (fmax is not None and fmax>500.0)) else 'main'
branch='staged'
json.dump({'fmax':fmax,'mindist_nm':mind,
           'equilibration':branch,'declared':True,
           'legacy_gate_branch':legacy,'legacy_gate_agrees':(legacy=='repair'),
           'thresholds':{'fmax':500.0,'mindist_nm':0.12},
           'note':'System-SOL mindist 恒等于水分子自身 O-H 键长（TIP3P 0.09572 nm，EM 段无约束），必 < 0.12 nm 阈值；V2 因此 138/138 走分段路径。本值仅作诊断记录，不再决定分支。'}, open('qc.json','w'), indent=1)
print('[QC] fmax=',fmax,' mindist=',mind,' legacy_gate=',legacy,' -> equilibration=',branch)
"
BRANCH=$(/mnt/hpc/home/wanghongmei/anaconda3/envs/md_env/bin/python -c "import json;print(json.load(open('qc.json'))['equilibration'])")

if [ "$BRANCH" = "staged" ]; then
  cat > em_strong.mdp <<'MDPEMS'
integrator  = steep
emtol       = 100.0
emstep      = 0.01
nsteps      = 200000
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
MDPEMS
  gmx grompp -f em_strong.mdp -c em.gro -r em.gro -p topol.top -o em_strong.tpr -n index.ndx -maxwarn 5 > grompp_ems.log 2>&1
  gmx mdrun -deffnm em_strong -ntmpi 1 -ntomp "${OMP_NUM_THREADS:-8}" >> mdrun_ems.log 2>&1
  EM_GRO=em_strong.gro
else
  EM_GRO=em.gro
fi
cd ..

# ---- 4) NVT -----------------------------------------------------------------
cd 04_nvt
for f in topol.top posre_prot.itp index.ndx; do cp -f ../03_min/$f .; done
mkdir -p ligand.acpype && cp -f ../03_min/ligand.acpype/* ligand.acpype/
cp -f "../03_min/$EM_GRO" start.gro
if [ "$BRANCH" = "staged" ]; then
  cat > nvt50.mdp  <<'MDPN50'
define = -DPOSRES -DPOSRESLIG
integrator  = md
dt          = 0.002
constraint_algorithm = lincs
constraints = h-bonds
lincs_iter  = 2
lincs_order = 4
cutoff-scheme = Verlet
nstlist     = 20
rcoulomb    = 1.0
rvdw        = 1.0
coulombtype = PME
tcoupl      = V-rescale
tc-grps     = Protein_MOL Water_and_ions
tau_t       = 0.2 0.2
ref_t       = 50 50
pbc         = xyz
nstxout     = 0
nstvout     = 0
nsteps      = 10000
continuation = no
pcoupl      = no
gen_vel     = yes
gen_temp    = 50
gen_seed    = 43
nstenergy   = 500
nstlog      = 500
MDPN50
  cat > nvt300.mdp <<'MDPN300'
define = -DPOSRES -DPOSRESLIG
integrator  = md
dt          = 0.002
constraint_algorithm = lincs
constraints = h-bonds
lincs_iter  = 2
lincs_order = 4
cutoff-scheme = Verlet
nstlist     = 20
rcoulomb    = 1.0
rvdw        = 1.0
coulombtype = PME
tcoupl      = V-rescale
tc-grps     = Protein_MOL Water_and_ions
tau_t       = 0.2 0.2
ref_t       = 300 300
pbc         = xyz
nstxout     = 0
nstvout     = 0
nsteps      = 40000
continuation = yes
pcoupl      = no
gen_vel     = no
nstenergy   = 500
nstlog      = 500
MDPN300
  gmx grompp -f nvt50.mdp -c start.gro -r start.gro -p topol.top -o nvt50.tpr -n index.ndx -maxwarn 5 > grompp_nvt50.log 2>&1
  gmx mdrun -deffnm nvt50 -ntmpi 1 -ntomp "${OMP_NUM_THREADS:-8}" >> mdrun_nvt50.log 2>&1
  gmx grompp -f nvt300.mdp -c nvt50.gro -r nvt50.gro -p topol.top -o nvt300.tpr -n index.ndx -maxwarn 5 > grompp_nvt300.log 2>&1
  gmx mdrun -deffnm nvt300 -ntmpi 1 -ntomp "${OMP_NUM_THREADS:-8}" >> mdrun_nvt300.log 2>&1
  NVT_GRO=nvt300.gro
else
  cat > nvt.mdp <<'MDPNVT'
define = -DPOSRES -DPOSRESLIG
integrator  = md
dt          = 0.002
constraint_algorithm = lincs
constraints = h-bonds
lincs_iter  = 2
lincs_order = 4
cutoff-scheme = Verlet
nstlist     = 20
rcoulomb    = 1.0
rvdw        = 1.0
coulombtype = PME
tcoupl      = V-rescale
tc-grps     = Protein_MOL Water_and_ions
tau_t       = 0.2 0.2
ref_t       = 300 300
pbc         = xyz
nstxout     = 0
nstvout     = 0
nsteps      = 50000
continuation = no
pcoupl      = no
gen_vel     = yes
gen_temp    = 300
gen_seed    = 43
nstenergy   = 500
nstlog      = 500
MDPNVT
  gmx grompp -f nvt.mdp -c start.gro -r start.gro -p topol.top -o nvt.tpr -n index.ndx -maxwarn 5 > grompp_nvt.log 2>&1
  gmx mdrun -deffnm nvt -ntmpi 1 -ntomp "${OMP_NUM_THREADS:-8}" >> mdrun_nvt.log 2>&1
  NVT_GRO=nvt.gro
fi
cd ..

# ---- 5) NPT -----------------------------------------------------------------
cd 05_npt
for f in topol.top posre_prot.itp index.ndx; do cp -f ../04_nvt/$f .; done
mkdir -p ligand.acpype && cp -f ../04_nvt/ligand.acpype/* ligand.acpype/
cp -f "../04_nvt/$NVT_GRO" nvt.gro
cat > npt.mdp <<'MDPNPT'
define = -DPOSRES -DPOSRESLIG
integrator  = md
dt          = 0.002
constraint_algorithm = lincs
constraints = h-bonds
lincs_iter  = 2
lincs_order = 4
cutoff-scheme = Verlet
nstlist     = 20
rcoulomb    = 1.0
rvdw        = 1.0
coulombtype = PME
tcoupl      = V-rescale
tc-grps     = Protein_MOL Water_and_ions
tau_t       = 0.2 0.2
ref_t       = 300 300
pbc         = xyz
nstxout     = 0
nstvout     = 0
nsteps      = 50000
continuation = yes
pcoupl      = C-rescale
pcoupltype  = isotropic
refcoord_scaling = com
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
gen_vel     = no
nstenergy   = 500
nstlog      = 500
MDPNPT
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr -n index.ndx -maxwarn 5 > grompp_npt.log 2>&1
gmx mdrun -deffnm npt -ntmpi 1 -ntomp "${OMP_NUM_THREADS:-8}" >> mdrun_npt.log 2>&1
cd ..

# ---- 6) 生产 MD -------------------------------------------------------------
cd 06_md
for f in topol.top index.ndx; do cp -f ../05_npt/$f .; done
mkdir -p ligand.acpype && cp -f ../05_npt/ligand.acpype/* ligand.acpype/
cp -f ../05_npt/npt.gro .
cat > md.mdp <<'MDPMD'
integrator  = md
dt          = 0.002
constraint_algorithm = lincs
constraints = h-bonds
lincs_iter  = 2
lincs_order = 4
cutoff-scheme = Verlet
nstlist     = 20
rcoulomb    = 1.0
rvdw        = 1.0
coulombtype = PME
tcoupl      = V-rescale
tc-grps     = Protein_MOL Water_and_ions
tau_t       = 0.5 0.5
ref_t       = 300 300
pbc         = xyz
nstxout     = 0
nstvout     = 0
nsteps      = 50000000
continuation = yes
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
nstxout-compressed = 5000
compressed-x-grps   = System
nstenergy   = 1000
nstlog      = 1000
MDPMD
gmx grompp -f md.mdp -c npt.gro -t ../05_npt/npt.cpt -p topol.top -o md.tpr -n index.ndx -maxwarn 5 > grompp_md.log 2>&1
gmx mdrun -deffnm md -ntmpi 1 -ntomp "${OMP_NUM_THREADS:-8}" -nb gpu -pme gpu >> mdrun_md.log 2>&1
cd ..

cp -f 03_min/qc.json 03_min/ndx_groups.json 06_md/ 2>/dev/null || true
echo "[MD2 DONE] $(pwd) branch=$BRANCH $(date -Is)"
