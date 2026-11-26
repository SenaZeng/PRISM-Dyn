#!/bin/bash
#SBATCH --job-name=prism_n04
#SBATCH --partition=cpuPartition
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/n04_%A_%a.out
#SBATCH --error=logs/n04_%A_%a.err
# 用法（--array 与 WM/NT 在提交时给出）：
#   sbatch --array=0-352%32 --export=WM=<work_manifest>,NT=353 scripts/node04_array.sbatch
# Slurm 环境陷阱（已实测）：srun/sbatch 进入 base 而非 md_env。
# 故一律绝对路径，不依赖 PATH / conda activate。
set -euo pipefail
ROOT=/mnt/hpc/home/wanghongmei/lang2
PY=/mnt/hpc/home/wanghongmei/anaconda3/envs/md_env/bin/python
VINA=/mnt/hpc/home/wanghongmei/apps/vina/vina_1.2.7_linux_x86_64
cd "$ROOT"
"$PY" -m core.node04_worker \
      --work-manifest "${WM}" \
      --task-id "${SLURM_ARRAY_TASK_ID}" --n-tasks "${NT}" \
      --vina "$VINA" --workers "${SLURM_CPUS_PER_TASK}"
