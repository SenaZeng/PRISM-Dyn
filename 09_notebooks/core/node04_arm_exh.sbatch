#!/bin/bash
#SBATCH --job-name=n04_exh
#SBATCH --partition=cpuPartition
#SBATCH --cpus-per-task=5
#SBATCH --output=logs/n04exh_%A_%a.out
#SBATCH --error=logs/n04exh_%A_%a.err
# exhaustiveness 对照臂专用 —— 与 core/node04_array.sbatch 分开，
# 使 Node 04 的冻结参数不被任何方式改动。
# 用法（465 = 93 × 5：每个 array 任务 = 一个"配体-靶点"的 5 个种子）：
#   sbatch --array=0-92%93 \
#          --export=ALL,WM=<manifest>,NT=93,EXH=128 core/node04_arm_exh.sbatch
set -euo pipefail
ROOT=/mnt/hpc/home/wanghongmei/lang2
PY=/mnt/hpc/home/wanghongmei/anaconda3/envs/md_env/bin/python
VINA=/mnt/hpc/home/wanghongmei/apps/vina/vina_1.2.7_linux_x86_64
cd "$ROOT"
"$PY" -m core.node04_worker \
      --work-manifest "${WM}" \
      --task-id "${SLURM_ARRAY_TASK_ID}" --n-tasks "${NT}" \
      --vina "$VINA" --workers "${SLURM_CPUS_PER_TASK}" \
      --exhaustiveness "${EXH}"
