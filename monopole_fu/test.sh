#!/bin/bash
#
#SBATCH --job-name=test               # 作业名
#SBATCH --nodes=1                     # 使用 1 个节点
#SBATCH --ntasks-per-node=1           # 每个节点一个任务
#SBATCH --cpus-per-task=8             # 每个任务使用 8 个 CPU 核
#SBATCH --time=00:05:00               # 最长运行时间 5 分钟
#SBATCH --output=slurm-%j.out         # 标准输出重定向
#SBATCH --error=slurm-%j.err          # 标准错误重定向

# —— 加载运行环境 —— 
module load matlab/2019b

# 切换到脚本所在目录
cd /gs/home/by2430102

# —— 运行 MATLAB 程序，并添加所需路径 —— 
matlab -nodesktop -nosplash -r "\
addpath(genpath('/gs/home/by2430102/fwd_matlab_patch')); \
addpath(genpath('/gs/home/by2430102/irfu-matlab-master')); \
Monopole_search_MMS_HPC; \
exit"
