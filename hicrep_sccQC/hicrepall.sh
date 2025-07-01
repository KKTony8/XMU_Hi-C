# 11:27 2025/7/1 tested by wky
# !/bin/bash

# 分析指定染色体的相关性脚本

# 输入和输出目录
DIRIN="/data5/Wuky/2XMU_T7-HiC/MME12scc/mcool"
OUTDIR="/data5/Wuky/2XMU_T7-HiC/MME12scc/scc_output_chr_selected"
mkdir -p "$OUTDIR"

# 收集所有 .mcool 文件的样本名（不包含路径和扩展名）
samples=($(ls "$DIRIN"/*.mcool | xargs -n1 basename | sed 's/\.mcool$//'))

total=$((${#samples[@]} - 1))

# 分析以下染色体：chr1 - chr22 和 chrX
CHR_LIST=(
  chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10
  chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19
  chr20 chr21 chr22 chrX
)
CHR_ARGS="${CHR_LIST[*]}"

# 两两配对计算相关性（去除重复组合）
for i in $(seq 0 $total); do
  for j in $(seq 0 $total); do
    if [ $j -gt $i ]; then
      F1="${DIRIN}/${samples[$i]}.mcool"
      F2="${DIRIN}/${samples[$j]}.mcool"
      OUT="${OUTDIR}/${samples[$i]}_${samples[$j]}_selected_chr_scc.txt"

      echo "Running: ${samples[$i]} vs ${samples[$j]} on selected chromosomes"
      hicrep "$F1" "$F2" "$OUT" \
        --binSize 500000 \
        --h 10 \
        --dBPMax 5000000 \
        --chrNames $CHR_ARGS
    fi
  done
 done
