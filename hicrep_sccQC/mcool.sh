# 11:27 2025/7/1 tested by wky
# 输入hicPro处理得到类似SRR5005048.allValidPairs.hic文件
# !/bin/bash

# 输入目录
IN_DIR="/data5/Wuky/2XMU_T7-HiC/MME12scc"
# 输出目录
OUT_DIR="${IN_DIR}/mcool"
mkdir -p "$OUT_DIR"

# 遍历所有 .hic 文件
for infile in "$IN_DIR"/*.hic; do
    basename=$(basename "$infile" .hic)
    outfile="${OUT_DIR}/${basename}.mcool"

    echo "Converting $infile to $outfile (all resolutions)"

    hic2cool convert "$infile" "$outfile" -r 0
done

