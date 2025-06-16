####主要用于统计四种模式的序列数量，根据设计的barcode和umi来做区分
##22:55 2025/6/16 wky
#!/bin/bash

DEFAULT_READS=1000000

# 定义四个模式及其按照NNN分割的片段
declare -A PATTERNS
PATTERNS[1]="TAATACGACTCACTATAGGGNNNNNNNNCGAAACATCGGCCACCCCNNNNNNNNCAACCCTGCGAC"
FRAGMENTS1=(
    "TAATACGACTCACTATAGGG"    # 第一段
    "CGAAACATCGGCCACCCC"      # NNNNNNNN后的第二段
    "CAACCCTGCGAC"            # NNNNNNNN后的第三段
)

PATTERNS[2]="CGCAATGAAGTCGCAGGGTTGNNNNNNNNGGGGTGGCCGATGTTTCGNNNNNNNNCCCTATAGTGAGTCGTATTA"
FRAGMENTS2=(
    "CGCAATGAAGTCGCAGGGTTG"    # 第一段
    "GGGGTGGCCGATGTTTCG"       # NNNNNNNN后的第二段
    "CCCTATAGTGAGTCGTATTA"     # NNNNNNNN后的第三段
)

# PATTERNS[3]="TAATACGACTCACTATAGGGCGAAACATCGGCCACNNNNNNNNCCCNNNNNNNNCAACCCTGCGA"
# FRAGMENTS3=(
#     "TAATACGACTCACTATAGGGCGAAACATCGGCCAC"    # 第一段
#     "CCC"                                     # NNNNNNNN后的第二段
#     "CAACCCTGCGA"                            # NNNNNNNN后的第三段
# )

PATTERNS[3]="TAATACGACTCACTATAGGGGTGGCCGATGTTTCGNNNNNNNNCCCNNNNNNNNCAACCCTGCGACTTCA"
FRAGMENTS3=(
    "TAATACGACTCACTATAGGGGTGGCCGATGTTTCG"    # 第一段
    "CCC"                                     # NNNNNNNN后的第二段
    "CAACCCTGCGACTTCA"                            # NNNNNNNN后的第三段
)

PATTERNS[4]="CAATGAAGTCGCAGGGTTGNNNNNNNNGGGNNNNNNNNCGAAACATCGGCCACCCCTATAGTGAGTCGTATTA"
FRAGMENTS4=(
    "CAATGAAGTCGCAGGGTTG"                    # 第一段
    "GGG"                                     # NNNNNNNN后的第二段
    "CGAAACATCGGCCACCCCTATAGTGAGTCGTATTA"    # NNNNNNNN后的第三段
)

mkdir -p head${DEFAULT_READS}_merged

all_samples_file="head${DEFAULT_READS}_merged/all_samples_analysis.txt"

# 写入文件头
{
    echo "===== Hi-C Samples Analysis Report ====="
    echo "Analysis Date: $(date)"
    echo "Reads analyzed per sample: $DEFAULT_READS"
    echo -e "\n=== Sample Statistics ===\n"
    echo -e "Sample\tTotal_Reads\tMatched_Reads\tPercentage\tPattern"
    echo "--------------------------------------------------------------------------------"
} > $all_samples_file

# 函数：使用多个片段进行grep
grep_fragments() {
    local file=$1
    shift
    local fragments=("$@")
    local temp_file=$(mktemp)
    
    cp "$file" "$temp_file"
    
    for fragment in "${fragments[@]}"; do
        grep "$fragment" "$temp_file" > "${temp_file}.new"
        mv "${temp_file}.new" "$temp_file"
    done
    
    local count=$(wc -l < "$temp_file")
    rm "$temp_file"
    echo "$count"
}

# 处理每个样本
for r1 in *_1.fq.gz; do
    base_name=${r1%_1.fq.gz}
    r2="${base_name}_2.fq.gz"
    merged_file="head${DEFAULT_READS}_merged/${base_name}_merged.fq"

    echo "处理样本: $base_name"

    # 提取并合并reads
    paste <(zcat $r1 | head -n $((DEFAULT_READS*4))) \
          <(zcat $r2 | head -n $((DEFAULT_READS*4))) > $merged_file

    total_reads=$DEFAULT_READS

    # 处理所有模式
    for i in {1..4}; do
        echo "分析模式$i..."
        pattern=${PATTERNS[$i]}
        fragments_var="FRAGMENTS${i}[@]"
        matched_count=$(grep_fragments "$merged_file" "${!fragments_var}")
        percentage=$(awk "BEGIN {printf \"%.4f\", ($matched_count/$total_reads)*100}")
        echo -e "${base_name}\t${total_reads}\t${matched_count}\t${percentage}%\t${pattern}" >> $all_samples_file
    done

    echo "--------------------------------------------------------------------------------" >> $all_samples_file
done

echo "分析完成！结果已保存到：$all_samples_file"
