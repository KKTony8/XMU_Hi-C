####主要是用于统计占前10份额的reads是些什么？
##22:55 2025/6/16 wky
export TMPDIR="/data5/Wuky/tmp"

# 创建临时目录（如果不存在）
mkdir -p "$TMPDIR"

# 执行原命令，sort 使用 -T 参数指定临时目录
zcat XM-QGY-MME-2_L5_1.fq.gz | \
awk 'NR%4==2 {count[$0]++} END {for (seq in count) print count[seq], seq}' | \
sort -nr -T "$TMPDIR" | \
head -10
