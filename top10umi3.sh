export TMPDIR="/data5/Wuky/tmp"

# ������ʱĿ¼����������ڣ�
mkdir -p "$TMPDIR"

# ִ��ԭ���sort ʹ�� -T ����ָ����ʱĿ¼
zcat XM-QGY-MME-2_L5_1.fq.gz | \
awk 'NR%4==2 {count[$0]++} END {for (seq in count) print count[seq], seq}' | \
sort -nr -T "$TMPDIR" | \
head -10

