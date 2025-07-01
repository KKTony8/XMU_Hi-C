# mds图
# 11:27 2025/7/1 tested by wky
import os 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
import re

# 输入目录
input_dir = "/data5/Wuky/2XMU_T7-HiC/MME12scc/scc_output_chr_selected"

# 获取所有 *_selected_chr_scc.txt 文件
files = [f for f in os.listdir(input_dir) if f.endswith("_selected_chr_scc.txt")]

# 使用正则提取样本名
sample_set = set()
pattern = r"(.+)\.allValidPairs_(.+)\.allValidPairs"

for f in files:
    match = re.match(pattern, f.replace("_selected_chr_scc.txt", ""))
    if match:
        sample_set.update([match.group(1), match.group(2)])

samples = sorted(sample_set)

# 初始化空的 SCC 相似性矩阵
scc_matrix = pd.DataFrame(index=samples, columns=samples, dtype=float)

# 读取每个文件并填入矩阵
for f in files:
    path = os.path.join(input_dir, f)
    match = re.match(pattern, f.replace("_selected_chr_scc.txt", ""))
    if not match:
        print(f"⚠️ 跳过不匹配的文件名: {f}")
        continue
    s1, s2 = match.group(1), match.group(2)

    with open(path) as fh:
        scores = [float(line.strip()) for line in fh if line.strip() and not line.startswith("#")]
    avg_scc = np.mean(scores)

    scc_matrix.loc[s1, s2] = avg_scc
    scc_matrix.loc[s2, s1] = avg_scc

# 自己与自己设为 1
np.fill_diagonal(scc_matrix.values, 1.0)

# 转换为距离矩阵
distance_matrix = 1 - scc_matrix

# 多维尺度分析 (MDS)
mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
coords = mds.fit_transform(distance_matrix)

# 可视化 MDS 散点图
plt.figure(figsize=(10, 8))
for i, sample in enumerate(scc_matrix.index):
    plt.scatter(coords[i, 0], coords[i, 1], s=100)
    plt.text(coords[i, 0] + 0.01, coords[i, 1] + 0.01, sample, fontsize=9)

plt.title("MDS Plot of Samples Based on SCC (Selected Chromosomes)")
plt.xlabel("MDS1")
plt.ylabel("MDS2")
plt.tight_layout()
plt.savefig("scc_mds_selected_chr.png", dpi=300)
plt.show()

