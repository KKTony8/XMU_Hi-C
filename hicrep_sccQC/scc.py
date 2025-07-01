# 画图1
# 11:27 2025/7/1 tested by wky
import os
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import re

# 输入目录（你的文件目录）
input_dir = "/data5/Wuky/2XMU_T7-HiC/MME12scc/scc_output_chr_selected"

# 获取所有 *_selected_chr_scc.txt 文件
files = [f for f in os.listdir(input_dir) if f.endswith("_selected_chr_scc.txt")]

# 提取所有样本名
sample_set = set()
pattern = r"(.+)\.allValidPairs_(.+)\.allValidPairs"

for f in files:
    match = re.match(pattern, f.replace("_selected_chr_scc.txt", ""))
    if match:
        sample_set.add(match.group(1))
        sample_set.add(match.group(2))

samples = sorted(sample_set)

# 初始化空矩阵
scc_matrix = pd.DataFrame(index=samples, columns=samples, dtype=float)

# 读取每个文件并写入矩阵
for f in files:
    path = os.path.join(input_dir, f)
    match = re.match(pattern, f.replace("_selected_chr_scc.txt", ""))
    if not match:
        print(f"⚠️ 跳过不匹配的文件名: {f}")
        continue

    s1, s2 = match.group(1), match.group(2)

    # 读取数值
    with open(path) as fh:
        scores = [float(line.strip()) for line in fh if line.strip() and not line.startswith("#")]

    avg_scc = np.mean(scores)
    scc_matrix.loc[s1, s2] = avg_scc
    scc_matrix.loc[s2, s1] = avg_scc

# 自己与自己为 1
np.fill_diagonal(scc_matrix.values, 1.0)

# 绘制热图
plt.figure(figsize=(12, 10))
sns.heatmap(scc_matrix, annot=True, fmt=".2f", cmap="YlGnBu", linewidths=0.5, cbar_kws={'label': 'SCC'})
plt.title("SCC Matrix Heatmap (Selected Chromosomes)")
plt.tight_layout()
plt.savefig("scc_heatmap_selected_chr.png", dpi=300)
plt.show()


