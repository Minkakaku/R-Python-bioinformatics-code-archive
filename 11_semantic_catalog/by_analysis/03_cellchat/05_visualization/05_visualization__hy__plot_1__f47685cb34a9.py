import matplotlib.pyplot as plt
import numpy as np

# ========================
# 1. 整理数据（完全基于您的 CellChat 结果）
# ========================
# 唯一的相互作用名称（取并集）
interactions = [
    'JAG1_NOTCH3',
    'SEMA4C_PLXNB2',
    'PDGFB_PDGFRB',
    'COL4A1_CD44',
    'COL4A2_CD44',
    'LAMC1_CD44',
    'CD46_JAG1',
    'CD99_CD99',
    'DLL4_NOTCH3'
]

# 正常组 (out_a) 的概率（没有的填 0）
out_a_probs = [0.000310, 0.000286, 0, 0, 0, 0, 0, 0, 0]

# 肿瘤组 (in_d) 的概率（没有的填 0）
in_d_probs  = [0.001123, 0.001125, 0.002527, 0.012314, 0.006875, 0.000863, 0.000723, 0.001864, 0.001663]

# Y 轴位置（0,1,2... 对应每个相互作用）
y_pos = np.arange(len(interactions))

# 气泡大小缩放系数（让小的概率也能看见，大的不会太大）
scale_factor = 50000  

# ========================
# 2. 开始绘图
# ========================
fig, ax = plt.subplots(figsize=(12, 6))

# 绘制正常组 (X=0)
ax.scatter(
    [0] * len(interactions), 
    y_pos, 
    s=[p * scale_factor for p in out_a_probs], 
    c='#74b9ff', 
    label='Normal_Endo (out_a)', 
    alpha=0.8, 
    edgecolors='black'
)

# 绘制肿瘤组 (X=1)
ax.scatter(
    [1] * len(interactions), 
    y_pos, 
    s=[p * scale_factor for p in in_d_probs], 
    c='#ff7675', 
    label='Tumor_Endo (in_d)', 
    alpha=0.8, 
    edgecolors='black'
)

# ========================
# 3. 美化坐标轴和标签
# ========================
ax.set_xticks([0, 1])
ax.set_xticklabels(['Normal Endo (out_a)', 'Tumor Endo (in_d)'], fontsize=12, fontweight='bold')
ax.set_yticks(y_pos)
ax.set_yticklabels(interactions, fontsize=11)

ax.set_xlabel('Group', fontsize=12)
ax.set_ylabel('Ligand-Receptor Interaction', fontsize=12)
ax.set_title('Bubble Plot: Signaling Strength to MCAM+ CAF', fontsize=14, fontweight='bold')

# 添加网格线（仅水平，方便读数）
ax.grid(axis='y', linestyle='--', alpha=0.5)

# 添加图例
ax.legend(fontsize=11, loc='upper right')

# 自动调整布局防止文字被切
plt.tight_layout()

# ========================
# 4. 保存为 PDF 和 PNG
# ========================
plt.savefig('bubble_plot.pdf')
plt.savefig('bubble_plot.png', dpi=300)

