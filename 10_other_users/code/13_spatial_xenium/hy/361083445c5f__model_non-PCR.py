#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
多模型预测、可视化及模型保存脚本
输出: results/model_*.pkl, 模型评估表, 高质量图片
"""

import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import (accuracy_score, roc_auc_score, classification_report,
                             confusion_matrix, roc_curve, precision_score,
                             recall_score, f1_score)
import os
import joblib

# ---------- 可选模型 ----------
try:
    from xgboost import XGBClassifier
    HAS_XGB = True
except ImportError:
    HAS_XGB = False
    print("XGBoost 未安装，跳过 XGBoost 模型。可用 pip install xgboost 安装。")

try:
    from lightgbm import LGBMClassifier
    HAS_LGB = True
except ImportError:
    HAS_LGB = False
    print("LightGBM 未安装，跳过 LightGBM 模型。可用 pip install lightgbm 安装。")

# ================= 配置参数 =================
INPUT_FILE = "integrated.tsv"
TARGET_COL = "state_state"
ID_COL = "sample"
POS_LABEL = "non-pCR"
TEST_SIZE = 0.2
RANDOM_STATE = 42
N_FEATURES_TO_SELECT = 50
SAVE_DIR = "results-npCR"

# ================= 准备输出目录 =================
os.makedirs(SAVE_DIR, exist_ok=True)

# 全局绘图风格 (使用 seaborn 设置，兼容性更好)
sns.set_style("whitegrid")
sns.set_context("notebook", font_scale=1.2)
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']  # 可选中文支持
plt.rcParams['axes.unicode_minus'] = False

# ================= 1. 数据加载与预处理 =================
print(f"读取 {INPUT_FILE} ...")
df = pd.read_csv(INPUT_FILE, sep='\t', dtype=str, keep_default_na=False)
df.set_index(ID_COL, inplace=True)
print(f"数据形状: {df.shape}")

if TARGET_COL not in df.columns:
    raise ValueError(f"未找到目标列 '{TARGET_COL}'")

y_raw = df[TARGET_COL]
X_raw = df.drop(columns=[TARGET_COL])
print(f"目标分布:\n{y_raw.value_counts()}")

# 转换为数值
X = X_raw.apply(pd.to_numeric, errors='coerce')

# 编码目标
le = LabelEncoder()
y = le.fit_transform(y_raw)
print(f"类别映射: {dict(zip(le.classes_, le.transform(le.classes_)))}")
if le.classes_[1] != POS_LABEL:
    y = (y_raw == POS_LABEL).astype(int)

# 删除高缺失列 (>50%)
missing_ratio = X.isnull().mean()
X = X.loc[:, missing_ratio < 0.5]

# 删除常量列
constant_cols = X.columns[X.nunique(dropna=True) <= 1]
X.drop(columns=constant_cols, inplace=True)
print(f"预处理后特征数: {X.shape[1]}")

# 填补缺失值
imputer = SimpleImputer(strategy='median')
X_imputed = pd.DataFrame(imputer.fit_transform(X), columns=X.columns, index=X.index)

# 标准化
scaler = StandardScaler()
X_scaled = pd.DataFrame(scaler.fit_transform(X_imputed), columns=X.columns, index=X.index)

# ================= 2. 拆分数据 =================
X_train, X_test, y_train, y_test = train_test_split(
    X_scaled, y, test_size=TEST_SIZE, stratify=y, random_state=RANDOM_STATE
)
print(f"训练集: {X_train.shape[0]} 样本, 测试集: {X_test.shape[0]} 样本")

# ================= 3. 特征选择 =================
print(f"ANOVA 选择前 {N_FEATURES_TO_SELECT} 个特征...")
selector = SelectKBest(f_classif, k=min(N_FEATURES_TO_SELECT, X_train.shape[1]))
selector.fit(X_train, y_train)
selected_features = X_train.columns[selector.get_support()]
print(f"选中特征: {len(selected_features)} 个")

X_train_sel = selector.transform(X_train)
X_test_sel = selector.transform(X_test)

# 保存选中特征
pd.Series(selected_features).to_csv(os.path.join(SAVE_DIR, "selected_features.txt"), index=False, header=False)

# ================= 4. 模型定义 =================
models = {
    'Logistic Regression': LogisticRegression(penalty='l2', C=1.0, solver='lbfgs',
                                              max_iter=1000, random_state=RANDOM_STATE,
                                              class_weight='balanced'),
    'Random Forest': RandomForestClassifier(n_estimators=200, max_depth=10,
                                            random_state=RANDOM_STATE, class_weight='balanced'),
    'SVM (RBF)': SVC(kernel='rbf', probability=True, random_state=RANDOM_STATE,
                     class_weight='balanced'),
    'Gradient Boosting': GradientBoostingClassifier(n_estimators=200, max_depth=5,
                                                    random_state=RANDOM_STATE),
    'KNN (k=5)': KNeighborsClassifier(n_neighbors=5)
}

if HAS_XGB:
    models['XGBoost'] = XGBClassifier(n_estimators=200, max_depth=5, random_state=RANDOM_STATE,
                                      eval_metric='logloss', verbosity=0)
if HAS_LGB:
    models['LightGBM'] = LGBMClassifier(n_estimators=200, max_depth=5, random_state=RANDOM_STATE,
                                        verbose=-1)

# ================= 5. 训练与评估 =================
results = {}
roc_data = {}
best_auc = 0.0
best_model = None
best_model_name = ""

for name, model in models.items():
    print(f"\n训练模型: {name}")
    model.fit(X_train_sel, y_train)
    y_pred = model.predict(X_test_sel)
    y_prob = model.predict_proba(X_test_sel)[:, 1]

    # ---------- 基础指标 ----------
    acc = accuracy_score(y_test, y_pred)
    auc_val = roc_auc_score(y_test, y_prob)
    # 针对阳性类 (pCR) 的指标
    prec_pos = precision_score(y_test, y_pred, pos_label=1)
    rec_pos = recall_score(y_test, y_pred, pos_label=1)
    f1_pos = f1_score(y_test, y_pred, pos_label=1)
    # 宏平均指标（适合不平衡数据）
    prec_macro = precision_score(y_test, y_pred, average='macro')
    rec_macro = recall_score(y_test, y_pred, average='macro')
    f1_macro = f1_score(y_test, y_pred, average='macro')

    results[name] = {
        'Accuracy': acc,
        'AUC': auc_val,
        'Precision (pCR)': prec_pos,
        'Recall (pCR)': rec_pos,
        'F1 (pCR)': f1_pos,
        'Precision (macro)': prec_macro,
        'Recall (macro)': rec_macro,
        'F1 (macro)': f1_macro
    }

    print(f"{name} | Acc={acc:.4f} AUC={auc_val:.4f} | "
          f"Prec(pCR)={prec_pos:.3f} Rec(pCR)={rec_pos:.3f} F1(pCR)={f1_pos:.3f}")
    print(classification_report(y_test, y_pred, target_names=le.classes_))

    # ROC数据
    fpr, tpr, _ = roc_curve(y_test, y_prob)
    roc_data[name] = (fpr, tpr, auc_val)

    # ---------- 保存模型 ----------
    pkl_name = name.replace(' ', '_').replace('(', '').replace(')', '') + ".pkl"
    joblib.dump(model, os.path.join(SAVE_DIR, pkl_name))
    print(f"模型已保存: {pkl_name}")

    # ---------- 记录最佳模型 ----------
    if auc_val > best_auc:
        best_auc = auc_val
        best_model = model
        best_model_name = name
        best_y_pred = y_pred
        best_y_prob = y_prob

# ---------- 交叉验证（最佳模型） ----------
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)
cv_scores = cross_val_score(best_model, X_train_sel, y_train, cv=cv, scoring='roc_auc')
print(f"\n最佳模型 ({best_model_name}) 5折CV AUC: {cv_scores.mean():.4f} ± {cv_scores.std():.4f}")

# ---------- 保存详细评估表 ----------
res_df = pd.DataFrame(results).T
res_df.index.name = 'Model'
res_df = res_df.round(4)
res_df.to_csv(os.path.join(SAVE_DIR, "model_metrics.csv"))
print("\n详细评估表已保存到 model_metrics.csv")
print(res_df)

# ================= 6. 可视化 =================
# 统一调色板（美观兼顾色盲友好）
palette = sns.color_palette("Set2", len(models))
model_names = list(models.keys())

# ---------- 6.1 特征重要性（树模型） ----------
tree_models = {k: v for k, v in models.items() if hasattr(v, 'feature_importances_')}
if tree_models:
    n_tree = len(tree_models)
    fig, axes = plt.subplots(1, n_tree, figsize=(6 * n_tree, 8))
    if n_tree == 1:
        axes = [axes]
    for ax, (name, model) in zip(axes, tree_models.items()):
        importances = model.feature_importances_
        idx = np.argsort(importances)[-20:]
        feat_names = selected_features[idx]
        colors = sns.color_palette("viridis", len(idx))
        ax.barh(range(len(idx)), importances[idx], color=colors)
        ax.set_yticks(range(len(idx)))
        ax.set_yticklabels(feat_names, fontsize=9)
        ax.set_title(f'{name} - Top 20 Features', fontweight='bold')
        ax.set_xlabel('Importance')
        ax.xaxis.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(os.path.join(SAVE_DIR, "feature_importance.png"), dpi=200, bbox_inches='tight')
    plt.close()
    print("特征重要性图已保存")

# ---------- 6.2 模型性能对比（柱状图） ----------
metrics_plot = ['Accuracy', 'AUC', 'F1 (pCR)', 'Recall (pCR)', 'Precision (pCR)']
fig, axes = plt.subplots(2, 3, figsize=(18, 10))
axes = axes.flatten()
for i, metric in enumerate(metrics_plot):
    if i >= len(axes):
        break
    ax = axes[i]
    values = [results[name][metric] for name in model_names]
    bars = ax.bar(model_names, values, color=palette, edgecolor='white', linewidth=0.8)
    ax.set_title(metric, fontsize=13, fontweight='bold')
    ax.set_ylim(0, 1.1)
    ax.set_ylabel(metric)
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    # 在柱子上标注数值
    for bar, v in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, v + 0.02, f'{v:.3f}',
                ha='center', va='bottom', fontsize=9, fontweight='bold')
    ax.tick_params(axis='x', rotation=45)
# 隐藏多余子图
for j in range(i+1, len(axes)):
    axes[j].set_visible(False)
plt.suptitle('Model Performance Comparison', fontsize=16, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(os.path.join(SAVE_DIR, "model_comparison.png"), dpi=200, bbox_inches='tight')
plt.close()
print("模型对比图已保存")

# ---------- 6.3 精美 ROC 曲线 ----------
fig, ax = plt.subplots(figsize=(9, 7))
# 使用更协调的颜色循环
colors_roc = sns.color_palette("Dark2", len(roc_data))
for idx, (name, (fpr, tpr, roc_auc_val)) in enumerate(roc_data.items()):
    ax.plot(fpr, tpr, color=colors_roc[idx], lw=2.5,
            label=f'{name} (AUC = {roc_auc_val:.3f})')

# 对角线参考线
ax.plot([0, 1], [0, 1], 'k--', lw=1.5, alpha=0.6, label='Random Classifier')

# 美化细节
ax.set_xlim([-0.02, 1.02])
ax.set_ylim([0.0, 1.05])
ax.set_xlabel('False Positive Rate (1 - Specificity)', fontsize=12, fontweight='bold')
ax.set_ylabel('True Positive Rate (Sensitivity)', fontsize=12, fontweight='bold')
ax.set_title('ROC Curves (Test Set)', fontsize=15, fontweight='bold', pad=20)
ax.legend(loc='lower right', fontsize=9, frameon=True, fancybox=True, framealpha=0.9)
ax.grid(True, linestyle='--', alpha=0.5)

# 添加注释突出最佳模型
ax.annotate(f'Best: {best_model_name}\nAUC = {best_auc:.3f}',
            xy=(0.6, 0.2), xycoords='axes fraction',
            fontsize=10, fontweight='bold', color='darkred',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray', alpha=0.8))

# 紧凑布局
plt.tight_layout()
plt.savefig(os.path.join(SAVE_DIR, "roc_curves.png"), dpi=300, bbox_inches='tight')
plt.close()
print("精美 ROC 曲线图已保存")

# ---------- 6.4 混淆矩阵（最佳模型） ----------
cm = confusion_matrix(y_test, best_y_pred)
fig, ax = plt.subplots(figsize=(7, 6))
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', cbar=False,
            xticklabels=le.classes_, yticklabels=le.classes_,
            annot_kws={'fontsize': 14, 'fontweight': 'bold'})
ax.set_xlabel('Predicted Label', fontsize=12, fontweight='bold')
ax.set_ylabel('True Label', fontsize=12, fontweight='bold')
ax.set_title(f'Confusion Matrix - {best_model_name}', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(SAVE_DIR, "confusion_matrix.png"), dpi=200, bbox_inches='tight')
plt.close()
print("混淆矩阵已保存")

print(f"\n所有结果已保存至 '{SAVE_DIR}' 目录。")
print("脚本执行完毕。")
