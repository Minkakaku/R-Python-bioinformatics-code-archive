import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats as stats
from scipy.stats import mannwhitneyu
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import itertools

# Read data
df = pd.read_csv('aaa', sep='\t')

# Display basic information about the data
print("Data preview:")
print(df.head())
print("\nData shape:", df.shape)
print("\nGroups in data:", df['group'].unique())

# Calculate CN11 percentage for each sample
def calculate_cn11_percentage(group):
    total_cells = len(group)
    cn11_cells = len(group[group['CNs_15'] == 'CN11'])
    percentage = (cn11_cells / total_cells) * 100 if total_cells > 0 else 0
    return percentage

# Calculate CN11 percentage for each sample
sample_stats = df.groupby(['sample', 'group']).apply(calculate_cn11_percentage).reset_index()
sample_stats.columns = ['sample', 'group', 'cn11_percentage']

print("\nCN11 percentage by sample:")
print(sample_stats)

# Create boxplot
plt.figure(figsize=(10, 6))

# Create boxplot using seaborn
sns.boxplot(x='group', y='cn11_percentage', data=sample_stats, palette='Set2')

# Add individual data points
sns.stripplot(x='group', y='cn11_percentage', data=sample_stats, 
              color='black', alpha=0.6, jitter=True, size=6)

# Set titles and labels in English
plt.title('CN11 Percentage Distribution by Group', fontsize=16, fontweight='bold')
plt.xlabel('Group', fontsize=12)
plt.ylabel('CN11 Cell Percentage (%)', fontsize=12)

# Add grid for better readability
plt.grid(axis='y', alpha=0.3)

# Adjust layout
plt.tight_layout()

# Save the plot
plt.savefig('cn11_percentage_boxplot.pdf', dpi=300, bbox_inches='tight')

# Display the plot
plt.show()

# Print statistical summary
print("\nStatistical summary by group:")
group_summary = sample_stats.groupby('group')['cn11_percentage'].describe()
print(group_summary)

# Statistical tests for pairwise group comparisons
print("\n" + "="*60)
print("PAIRWISE GROUP COMPARISONS")
print("="*60)

# Get unique groups
groups = sample_stats['group'].unique()
print(f"Groups to compare: {groups}")

# Perform pairwise comparisons using Mann-Whitney U test (non-parametric)
print("\nMann-Whitney U Test Results (non-parametric):")
print("-" * 50)

# Generate all possible pairwise combinations
pairwise_combinations = list(itertools.combinations(groups, 2))

# Store results
results = []

for group1, group2 in pairwise_combinations:
    data1 = sample_stats[sample_stats['group'] == group1]['cn11_percentage']
    data2 = sample_stats[sample_stats['group'] == group2]['cn11_percentage']
    
    # Perform Mann-Whitney U test
    stat, p_value = mannwhitneyu(data1, data2, alternative='two-sided')
    
    # Apply Bonferroni correction for multiple comparisons
    n_comparisons = len(pairwise_combinations)
    corrected_alpha = 0.05 / n_comparisons
    significant = p_value < corrected_alpha
    
    results.append({
        'group1': group1,
        'group2': group2,
        'statistic': stat,
        'p_value': p_value,
        'significant': significant,
        'corrected_alpha': corrected_alpha
    })
    
    print(f"{group1} vs {group2}:")
    print(f"  Statistic = {stat:.4f}, p-value = {p_value:.4f}")
    print(f"  Significant (after Bonferroni correction, α={corrected_alpha:.4f}): {significant}")
    print(f"  Interpretation: {'Significant difference' if significant else 'No significant difference'}")
    print()

# If data meets assumptions, also perform parametric test (T-test)
print("\nT-test Results (parametric):")
print("-" * 30)

for group1, group2 in pairwise_combinations:
    data1 = sample_stats[sample_stats['group'] == group1]['cn11_percentage']
    data2 = sample_stats[sample_stats['group'] == group2]['cn11_percentage']
    
    # Perform T-test
    t_stat, p_value = stats.ttest_ind(data1, data2)
    
    # Apply Bonferroni correction
    corrected_alpha = 0.05 / len(pairwise_combinations)
    significant = p_value < corrected_alpha
    
    print(f"{group1} vs {group2}:")
    print(f"  t-statistic = {t_stat:.4f}, p-value = {p_value:.4f}")
    print(f"  Significant (after Bonferroni correction): {significant}")
    print()

# Create a summary table of significant results
significant_pairs = [f"{res['group1']}-{res['group2']}" for res in results if res['significant']]
non_significant_pairs = [f"{res['group1']}-{res['group2']}" for res in results if not res['significant']]

print("\nSUMMARY OF SIGNIFICANT DIFFERENCES:")
print("-" * 40)
if significant_pairs:
    print("Significantly different pairs:")
    for pair in significant_pairs:
        print(f"  - {pair}")
else:
    print("No significant differences found after multiple comparison correction.")

if non_significant_pairs:
    print("\nNon-significant pairs:")
    for pair in non_significant_pairs:
        print(f"  - {pair}")

# Kruskal-Wallis test for overall group differences (non-parametric ANOVA)
print("\n" + "="*60)
print("OVERALL GROUP DIFFERENCES (Kruskal-Wallis Test)")
print("="*60)

# Prepare data for Kruskal-Wallis test
group_data = []
for group in groups:
    group_data.append(sample_stats[sample_stats['group'] == group]['cn11_percentage'])

# Perform Kruskal-Wallis test
h_stat, p_value_kw = stats.kruskal(*group_data)

print(f"Kruskal-Wallis H-statistic: {h_stat:.4f}")
print(f"P-value: {p_value_kw:.4f}")
print(f"Overall significant difference: {p_value_kw < 0.05}")

# Save results to Excel
with pd.ExcelWriter('cn11_analysis_results.xlsx') as writer:
    sample_stats.to_excel(writer, sheet_name='Sample_Stats', index=False)
    group_summary.to_excel(writer, sheet_name='Group_Summary')
    
    # Create results dataframe
    results_df = pd.DataFrame(results)
    results_df.to_excel(writer, sheet_name='Statistical_Tests', index=False)

print("\nAnalysis completed! Results saved to current directory.")
print(f"Boxplot saved as: cn11_percentage_boxplot.png")
print(f"Statistical results saved as: cn11_analysis_results.xlsx")
