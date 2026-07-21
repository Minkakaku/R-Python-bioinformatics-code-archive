import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import matplotlib.colors as mcolors
from io import StringIO
import matplotlib.font_manager as font_manager

# 添加字体文件
font_path = '/home/hy/software/miniconda3/envs/sprna/lib/python3.11/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf'
try:
    font_manager.fontManager.addfont(font_path)
    # 设置全局字体（使用字体的 family name，通常是 'Arial'）
    plt.rcParams['font.family'] = 'Arial'
except Exception as e:
    print(f"警告: 未能加载 Arial 字体: {e}")
    plt.rcParams['font.family'] = 'sans-serif'  # 回退

# Data preparation
data_text = """source	prob	group
M2	0.00555352892859266	MPC
M2	0.00182258313372215	MPC
M2	0.00217631560383282	MPC
TAM	0.000454113381772543	MPC
TAM	0.00286061565261657	MPC
TAM	0.00135905872945475	MPC
Plasma	0.00198849039347941	TNBC
EMT_M2	0.000441411744114992	TNBC
M2	0.000419089154014708	TNBC
TAM	0.0014738183083879	TNBC
TAM	0.000587025147818693	TNBC
TAM	0.000962077647263893	TNBC
TAM	0.000817232202864425	TNBC
TAM	0.000286720990373654	TNBC
TAM	0.000507651577712132	TNBC
TAM	0.000786641351321821	TNBC
TAM	0.000972033716522412	TNBC
TAM	0.00128951470347535	TNBC
TAM	0.00058623657283379	TNBC
TAM	0.000461635509617658	TNBC
TAM	0.000823164524002055	TNBC
TAM	0.00108317615190638	TNBC
TAM	0.000930892582869553	TNBC
TAM	0.000459543367815981	TNBC
TAM	0.000471855273288566	TNBC
TAM	0.00115860320487449	TNBC
TAM	0.00118034434198157	TNBC
TAM	0.000998150822304758	TNBC
TAM	0.000978861717540318	TNBC
TAM	0.00127646946806862	TNBC
TAM	0.00129702382902558	TNBC
TAM	0.00302450687829117	TNBC
TAM	0.00140679520059449	TNBC
TAM	0.00323255932210229	TNBC
TAM	0.00249924206098873	TNBC
TAM	0.00209719855827243	TNBC
TAM	0.000903941618911605	TNBC
TAM	0.00348500002624691	TNBC
TAM	0.000921590840250826	TNBC
TAM	0.00120891764621607	TNBC
TAM	0.000857461066653748	TNBC
TAM	0.00163361455548266	TNBC
TAM	0.00117463211327686	TNBC
TAM	0.00102021584504643	TNBC
TAM	0.00103327166250966	TNBC
TAM	0.00151608188498444	TNBC
TAM	0.00319269032838356	TNBC
TAM	0.00126977703370197	TNBC
TAM	0.00359652650512519	TNBC
TAM	0.000910588225698885	TNBC
TAM	0.000894525029494343	TNBC
TAM	0.000888097100753974	TNBC
TAM	0.00122489750523733	TNBC
TAM	0.00129352627693648	TNBC
TAM	0.0015686112430082	TNBC
TAM	0.00156614291099353	TNBC
TAM	0.000766121409646375	TNBC
TAM	0.00342301267457813	TNBC
TAM	0.00189348479586798	TNBC
EMT_M2	0.000813349175028155	TNBC
TAM	0.000989266045827652	TNBC
EMT_M2	0.000458612894326154	TNBC
TAM	0.000557847747053862	TNBC
EMT_M2	0.000451758556740567	TNBC
TAM	0.00054951108037052	TNBC
EMT_M2	0.000500950054331834	TNBC
TAM	0.000609340255024886	TNBC
TAM	0.000268388461441982	TNBC
EMT_M2	0.000557847747053862	TNBC
TAM	0.000678540502701879	TNBC
EMT_M2	0.000820889320734857	TNBC
TAM	0.000998435397584158	TNBC
EMT_M2	0.000425590752036504	TNBC
TAM	0.000517683960327461	TNBC
EMT_M2	0.000488010086109336	TNBC
TAM	0.000593602138050502	TNBC
EMT_M2	0.000913189463467507	TNBC
MKI67+_Plasma	0.00070755201514371	TNBC
MKI67+_T	0.000694769195728792	TNBC
pDC	0.000820693753900979	TNBC
EMT_M2	0.00115921359769468	TNBC
EMT_M2	0.00276700587579141	TNBC
MKI67+_Plasma	0.00214481117142934	TNBC
MKI67+_T	0.00210611717884227	TNBC
pDC	0.00248720709284128	TNBC
EMT_M2	0.00231598916987919	TNBC
MKI67+_Plasma	0.00179502847066394	TNBC
MKI67+_T	0.001762633672051	TNBC
pDC	0.00208170185456759	TNBC
Plasma	0.00134956581367739	TNBC
Plasma	0.00608255855902309	TNBC
CD8+_Tex	0.000289329426838063	TNBC
CD8+_Tex	0.000851865429364284	TNBC
CD8+_Tex	0.000956324632800727	TNBC
CD8+_Tex	0.00101147111790347	TNBC
CD8+_Tex	0.000795290342142928	TNBC
CD8+_Tex	0.000813984458731083	TNBC
CD8+_Tex	0.00119060534656959	TNBC
CD8+_Tex	0.00124042725562389	TNBC
CD8+_Tex	0.00292026149440131	TNBC
CD8+_Tex	0.000442288686810524	TNBC
CD8+_Tex	0.000817148405123764	TNBC
CD8+_Tex	0.000965151190116266	TNBC
CD8+_Tex	0.0007715523054052	TNBC
CD8+_Tex	0.000632637132515376	TNBC
CD8+_Tex	0.000470073851218137	TNBC
CD8+_Tex	0.00181725593525824	TNBC
CD8+_Tex	0.00229140492851075	TNBC
CD8+_Tex	0.00111363576204303	TNBC
CD8+_Tex	0.000627638357929477	TNBC
CD8+_Tex	0.000303396683597245	TNBC
CD8+_Tex	0.00056553826979568	TNBC
CD8+_Tex	0.000772744826170499	TNBC
CD8+_Tex	0.000463850233113275	TNBC
Plasma	0.00172526644688183	TNBC"""

# Read data
df = pd.read_csv(StringIO(data_text), sep='\t', header=0, names=['cell_type', 'intensity', 'group'])

# Prepare bubble data
def prepare_bubble_data(df, group_name):
    group_df = df[df['group'] == group_name]
    bubble_data = group_df.groupby('cell_type').agg(
        count=('intensity', 'size'),
        mean_intensity=('intensity', 'mean')
    ).reset_index()
    bubble_data = bubble_data.sort_values('count', ascending=False)
    return bubble_data

mpc_data = prepare_bubble_data(df, 'MPC')
tnbc_data = prepare_bubble_data(df, 'TNBC')

# Get all unique cell types for x-axis
all_cell_types = sorted(set(mpc_data['cell_type'].tolist() + tnbc_data['cell_type'].tolist()))

# Color mapping range
all_intensities = pd.concat([mpc_data['mean_intensity'], tnbc_data['mean_intensity']])
vmin, vmax = all_intensities.min(), all_intensities.max()

# Custom red colormap (CNS style)
colors = ['#FFF5F5', '#FFE5E5', '#FFCCCC', '#FF9999', '#FF6666', '#FF3333', '#FF0000', '#CC0000', '#990000']
red_cmap = mcolors.LinearSegmentedColormap.from_list('custom_red', colors, N=256)
norm = Normalize(vmin=vmin, vmax=vmax)

# Create figure with adjusted layout for external legends
fig, ax = plt.subplots(1, 1, figsize=(16, 6))

# Y-axis positions for two groups
y_positions = {'MPC': 0.33, 'TNBC': -0.33}

# Plot MPC group
for idx, cell_type in enumerate(all_cell_types):
    if cell_type in mpc_data['cell_type'].values:
        cell_data = mpc_data[mpc_data['cell_type'] == cell_type].iloc[0]
        size = cell_data['count'] * 150  # Adjusted scaling factor
        color_val = cell_data['mean_intensity']
        ax.scatter(idx, y_positions['MPC'], s=size, c=[color_val], 
                   cmap=red_cmap, norm=norm, alpha=0.85, edgecolors='black', 
                   linewidth=1.2, zorder=3)

# Plot TNBC group
for idx, cell_type in enumerate(all_cell_types):
    if cell_type in tnbc_data['cell_type'].values:
        cell_data = tnbc_data[tnbc_data['cell_type'] == cell_type].iloc[0]
        size = cell_data['count'] * 150  # Adjusted scaling factor
        color_val = cell_data['mean_intensity']
        ax.scatter(idx, y_positions['TNBC'], s=size, c=[color_val], 
                   cmap=red_cmap, norm=norm, alpha=0.85, edgecolors='black', 
                   linewidth=1.2, zorder=3)

# Formatting
ax.set_xlim(-0.8, len(all_cell_types)-0.2)
ax.set_ylim(-0.8, 0.8)
ax.set_xticks(range(len(all_cell_types)))
ax.set_xticklabels(all_cell_types, rotation=45, ha='right', fontsize=14, fontname='Arial')
ax.set_yticks(list(y_positions.values()))
ax.set_yticklabels(['MPC', 'TNBC'], fontsize=16, fontweight='bold', fontname='Arial')
ax.tick_params(axis='both', which='major', length=4, width=1)

# Add horizontal lines for grouping
ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.8, alpha=0.5, zorder=1)
ax.axhline(y=y_positions['MPC']+0.2, color='lightgray', linestyle='--', linewidth=0.5, alpha=0.5)
ax.axhline(y=y_positions['TNBC']-0.2, color='lightgray', linestyle='--', linewidth=0.5, alpha=0.5)

# Remove spines except bottom
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_linewidth(1.2)

# Set background
ax.set_facecolor('white')
fig.patch.set_facecolor('white')

# Add colorbar only
cbar_ax = fig.add_axes([0.92, 0.25, 0.012, 0.5])
sm = ScalarMappable(cmap=red_cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, cax=cbar_ax)
cbar.set_label('Mean intensity', fontsize=14, fontname='Arial', rotation=270, labelpad=15)
cbar.ax.tick_params(labelsize=13)

# Adjust layout to accommodate external colorbar
plt.subplots_adjust(right=0.88)

# Save as PDF and PNG
plt.savefig('bubble_plot.pdf', format='pdf', dpi=300, bbox_inches='tight')
plt.savefig('bubble_plot.png', format='png', dpi=300, bbox_inches='tight')
print("Figures saved as bubble_plot.pdf and bubble_plot.png")

# Display
plt.show()

# Output statistics
print("\n=== MPC Group Statistics ===")
print(mpc_data[['cell_type', 'count', 'mean_intensity']].to_string(index=False))
print("\n=== TNBC Group Statistics ===")
print(tnbc_data[['cell_type', 'count', 'mean_intensity']].to_string(index=False))
