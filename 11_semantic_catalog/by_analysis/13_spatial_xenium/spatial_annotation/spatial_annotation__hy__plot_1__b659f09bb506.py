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
data_text = """M2	0.00545166579787127	MPC
M2	0.00182910354727164	MPC
M2	0.00170923145620583	MPC
TAM	0.00154693515797679	MPC
TAM	0.000641502889653808	MPC
EMT_M2	0.000915339984513329	MPC
M2	0.000835044557623978	MPC
M2	0.0054549741564803	TNBC
M2	0.00280095857432525	TNBC
M2	0.00927320974838455	TNBC
Plasma	0.0024198072346738	TNBC
EMT_M2	0.000651092694431434	TNBC
M2	0.000628992631748082	TNBC
Plasma	0.000378826417261969	TNBC
TAM	0.00194029669792495	TNBC
TAM	0.00188938162545932	TNBC
TAM	0.000844219565255443	TNBC
TAM	0.0017560635398065	TNBC
TAM	0.00166220368571021	TNBC
TAM	0.00130019845357917	TNBC
TAM	0.000649714183941881	TNBC
TAM	0.00128115193066788	TNBC
TAM	0.000887040813600131	TNBC
TAM	0.000909088286127572	TNBC
TAM	0.00287072029577308	TNBC
TAM	0.00175489587082543	TNBC
TAM	0.000914066493357609	TNBC
TAM	0.00119490592738643	TNBC
TAM	0.00126162335912151	TNBC
TAM	0.00130023284893752	TNBC
TAM	0.000705788864464765	TNBC
TAM	0.00139772567360608	TNBC
TAM	0.00140674757310219	TNBC
TAM	0.0016976492415185	TNBC
TAM	0.00143863994504641	TNBC
TAM	0.00121082583658931	TNBC
TAM	0.000780841823740699	TNBC
TAM	0.00108202134361006	TNBC
TAM	0.00127077385566119	TNBC
TAM	0.00151704666400672	TNBC
TAM	0.00151748986327856	TNBC
TAM	0.00180375649136567	TNBC
TAM	0.00200770844596778	TNBC
TAM	0.00456694760933359	TNBC
TAM	0.00129271532476743	TNBC
TAM	0.00490385531365691	TNBC
TAM	0.00176983743434707	TNBC
TAM	0.00343559539966315	TNBC
TAM	0.00117425939329746	TNBC
TAM	0.00652290974924774	TNBC
TAM	0.0018607296581162	TNBC
TAM	0.00182508961012983	TNBC
TAM	0.00134519442598953	TNBC
TAM	0.00219719869350334	TNBC
TAM	0.00194835624535158	TNBC
TAM	0.00185898386575851	TNBC
TAM	0.00229499086647696	TNBC
TAM	0.0048751232174389	TNBC
TAM	0.00556417975467369	TNBC
TAM	0.00140545148441842	TNBC
TAM	0.00131723556698685	TNBC
TAM	0.0011117427871007	TNBC
TAM	0.00107672094908343	TNBC
TAM	0.00162124194095223	TNBC
TAM	0.00134095516636341	TNBC
TAM	0.00179999502236863	TNBC
TAM	0.00105903646604981	TNBC
TAM	0.000830052506958743	TNBC
M2	0.00046009436409258	TNBC
TAM	0.0010743544346271	TNBC
M2	0.000392306209529828	TNBC
TAM	0.000916147021898349	TNBC
M2	0.00197614334023278	TNBC
TAM	0.00460511591917613	TNBC
M2	0.00114950514419169	TNBC
TAM	0.00268171027258272	TNBC
EMT_M2	0.000666358594591945	TNBC
CD8+_Tem	0.00123989525574929	TNBC
EMT_M2	0.000841504540206756	TNBC
TAM	0.000933089859696571	TNBC
Treg	0.00107168117003257	TNBC
CD8+_Tem	0.00160084407707539	TNBC
EMT_M2	0.00108660308104066	TNBC
TAM	0.00120483158419629	TNBC
Treg	0.00138372863780819	TNBC
CD8+_Tem	0.000841504540206756	TNBC
EMT_M2	0.000571047547351511	TNBC
TAM	0.000633216279040957	TNBC
Treg	0.000727299959889105	TNBC
CD8+_Tem	0.000754117895210175	TNBC
EMT_M2	0.000511732319033766	TNBC
TAM	0.000567447194394516	TNBC
Treg	0.000651765235559171	TNBC
CD8+_Tem	0.000717946734301477	TNBC
EMT_M2	0.000487181480638064	TNBC
TAM	0.000540224827296839	TNBC
Treg	0.000620500357323316	TNBC
CD8+_Tem	0.000933089859696571	TNBC
EMT_M2	0.000633216279040957	TNBC
TAM	0.000702148433317156	TNBC
Treg	0.000806465832813286	TNBC
CD8+_Tem	0.00107168117003257	TNBC
EMT_M2	0.000727299959889105	TNBC
TAM	0.000806465832813286	TNBC
Treg	0.000926267184831123	TNBC
CD8+_Tem	0.00112177138239964	TNBC
EMT_M2	0.0007613061150605	TNBC
TAM	0.000844170396948419	TNBC
Treg	0.000969567366846619	TNBC
CD8+_Tem	0.00136768783940531	TNBC
EMT_M2	0.000928274255023023	TNBC
TAM	0.00102929344582671	TNBC
Treg	0.00118215691283218	TNBC
CD8+_Tem	0.000972421254676811	TNBC
EMT_M2	0.000659915819209548	TNBC
TAM	0.00073175236693063	TNBC
Treg	0.000840464290252409	TNBC
CD8+_Tem	0.00159690571903837	TNBC
EMT_M2	0.00108392847281302	TNBC
TAM	0.0012018663137468	TNBC
Treg	0.00138032368347368	TNBC
EMT_M2	0.00138280219332649	TNBC
MKI67+_Plasma	0.000828025927258347	TNBC
EMT_M2	0.00138842739034965	TNBC
EMT_M2	0.00361767688210125	TNBC
MKI67+_Plasma	0.00216822208826587	TNBC
EMT_M2	0.00276483530662263	TNBC
MKI67+_Plasma	0.00165651108093306	TNBC
Plasma	0.00263337852755136	TNBC
Plasma	0.00849166079014639	TNBC
Plasma	0.0136188250020638	TNBC
CD8+_Tem	0.00107129044061207	TNBC
CD8+_Tex	0.000466249811779121	TNBC
CD8+_Tem	0.00159857562462838	TNBC
CD8+_Tex	0.000695943714367758	TNBC
CD8+_Tem	0.00106860172576315	TNBC
CD8+_Tex	0.000465078915318225	TNBC
CD8+_Tem	0.00157355114688895	TNBC
CD8+_Tex	0.000685039553129211	TNBC
CD8+_Tem	0.00158555472722809	TNBC
CD8+_Tex	0.000690269952459127	TNBC
CD8+_Tem	0.00109447745670193	TNBC
CD8+_Tex	0.000476347569927373	TNBC
CD8+_Tem	0.00205723794195742	TNBC
CD8+_Tex	0.000895855831692534	TNBC
CD8+_Tem	0.00599879661782808	TNBC
CD8+_Tex	0.00261810572147212	TNBC
CD8+_Tem	0.00141613619065777	TNBC
CD8+_Tex	0.00061645466612874	TNBC
CD8+_Tem	0.001197055054652	TNBC
CD8+_Tex	0.000521022461589637	TNBC
CD8+_Tem	0.0014397766738899	TNBC
CD8+_Tex	0.000626753923919822	TNBC
CD8+_Tem	0.00193296740399901	TNBC
CD8+_Tex	0.000841681149024708	TNBC
CD8+_Tem	0.000804601179153954	TNBC
CD8+_Tex	0.000350127824967873	TNBC
CD8+_Tem	0.000474167030442112	TNBC
CD8+_Tex	0.000206298561092566	TNBC
CD8+_Tem	0.0028927341735248	TNBC
CD8+_Tex	0.00126028114495084	TNBC
CD8+_Tem	0.00369755798038671	TNBC
CD8+_Tex	0.00161165386065406	TNBC
CD8+_Tem	0.00174976647042081	TNBC
CD8+_Tex	0.000761830151243012	TNBC
CD8+_Tem	0.00114156495431725	TNBC
CD8+_Tex	0.000496854609442175	TNBC
CD8+_Tem	0.00303800705244729	TNBC
CD8+_Tex	0.00132368119569747	TNBC
CD8+_Tem	0.000779335338577501	TNBC
CD8+_Tex	0.00033912837422768	TNBC
CD8+_Tem	0.000879770612784026	TNBC
CD8+_Tex	0.000382854597917298	TNBC
CD8+_Tem	0.000852797168521802	TNBC
CD8+_Tex	0.000371110757439011	TNBC
CD8+_Tem	0.00135836058413037	TNBC
CD8+_Tex	0.000591285194563794	TNBC
CD8+_Tem	0.00135016970630721	TNBC
CD8+_Tex	0.0005877170385785	TNBC
Plasma	0.00307754716910107	TNBC"""

# Read data
df = pd.read_csv(StringIO(data_text), sep='\t', header=None, names=['cell_type', 'intensity', 'group'])

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
        size = cell_data['count'] * 120
        color_val = cell_data['mean_intensity']
        ax.scatter(idx, y_positions['MPC'], s=size, c=[color_val], 
                   cmap=red_cmap, norm=norm, alpha=0.85, edgecolors='black', 
                   linewidth=1.2, zorder=3)

# Plot TNBC group
for idx, cell_type in enumerate(all_cell_types):
    if cell_type in tnbc_data['cell_type'].values:
        cell_data = tnbc_data[tnbc_data['cell_type'] == cell_type].iloc[0]
        size = cell_data['count'] * 120
        color_val = cell_data['mean_intensity']
        ax.scatter(idx, y_positions['TNBC'], s=size, c=[color_val], 
                   cmap=red_cmap, norm=norm, alpha=0.85, edgecolors='black', 
                   linewidth=1.2, zorder=3)

# Formatting
ax.set_xlim(-0.8, len(all_cell_types)-0.2)
ax.set_ylim(-0.8, 0.8)
ax.set_xticks(range(len(all_cell_types)))
ax.set_xticklabels(all_cell_types, rotation=45, ha='right', fontsize=10, fontname='Arial')
ax.set_yticks(list(y_positions.values()))
ax.set_yticklabels(['MPC', 'TNBC'], fontsize=12, fontweight='bold', fontname='Arial')
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

# Add colorbar only (no size legend)
cbar_ax = fig.add_axes([0.92, 0.25, 0.012, 0.5])  # [left, bottom, width, height]
sm = ScalarMappable(cmap=red_cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, cax=cbar_ax)
cbar.set_label('Mean intensity', fontsize=10, fontname='Arial', rotation=270, labelpad=15)
cbar.ax.tick_params(labelsize=9)

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
