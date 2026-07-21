# menu/tools.py
import streamlit as st
from moudle.tools.run_cellchat import *
from moudle.tools.run_Proption import *
from moudle.tools.run_infercnv import *
from moudle.tools.run_monocle2 import *
from moudle.tools.run_SCENIC import *
from moudle.tools.run_scscore import *
# from moudle.tools.sc_score import *
# from moudle.tools.sc_paga import *
from moudle.tools.styles import *
import streamlit_antd_components as sac
if 'tool' not in st.session_state:
    st.session_state['tool'] = None

# 工具卡片列表
tool_list = [
    {'name': 'cellchat细胞通讯分析', 'img': img_to_html('moudle/template/tools/images/CellChat_Logo.png'), 'description': ''},
    {'name': '细胞比例分析', 'img': img_to_html('moudle/template/tools/images/proportion.png'), 'description': ''},
    {'name': 'inferCNV分析', 'img': img_to_html('moudle/template/tools/images/inferCNV.png'), 'description': ''},
    {'name': 'monocle2拟时序分析', 'img': img_to_html('moudle/template/tools/images/pseudotime.png'), 'description': ''},
    {'name': 'SCENIC转录因子分析', 'img': img_to_html('moudle/template/tools/images/SCENIC.jpg'), 'description': ''},
    {'name': '基因集评分分析', 'img': img_to_html('moudle/template/tools/images/score.png'), 'description': '基因集评分和可视化'},
    # {'name': 'paga轨迹分析', 'img': img_to_html('moudle/template/tools/images/paga.png'), 'description': '拟时序分析与轨迹推断'},
    # {'name': '差异基因火山图', 'img': img_to_html('moudle/template/tools/images/volcano.png'), 'description': '可视化差异基因'},
    # {'name': '组间细胞类型箱型图', 'img': img_to_html('moudle/template/tools/images/boxplot.png'), 'description': '细胞类型在不同组之间的表达差异'}
]

def main():
    st.subheader('🍒 scRNA-seq data visualization tools')
    sac.alert(
        label='**使用说明**',
        description='**点击下面任意模块卡片，进入相应的功能分析页面。**',
        color='cyan',
        icon=True,
        size=15
    )
    
    horizontal_bar = "<hr style='margin-top: 2; margin-bottom: 0; height: 1px; border: 1px solid #635985;'><br>"
    N_cards_per_col = 3
    
    # 功能页面加载
    if st.session_state['tool']:
        if st.session_state['tool'] == 'cellchat细胞通讯分析':
            cellchat()
        elif st.session_state['tool'] == '细胞比例分析':
             Proption()
        elif st.session_state['tool'] == 'inferCNV分析':
             inferCNV()
             inferCNV_anno()
        elif st.session_state['tool'] == 'monocle2拟时序分析':
             monocle2()
        elif st.session_state['tool'] == 'SCENIC转录因子分析':
             SCENIC()
        elif st.session_state['tool'] == '基因集评分分析':
             scscore()

    else:
        for index, item in enumerate(tool_list):
            i = index % N_cards_per_col
            if i == 0:
                st.markdown(horizontal_bar, unsafe_allow_html=True)
                cols = st.columns(N_cards_per_col, gap="large")
            with cols[i]:
                st.button(label=item["name"],
                          key=item["name"],
                          on_click=lambda name=item['name']: st.session_state.update({'tool': name}),
                          type='primary',
                          use_container_width=True)
                st.markdown(item["img"], unsafe_allow_html=True)
                st.markdown(f"**:gray[{item['description']}]**")

        st.markdown(horizontal_bar, unsafe_allow_html=True)