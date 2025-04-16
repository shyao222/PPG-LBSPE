import streamlit as st
import pandas as pd
import plotly.express as px
from rdkit import Chem
from rdkit.Chem import Draw
from streamlit_plotly_events import plotly_events
from st_aggrid import AgGrid, GridOptionsBuilder

df = pd.read_csv("gen_pred_final.csv")
df['t+ class'] = df['t+ class'].replace({0: '<0.5', 1: '>0.5'})
df = df.reset_index().rename(columns={'index': 'Sample ID'})

# 页面配置
st.set_page_config(
    page_title="GP-LBSPE",
    page_icon=":battery:",
    layout="wide",
    initial_sidebar_state="expanded"
)

# 初始化session state
if 'selected_index' not in st.session_state:
    st.session_state.selected_index = None

# ===== 侧边栏筛选 =====
with st.sidebar:
    st.header("Filter Samples")
    
    # 连续变量筛选
    filters = {
        "conductivity": st.slider(
            "Log Conductivity (log S/cm)",
            float(df["conductivity"].min()),
            float(df["conductivity"].max()),
            (float(df["conductivity"].min()), float(df["conductivity"].max()))
        ),
        "Tg": st.slider(
            "Glass Transition Temperature (K)",
            float(df["Tg"].min()),
            float(df["Tg"].max()),
            (float(df["Tg"].min()), float(df["Tg"].max()))
        ),
        "Td": st.slider(
            "Thermal Decomposition Temperature (K)",
            float(df["Td"].min()),
            float(df["Td"].max()),
            (float(df["Td"].min()), float(df["Td"].max()))
        )
    }

    # 类别变量筛选（复选框）
    st.markdown("Lithium ion transference number")
    col1, col2 = st.columns(2)
    with col1:
        select_less = st.checkbox("<0.5", value=True)
    with col2:
        select_greater = st.checkbox("&gt;0.5", value=True) #用 “>” 不会显示，将 > 替换为 HTML 实体 &gt;
    
    t_class_selected = []
    if select_less:
        t_class_selected.append("<0.5")
    if select_greater:
        t_class_selected.append(">0.5")
    
    # 如果 t_class_selected 为空（两个复选框都未选中），直接返回空的 DataFrame
    if not t_class_selected:
        filtered_df = pd.DataFrame(columns=df.columns)
    else:
        # 应用筛选条件
        filter_conditions = ((df["conductivity"].between(*filters["conductivity"])) &
                            (df["Tg"].between(*filters["Tg"])) &
                            (df["Td"].between(*filters["Td"])) &
                            (df["t+ class"].isin(t_class_selected)))
        filtered_df = df[filter_conditions]

# ===== 主界面 =====
st.header("Generated Polymers for Lithium Battery Solid Ploymer Electrolyte")

# 创建两列布局
col_table, col_struct = st.columns([5, 2])

with col_table:
    #st.subheader("📊 Filtered Molecules Data")
    # 搜索栏（带功能按钮）
    st.markdown("""
                <style>
                    div[data-testid="column"] {
                        align-self: flex-end !important;
                    }
                    div[data-testid="stVerticalBlock"] > div[data-testid="stHorizontalBlock"] {
                        align-items: end;
                    }
                    .stTextInput input {
                        height: 46px;
                    }
                    .stButton button {
                        height: 46px !important;
                        padding: 0 20px !important;
                    }
                </style>
                """, unsafe_allow_html=True)
    
    cols = st.columns([5, 1])  # 调整列宽比例
    with cols[0]:
        search_input = st.text_input("Search",
                                    value=st.session_state.get("search_term", ""),
                                    placeholder="Search by Sample ID or SMILES",
                                    label_visibility="collapsed")
    with cols[1]:
        search_clicked = st.button("🔍 Search", use_container_width=True)
 
    # 处理搜索逻辑
    if search_clicked or (search_input and search_input != st.session_state.get("prev_search", "")):
        st.session_state.search_term = search_input
        st.session_state.prev_search = search_input
        st.rerun()

    # 应用搜索条件
    if "search_term" in st.session_state and st.session_state.search_term:
        search_term = st.session_state.search_term
        try:
            # 尝试精确匹配 Sample ID
            sample_id = int(search_term)
            search_condition = (filtered_df["Sample ID"] == sample_id)
        except ValueError:
            # 模糊匹配 SMILES
            search_condition = (filtered_df["SMILES"].str.contains(search_term, case=False, regex=False, na=False))
        filtered_df = filtered_df[search_condition]
    
    # 数据展示（使用AgGrid）
    # 如果 filtered_df 为空，创建一个空的 DataFrame，列名与原始数据一致
    if filtered_df.empty:
            filtered_df = pd.DataFrame(columns=df.columns)
    
    # 更新 Grid
    gb = GridOptionsBuilder.from_dataframe(filtered_df)
    gb.configure_column("Sample ID", width=75)
    gb.configure_column("SMILES", width=250)
    gb.configure_column("conductivity", width=90)
    gb.configure_column("t+ class", width=60)
    gb.configure_column("Td", width=70)
    gb.configure_column("Tg", width=70)
    
    gb.configure_selection('single', use_checkbox=False)
    grid_options = gb.build()

    # 显示表格
    grid_response = AgGrid(filtered_df,
                           gridOptions=grid_options,
                           height=500,
                           width='100%', 
                           fit_columns_on_grid_load=True,
                           theme='streamlit',
                           update_mode='SELECTION_CHANGED')
    st.markdown(f'<p style="text-align:right; color:#808080; font-size:0.9em;">{len(filtered_df)} samples shown</p>', 
                unsafe_allow_html=True)
    
    # 获取选中行
    selected_rows = pd.DataFrame(grid_response['selected_rows'])
    if not selected_rows.empty:  # 正确判断是否为空
        st.session_state.selected_index = selected_rows.iloc[0]['Sample ID']

    
# 分子结构展示
with col_struct:
    selected_sample_id = st.session_state.get('selected_index', None)
    if selected_sample_id is not None:
        # 检查 selected_sample_id 是否在原始数据中
        if selected_sample_id not in df["Sample ID"].values:
            st.warning(f"Sample ID {selected_sample_id} not found in original data.")
        else:
            # 获取选中行的数据
            selected_row = df[df["Sample ID"] == selected_sample_id].iloc[0]
            smiles = selected_row['SMILES']
            # 绘制分子结构
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                img = Draw.MolToImage(mol, size=(400, 300))
                st.image(img, caption=f"Sample ID: {selected_sample_id}")
                st.write(f"**Polymer SMILES:** {smiles}")
                st.write(f"**Log Conductivity:** {selected_row['conductivity']:.4f} log S/cm")
                st.write(f"**Glass Transition Temperature:** {selected_row['Tg']:.3f} K")
                st.write(f"**Thermal Decomposition Temperature:** {selected_row['Td']:.3f} K")
                st.write(f"**Lithium ion transference number:** {selected_row['t+ class']}")
            else:
                st.warning("Invalid SMILES structure.")
    else:
        st.info("Select a sample to view structure.")

# 性能指标全称映射
performance_metrics_full = {
    "conductivity": "Log Conductivity (log S/cm)",
    "Tg": "Glass Transition Temperature (K)",
    "Td": "Thermal Decomposition Temperature (K)"
}
performance_metrics = list(performance_metrics_full.keys())

COLOR_DISCRETE_MAP = {"<0.5": "#FF5252", ">0.5": "#4CAF50"}
CUSTOM_DATA_COLS = ["Sample ID", "t+ class"]
HOVER_TEMPLATE = """
<b>Sample ID:</b> %{{customdata[0]}}<br>
<b>X ({x_label}):</b> %{{x}}<br>
<b>Y ({y_label}):</b> %{{y}}<br>
<b>Li<sup>+</sup> transference number:</b> %{{customdata[1]}}
"""

# 创建并排布局
col1, col2 = st.columns(2)

with col1:
    # X 轴选择器
    x_axis = st.selectbox(
        "Select X Axis for Scatter Plot:",
        options=performance_metrics,
        format_func=lambda x: performance_metrics_full[x],
        key="x_axis_selectbox"
    )

with col2:
    # Y 轴选择器（动态排除已选X轴）
    y_axis = st.selectbox(
        "Select Y Axis for Scatter Plot:",
        options=[m for m in performance_metrics if m != x_axis],
        format_func=lambda x: performance_metrics_full[x],
        key="y_axis_selectbox"
    )
# 数据验证和可视化
if filtered_df.empty:
    st.warning("⚠️ No data available for visualization with current filters")
else:
    # 创建交互式散点图
    fig = px.scatter(
        filtered_df,
        x=x_axis,
        y=y_axis,
        color="t+ class",
        color_discrete_map=COLOR_DISCRETE_MAP,
        labels=performance_metrics_full,
        custom_data=CUSTOM_DATA_COLS,
        # hover_name="Sample ID"
    )

    # 配置悬停信息
    fig.update_traces(
        hovertemplate=HOVER_TEMPLATE.format(
            x_label=performance_metrics_full[x_axis],
            y_label=performance_metrics_full[y_axis]
        )
    )
    
    # 优化图表布局
    fig.update_layout(
        hovermode="closest",
        plot_bgcolor="rgba(245,245,245,1)",
        legend=dict(
            title="Li<sup>+</sup> transference number",
            orientation="h",
            yanchor="bottom",
            y=1.05
        )
    )
    
    st.plotly_chart(fig)
