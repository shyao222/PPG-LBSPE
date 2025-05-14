# Libraries
# ---------------------------------------------
# General
import os
import glob

# Data manipulation 
import pandas as pd
import plotly.express as px
from st_aggrid import AgGrid, GridOptionsBuilder

# Streamlit Libraries
import streamlit as st
from streamlit_plotly_events import plotly_events
import streamlit_antd_components as sac
from streamlit_ketcher import st_ketcher

# RDKit - To handle molecules
from rdkit import Chem
from rdkit.Chem import Draw

# Utilities
from scripts import project_overview_page
# from scripts.chat_page import chat_page_fn
from scripts import utils
from scripts.utils import display_molecule_in_dataframe_as_html
# from scripts.utils import get_smiles_from_name, get_iupac_from_smiles, display_molecule_in_dataframe_as_html
from scripts.predict_property import prediction_dataframe
from scripts.molecule_generation import generate


# Set page config
st.set_page_config(
    page_title="PPG-LBSPE",
    page_icon=":battery:",
    layout="wide",
    initial_sidebar_state="expanded")


# Loading the dataset
df_gp = pd.read_csv("utilities/dataset/gen_pred_final.csv")
df_gp['t+ class'] = df_gp['t+ class'].replace({0: '<0.5', 1: '>0.5'})
df_gp = df_gp.reset_index().rename(columns={'index': 'Sample ID', 'conductivity': 'Log Conductivity'})
filtered_df = pd.DataFrame(columns=df_gp.columns)  # initialization


# ===== SIDEBAR =====
with st.sidebar:
    # LOGO
    with st.container(border=True):
        #st.markdown('<h2 style="text-align: center;font-size: 1.5em;">SPE Prediction & Gerneration</h2>', unsafe_allow_html=True)
        st.image('utilities/logo/PPG-SPE.png')

    # with st.expander("Property Selection", expanded=False):
    #     property_selection = sac.chip(items=[
    #         sac.ChipItem(label='Log Conductivity (log S/cm)'),
    #         sac.ChipItem(label='Glass Transition Temperature (°C)'),
    #         sac.ChipItem(label='Thermal Decomposition Temperature (°C)'),
    #         sac.ChipItem(label='Li⁺ transference number classes'),
    #         ], 
    #         label='### Select Properties',
    #         index=[0,1,2,3],
    #         align='center',
    #         multiple=True, radius='md',
    #         variant='filled',
    #         )
    
    with st.expander("Input Selection", expanded=True):
        input_selection = sac.chip(items=[
            sac.ChipItem(label='SMILES input'),
            sac.ChipItem(label='Upload SMILES as file input'),
            sac.ChipItem(label='Draw molecule')],
                        label='Select Input Type',
                        index=[0], 
                        align='center',
                        multiple=False,
                        radius='md',
                        variant='filled')
    
    #Filter Samples for Visualization
    with st.expander("Filter Samples for Visualization", expanded=True):
        # Initialize session_state
        if "filtered_df" not in st.session_state:
            st.session_state.filtered_df = pd.DataFrame()
        # Continuous Variable Selection
        filters = {
            "Log Conductivity": st.slider(
                "Log Conductivity (log S/cm)",
                float(df_gp["Log Conductivity"].min()),
                float(df_gp["Log Conductivity"].max()),
                (float(df_gp["Log Conductivity"].min()), float(df_gp["Log Conductivity"].max()))
            ),
            "Tg": st.slider(
                "Glass Transition Temperature (°C)",
                float(df_gp["Tg"].min()),
                float(df_gp["Tg"].max()),
                (float(df_gp["Tg"].min()), float(df_gp["Tg"].max()))
            ),
            "Td": st.slider(
                "Thermal Decomposition Temperature (°C)",
                float(df_gp["Td"].min()),
                float(df_gp["Td"].max()),
                (float(df_gp["Td"].min()), float(df_gp["Td"].max()))
            )}

        # Categorical Variable Filtering（checkbox）
        st.markdown("Li⁺ transference number")
        col1, col2 = st.columns(2)
        with col1:
            select_less = st.checkbox("<0.5", value=True)
        with col2:
            # When the ">" symbol fails to display properly, replace it with HTML entity &gt;
            select_greater = st.checkbox("&gt;0.5", value=True) 

        t_class_selected = []
        if select_less:
            t_class_selected.append("<0.5")
        if select_greater:
            t_class_selected.append(">0.5")
        
        # If t_class_selected is empty (meaning neither checkbox is selected)
        # return an empty DataFrame directly.
        if not t_class_selected:
            filtered_df = pd.DataFrame(columns=df_gp.columns)
        else:
            # Apply filter conditions
            filter_conditions = ((df_gp["Log Conductivity"].between(*filters["Log Conductivity"])) &
                                (df_gp["Tg"].between(*filters["Tg"])) &
                                (df_gp["Td"].between(*filters["Td"])) &
                                (df_gp["t+ class"].isin(t_class_selected)))
            filtered_df = df_gp[filter_conditions]
        st.session_state.filtered_df = filtered_df
    st.divider()
    st.markdown(
    """
    <div style='font-size: 12px; color: gray;'>
    Copyright © 2009–2025 Changchun Institute of Applied Chemistry, Chinese Academy of Sciences<br>
    Address: 5625 Renmin Street, Changchun, Jilin Province, China 130022<br>
    Tel: 86-0431-85687300<br>
    ICP License: Ji ICP No. 12000082
    </div>
    """,
    unsafe_allow_html=True
    )

# ===== Main Interface =====
st.markdown("""
<style>
h1 {
    font-size: 26px !important;
}
</style>
""", unsafe_allow_html=True)
st.title("Prediction and Generation of Polymer Electrolytes for Solid-State Lithium Batteries")
st.divider()

def Generate():
    with st.expander("How to Make Generations", expanded=True):
        st.info("1. Use the slider or text box below to specify how many polymers to generate")
        st.info("2. Click the 'Start Generation' button")
    
    # Initialize session state
    if 'gen_num' not in st.session_state:
        st.session_state.gen_num = 1
    if 'text' not in st.session_state:
        st.session_state.text = str(st.session_state.gen_num)
    # Text Input Callback Function
    def update_from_text():
        try:
            val = int(st.session_state.text)
            if 1 <= val <= 100:
                st.session_state.gen_num = val
            else:
                st.error("Please enter an integer between 1-100")
                st.session_state.text = str(st.session_state.gen_num)
        except ValueError:
            st.error("Please enter a valid integer")
            st.session_state.text = str(st.session_state.gen_num)

    # Slider Update Synchronizes Text
    def update_from_slider():
        st.session_state.text = str(st.session_state.gen_num)
    col1, col2 = st.columns([1, 1])
    with col1:
        st.slider(
            'Select number of polymers to generate:',
            min_value=1,
            max_value=100,
            key='gen_num',
            on_change=update_from_slider
        )
    with col2:
        st.text_input(
            'Enter the number of polymers to generate (1–100):',
            key='text',
            on_change=update_from_text
        )
    if st.button('Start Generation', type='primary', use_container_width=True):
        with st.spinner("Generating..."):
            df_generated = generate(st.session_state.gen_num)

        if df_generated.empty:
            st.warning("No valid polymers were generated. Please try again or adjust your settings.")
            return
        # Set Display Count
        display_n = min(len(df_generated), 10)
        col_table, col_struct = st.columns([2,1])
        with col_table:
            st.write(f"Displaying {display_n} of {st.session_state.gen_num} generated molecules:")
            st.dataframe(df_generated.head(display_n))
        with col_struct:
            st.markdown("Molecular Structures:")
            images = []
            captions = []
            for i, smi in enumerate(df_generated['SMILES'].head(display_n), start=1):
                mol = Chem.MolFromSmiles(smi)
                if mol:
                    img = Draw.MolToImage(mol, size=(300,300), fitImage=True)
                    images.append(img)
                    captions.append(f"{i-1}: {smi}")
            if images:
                st.image(images, caption=captions, use_container_width=True)
            else:
                st.info("No valid molecular structures to display.")


def Predict():
    # Initialize the session state
    if "prediction_df_html" not in st.session_state:
        st.session_state.prediction_df_html = []
    if "prediction_df" not in st.session_state:
        st.session_state.prediction_df = []
    
    if input_selection == 'SMILES input':
        with st.expander("How to Make Predictions", expanded=True):
            st.info("1. Input SMILES of the polymer monomer in the box below")
            st.info("Note that use \'K\' to replace the polymer connection sites, For example: [K]CCO[K]")
            st.info("2. Click the button below to get the prediction")

        smiles_input = st.text_input('Please input polymer monomer SMILES in the box below: (Or you can choose other input type in the sidebar)',"")
        prediction = st.button('Predict property of SPE', use_container_width=True, type='primary')
        if prediction:
            with st.spinner("Predicting..."):
                files = glob.glob('molecule_images/*.png')
                for f in files:
                    os.remove(f)
                smiles_list = smiles_input.split(",")
                output_df = prediction_dataframe(smiles_list)
                # Convert each element in the array to a numerical value
                for column in output_df.columns:
                    try:
                        output_df[column] = output_df[column].apply(lambda x: x[0]).astype(float)
                    except:
                        continue
                st.session_state.prediction_df.append(output_df)
                html_df = display_molecule_in_dataframe_as_html(output_df)
                st.session_state.prediction_df_html.append(html_df)
        if st.session_state.prediction_df != []:
            molecular_weight = str(round(st.session_state.prediction_df[-1]['monomer Molecular Weight'].tolist()[0], 3)) + ' g/mol'
            predicted_cond = str(round(st.session_state.prediction_df[-1]['log conductivity'].tolist()[0], 3)) + ' log S/cm'
            predicted_t = str(st.session_state.prediction_df[-1]['t+ class'].tolist()[0])
            predicted_Tg = str(round(st.session_state.prediction_df[-1]['Tg'].tolist()[0], 3)) + ' °C'
            predicted_Td = str(round(st.session_state.prediction_df[-1]['Td'].tolist()[0], 3)) + ' °C'
            predicted_property_names = ['Log Conductivity',
                                        'Li⁺ transference number',
                                        'Glass Transition Temperature',
                                        'Thermal Decomposition Temperature']
            predicted_property_values = [predicted_cond, predicted_t, predicted_Tg, predicted_Td]

            with st.container(border=True):
                    st.markdown("## Predicted Properties")
                    structure, properties = st.columns([1,2])
                    with structure:
                        with st.container(border=True):
                            st.info('Structure')
                            images = glob.glob('molecule_images/*.png')
                            st.image(images)
                        with st.container(border=True):
                            st.info('monomer Molecular Weight')
                            st.markdown(f'<div id="" style="display: flex; justify-content: center; align-items: center; font-size: 20px; height:50px; ">{molecular_weight}</div>', unsafe_allow_html=True)
                    with properties:
                        with st.container(border=True):
                            col1, col2 = st.columns(2)
                            for i, col in enumerate([col1, col2]):
                                with col:
                                    for j in range(2):
                                        with st.container(border=True):
                                            st.markdown(f'<div id="" style=" height:50px; background-color: #ff4b4b; border-radius: 10px; display: flex; justify-content: center; align-items: center; font-weight: bold">{predicted_property_names[2*i+j]}</div>',
                                                        unsafe_allow_html=True)
                                            st.markdown(
                                                f'<div id="" style="display: flex; justify-content: center; align-items: center; font-size: 20px; height:100px; ">{predicted_property_values[2 * i + j]}</div>',
                                                unsafe_allow_html=True)

    elif input_selection == 'Upload SMILES as file input':
        with st.expander("How to Make Predictions", expanded=True):
            st.info("1. Upload SMILES of polymer monomers in CSV format")
            st.info("Note that SMILEs of the polymers should be in \'SMILES\' column, and use \'K\' to replace the polymer connection sites")
            st.info("2. Click the button below to get the prediction")
        many_SMILES = st.file_uploader('Please upload a CSV file:')
        prediction = st.button(f'Predict property of SPEs', type='primary', use_container_width=True)
        if prediction:
            with st.spinner("Predicting..."):
                files = glob.glob('molecule_images/*.png')
                for f in files:
                    os.remove(f)
                read_df = pd.read_csv(many_SMILES)
                read_smiles = read_df["SMILES"].tolist()
                output_df = prediction_dataframe(read_smiles)
                html_df = display_molecule_in_dataframe_as_html(output_df)
                st.markdown(html_df,unsafe_allow_html=True)

    elif input_selection == 'Draw molecule':
            st.session_state.molfile = Chem.MolToSmiles(Chem.MolFromSmiles("c1ccccc1"))
            files = glob.glob('molecule_images/*.png')
            for f in files:
                os.remove(f)
            with st.expander("How to Make Predictions", expanded=True):
                st.info("1. Draw your polymer monomer in the box below, then click the \'Apply\' button at the lower right corner of the box")
                st.info("Note that use \'K\' to replace the polymer connection sites")
                st.info("2. The SMILES of the structure you drew will be generated in the \'SMILES Representation\' box below")
                st.info("3. Click the predict button to get the predicted SPE properties of your polymer")
            st.markdown("Draw your polymer here:")
            smiles = st_ketcher(st.session_state.molfile)
            smile_code = smiles.replace('.',', ') # Replace the ". "in the SMILES with ", "
            molecule = st.text_input("SMILES Representation", smile_code)
            moleculesList = molecule.split(",")
            prediction2 = st.button('Predict property of SPE', type='primary', use_container_width=True)

            if prediction2:
                with st.spinner("Predicting..."):
                    output_df = prediction_dataframe(moleculesList)
                    # Convert each element in the array to a numerical value
                    for column in output_df.columns:
                        try:
                            output_df[column] = output_df[column].apply(lambda x: x[0]).astype(float)
                        except:
                            continue
                    st.session_state.prediction_df.append(output_df)
                    html_df = display_molecule_in_dataframe_as_html(output_df)
                    st.session_state.prediction_df_html.append(html_df)

            if st.session_state.prediction_df != []:
                molecular_weight = str(round(st.session_state.prediction_df[-1]['monomer Molecular Weight'].tolist()[0], 3)) + ' g/mol'
                predicted_cond = str(round(st.session_state.prediction_df[-1]['log conductivity'].tolist()[0], 3)) + ' log S/cm'
                predicted_t = str(st.session_state.prediction_df[-1]['t+ class'].tolist()[0])
                predicted_Tg = str(round(st.session_state.prediction_df[-1]['Tg'].tolist()[0], 3)) + ' °C'
                predicted_Td = str(round(st.session_state.prediction_df[-1]['Td'].tolist()[0], 3)) + ' °C'
                predicted_property_names = ['Log Conductivity',
                                            'Li⁺ transference number',
                                            'Glass Transition Temperature',
                                            'Thermal Decomposition Temperature']
                predicted_property_values = [predicted_cond, predicted_t, predicted_Tg, predicted_Td]

                with st.container(border=True):
                    st.markdown("## Predicted Properties")
                    structure, properties = st.columns([1,2])

                    with structure:
                        with st.container(border=True):
                            st.info('Structure')
                            images = glob.glob('molecule_images/*.png')
                            st.image(images)
                        with st.container(border=True):
                            st.info('monomer Molecular Weight')
                            st.markdown(f'<div id="" style="display: flex; justify-content: center; align-items: center; font-size: 20px; height:50px; ">{molecular_weight}</div>', unsafe_allow_html=True)
                    with properties:
                        with st.container(border=True):
                            col1, col2 = st.columns(2)
                            for i, col in enumerate([col1, col2]):
                                with col:
                                    for j in range(2):
                                        with st.container(border=True):
                                            st.markdown(f'<div id="" style=" height:50px; background-color: #ff4b4b; border-radius: 10px; display: flex; justify-content: center; align-items: center; font-weight: bold">{predicted_property_names[2 * i + j]}</div>',
                                                        unsafe_allow_html=True)
                                            st.markdown(
                                                f'<div id="" style="display: flex; justify-content: center; align-items: center; font-size: 20px; height:100px; ">{predicted_property_values[2 * i + j]}</div>',
                                                unsafe_allow_html=True)


def Visualize():
    filtered_df = st.session_state.filtered_df
    if filtered_df.empty:
        st.warning("No data found.")
    # Initialize session state
    if 'selected_index' not in st.session_state:
        st.session_state.selected_index = None
    
    st.info('These 26,227 samples were pre-generated using the generative model and can be used for polymer screening, structural analysis and performance visualization. You can use the \'Filter Samples for Visualization\' box in the sidebar to filter the samples.')
    col_table, col_struct = st.columns([5, 2])
    with col_table:
        # Search Bar (with Function Buttons)
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
        cols = st.columns([5, 1])
        with cols[0]:
            search_input = st.text_input("Search",
                                        value=st.session_state.get("search_term", ""),
                                        placeholder="Search by Sample ID or SMILES",
                                        label_visibility="collapsed")
        with cols[1]:
            search_clicked = st.button("🔍 Search", use_container_width=True)
        # Process Search Logic
        if search_clicked or (search_input and search_input != st.session_state.get("prev_search", "")):
            st.session_state.search_term = search_input
            st.session_state.prev_search = search_input
            st.rerun()
        # apply search condition
        if "search_term" in st.session_state and st.session_state.search_term:
            search_term = st.session_state.search_term
            try:
                # Perform exact matching on Sample ID
                sample_id = int(search_term)
                search_condition = (filtered_df["Sample ID"] == sample_id)
            except ValueError:
                # Fuzzy Matching SMILES
                search_condition = (filtered_df["SMILES"].str.contains(search_term, case=False, regex=False, na=False))
            filtered_df = filtered_df[search_condition]
        
        # Data Display (using AgGrid)
        # If filtered_df is empty, create an empty DataFrame with columns matching the original data
        if filtered_df.empty:
                filtered_df = pd.DataFrame(columns=df_gp.columns)
        
        # update Grid
        gb = GridOptionsBuilder.from_dataframe(filtered_df)
        gb.configure_column("Sample ID", width=75)
        gb.configure_column("SMILES", width=250)
        gb.configure_column("Log Conductivity", width=90)
        gb.configure_column("t+ class", width=60)
        gb.configure_column("Td", width=70)
        gb.configure_column("Tg", width=70)
        
        gb.configure_selection('single', use_checkbox=False)
        grid_options = gb.build()

        # Display Table
        grid_response = AgGrid(filtered_df,
                            gridOptions=grid_options,
                            height=500,
                            width='100%', 
                            fit_columns_on_grid_load=True,
                            theme='streamlit',
                            update_mode='SELECTION_CHANGED')
        st.markdown(f'<p style="text-align:right; color:#808080; font-size:0.9em;">{len(filtered_df)} samples shown</p>', 
                    unsafe_allow_html=True)
        
        # Get Selected Rows
        selected_rows = pd.DataFrame(grid_response['selected_rows'])
        if not selected_rows.empty:  # Correctly Determine if Empty
            st.session_state.selected_index = selected_rows.iloc[0]['Sample ID']

        
    # Molecular Structure Visualization
    with col_struct:
        selected_sample_id = st.session_state.get('selected_index', None)
        if selected_sample_id is not None:
            # Check if selected_sample_id Exists in Original Data
            if selected_sample_id not in df_gp["Sample ID"].values:
                st.warning(f"Sample ID {selected_sample_id} not found in original data.")
            else:
                # Get Selected Rows
                selected_row = df_gp[df_gp["Sample ID"] == selected_sample_id].iloc[0]
                smiles = selected_row['SMILES']
                # Generate Molecular Structure
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    img = Draw.MolToImage(mol, size=(400, 300))
                    st.image(img, caption=f"Sample ID: {selected_sample_id}")
                    st.write(f"**Polymer SMILES:** {smiles}")
                    st.write(f"**Log Conductivity:** {selected_row['Log Conductivity']:.4f} log S/cm")
                    st.write(f"**Glass Transition Temperature:** {selected_row['Tg']:.3f} °C")
                    st.write(f"**Thermal Decomposition Temperature:** {selected_row['Td']:.3f} °C")
                    st.write(f"**Li⁺ transference number:** {selected_row['t+ class']}")
                else:
                    st.warning("Invalid SMILES structure.")
        else:
            st.info("Select a sample to view structure.")

    # Performance Metrics Full Name Mapping
    performance_metrics_full = {
        "Log Conductivity": "Log Conductivity (log S/cm)",
        "Tg": "Glass Transition Temperature (°C)",
        "Td": "Thermal Decomposition Temperature (°C)"
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

    col1, col2 = st.columns(2)
    with col1:
        x_axis = st.selectbox(
            "Select X Axis for Scatter Plot:",
            options=performance_metrics,
            format_func=lambda x: performance_metrics_full[x],
            key="x_axis_selectbox"
        )
    with col2:
        y_axis = st.selectbox(
            "Select Y Axis for Scatter Plot:",
            options=[m for m in performance_metrics if m != x_axis],
            format_func=lambda x: performance_metrics_full[x],
            key="y_axis_selectbox"
        )
    if filtered_df.empty:
        st.warning("⚠️ No data available for visualization with current filters")
    else:
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
        # Configure Hover Information
        fig.update_traces(
            hovertemplate=HOVER_TEMPLATE.format(
                x_label=performance_metrics_full[x_axis],
                y_label=performance_metrics_full[y_axis]
            )
        )
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


def Project_Overview():
    st.title("Project Overview")
    st.write("This is the project overview page.")


def page_selection():
    selected = sac.segmented(items=[
        sac.SegmentedItem(label='Generate'),
        sac.SegmentedItem(label='Predict'),
        sac.SegmentedItem(label='Visualize'),
        sac.SegmentedItem(label='About')],
            format_func='title',
            align='center', 
            use_container_width=True)

    if selected == "About":
        project_overview_page.project_overview()
    if selected == "Predict":
        Predict()
    if selected == "Generate":
        Generate()
    if selected == "Visualize":
        Visualize()

page_selection()