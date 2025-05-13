import streamlit as st
def project_overview():
    st.title("Project Overview")
    st.markdown(""" 
    Welcome to PPG-LBSPE, a Streamlit-based web application for Generation, Prediction, and
    Visualization of polymer candidates for lithium battery solid polymer electrolytes (SPEs). This 
    tool accelerates the discovery of high-performance polymers by integrating machine learning
    models with interactive analysis interfaces. 
    """)
    st.markdown("---")

    st.markdown(""" 
    ##  1.Background & Motivation
    Solid polymer electrolytes (SPEs) are key materials for next-generation lithium batteries, offering 
    enhanced safety and mechanical flexibility compared to liquid electrolytes. However, identifying
    polymers with high ionic conductivity, thermal stability, and lithium-ion transference number
    remains a major challenge. 
    
    To address this, PPG-LBSPE combines predictive modeling and generative design to:
    - Generate novel polymer structures using a trained molecular generation model
    - Predict key SPE properties from polymer SMILES using machine learning models
    - Visualize 26k+ pre-generated samples to explore structure–property relationships
    """)

    st.markdown("""
    ## 2.Objective
    - **Accelerate Polymer Discovery:** Provide a streamlined interface that allows users to screen thousands of candidate polymers in minutes instead of weeks.
    - **Integrate Predictive & Generative Models:** Use in-house models to predict key SPE properties and generate promising new polymer structures automatically.
    - **Enable Interactive Visualization:** Empower users to explore large datasets using dynamic filters, search functions, and scatter plots to uncover insights.
    """)

    st.markdown("""
                ## 3.Key Features
                ### 3.1. Polymer Generation(Generate)
                - Choose how many polymers to generate (1–100)
                - Generate candidates using the trained generative model
                - Display up to 10 valid SMILES results in a table
                - Preview each generated structure as an image
                - Download generated SMILES in a CSV file
                
                ### 3.2. Polymer Prediction(Predict)
                - Accept multiple input modes:
                    - Direct SMILES input via text box
                    - CSV file upload with a “SMILES” column
                    - Molecule sketching using the embedded editor
                - Predict the following SPE properties:
                    - Log ionic conductivity (σ, log S/cm)
                    - Glass transition temperature ($T_g$, °C)
                    - Thermal decomposition temperature ($T_d$, °C)
                    - Li$^+$ transference number ($t^+$, classification, <0.5 or >0.5)
                - Results are presented with structure images and metric cards

                ### 3.3. Data Visualization(Visualize)
                - Load and explore a dataset of 26,227 pre-generated polymers
                - Filter samples by value ranges for conductivity, $T_g$, $T_d$, or $t^+$ class
                - Search by Sample ID or SMILES substring
                - View interactive scatter plots of any two properties
                    - Click on any point to view detailed structure and property values
                    - Color-coded by transference number class
                """)
    
    st.markdown("""
                ## 4. Requirements
                - **Data Analysis:** Pandas, NumPy
                - **Data Preparation:** RDKit, Mordred Molecular Descriptors
                - **Machine Learning:** Scikit-learn, XGBoost, PyTorch
                - **Frontend:** Streamlit
                """)

    st.markdown("""
                ## 5. Reference
                - **Generative Model:** https://openreview.net/forum?id=l4IHywGq6a
                - **Predictive Models Data:**
                    - https://pubs.acs.org/doi/10.1021/acscentsci.2c01123
                    - https://polymer.nims.go.jp/
                - **Frontend:** https://github.com/CubeStar1/ChemPredictor
                """)
    
    st.markdown("""

                ## 6. Contact
                If you have any questions or need more information, feel free to contact:
                - **Email:** ysh222@ciac.ac.cn, lyliu@ciac.ac.cn, hfli@ciac.ac.cn
                - **Github:** https://github.com/shyao222/PPG-LBSPE
                """)