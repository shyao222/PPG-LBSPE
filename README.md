# PPG-LBSPE: A Streamlit App for Polymer Electrolyte Discovery

**PPG-LBSPE** (Polymers Prediction and Generation for Lithium Battery Solid Polymer Electrolyte) is a Streamlit-based web application for **Polymer Generation, Property Prediction, and Data Visualization** aimed at accelerating the discovery of high-performance solid polymer electrolytes (SPEs) for lithium batteries.

## Project Overview

Solid polymer electrolytes offer safety and flexibility over liquid electrolytes but are hard to design due to complex structure–property relationships. This tool integrates predictive and generative machine learning models to streamline the design of novel polymer candidates.

## Key Features

### 1. Polymer Generation (Generate)
- Generate 1–100 polymer candidates using a pretrained generative model.
- Display valid SMILES in a table with molecular structure previews.
- Download generated molecules in CSV format.

### 2. Property Prediction (Predict)
- Input SMILES via text box, CSV upload, or molecular sketcher.
- Predict 4 SPE-relevant properties:
  - **Log ionic conductivity** (σ, log S/cm)
  - **Glass transition temperature** (T<sub>g</sub>, °C)
  - **Thermal decomposition temperature** (T<sub>d</sub>, °C)
  - **Li⁺ transference number** (t⁺ classification: <0.5 or >0.5)
- Visual results with structure images and metrics.

### 3. Data Visualization (Visualize)
- Explore 26,000+ pre-generated polymers.
- Filter samples by conductivity, T<sub>g</sub>, T<sub>d</sub>, or t⁺.
- Search by Sample ID or SMILES substring.
- Interactive scatter plots to investigate structure–property relationships.

## Installation

### 1. Clone this repository

```bash
git clone https://github.com/shyao222/PPG-LBSPE.git
cd PPG-LBSPE
````

### 2. Install required dependencies

```bash
pip install -r requirements.txt
```

### 3. Run the app

```bash
streamlit run speapp.py
```

## Hosted Version
For a quick demo, you can also access the hosted version of PPG-LBSPE at https://ppg-lbspe.streamlit.app/


## Requirements

* **Frontend:** Streamlit==1.41.1
* **Data Handling:** Pandas==2.0.3, NumPy==1.24.3
* **Chemical Descriptors:** RDKit==2023.09.4, Mordred==1.2.0
* **Modeling:** Scikit-learn==1.3.0, XGBoost==2.0.3, PyTorch==2.1.2+cpu

## References

* **Generative Model:** [OpenReview Paper](https://openreview.net/forum?id=l4IHywGq6a)
* **Predictive Model Datasets:**

  * [https://pubs.acs.org/doi/10.1021/acscentsci.2c01123](https://pubs.acs.org/doi/10.1021/acscentsci.2c01123)
  * [https://polymer.nims.go.jp/](https://polymer.nims.go.jp/)
* **Frontend Template:** [ChemPredictor](https://github.com/CubeStar1/ChemPredictor)

## Contact

If you have any questions or need more information, feel free to contact:

* **Email:** [ysh222@ciac.ac.cn](mailto:ysh222@ciac.ac.cn) | [lyliu@ciac.ac.cn](mailto:lyliu@ciac.ac.cn) | [hfli@ciac.ac.cn](mailto:hfli@ciac.ac.cn)
* **GitHub:** [shyao222/PPG-LBSPE](https://github.com/shyao222/PPG-LBSPE)
