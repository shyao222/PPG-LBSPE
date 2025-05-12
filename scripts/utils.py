# Standard Libraries
import glob
from typing import Optional, Tuple
import requests

# Data Manipulation Libraries
import pandas as pd
import numpy as np

# Streamlit Libraries
import streamlit as st

# Chemistry Libraries
from chembl_webresource_client.new_client import new_client as ch
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from mordred import Calculator, descriptors


from scripts.image_handling import convert_df



def display_molecule_in_dataframe_as_html(dataframe):
    df = dataframe
    for index, i in enumerate(df['SMILES']):
        Draw.MolToFile(Chem.MolFromSmiles(str(i)), f'molecule_images/{i}.png', size=(300, 300), fitImage=True, imageType='png')

    images = glob.glob('molecule_images/*.png')
    df['Structure'] = images
    html_df = convert_df(df)
    return html_df


def display_molecule(molecule):
    img = Draw.MolToImage(molecule, size=(1000, 1000), fitImage=True)
    st.image(img)
