import pandas as pd
import numpy as np
import math
from rdkit import Chem
from rdkit.Chem.Lipinski import HeavyAtomCount
from rdkit.Chem import MolFromSmiles, MolToSmiles, AllChem
from rdkit.Chem import Descriptors
from mordred import Calculator, descriptors, Polarizability, HydrogenBond, RotatableBond, VdwVolumeABC, Weight, AcidBase, AtomCount, KappaShapeIndex
import pickle
import re
import itertools



def replace_K_with_Au_Cu(smiles):
    metals = itertools.cycle(['Au', 'Cu'])
    new_smiles = re.sub(r'\[K\]', lambda m: f'[{next(metals)}]', smiles)
    return new_smiles


def create_long_smiles(smile,req_length):
    #check if smiles is a polymer
    if 'Cu' in smile:
        #calculate required repeats so smiles > 30 atoms long
        num_heavy = HeavyAtomCount(MolFromSmiles(smile))-2
        repeats = math.ceil(req_length/num_heavy)-1
        
        #if polymer is less than 30 long, repeat until 30 long
        if repeats > 0:
            try:
                #code to increase length of monomer
                mol = MolFromSmiles(smile)
                new_mol = mol

                #join repeats number of monomers into polymer
                for i in range(repeats):
                    #join two polymers together at Cu and Au sites
                    rxn = AllChem.ReactionFromSmarts('[Cu][*:1].[*:2][Au]>>[*:1]-[*:2]')
                    results = rxn.RunReactants((mol,new_mol))
                    assert len(results) == 1 and len(results[0]) == 1, smile
                    new_mol = results[0][0]

                new_smile = MolToSmiles(new_mol)

            except:
                #make smile none if reaction fails
                return 'None'
                
        #if monomer already long enough use 1 monomer unit     
        else:
            new_smile = smile

        #caps ends of polymers with carbons
        new_smile = new_smile.replace('[Cu]','C').replace('[Au]','C').replace('[Ca]','C')

    else:
        new_smile = smile

    #make sure new smile in cannonical
    long_smile = MolToSmiles(MolFromSmiles(new_smile))
    return long_smile


def is_numeric(value):    
    if isinstance(value, (int, float, np.number)):
        return True
    else:
        return False


def conductivity_prediction(input_smiles):
    """
    This function takes a trained model and a SMILES string as input and returns the predicted conductivity of the molecule.
    """
    try:
        req_length = 30
        df = pd.DataFrame({'monomer_K': [input_smiles]})
        
        #除了SMILES之外的特征设置为数据集中的众数
        df['salt_smiles'] = 'O=S(=O)([N-]S(=O)(=O)C(F)(F)F)C(F)(F)F.[Li+]'
        df['mw'] = 5000000
        df['molality'] = 1.135074
        df['temperature'] = 60
        
        # Replace K with Au and Cu
        df['monomer_AuCu'] = df['monomer_K'].apply(replace_K_with_Au_Cu)
        
        #create long smiles
        df['smiles'] = create_long_smiles(df['monomer_AuCu'].iloc[0], req_length)
        
        #mordred for polymer
        mol_list_poly = [Chem.MolFromSmiles(smiles) for smiles in df['smiles']]
        MRD_descrip_poly = [Polarizability, HydrogenBond, RotatableBond, VdwVolumeABC, Weight, AcidBase, AtomCount, KappaShapeIndex]
        calc_poly = Calculator(MRD_descrip_poly)
        features_poly = calc_poly.pandas(mol_list_poly)
        newcolumns = {}
        for column in features_poly:
            newcolumns[column] = 'polymer ' + column
        poly_correct = features_poly.rename(columns = newcolumns)
        df_Mord_poly = pd.concat([df, poly_correct], axis=1)
    
        #mordred for salt
        mol_list_salt = [Chem.MolFromSmiles(smiles) for smiles in df['salt_smiles']]
        MRD_descrip_salt = [Polarizability, HydrogenBond, RotatableBond, Weight, AcidBase, AtomCount]
        calc_salt = Calculator(MRD_descrip_salt)
        features_salt = calc_salt.pandas(mol_list_salt)
        newcolumns = {}
        for column in features_salt:
            newcolumns[column] = 'salt ' + column
        salt_correct = features_salt.rename(columns = newcolumns)
        df_Mord_salt = pd.concat([df_Mord_poly, salt_correct], axis=1)
    
        #drop columns
        dropCol = ['smiles','salt_smiles','monomer_K','monomer_AuCu',
               'polymer Kier3',
               'salt RotRatio',
              ]
        dfml = df_Mord_salt.drop(columns = dropCol)

        #For some molecules the Vabc could not be calculated
        dfml['polymer Vabc'] = dfml['polymer Vabc'].apply(lambda x: x if is_numeric(x) else 0)
    
        #StandardScaler
        with open("utilities/models/conductivity_standardscaler.pkl", "rb") as f:
            conductivity_standardscaler = pickle.load(f)
        dfml_std = pd.DataFrame(conductivity_standardscaler.transform(dfml), columns=dfml.columns)
    
        #Normalizer
        with open("utilities/models/conductivity_normalizer.pkl", "rb") as f:
            conductivity_normalizer = pickle.load(f)
        dfml_norm = pd.DataFrame(conductivity_normalizer.transform(dfml_std), columns=dfml_std.columns)
    
        #features selection
        selected_features = ['mw', 'salt nC', 'salt nAcid', 'salt AMW', 'salt nHBAcc',
                         'salt bpol', 'salt apol', 'polymer nP', 'salt nN', 'polymer nN',
                         'polymer nS', 'polymer nBase', 'polymer AMW', 'polymer RotRatio',
                         'polymer nRot', 'polymer nHBDon', 'temperature', 'molality',
                         'polymer nB', 'salt nP']
        dfml_selected = dfml_norm[selected_features]
        
        # Load the trained model
        with open("utilities/models/conductivity_xgb.pkl", "rb") as f:
            conductivity_xgb = pickle.load(f)
    
        # Make the prediction
        prediction = conductivity_xgb.predict(dfml_selected)
    
        return prediction
    
    except Exception as e:
        print(f'Error predicting {input_smiles}: {e}')
        return None


def transference_number_prediction(input_smiles):
    try:
        req_length = 40
        df = pd.DataFrame({'monomer_K': [input_smiles]})
        
        #除了SMILES之外的特征设置为数据集中的众数
        df['salt_smiles'] = 'O=S(=O)([N-]S(=O)(=O)C(F)(F)F)C(F)(F)F.[Li+]'
        df['Polymer1_molecular_weight'] = 600000
        df['Electrolyte_Salt_molality'] = 1.135074
        df['T for Li+ transference number'] = 60
        df['Electrolyte_Filler_Amount (wt%)'] = 0
        df['Electrolyte_Plasticizer_Amount (wt%)'] = 0
        
        # Replace K with Au and Cu
        df['monomer_AuCu'] = df['monomer_K'].apply(replace_K_with_Au_Cu)
        
        #create long smiles
        df['smiles'] = create_long_smiles(df['monomer_AuCu'].iloc[0], req_length)
        
        #mordred for polymer
        mol_list_poly = [Chem.MolFromSmiles(smiles) for smiles in df['smiles']]
        MRD_descrip_poly = [Polarizability, HydrogenBond, RotatableBond, VdwVolumeABC, Weight, AcidBase, AtomCount, KappaShapeIndex]
        calc_poly = Calculator(MRD_descrip_poly)
        features_poly = calc_poly.pandas(mol_list_poly)
        newcolumns = {}
        for column in features_poly:
            newcolumns[column] = 'polymer ' + column
        poly_correct = features_poly.rename(columns = newcolumns)
        df_Mord_poly = pd.concat([df, poly_correct], axis=1)
    
        #mordred for salt
        mol_list_salt = [Chem.MolFromSmiles(smiles) for smiles in df['salt_smiles']]
        MRD_descrip_salt = [Polarizability, HydrogenBond, RotatableBond, Weight, AcidBase, AtomCount]
        calc_salt = Calculator(MRD_descrip_salt)
        features_salt = calc_salt.pandas(mol_list_salt)
        newcolumns = {}
        for column in features_salt:
            newcolumns[column] = 'salt ' + column
        salt_correct = features_salt.rename(columns = newcolumns)
        df_Mord_salt = pd.concat([df_Mord_poly, salt_correct], axis=1)
    
        #drop columns
        dropCol = ['smiles', 'salt_smiles', 'monomer_K','monomer_AuCu', 'polymer Vabc', 'polymer Kier3', 'salt RotRatio',]
        dfml = df_Mord_salt.drop(columns = dropCol)
    
        #StandardScaler
        with open("utilities/models/t+_standardscaler.pkl", "rb") as f:
            t_standardscaler = pickle.load(f)
        dfml_std = pd.DataFrame(t_standardscaler.transform(dfml), columns=dfml.columns)
    
        #Normalizer
        with open("utilities/models/t+_normalizer.pkl", "rb") as f:
            t_normalizer = pickle.load(f)
        dfml_norm = pd.DataFrame(t_normalizer.transform(dfml_std), columns=dfml_std.columns)
    
        #features selection
        selected_features = ['Polymer1_molecular_weight', 'Electrolyte_Salt_molality',
                             'T for Li+ transference number', 'Electrolyte_Filler_Amount (wt%)',
                             'Electrolyte_Plasticizer_Amount (wt%)', 'polymer apol',
                             'polymer nHBDon', 'polymer MW', 'salt apol', 'salt nHBAcc',
                             'salt nAcid']
        dfml_selected = dfml_norm[selected_features]
        
        # Load the trained model
        with open("utilities/models/t+_gbdt.pkl", "rb") as f:
            t_gbdt = pickle.load(f)
    
        # Make the prediction
        prediction = t_gbdt.predict(dfml_selected)
        
        return prediction[0]
    
    except Exception as e:
        print(f'Error predicting {input_smiles}: {e}')
        return None


def glass_transition_temperature_prediction(input_smiles):
    try:
        req_length = 30

        df = pd.DataFrame({'monomer_K': [input_smiles]})

        # Replace K with Au and Cu
        df['monomer_AuCu'] = df['monomer_K'].apply(replace_K_with_Au_Cu)
    
        #create long smiles
        df['smiles'] = create_long_smiles(df['monomer_AuCu'].iloc[0], req_length)
        
        #mordred for polymer
        mol_list_poly = [Chem.MolFromSmiles(smiles) for smiles in df['smiles']]
        MRD_descrip_poly = [Polarizability,HydrogenBond,RotatableBond,Weight,AcidBase,AtomCount,KappaShapeIndex,VdwVolumeABC]
        calc_poly = Calculator(MRD_descrip_poly)
        features_poly = calc_poly.pandas(mol_list_poly)
        newcolumns = {}
        for column in features_poly:
            newcolumns[column] = 'polymer ' + column
        poly_correct = features_poly.rename(columns = newcolumns)
        df_Mord_poly = pd.concat([df, poly_correct], axis=1)
    
        #drop columns
        dropCol = ['smiles','monomer_K','monomer_AuCu']
        dfml = df_Mord_poly.drop(columns = dropCol)

        #For some molecules the Vabc could not be calculated
        dfml['polymer Vabc'] = dfml['polymer Vabc'].apply(lambda x: x if is_numeric(x) else 0)

        #StandardScaler
        with open("utilities/models/tg_standardscaler.pkl", "rb") as f:
            tg_standardscaler = pickle.load(f)
        dfml_std = pd.DataFrame(tg_standardscaler.transform(dfml), columns=dfml.columns)
    
        #Normalizer
        with open("utilities/models/tg_normalizer.pkl", "rb") as f:
            tg_normalizer = pickle.load(f)
        dfml_norm = pd.DataFrame(tg_normalizer.transform(dfml_std), columns=dfml_std.columns)
    
        #features selection
        selected_features = ['polymer apol', 'polymer bpol', 'polymer nHBAcc', 'polymer nHBDon',
                            'polymer nRot', 'polymer AMW', 'polymer nBase', 'polymer nBridgehead',
                            'polymer nHetero', 'polymer nH', 'polymer nN', 'polymer nO',
                            'polymer nS', 'polymer nP', 'polymer nF', 'polymer nCl',
                            'polymer Kier2']
        dfml_selected = dfml_norm[selected_features]
        
        # Load the trained model
        with open("utilities/models/tg_xgb.pkl", "rb") as f:
            tg_xgb = pickle.load(f)
    
        # Make the prediction
        prediction = tg_xgb.predict(dfml_selected)
    
        return prediction
        
    except Exception as e:
        print(f'Error predicting {input_smiles}: {e}')
        return None


def thermal_decomposition_temperature_prediction(input_smiles):
    try:
        req_length = 30
        # 统一输入格式为列表

        df = pd.DataFrame({'monomer_K': [input_smiles]})
    
        # Replace K with Au and Cu
        df['monomer_AuCu'] = df['monomer_K'].apply(replace_K_with_Au_Cu)
    
        #create long smiles
        df['smiles'] = create_long_smiles(df['monomer_AuCu'].iloc[0], req_length)
        
        #mordred for polymer
        mol_list_poly = [Chem.MolFromSmiles(smiles) for smiles in df['smiles']]
        MRD_descrip_poly = [Polarizability,HydrogenBond,RotatableBond,Weight,AcidBase,AtomCount,KappaShapeIndex,VdwVolumeABC]
        calc_poly = Calculator(MRD_descrip_poly)
        features_poly = calc_poly.pandas(mol_list_poly)
        newcolumns = {}
        for column in features_poly:
            newcolumns[column] = 'polymer ' + column
        poly_correct = features_poly.rename(columns = newcolumns)
        df_Mord_poly = pd.concat([df, poly_correct], axis=1)
    
        #drop columns
        dropCol = ['smiles','monomer_K','monomer_AuCu',]
        dfml = df_Mord_poly.drop(columns = dropCol)
        
        #For some molecules the Vabc could not be calculated
        dfml['polymer Vabc'] = dfml['polymer Vabc'].apply(lambda x: x if is_numeric(x) else 0)
    
        #StandardScaler
        with open("utilities/models/td_standardscaler.pkl", "rb") as f:
            td_standardscaler = pickle.load(f)
        dfml_std = pd.DataFrame(td_standardscaler.transform(dfml), columns=dfml.columns)
    
        #Normalizer
        with open("utilities/models/td_normalizer.pkl", "rb") as f:
            td_normalizer = pickle.load(f)
        dfml_norm = pd.DataFrame(td_normalizer.transform(dfml_std), columns=dfml_std.columns)
    
        #features selection
        selected_features = ['polymer apol', 'polymer bpol', 'polymer nHBAcc', 'polymer nHBDon',
                             'polymer nRot', 'polymer MW', 'polymer AMW', 'polymer nBase',
                             'polymer nHeavyAtom', 'polymer nBridgehead', 'polymer nHetero',
                             'polymer nC', 'polymer nN', 'polymer nO', 'polymer nS', 'polymer nP',
                             'polymer nF', 'polymer nCl', 'polymer Kier1', 'polymer Kier2']
        dfml_selected = dfml_norm[selected_features]
        
        # Load the trained model
        with open("utilities/models/td_xgb.pkl", "rb") as f:
            td_xgb = pickle.load(f)
    
        # Make the prediction
        prediction = td_xgb.predict(dfml_selected)
    
        return prediction
        
    except Exception as e:
        print(f'Error predicting {input_smiles}: {e}')
        return None


def prediction_dataframe(input_smiles):
    df = pd.DataFrame({'SMILES': input_smiles})
    df['log conductivity'] = df['SMILES'].apply(conductivity_prediction)
    df['t+ class'] = df['SMILES'].apply(transference_number_prediction).replace({0: '<0.5', 1: '>0.5'})
    df['Tg'] = df['SMILES'].apply(glass_transition_temperature_prediction)
    df['Td'] = df['SMILES'].apply(thermal_decomposition_temperature_prediction)
    df['monomer Molecular Weight'] = df['SMILES'].apply(lambda x: Descriptors.MolWt(Chem.MolFromSmiles(x)))-39.098*2
    return df