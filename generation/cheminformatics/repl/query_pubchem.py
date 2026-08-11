import pubchempy as pcp 
from rdkit import Chem 
from rdkit.Chem.MolStandardize import rdMolStandardize
import os 
import pandas as pd
from pathlib import Path 
#|%%--%%| <ONX6bFahdF|QjdcodofiS>
os.chdir('/home/yang/projects/t99_calc/data/')
base = Path.cwd()

qchem_data = base / "./storage"
json = qchem_data / "./json/df1_fbd.json"
json_pc = qchem_data / "./json/df1_fbd_pc.json"

df = pd.read_json(json)
df_pc = pd.read_json(json_pc)

smiles = df["SMILES"].to_list()

results = [];
i = 0
for smile in smiles:
    i += 1
    try: 
        result = pcp.get_compounds(f"{smile}", "smiles")
        results.append(result)
        print(result)
    except Exception as e: 
        results.append("NONE")
        continue 

#|%%--%%| <QjdcodofiS|YLIzzhocY8>
weights = []
names = [] 
formulas = []
for cid_list in results:
    for cid in cid_list: 
        name = cid.iupac_name
        weight = cid.molecular_weight
        formula = cid.molecular_formula
        weights.append(weight)
        names.append(name)
        formulas.append(formula)

#|%%--%%| <YLIzzhocY8|shOTBbcpyR>
weights
names
formulas

#|%%--%%| <shOTBbcpyR|gvMWo9ZYia>
df["mw"] = weights
df["iupac"] = names
df["mol_formula"] = formulas 

#|%%--%%| <gvMWo9ZYia|QzWYjWhv8S>
df.to_json("./storage/json/df1_fbd_pc.json", orient="records")


#|%%--%%| <QzWYjWhv8S|q9RNbWb9ud>
class InternalValid: 
     @staticmethod
     def validator(non_canon): 
           custom_SO3 = (
               "SULFONIC_ACID\t"
               "[S:1]([O:2])([O:3])([O:4])[#6:5]>>"
               "[S:1](=[O:2])(=[O:3])([O:4])[#6:5]\n"
           )
           params = rdMolStandardize.CleanupParameters()
           norm_SO3 = rdMolStandardize.NormalizerFromData(custom_SO3,
           params)
           if not isinstance(non_canon, str):
               raise TypeError("SMILES must be a string")
           smi = non_canon.strip()
           mol = Chem.MolFromSmiles(smi) 
           if mol is None:
               raise ValueError(f"Bad SMILES: {smi}")
           mol = norm_SO3.normalize(mol)
           mol = rdMolStandardize.Cleanup(mol)
           return Chem.MolToSmiles(mol, canonical=True)

smiles = 'O=S(=O)(O)CCC(F)(F)C(F)(F)F'
canonsmiles = InternalValid.validator(smiles)
