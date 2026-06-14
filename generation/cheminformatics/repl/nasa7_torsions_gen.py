import os 
os.chdir('/home/yang/projects/t99_calc/data/generation/cheminformatics')
import inspect
import pandas as pd
import numpy as np 
import os
from pathlib import Path 
from rdkit import Chem 
from rdkit.Chem import MolFromSmiles 
from modules.data import bridge
from modules.data import sorting
#|%%--%%| <Ro7OW3k33G|a0RkJocbnJ>
os.chdir('/home/yang/projects/t99_calc/data/')
base = Path.cwd()

base = Path.cwd()

# parent directories 
qchem_data = base / "./storage"
csv = qchem_data / "./csv"
sdf = qchem_data / "./sdf"
combined = qchem_data / "./combined" 

# child directories 
qchem_out = combined / "./log_fchk"
smiles_out = combined / "./log_fchk_to_smiles"

# files 
name_nasa7_202_csv = csv / "./nasa7_202_concat.csv"

""" mols.sdf is a source file for all 289 Mol objects produced through the function FileToMol() 
f2mol_289.csv (also a product of FileToMol()) lists the input files (.log and .fchk) in the column ["Input Files"], and the 
singluar .sdf files (which are not sourced here) that are stored in ./qchem_data/sdf 
"""
combined_289_sdf = combined / "./mols.sdf"
file_289_csv = csv / "./f2mol_289.csv"

# imported module classes and functions
BytesPDB = sorting.BytesPDB
MoleculeSorter = sorting.MoleculeSorter
csv_generator = bridge.csv_generator 

# functions
def merger(left, right):
    new = pd.merge(left, right, on="Molecule", how="inner") 
    csv_generator(new, "{left}_{right}_merged")
    return new 

load_mols_sdf = bridge.load_mols_sdf
mols_289 = load_mols_sdf(combined_289_sdf)

file_289_df = pd.read_csv(file_289_csv)
name_nasa7_202_df = pd.read_csv(name_nasa7_202_csv)

# Names, SMILES, and Mol objects for all 289 .log and .fchk 
file_mol_289_df = file_289_df.assign(mol=mols_289)
file_mol_289_df.keys()

smiles_289_dict = {
        "File": [],
        "SMILES": []
    }

for row_index, row in file_mol_289_df.iterrows():
    file_289 = row["Input File"]
    mol_289 = row["mol"]
    name_289 = Path(file_289).stem
    for smiles in smiles_out.iterdir():
        smile_289 = Path(smiles).stem
        if smile_289 == name_289:
            with open(smiles) as f:
                for line in f:
                    line = line.strip()
                    parts = line.split(maxsplit=1)
                    smiles = parts[0]
                    name = parts[1]
                    smiles_289_dict['File'].append(Path(name).stem) 
                    smiles_289_dict['SMILES'].append(smiles)

""" Names, SMILES, and Mol objects correlating to the 202 molecules with 
available NASA7 parameters: the three required arguments for 
BytesPDB() 
""" 

smiles_202_list = []
mol_202_list = [] 
name_202_list = []

for row_index, row in file_mol_289_df.iterrows():
    file_289 = row["Input File"]
    mol_289 = row["mol"]
    name_289 = Path(file_289).stem
    for row_index, row in name_nasa7_202_df.iterrows():
        name_202 = row["Molecule"]
        if name_289 == name_202: 
            name_202_list.append(name_202)
            mol_202_list.append(mol_289)
            
for smiles_name, smiles in zip(smiles_289_dict["File"], smiles_289_dict["SMILES"]):
    for name in name_202_list:
        if name == smiles_name:
            smiles_202_list.append(smiles)
            
# Functions to reverse asterik-containing labels, normalize/combine them. 
def rev_tor_label(label: str, sep='-') -> str:
    reverse_label = sep.join(label.split(sep)[::-1])
    return reverse_label

def normalize_columns(input_dict: dict):
    seen = set() 
    matches = {} 
    idxs = [] 
    idys = []
    for x, (key, val) in enumerate(input_dict.items()):
        if "*" in key:
            idxs.append(x)
    for i in idxs: 
        ikey = list(input_dict)[i] 
        istrip = ikey.lstrip("*").split()[0]
        irev = rev_tor_label(istrip)
       # print(irev)
        ival = input_dict[ikey]
        seen.add(ikey) 
        for j in idxs:
            jkey = list(input_dict)[j] 
            jstrip = jkey.lstrip("*").split()[0]
            jval = input_dict[jkey]
            if irev == jstrip and jkey not in seen:
                idxs.append(i) 
                idys.append(j)
                test = [k for k in zip(idxs, idys)]
                print(test)
#                print(ival)
#                print(type(jval))
                matches[f"{ikey}"] = ival + jval 
    return matches, test

# function to remove rows in a pandas dataframe based on a column identifier 
def drop_df_rows(input_df: pd.DataFrame, col_name: str) -> pd.DataFrame:
    index_list = []
    df = input_df.copy()
    for row_index, row in df.iterrows():
        if row[f"{col_name}"] == 0:
            index_list.append(row_index)
    dropped_df = df.drop(index=index_list)
    return dropped_df 

#|%%--%%| <a0RkJocbnJ|nFRoV1sbna>
bpdb_inst = BytesPDB(name=name_202_list, mol=mol_202_list, smiles=smiles_202_list, canon_tor=True)
ms_inst = MoleculeSorter(bpdb_inst)

# Order: 1) mol_sorted or yes_sn 2) forbidden or no_sn 3) no_rot or no_rot
torsions_202_dict_no_sn = ms_inst.analyze_all()[0]
torsions_202_dict_yes_sn = ms_inst.analyze_all()[1]
torsions_202_dict_no_rot = ms_inst.analyze_all()[2]
torsions_202_dict_custom = ms_inst.analyze_all()[3]

len(torsions_202_dict_no_rot)
len(torsions_202_dict_custom)

torsions_202_df_no_rot = pd.DataFrame(torsions_202_dict_no_rot)
torsions_202_df_yes_sn = pd.DataFrame(torsions_202_dict_yes_sn)
torsions_202_df_no_sn = pd.DataFrame(torsions_202_dict_no_sn)
torsions_202_df_custom = pd.DataFrame(torsions_202_dict_custom)

no_rot_202_T = torsions_202_df_no_rot.T
yes_sn_202_T = torsions_202_df_yes_sn.T
no_sn_202_T = torsions_202_df_no_sn.T
custom_202_T = torsions_202_df_custom.T

no_sn_202_T_df = pd.DataFrame(no_sn_202_T)
yes_sn_202_T_df = pd.DataFrame(yes_sn_202_T)
no_rot_202_T_df = pd.DataFrame(no_rot_202_T)
custom_202_T_df = pd.DataFrame(custom_202_T)

merged_no_sn_df = merger(no_sn_202_T_df, name_nasa7_202_df)
no_sn_df = merged_no_sn_df.fillna(0) # 126 rows 

merged_yes_sn_df = merger(yes_sn_202_T_df, name_nasa7_202_df)
yes_sn_df = merged_yes_sn_df.fillna(0) # 20 rows 

merged_no_rot_df = merger(no_rot_202_T_df, name_nasa7_202_df)
no_rot_df = merged_no_rot_df.fillna(0) # 55 rows 

merged_custom_df = merger(custom_202_T_df, name_nasa7_202_df)
custom_202_df = merged_custom_df.fillna(0) 
custom_202_df_dedup = custom_202_df.drop_duplicates(subset=['Molecule'])

len(custom_202_df.columns)
custom_202_df.columns
csv_generator(no_sn_df, "no_sn_202")
csv_generator(yes_sn_df, "yes_sn_202")
csv_generator(no_rot_df, "no_rot_202")
csv_generator(custom_202_df_dedup, "custom_202_df")

filtered_normalized_145_df = drop_df_rows(input_df=custom_202_df_dedup, col_name="Global: Torsional Axes")
csv_generator(filtered_normalized_145_df, "filtered_normalized_145_all")

#|%%--%%| <nFRoV1sbna|9DkIEAdTIa>
"""
    1. torsions normalized; xhyb-xhyb-xhyb-xhyb
    2. torsions unnormalized; xhyb-xhyb-xhyb-xhyb
    3. torsions normalized; x-xhyb-xhyb-x 
    4. torsions unnormalized; x-xhyb-xhyb-x 
"""

#|%%--%%| <9DkIEAdTIa|Kd2r7r9bev>
amalgam_list = []
reverse_list = []
for row in custom_202_df.columns:
    sep = "-"
    if "Torsional Axes" in row and "Global" not in row:
        strip = row.split()[0]
        rev = sep.join(strip.split(sep)[::-1])
        amalgam_list.append(strip)
        reverse_list.append(rev)

# Check to ensure that torsions have been normalized succesfully!
rev_check = []
for r1 in amalgam_list:
    for r2 in reverse_list:
        if r1 == r2: 
            rev_check.append(r2)

# No Biggie!! IF rev is simply the identical torsions i.e., Y-X-X-Y  
rev_check

