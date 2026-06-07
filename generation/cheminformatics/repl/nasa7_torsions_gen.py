import inspect
import pandas as pd
import os
from pathlib import Path 
from rdkit import Chem 
from rdkit.Chem import MolFromSmiles 
from modules.data import bridge
from modules.data import sorting

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
            

bpdb_inst = BytesPDB(name=name_202_list, mol=mol_202_list, smiles=smiles_202_list)
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
#merged_custom_df.drop_duplicates()
custom_202_df = merged_custom_df.fillna(0) 
custom_202_df_dedup = custom_202_df.drop_duplicates(subset=['Molecule'])

csv_generator(no_sn_df, "no_sn_202")
csv_generator(yes_sn_df, "yes_sn_202")
csv_generator(no_rot_df, "no_rot_202")
csv_generator(custom_202_df_dedup, "custom_202_df")

custom_202_df_dedup.columns
custom_202_df_dedup.filter(regex=r"H-C-C-C|C-C-C-H").head(20)
custom_202_df_dedup.filter(regex=r"F-C-C-C|C-C-C-F").head(20)

custom_202_df_dedup.filter(regex=r"SMILES").head(20)

no_rot_df.filter(regex=r"H-C-C-C|C-C-C-H").head(20)
no_rot_df.filter(regex=r"F-C-C-C|C-C-C-F").head(30)

#|%%--%%| <Ro7OW3k33G|0FmkgmdXkK>
custom_202_df_dedup.columns
index_list = []
for row_index, row in custom_202_df_dedup.iterrows():
    if row['Global Sum: Torsional Axes'] == 0:
        index_list.append(row_index)

custom_202_df_dedup.drop(index=index_list)

print(custom_202_df_dedup['Global Sum: Torsional Axes'].head(100).to_string(index=False))
#|%%--%%| <0FmkgmdXkK|daKd1H6cUe>
custom_202_df_dedup.filter(like="torsional")

rev_tor_label('H-C-C-C')
mirror_labels_templates = {'F-C-C-C': [],
                           'C-C-C-F': [],
                           'H-C-C-C': [], 
                           'C-C-C-H': []}
for key in mirror_labels_templates:
    print(key)

#|%%--%%| <daKd1H6cUe|mADerLloVt>
i = 0
for key in custom_202_df_dedup["forward"]:
    i += 1 
    print(key)
print(i)


    for label, val in key['F-C-C-C'].items():
        print(key)

#|%%--%%| <mADerLloVt|d0QQhZWxhU>

# columns that contain "Global" and "torsional" at once
custom_202_df_dedup.filter(regex=r"(?=.*Global)(?=.*Torsional)")

# columns that contain "Local" and "torsional" at once
torsions_labels_subset = custom_202_df_dedup.filter(regex=r"(?=.*Local)(?=.*Torsional)")
#custom_202_df.columns
#|%%--%%| <d0QQhZWxhU|hdJhochXjC>
def check_reverse(ref1, ref2, sep='-') -> bool:
    reverse_string = sep.join(ref1.split(sep)[::-1])
    if reverse_string == ref2:
        return None

def rev_tor_label(label: str, sep='-') -> str:
    reverse_label = sep.join(label.split(sep)[::-1])
    return reverse_label

#|%%--%%| <hdJhochXjC|P6YYY6LvII>
new_torsions_dict = {
        "forward" : {},
        "reverse" : {}, 
        "combined": {}, 
} 

for label_out, val_out in torsions_labels_subset.to_dict().items():
    reverse_label_out = rev_tor_label(label_out) 
    for label_in, val_in in torsions_labels_subset.to_dict().items():
        if label_in == reverse_label_out:
            new_torsions_dict["forward"][label_out] = val_out
            new_torsions_dict["reverse"][label_in] = val_in 
            new_torsions_dict["combined"][f"{label_in}*"] = val_in + val_out

#|%%--%%| <P6YYY6LvII|HIfi43pLkT>
            
rev_torsions_dict = {} 
for label_out, val_out in torsions_labels_subset.to_dict().items():
    reverse_label = rev_tor_label(label) 
    rev_torsions_dict[reverse_label] = val_out 

    for label_in, val_in in torsion_labels_subset.to_dict().items():
        if label_in == reverse_label_out:
            new_torsions_dict["forward"][label_out] = val_out
            new_torsions_dict["reverse"][label_in] = val_in 
            new_torsions_dict["combined"][f"{label_in}*"] = val_in + val_out
            

    


