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

# function to check if dicts have been normalized succesfully 
def normal_check(input_df: pd.DataFrame): 
    rev_check = []
    for row in input_df.columns:
        sep = "-"
        if "Torsional Axes" in row and "Global" not in row:
            strip = row.split()[0]
            rev = sep.join(strip.split(sep)[::-1])
            amalgam_list.append(strip)
            reverse_list.append(rev)

    # Check to ensure that torsions have been normalized succesfully!
    for r1 in amalgam_list:
        for r2 in reverse_list:
            if r1 == r2: 
                rev_check.append(r2)

    return rev_check 

#|%%--%%| <a0RkJocbnJ|9DkIEAdTIa>
""" "df" + numeric suffix 1-4; descriptions below:
        1. torsions normalized; xhyb-xhyb-xhyb-xhyb
        2. torsions unnormalized; xhyb-xhyb-xhyb-xhyb
        3. torsions normalized; x-xhyb-xhyb-x 
        4. torsions unnormalized; x-xhyb-xhyb-x 
    "drop" meaning 0 torsion mols. removed """
#1 ---------------------------------
bpdb_inst1 = BytesPDB(name=name_202_list, mol=mol_202_list, smiles=smiles_202_list, canon_tor=True, hybs="all", many_one="one")
ms_inst1 = MoleculeSorter(bpdb_inst1)
dict1 = ms_inst1.analyze_all()
df1 = pd.DataFrame(dict1).T
df1_merged = merger(df1, name_nasa7_202_df).fillna(0)
len(df1_merged.columns) 
df1_dedup = df1_merged.drop_duplicates(subset=["Molecule"])
df1_dedup_drop = drop_df_rows(input_df=df1_dedup, col_name="Global: Torsional Axes")
csv_generator(df1_dedup, "df1")
csv_generator(df1_dedup_drop, "df1_drop")

#2 ---------------------------------
bpdb_inst2 = BytesPDB(name=name_202_list, mol=mol_202_list, smiles=smiles_202_list, canon_tor=False, hybs="all", many_one="one")
ms_inst2 = MoleculeSorter(bpdb_inst2)
dict2 = ms_inst2.analyze_all()
df2 = pd.DataFrame(dict2).T
df2_merged = merger(df2, name_nasa7_202_df).fillna(0)
len(df2_merged.columns)
df2_dedup = df2_merged.drop_duplicates(subset=["Molecule"])
df2_dedup_drop = drop_df_rows(input_df=df2_dedup, col_name="Global: Torsional Axes")
csv_generator(df2_dedup, "df2")
csv_generator(df2_dedup_drop, "df2_drop")

#3 ---------------------------------
bpdb_inst3 = BytesPDB(name=name_202_list, mol=mol_202_list, smiles=smiles_202_list, canon_tor=True, hybs="central", many_one="one")
ms_inst3 = MoleculeSorter(bpdb_inst3)
dict3 = ms_inst3.analyze_all()
df3 = pd.DataFrame(dict3).T
df3_merged = merger(df3, name_nasa7_202_df).fillna(0)
len(df3_merged.columns)
df3_dedup = df3_merged.drop_duplicates(subset=["Molecule"])
df3_dedup_drop = drop_df_rows(input_df=df3_dedup, col_name="Global: Torsional Axes")
csv_generator(df3_dedup, "df3")
csv_generator(df3_dedup_drop, "df3_drop")

#4 ---------------------------------
bpdb_inst4 = BytesPDB(name=name_202_list, mol=mol_202_list, smiles=smiles_202_list, canon_tor=False, hybs="central", many_one="one")
ms_inst4 = MoleculeSorter(bpdb_inst4)
dict4 = ms_inst4.analyze_all()
df4 = pd.DataFrame(dict4).T
df4_merged = merger(df4, name_nasa7_202_df).fillna(0)
len(df4_merged.columns)
df4_dedup = df4_merged.drop_duplicates(subset=["Molecule"])
df4_dedup_drop = drop_df_rows(input_df=df4_dedup, col_name="Global: Torsional Axes")
csv_generator(df4_dedup, "df4")
csv_generator(df4_dedup_drop, "df4_drop")

