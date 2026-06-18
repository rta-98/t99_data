import pandas as pd
import os
from pathlib import Path 
from rdkit import Chem 
from rdkit.Chem import MolFromSmiles 
from modules.bridge import *
from modules.sorting import *

# File paths ---------------------------------
base = Path.cwd() 

qdata = base / "./qchem_data" 
combined = qdata / "./combined" 
csv = base / "./qchem_data/csv" 

dup_290_csv = base / "./nasa7_290.csv" 

no_rot = csv / "./no_rot.csv" 
no_sn = csv / "./no_sn.csv" 
yes_sn = csv / "./yes_sn.csv" 

original_nasa7_csv = csv / "./nasa7_parms_final.csv"
rdkit_mol_obab_smi_csv = csv / "./merged_smi_221_mol_221.csv" 

unmatch_molecules_csv = csv / "./unmatch_molecules.csv" 

pfas_data_130_csv = csv / "./PFAS_data_130_strip.csv"

Molecule_unmatch_log_csv = csv / "./Molecule_unmatch_log.csv" 

#nasagen_fit_csv = csv / "./nasagen_fit_results.csv" 
nasagen_fit_csv = csv / "./nasagen_fit_130.csv" 

nasa7_242_parms_csv = csv / "./nasa7_242_parms_df.csv"

nasagen_112_clean_csv = csv / "./nasagen_plus_112_clean.csv"

nasa7_202_clean_csv = csv / "./nasa7_202_clean.csv"

pfas_data_130_personal_csv  = csv / "./PFAS_data_130_personal.csv"

sdf_dir = base / "./qchem_data/sdf" 
smi_dir = base / "./qchem_data/combined/log_fchk_to_smiles"

f2mol_csv = csv / "./f2mol_289.csv"
combined_sdf = base / "./qchem_data/combined/mols.sdf"

log_fchk_path = combined / "./log_fchk"

nasa7_parms_final_csv = qdata / "./nasa7_parms_final.csv" 

# Junk ---------------------------------
#nasagen_fit_csv = csv / "./nasagen_fit_results.csv" 
#|%%--%%| <YtzL7VqY3s|JSnmY8W2m6>
# Problem: There are duplictes (obviously) in the return dictionary from FileToMol(); solve this by creating a df from it 
# and dropping them
# Generate dict and store mols in RAM

f2mol_inst = FileToMol(log_fchk_path)
f2mol_dict = f2mol_inst.file_mol_list_gen()
#|%%--%%| <JSnmY8W2m6|G6Jg1VAZXV>
# Dictionary from instance; duplicates and NaNs dropped 
f2mol_df_org0 = pd.DataFrame(f2mol_dict)
f2mol_df_org1 = f2mol_df_org0.drop_duplicates(subset=['Output File'], keep="first")
f2mol_df = f2mol_df_org1.dropna(subset=['mol']) 
f2mol_df # 289 x 3 

# Extracting mol column from the dict
mols_289 = f2mol_df['mol']

# Saving these 289 mols to their own sdf
# save_mols_sdf(mols_289, "mols.sdf") # saves in home
#|%%--%%| <G6Jg1VAZXV|sFlDXBpJvF>
mols = load_mols_sdf(combined / "mols.sdf")
len(mols) # verified: 289 mols long from combined sdf file! 

f2mol_df_to_csv = f2mol_df.drop(columns=['mol'])
csv_generator(f2mol_df_to_csv, "f2mol_289")

#Problem solved! 

# some fchk files are likely corrupt; revist this later, for now use SDFtoMol
sdf2mol_inst = SDFtoMol(sdf_dir)
sdf2mol_dict = sdf2mol_inst.sdf_mol_dict_gen() 
sdf2mol_dict 

# Dictionary from instance; duplicates and NaNs dropped 
sdf2mol_df_org0 = pd.DataFrame(sdf2mol_dict)
sdf2mol_df_org1 = sdf2mol_df_org0.drop_duplicates(subset=['sdf'], keep="first")
sdf2mol_df = sdf2mol_df_org1.dropna(subset=['mol']) 
sdf2mol_df # 289 x 3 

# Extracting mol column from the dict
mols_289 = sdf2mol_df['mol']

# Saving these 289 mols to their own sdf
save_mols_sdf(mols_289, "mols.sdf") 
mols = load_mols_sdf(combined / "mols.sdf")
len(mols) # verified: 289 mols long from combined sdf file! 


sdf2mol_df_to_csv = sdf2mol_df.drop(columns=['mol'])
csv_generator(sdf2mol_df_to_csv, "sdf2mol_289")

# Junk ---------------------------------
#|%%--%%| <sFlDXBpJvF|RftDWjuiQ4>
# Creating a dataframe that houses all 290 smiles, mols, their file names ('Molecules'), nasa7 parameters, and torsional data 

# f2mol.csv contains the input and output files, which correlate to the mol objects that are loaded in the next line:
f2mol_df_org = pd.read_csv(f2mol_csv)

# Loading the mol objects from bridge.py into a new list, to be combined into the master df
f2mol_mol_list = load_mols_sdf(combined_sdf)

len(f2mol_mol_list)

f2mol_df = f2mol_df_org.assign(mol=f2mol_mol_list)
output_files = f2mol_df['Output File']

# Generating a smiles column 
file_name_smi = []
for sdf_file in output_files:
    sdf_file_name = Path(sdf_file).stem 
    #print(file_name)
    for smile_file in smi_dir.iterdir():
        smile_file_name = Path(smile_file).stem 
        if smile_file_name == sdf_file_name:
            with open(smile_file) as f:
                for line in f: 
                    line = line.strip() 
                    parts = line.split(maxsplit=1)
                    file = Path(parts[1]).stem
                    name = parts[1] 
                    smi = parts[0]
                    file_name_smi.append((file, name, smi))

len(file_name_smi) # 289

smiles_289 = [smiles[2] for smiles in file_name_smi]

f2mol_smile_df = f2mol_df.assign(SMILES=smiles_289)
f2mol_smile_df # smiles column added
f2mol_smile_df.keys()

# Generating a mol column based off nasa7_202_clean.csv 

# INSERT HERE ---------------------------------

name_smile_mol_202_dict = {
        "Molecule": [],
        "SMILES"  : [], 
        "mol"     : []
    }

i = 0
for name_202 in name_202_list: 
    for row_index, row in f2mol_smile_df.iterrows():
        path_289 = row["Output File"]
        name_289 = Path(path_289).stem
        smiles_289 = row["SMILES"]
        mol_289 = row["mol"]
        if name_289 == name_202:
            name_smile_mol_202_dict["Molecule"].append(name_289)
            name_smile_mol_202_dict["SMILES"].append(smiles_289)
            name_smile_mol_202_dict["mol"].append(mol_289)
            #print(name_289) 
            i += 1
            # i = 72: Success! 
         
len(name_smile_mol_202_dict['mol'])

name_smile_mol_202_df = pd.DataFrame(name_smile_mol_202_dict)
name_smile_mol_202_df.keys()
name_202_df = name_smile_mol_202_df.drop(columns={"mol", "SMILES"})
name_202_df # name_202_df.drop_duplicates produces 201 long 

nasa7_202_updated_names = pd.merge(name_202_df, nasa7_202_clean_df, on="Molecule", how="left", sort=True) 
print(nasa7_202_updated_names.head(100).to_string(index=False))
print(nasa7_202_clean_df['Molecule'].head(100).to_string(index=False))

nasa7_202_clean_df
name_202_df = pd.DataFrame(name_202_list, columns=["Molecule"])
name_202_df

nasa7_parms_final_df = pd.read_csv(original_nasa7_csv)
nasa7_parms_final_df.keys() 
nasa7_parms_final_df['Log_file']
print(nasa7_parms_final_df.head(100).to_string(index=False))

#Junk--------------------------------- 

#pfas_data_130_df = pd.read_csv(pfas_data_130_csv)
#
#nasa7_129_log = pfas_data_130_df.drop(columns={"Molecule"}) 
#nasa7_log_list = nasa7_129_log.T.values.tolist()
#stem_list_130 = []
#for i in nasa7_log_list[0]:
#    name = Path(i).stem
#    stem_list_130.append(name)
#
#stem_list_130
#
#name_202_to_swap = nasa7_202_clean_df["Molecule"].to_list()
#

#|%%--%%| <RftDWjuiQ4|c7gTYJMbO6>
files_130 = nasa7_parms_final_df['Log_file']
names_130 = []
for file in files_130:
    name = Path(file).stem 
    names_130.append(name)

pop_col_abv = nasa7_parms_final_df.pop('Abbreviation')
pop_col_log = nasa7_parms_final_df.pop('Log_file')
pop_col_smiles = nasa7_parms_final_df.pop('SMILES')

nasa7_parms_final_df.insert(0, 'Molecule', names_130)

nasa7_parms_final_df.keys()
nasa7_202_clean_df.keys()

nasa7_130 = nasa7_parms_final_df
nasa7_130
nasa7_72 = nasa7_202_clean_df.iloc[130:202, :]
nasa7_72

nasa7_202_concat = pd.concat([nasa7_130, nasa7_72], axis=0)
csv_generator(nasa7_202_concat, "nasa7_202_concat")

#|%%--%%| <c7gTYJMbO6|p6D3wDefs2>
# Now to run BytesPDB on the mol objects; extract dictonary values to lists to enter; utilize 
# bridge.py methods to flatten the nested dictionary 

"""
            name: list, # Usually derived from Path(<file>).stem and contained in "Molecules"  
            mol: list, # Mol objects stored in RAM  
            smiles: list): # Accepts non-canonical SMILES
            
             inst = BytesPDB(abbrv=abbrv_col, smiles=smiles_col)
             tmp = MoleculeSorter(inst)
             tmp.analyze_all()
"""

# Converting dictionary values to lists

name_smile_mol_202_dict.keys()
name_smile_mol_202_dict 

name_col = name_smile_mol_202_dict['Molecule']
smiles_col = name_smile_mol_202_dict['SMILES'] 
mol_col = name_smile_mol_202_dict['mol']

bpdb_inst = BytesPDB(name=name_col, mol=mol_col, smiles=smiles_col)
ms_inst = MoleculeSorter(bpdb_inst)

# Order: 1) mol_sorted or yes_sn 2) forbidden or no_sn 3) no_rot or no_rot
torsions_202_dict_no_sn = ms_inst.analyze_all()[0]
torsions_202_dict_yes_sn = ms_inst.analyze_all()[1]
torsions_202_dict_no_rot = ms_inst.analyze_all()[2]

torsions_202_df_no_rot = pd.DataFrame(torsions_202_dict_no_rot)
torsions_202_df_yes_sn = pd.DataFrame(torsions_202_dict_yes_sn)
torsions_202_df_no_sn = pd.DataFrame(torsions_202_dict_no_sn)

no_rot_202_T = torsions_202_df_no_rot.T
yes_sn_202_T = torsions_202_df_yes_sn.T
no_sn_202_T = torsions_202_df_no_sn.T

csv_generator(no_rot_202_T, "no_rot_202_T")
csv_generator(yes_sn_202_T, "yes_sn_202_T")
csv_generator(no_sn_202_T, "no_sn_202_T")

no_sn_202_T_df = pd.DataFrame(no_sn_202_T)
yes_sn_202_T_df = pd.DataFrame(yes_sn_202_T)
no_rot_202_T_df = pd.DataFrame(no_rot_202_T)

no_sn_202_T.keys()

# Merging these data sets with parameters 

nasa7_242_df = pd.read_csv(nasa7_242_parms_csv)

nasa7_242_df.keys()

def merger(left, right):
    new = pd.merge(left, right, on="Molecule", how="inner") 
    csv_generator(new, "{left}_{right}_merged")
    return new 

merged_no_sn_df = merger(no_sn_202_T_df, nasa7_202_concat)
no_sn_df = merged_no_sn_df.fillna(0) # 126 rows 

merged_yes_sn_df = merger(yes_sn_202_T_df, nasa7_202_concat)
yes_sn_df = merged_yes_sn_df.fillna(0) # 20 rows 

merged_no_rot_df = merger(no_rot_202_T_df, nasa7_202_concat)
no_rot_df = merged_no_rot_df.fillna(0) # 55 rows 

csv_generator(no_sn_df, "no_sn_202")
csv_generator(yes_sn_df, "yes_sn_202")
csv_generator(no_rot_df, "no_rot_202")


# df generation ---------------------------------
#nasa7_242_parms_df = pd.read_csv(nasa7_242_parms_csv, dtype={"big_id": "Int64"})
#|%%--%%| <p6D3wDefs2|Fte0BiOlyj>
# df generation ---------------------------------
#nasa7_242_parms_df = pd.read_csv(nasa7_242_parms_csv, dtype={"big_id": "Int64"})

arr1 = set() 
arr2 = set()

for sdf in sdf_dir.iterdir():
    file = sdf_dir / f"{sdf.stem}"
    arr1.add(sdf)
    print(file.stem)
    
len(arr1)

arr2 = set()

for smi in smi_dir.iterdir():
    file = smi_dir / f"{smi.stem}" 
    arr2.add(smi)

len(arr2)

file_name_smi = [] 
for smi in smi_dir.iterdir():
    with open(f"{smi}") as f:
        for line in f: 
            line = line.strip() 
            if not line:
                continue
            parts = line.split(maxsplit=1)
            file = Path(parts[1]).stem
            name = parts[1] 
            smi = parts[0]
            file_name_smi.append((file, name, smi))

nasa7_242_parms_df = pd.read_csv(nasa7_242_parms_csv)

nasa7_242_parms_df.keys() 

#listb = nasa7_242_parms_df['Molecule'].to_list()

listb = nasagen_112_clean_csv_df['Molecule'].to_list()

file_name_smi_df = pd.DataFrame(file_name_smi)

lista = file_name_smi_df.iloc[:, 0].to_list()

len(listb)
len(lista)

set(listb).issubset(set(lista)) # True! 

missing = list(set(listb) - set(lista))

missing # 0!  

# Generating a smiles column 

for sdf_file in output_files:
    sdf_file_name = Path(file).stem 
    #print(file_name)
    for smile_file in smiles_directory.iterdir():
        smile_file_name = Path(smile_file).stem 
        if smile_file_name == sdf_file_name:
            with open(smile_file) as f:

#no_rot = pd.read_csv(no_rot, dtype={"big_id": "Int64"}) 
#no_sn = pd.read_csv(no_sn, dtype={"big_id": "Int64"}) 
#yes_sn = pd.read_csv(yes_sn, dtype={"big_id": "Int64"}) 
#
#org_nasa_df = pd.read_csv(original_nasa7_csv, dtype={"big_id": "Int64"})
#org_nasa_df_clean = org_nasa_df.drop(columns=['Log_file', 'Veri. SMILES']) 
#rd_mol_df = pd.read_csv(rdkit_mol_obab_smi_csv, dtype={"big_id": "Int64"})
#
#unmatch_molecules_df = pd.read_csv(unmatch_molecules_csv, dtype={"big_id": "Int64"}) 
#
#pfas_data_130_df = pd.read_csv(pfas_data_130_csv, dtype={"big_id": "Int64"})
#
#Molecule_unmatch_log_df = pd.read_csv(Molecule_unmatch_log_csv, dtype={"big_id": "Int64"}) 
#
#nasagen_fit_df = pd.read_csv(nasagen_fit_csv, dtype={"big_id": "Int64"})
#
#nasagen_112_clean_df = pd.read_csv(nasagen_112_clean_csv, dtype={"big_id": "Int64"})
#
nasa7_202_clean_df = pd.read_csv(nasa7_202_clean_csv, dtype={"big_id": "Int64"}, index_col=0)
nasagen_112_clean_csv_df = pd.read_csv(nasagen_112_clean_csv)
nasagen_112_clean_csv_df 
# Junk ---------------------------------
#dup_290_df = pd.read_csv(dup_290_csv, dtype={"big_id": "Int64"}) 
#nasa7_242_parms_df = dedup_290_df


#|%%--%%| <Fte0BiOlyj|Onesi0p6GT>
# Df modification ---------------------------------


# Junk ---------------------------------
#''' Duplicate Removal ''' 
#df_diff_account = df_unmatch_smi_log_parm.drop_duplicates(subset=['Abbreviation'], keep='first') 
#
#''' 
#1. Renaming Log_file from org_nasa_df to Log Name; removing .log from file
#2. Canonicalizing the SMILES column from Tony's 130 PFAS data 
#''' 
#df_log_col = org_nasa_df["Log_file"].tolist()
#df_smile_col = org_nasa_df["Veri. SMILES"].tolist()
#
#new_list0 = [] 
#new_list1 = [] 
#
#for log, smile in zip(df_log_col, df_smile_col): 
#    fname = Path(log).stem
#    new_list0.append(fname) 
#
#    mol = Chem.MolFromSmiles(smile)
#    can_smile = Chem.MolToSmiles(mol
#    new_list1.append(can_smile) 
#
#mod_df = org_nasa_df.drop(columns=['Log_file', 'Veri. SMILES']) 
#mod_df.insert(int(1), "Log Name", new_list0)
#mod_df.insert(int(2), 'Veri. SMILES', new_list1) 
#
#org_nasa_df = mod_df

#nasa_Molecule = nasagen_fit_df['Molecule'].to_list()
#
#new_names = [] 
#
#for old_name in nasa_Molecule:
#    artifact = "_2000_" 
#    new_name = old_name.replace(artifact, "")  
#    new_names.append(new_name)
#   
#old_names = nasagen_fit_df.pop('Molecule') 
#nasagen_fit_df.insert(0, 'Molecule', new_names) 
#
#Log_Name_unmatch_strip = Log_Name_unmatch_df.drop(columns=['a0', 'a1', 'a2', 'a3', 'a4', 'H_f_0K', 'S(300K)'])
#
#nasagen_fit_df.rename(columns={"S": "S(300K)"}) 


#suffix = "_x" 
#unmatch_nasa7_112 = unmatch_nasa7_242.drop(columns=[c for c in unmatch_nasa7_242.columns if c.endswith(suffix)]) 
#unmatch_nasa7_112.rename(columns={c: c[:-2] for c in unmatch_nasa7_112.columns if c.endswith("_y")}, inplace=True)

#unmatch_nasa7_112.drop(columns=['Log Files', 'SMILES']) 
#
#col_molecule = unmatch_molecules_df.pop("Molecule")
#unmatch_molecules_df.insert(0, "Molecule", col_molecule)

#col_log_name = unmatch_molecules_df.pop("Log Name")
#unmatch_molecules_df.insert(0, "Log Name", col_log_name)
#
#i = 3
#key = left.columns[i]
#assert key == right.columns[i]
#pd.merge(left, right, on=key, how="inner")
#
#nasagen_fit_df.head(100)
#nasa7_242_parms_df.head(5)
#
#entropy_col_nasagen = nasagen_fit_df["S(300K)"].to_list() 
#entropy_col_242 = nasa7_242_parms_df["S(300K)"].to_list() 
#entropy_col_112 = unmatch_nasa7_112["S(300K)"].to_list() 
#
#ngen_rounded = [] 
#n7_rounded = [] 
#n7_112_rounded = []
#
#for nasagen in entropy_col_nasagen:
#    ngrnd = round(nasagen, 2)
#    ngen_rounded.append(ngrnd) 
#
#for nasa7 in entropy_col_242: 
#    n7rnd = round(nasa7, 2) 
#    n7_rounded.append(n7rnd) 

#for nasa7 in entropy_col_112: 
#    n7rnd = round(nasa7, 2) 
#    n7_112_rounded.append(n7rnd) 
#
##entropy_col_nasagen = nasagen_fit_df.pop("S(300K)") 
#nasagen_fit_df.insert(7,"S(300K)_rnd", ngen_rounded)
#
##entropy_col_242 = nasa7_242_parms_df.pop("S(300K)")
#nasa7_242_parms_df.insert(8,"S(300K)_rnd", n7_rounded) 
#
#col_112 = unmatch_nasa7_112.pop("S(300K)_rnd") 
#unmatch_nasa7_112.insert(8, "S(300K)_rnd", n7_112_rounded) 
#
#col = org_nasa_df.pop('Abbreviation') 
#org_nasa_df.insert(0, 'Molecule', col) 
#
#col_abbrv = org_nasa_df_clean.pop('Abbreviation')
#org_nasa_df_clean.insert(0, 'Molecule', col_abbrv) 
#
#entropy_unmatch_df["Molecule"] 
#
#nasagen_112_clean_df.drop(columns=['Unnamed: 0'])
#nasagen_112_clean_df.keys()
#|%%--%%| <Onesi0p6GT|sMVNkOfdM0>
# Merge 0 ---------------------------------
org_nasa_df.head(5)
nasa7_242_parms.head(5)

m2_Molecule = pd.merge(org_nasa_df, nasa7_242_parms, on="Molecule", how="right", indicator=True, sort=False) 
unmatch1_nasa7_242 = m2_Molecule.loc[m2_Molecule["_merge"] == "right_only"].drop(columns="_merge") 

unmatch1_nasa7_242.head(100)

# Junk ---------------------------------
# log_csv_df = pd.merge(csvdf, logfdf, on="Log Files", how="outer", sort=False)
#m_smi = pd.merge(org_nasa_df, rd_mol_df, on="Veri. SMILES", how="right", indicator=True, sort=False) 
#m_log = pd.merge(org_nasa_df, rd_mol_df, on="Log Name", how="right", indicator=True, sort=False) 
#df_unmatch_smi_parm = m_smi.loc[m_smi["_merge"] == "right_only"].drop(columns="_merge") # 97/221 
#df_unmatch_log_parm = m_log.loc[m_log["_merge"] == "right_only"].drop(columns="_merge") # 97/221 
#df_unmatch_smi_log_parm = pd.concat([df_unmatch_log_parm, df_unmatch_smi_parm], axis=0, ignore_index=True) 

#orgin_removed_df_m_log = pd.merge(org_nasa_df, rd_mol_df, on="Log Name", how="inner", sort=False) # 104/130 
#orgin_removed_df_m_smi = pd.merge(org_nasa_df, rd_mol_df, on="Veri. SMILES", how="inner", sort=False) # 127/130 
#
#m_smi = pd.merge(org_nasa_df, rd_mol_df, on="Veri. SMILES", how="right", indicator=True, sort=False) 
#df_unmatch_smi_log_parm = m_smi.loc[m_smi["_merge"] == "right_only"].drop(columns="_merge") # 97/221 

# m_entropy = pd.merge(pfas_data_130_df, nasa7_242_parms_df, on="Molecule", how="right", indicator=True, sort=False)
# unmatch_nasa7_242 = m_entropy.loc[m_entropy["_merge"] == "right_only"].drop(columns="_merge") 

#i = 0
#key = unmatch_nasa7_112.columns[i] 
#assert key == unmatch_molecules_df.columns[i]
#m_Molecule = pd.merge(unmatch_nasa7_112, unmatch_molecules_df, on=key, how='right', indicator=True, sort=False) 
#Molecule_unmatch = m_Molecule.loc[m_Molecule['_merge'] == "right_only"].drop(columns='_merge')

#i = 0
#key = unmatch_nasa7_112.columns[i] 
#assert key == Molecule_unmatch_log_df.columns[i]
#m_Log_Name = pd.merge(unmatch_nasa7_112, Molecule_unmatch_log_df, on=key, how='right', indicator=True, sort=False) 
#Log_Name_unmatch_df = m_Log_Name.loc[m_Log_Name['_merge'] == "right_only"].drop(columns='_merge')  
#|%%--%%| <sMVNkOfdM0|YjYSETzf3P>
# Correlating matches; merging nasagen fits, and the 60 remaining unmatches ---------------------------------
len(org_nasa_df)
len(nasagen_112_clean_df)
org_nasa_df.head(100)
nasagen_112_clean_df.head(100)

org_nasa_df_clean.keys()
nasagen_112_clean_df.keys()

nasa7_202_df = pd.concat([org_nasa_df_clean, nasagen_112_clean_df], axis=0, ignore_index=True) 

# Junk ---------------------------------
## 1. nasagen_fit_df + Log_Name_unmatch_strip merged on right (unmatched should be the df with log paths to see what remains) 
##    and key='Molecule'
#
#m1_Molecule = pd.merge(nasagen_fit_df, Log_Name_unmatch_strip, on="Molecule", how="right", indicator=True, sort=False) 
#Log_Name_unmatch_df1 = m1_Molecule.loc[m1_Molecule["_merge"] == "right_only"].drop(columns="_merge") 
#len(Log_Name_unmatch_df1)
#
#Log_Name_match_df1 = pd.merge(nasagen_fit_df, Log_Name_unmatch_strip, on="Molecule", how="inner", sort=False)
#
#m2_entropy = pd.merge(nasagen_fit_df, unmatch_nasa7_112, on="S(300K)_rnd", how="right", indicator=True, sort=False)
#entropy_unmatch_df = m2_entropy.loc[m2_entropy["_merge"] == "right_only"].drop(columns="_merge") 
#len(entropy_unmatch_df) # 0
#
#nasagen_unmatch_nasa7_112_merge_df = pd.merge(nasagen_fit_df, unmatch_nasa7_112, on="S(300K)_rnd", sort=False) 

|%%--%%| <YjYSETzf3P|6MLrIVQF6W>
# Conversion ---------------------------------
#df_unmatch_smi_log_parm.to_csv("df_unmatch_smi_log_parm.csv", index=False) 
#unmatch_molecules_df.to_csv("unmatch_molecules_df.csv", index=False) 
#unmatch_nasa7_112.to_csv("unmatch_nasa7_112.csv", index=False)
#Molecule_unmatch.to_csv("Molecule_unmatch.csv", index=False)
#Log_Name_unmatch_df.to_csv("Log_Name_unmatch.csv", index=False) 

Log_Name_match_df1.to_csv("./qchem_data/csv/Log_Name_match.csv", index=False)
nasa7_242_parms_df.to_csv("./qchem_data/csv/nasa7_242_parms_df.csv", index=False)

# Junk ---------------------------------
#|%%--%%| <6MLrIVQF6W|DIpRDtjC7d>
# Correlating unmatches; merge 
''' 
Merging unmatch_molecules_df + nasa7_242_parms_df on "Log Name" and "Molecule" 
Keep "Log Path", "SMILES" and "Veri. SMILES" in the final frame.

0. Cut molecules with known data from nasa7_242_parms_df; step 1 is complete, 
   "Molecule" column was used to make a perfect cut
1. Merge on "Log Name" first, using the same logic in # Merge cell. 
2. The unmatched molecules, then will be matched once more on "Molecule"  
'''

orgin_removed_df_m_log = pd.merge(org_nasa_df, rd_mol_df, on="Log Name", how="inner", sort=False) # 104/130 
orgin_removed_df_m_smi = pd.merge(org_nasa_df, rd_mol_df, on="Veri. SMILES", how="inner", sort=False) # 127/130 

m_smi = pd.merge(org_nasa_df, rd_mol_df, on="Veri. SMILES", how="right", indicator=True, sort=False) 
df_unmatch_smi_log_parm = m_smi.loc[m_smi["_merge"] == "right_only"].drop(columns="_merge") # 97/221 


#|%%--%%| <DIpRDtjC7d|0lRSAsXXvY>
# Merging PFAS_data_130_personal.csv with nasa7_202_clean.csv on the name column to provided a key to source .log files

name_log_130_df = pfas_130_df[['Molecule', 'Log Files']]
name_log_130_list = name_log_130_df.values.tolist()

name_130_list = [Path(mol[1]).stem for mol in name_log_130_list]
len(name_130_list)
len(name_72_list)
sorted(name_130_list)
sorted(name_72_list)
name_202_list = name_130_list + name_72_list
name_202_list
len(name_202_list) # 202 

#pfas_130_df = pd.read_csv(pfas_data_130_personal_csv)
#pfas_130_df.keys()
#
#name_log_130_df = pfas_130_df[['Molecule', 'Log Files']]
#
## In order to convert a pandas data frame to a list that is more than one column, one must first access the NumPy array using .values (or .to_numpy()) a
## name_log_130_list = name_log_130_df.to_list() # this will NOT work
#
## Now one can use .tolist(), NOT to_list()
#name_logs_130_list = name_log_130_df.values.tolist()
#
## molecules_202: copy and swap the first 130 indices with name_log_smile_130['Log Files']  
#molecules_202_list = molecules_202.to_list() 
#molecules_130_list = [] 
#for idx, i in enumerate(molecules_202_list):
#    if idx < 130:
#        molecules_130_list.append(i)
#    
#molecules_130_list     
#       
## also, can be done via list comprehension:
#molecules_130_list = [i for idx, i in enumerate(molecules_202_list) if idx < 130]
#len(molecules_130_list)
#
#name_130_arr = [Path(mol[1]).stem for mol in name_logs_130_list]
#
#name_logs_130_list
#
## This will be taken back to like 
#sorted_name_130_arr = sorted(name_130_arr)
#

#|%%--%%| <0lRSAsXXvY|6IguWlK26Q>
# Print ---------------------------------
print(unmatch_molecules_df.head(10).to_string(index=False))  
print(nasa7_242_parms_df.head(10).to_string(index=False))

print(pfas_data_130_df.head(10).to_string(index=False))

len(unmatch_nasa7_242) 
print(unmatch_nasa7_242.head(10).to_string(index=False))

nasagen_fit_df

# Junk ---------------------------------
#
#print(dedup_290_df.head(10).to_string(index=False))
#len(dedup_290_df)
#dup_290_df.keys()
#
#print(no_rot.head(10).to_string(index=False)) 
#no_rot.keys()
#
#org_nasa_df.keys()                 
#len(org_nasa_df)
#print(org_nasa_df.head(10).to_string(index=False))
#
#rd_mol_df.keys()
#len(rd_mol_df) # 221
#print(rd_mol_df.head(10).to_string(index=False))
#
#orgin_removed_df_m_log.keys()
#len(orgin_removed_df_m_log)
#print(orgin_removed_df_m_log.head(10).to_string(index=False))
#
#orgin_removed_df_m_smi.keys()
#len(orgin_removed_df_m_smi)
#print(orgin_removed_df_m_smi.head(10).to_string(index=False))
#
#df_unmatch_smi_log_parm.keys() 
#len(df_unmatch_smi_log_parm) 
#print(df_unmatch_smi_log_parm.head(10).to_string(index=False)) 
#
#m.keys() 
#len(m)
#print(m.head(10).to_string(index=False)) 
#
#df_diff_account # frd_903_cof_OH 


