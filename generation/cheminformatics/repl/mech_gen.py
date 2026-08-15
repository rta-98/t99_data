import pandas as pd 
import os 
os.chdir('/home/yang/projects/t99_calc/data/generation/cheminformatics')
from rdkit import Chem
from pathlib import Path
from modules.data import services
#|%%--%%| <iCm50vFZo8|DIpZAkcYji>
org_df_path = Path("/home/yang/projects/t99_calc/data/storage/json/df1_fbd_pc.json")
class_info_csv_path = Path("/home/yang/projects/t99_calc/data/storage/csv/pfas_class_case_info.csv")
class_case_pdbs_path = Path("/home/yang/projects/t99_calc/data/storage/pdb/pfas_class_case_pdbs") 

smiles_list = [

        ]
org_df = pd.read_json(org_df_path)
class_info_df = pd.read_csv(class_info_csv_path) #, sep=r"\s+", engine="python")
class_info_df

seen = set()
for class_smi in class_info_df["SMILES"]:
    veri_class_smi = services.InternalValid.validator(class_smi);
    seen.add(veri_class_smi)
    for df_smi in org_df["SMILES"]:
        veri_df_smi = services.InternalValid.validator(df_smi);
        if veri_smi == veri_df_smi:
            print("y")
        

#|%%--%%| <DIpZAkcYji|rQ1aqybf2E>
pdb_path_list = []
for name in class_case_pdbs_path.iterdir():
    pdb_path_list.append(Path(name).name)

new_class_info = class_info_df.assign(pdbs=pdb_path_list)

new_class_info.to_csv("pfas_class_cases_info.csv")

new_class_info.to_json("pfas_class_cases_info.json")


