import os 
import pandas as pd
from pathlib import Path 
#|%%--%%| <PEDqbsmODb|8OQnSJgwRK>
df1_fbd =  "/home/yang/exercises/js/img_gallery/static/storage/json/df1_fbd.json" 
df1 = pd.read_json(df1_fbd)

df1_nest_name = {}
mol_list = []

mol_names = df1["Molecule"].to_list()

#|%%--%%| <8OQnSJgwRK|Mvwlx462ZF>
df1_nest_name = {}
for name, (k, v) in zip(mol_names, df1.items()):
    df1_nest_name[name] = {k: v} 
df1 = pd.DataFrame(df1_nest_name)
df1.to_csv("df1_out.csv")

#|%%--%%| <Mvwlx462ZF|x9UXPt1UpN>
name = "52_heptene"
match = df1[df1["Molecule"] == f"{name}"]
nest = match.set_index("Molecule").to_dict(orient="index") 
nest[f"{name}"]["Motif"]
mol_names


#|%%--%%| <x9UXPt1UpN|XbIBSEtusu>
pdb = Path("/home/yang/projects/t99_calc/data/storage/pdb")
i = 0
for file in pdb.iterdir():
    for name in mol_names: 
        if Path(file).stem == name:
            print(file)
