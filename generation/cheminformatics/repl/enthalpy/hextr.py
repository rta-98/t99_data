from pathlib import Path
import pandas as pd 

step1_dir = Path('/home/yang/projects/62_ftab/KINETIC_NETWORK/STEP1')
step2_dir = Path('/home/yang/projects/62_ftab/KINETIC_NETWORK/STEP2')
step3_dir = Path('/home/yang/projects/62_ftab/KINETIC_NETWORK/STEP3')
step3a_dir = Path('/home/yang/projects/62_ftab/KINETIC_NETWORK/STEP3a')

step1_list = []
for file in step1_dir.iterdir():
    if file.is_file() and file.suffix == ".txt":
        print(file.stem)
        df = pd.read_csv(file, sep=r"\s+", engine="python")
        out = df['H'][3].copy()
        step1_list.append(out)

H_LIST = []
#|%%--%%| <g4g6ONeomd|lKwLnlIV0C>
total = sum(step1_list)
print(total)
#|%%--%%| <lKwLnlIV0C|67LiI7qZkI>
jens_data = Path('/home/yang/projects/62_ftab/KINETIC_NETWORK/data.csv')
jens_csv = pd.read_csv(jens_data, sep=r",", engine="python")
#|%%--%%| <67LiI7qZkI|5DdH4vsa2Q>
seen = set()
seen2 = set()
for row_index, row in jens_csv.iterrows():
    csv_name = row['file'].lower()
    if "step1" in csv_name:
        match_name = csv_name
    for file in step1_dir.iterdir():
        if file.suffix == ".txt":
            if file.stem.lower() in match_name and match_name not in seen and file.stem.lower() not in seen2:
                seen.add(match_name)
                seen2.add(file.stem.lower())
                df = pd.read_csv(file, sep=r"\s+", engine="python")
                out = df['H'][3].copy()
                H_LIST.append(out)
                print(out, match_name, file.stem.lower())

#|%%--%%| <5DdH4vsa2Q|88rCf4Ox12>
seen = set()
seen2 = set()
for row_index, row in jens_csv.iterrows():
    csv_name = row['file'].lower()
    if "step2" in csv_name:
        match_name = csv_name
    for file in step2_dir.iterdir():
        if file.suffix == ".txt":
            if file.stem.lower() in match_name and match_name not in seen and file.stem.lower() not in seen2:
                seen.add(match_name)
                seen2.add(file.stem.lower())
                df = pd.read_csv(file, sep=r"\s+", engine="python")
                out = df['H'][3].copy()
                H_LIST.append(out)
                print(out, match_name, file.stem.lower())

#|%%--%%| <88rCf4Ox12|m8xI9wkvUM>
seen = set()
seen2 = set()
for row_index, row in jens_csv.iterrows():
    csv_name = row['file'].lower()
    if "step3" in csv_name:
        match_name = csv_name
    for file in step3_dir.iterdir():
        if file.suffix == ".txt":
            if file.stem.lower() in match_name and match_name not in seen and file.stem.lower() not in seen2:
                seen.add(match_name)
                seen2.add(file.stem.lower())
                df = pd.read_csv(file, sep=r"\s+", engine="python")
                out = df['H'][3].copy()
                H_LIST.append(out)
                print(out, match_name, file.stem.lower())

#|%%--%%| <m8xI9wkvUM|JFdMA3ouiI>
seen = set()
seen2 = set()
for row_index, row in jens_csv.iterrows():
    csv_name = row['file'].lower()
    if "step3a" in csv_name:
        match_name = csv_name
    for file in step3a_dir.iterdir():
        if file.suffix == ".txt":
            if file.stem.lower() in match_name and match_name not in seen and file.stem.lower() not in seen2:
                seen.add(match_name)
                seen2.add(file.stem.lower())
                df = pd.read_csv(file, sep=r"\s+", engine="python")
                out = df['H'][3].copy()
                H_LIST.append(out)
                print(out, match_name, file.stem.lower())


#|%%--%%| <JFdMA3ouiI|45M3bnAap7>
H_LIST
len(H_LIST)
len(jens_csv)

new_jens_csv = jens_csv.assign(H_0K=H_LIST) 

    #csvdf = df.to_csv(filename, index=index) # include index positional argument for to_csv() 
new_jens_csv.to_csv("62_ftab_H_ext.csv")

#|%%--%%| <45M3bnAap7|DJqRFl4EGk>
print(-2689.205484)

hrt = 627.509
z = ((-2689.205484 + 2689.206428)) # thermal energies - thermal enthalpies
print((0.415081-z) * hrt) # proves thermal correction to enthalpy contains RT

#|%%--%%| <DJqRFl4EGk|eClUtIir04>
for file in step1_dir.iterdir():
    if file.name == "GS.txt":
        df = pd.read_csv(file, sep=r"\s+", engine="python")
        #out = df[['Enthalpy', 'H']].copy()
        out = df.copy()
        zed = out['Enthalpy'] - out['Zero_Point']         
        print(zed*hrt)
