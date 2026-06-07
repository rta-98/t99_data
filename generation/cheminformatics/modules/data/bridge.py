from ..app.services import *
from openbabel import openbabel as ob
from pathlib import Path
import pandas as pd
from typing import Optional, List 
from pathlib import Path 
import json 
import sqlite3 
import re 
import os
from typing import Optional, List 
from dataclasses import dataclass 
from rdkit import Chem
from rdkit.Chem import rdMolTransforms, rdMolDescriptors, rdchem, PandasTools 
from rdkit.Chem.rdchem import Mol
from pathlib import Path 
from itertools import zip_longest, product 
from collections import Counter 
import subprocess
from rdkit.Chem.MolStandardize import rdMolStandardize 

base = Path.cwd() 
qdata = base / "./qchem_data"
csv = qdata / "./csv" 
combined = qdata / "./combined"
log_path = qdata / "./log"
fchk_path = qdata / "./fchk" 
log_fchk_path = qdata / "./combined/log_fchk"
nasa7_202_clean_csv = csv / "nasa7_202_clean.csv" 
smiles_directory = qdata / "./file_to_smi"

class Flattener: 
    def __init__(self, dict_inst: Optional[dict] = None): 
        self.dict_inst = dict_inst 
        self.rows = [] 
        self.items = {}  
        self.final_df: pd.DataFrame = pd.DataFrame() 

    def amalgam(self): 
        for mol_id, payload in self.dict_inst.items():
            row = {"Molecule": mol_id, **self.flatten(payload, sep=" ")}
            self.rows.append(row) 

        final_df = pd.DataFrame(self.rows).fillna(0)
        cols = ["Molecule"] + sorted(c for c in final_df.columns if c != "Molecule")
        final_df = final_df[cols]
        final_df = final_df.rename(columns=self.canon_col) 
        final_df = final_df.rename(columns=self.rename_torsion_cols)
        final_df = final_df.rename(columns=self.rename_motif_cols) 
        final_df = final_df.T.groupby(final_df.columns, sort=False).sum().T
        final_df = final_df.T.groupby(final_df.columns, sort=False).first().T

        return final_df 

    def flatten(self, obj, parent_key="", sep="__"):
        self.items = {} 
        if isinstance(obj, dict):
            for k, v in obj.items():
                new_key = f"{parent_key}{sep}{k}" if parent_key else str(k) 
                self.items.update(self.flatten(v, new_key, sep=sep)) 
        elif isinstance(obj, (list, tuple)): 
            for i, v in enumerate(obj): 
                new_key = f"{parent_key}{sep}{i}" if parent_key else str(i) 
                self.items.update(self.flatten(v, new_key, sep=sep)) 
        else: 
            self.items[parent_key] = obj

        return self.items

    def canon_col(self, c):
        m = re.match(r"^(Motif\s+\w+\s+Pair Count)\s+([A-Za-z]+)-([A-Za-z]+)$", c)
        if not m:
            return c
        prefix, a, b = m.groups() 
        a, b = sorted([a, b])

        return f"{prefix} {a}-{b}"

    def rename_torsion_cols(self, c): 
        if "Torsion Counts" in c:
            tail = c.split()[-1]
            return f"{tail} Torsions" 

        return c 

    def rename_motif_cols(self, c):
        if c.startswith("Motif "):
            return c.replace("Motif ", "", 1) 
        return c 

class LogFchkToMol: 
    def __init__(self, 
                 source_dir: Path): 
        self.source_dir = source_dir
        self.file_mols_dict = {
                'Input File': [],
                'Output File': [],
                'mol': []
        }
        self.file_mols_df: pd.DataFrame = pd.DataFrame() 
    
    def file_mol_list_gen(self) -> dict: 
        sdf_out_dir = qdata / "./sdf"
        sdf_out_dir.mkdir(exist_ok=True)
        for file in self.source_dir.iterdir():
            sdf_file = sdf_out_dir / f"{file.stem}.sdf"
            fmt = "fch" if file.suffix.lower() == ".fchk" else "log"
            subprocess.run(
                ["obabel", f"-i{fmt}", str(file), "-osdf", "-O", str(sdf_file)],
                check=True
            )
            mol_obj = Chem.SDMolSupplier(str(sdf_file), removeHs=False)[0]
            self.file_mols_dict['Input File'].append(file) 
            self.file_mols_dict['Output File'].append(sdf_file)
            self.file_mols_dict['mol'].append(mol_obj)
        return self.file_mols_dict 

    def file_mol_df_gen(self):
        dict_inst = self.file_mols_dict
        file_mols_df = pd.DataFrame(dict_inst) 
        return self.file_mols_df

class SDFtoMol: 
    def __init__(self, 
                 source_dir: Path):
        self.source_dir = source_dir
        self.sdf_mol_dict = {
                'sdf': [],
                'mol': []
        }
        
    def sdf_mol_dict_gen(self) -> dict:
        for file in self.source_dir.iterdir():
            mol_obj = Chem.SDMolSupplier(str(file), removeHs=False)[0] 
            self.sdf_mol_dict['sdf'].append(file.stem)
            self.sdf_mol_dict['mol'].append(mol_obj) 
        return self.sdf_mol_dict 


class AppendToCSV: 
    def __init__(self, 
                 csv_path: Optional[Path] = None,
                 imp_df: Optional[pd.DataFrame] = None,
                 merge_key: Optional[str] = None):
        self.csv_path = csv_path
        self.imp_df = imp_df 
        self.merge_key = merge_key
        self.csv_df: pd.DataFrame = pd.DataFrame() 
        self.merged_df: pd.DataFrame = pd.DataFrame() 

    def converter(self):
        self.csv_df = pd.read_csv(self.csv_path) 
        df_left = df1.set_index('Log Files (Rel. Path)')
#        df_right = self.csv_df.set_index('Log Files (Rel. Path)') 
        self.merged_df = pd.merge(self.csv_df, self.imp_df, on=[f"{self.merge_key}"], how="left") 
        self.merged_df = self.csv_df.join(df_right, how="left", sort=False).reset_index() 
        return self.merged_df

# Mol objects generated in-situ saved to a single SDF file for reuse purposes 
# Identified by order 
def save_mols_sdf(mols, path): 
    writer = Chem.SDWriter(str(path)) 
    try: 
        for m in mols: 
            if m is not None:
                writer.write(m) 
    finally: 
        writer.close() 

def load_mols_sdf(path): 
    """Produces Mol objects from an .sdf
    Args:
        path (str): directory with the concatenated sdf
    Returns:
        list[Mol]: a list of mol objects; the order in which 
                    mol objects were fed in to save_mols_sdf(): 
                    is the order in which they are returned here. 
    out: Mol object
    """
    suppl = Chem.SDMolSupplier(str(path), sanitize=True, removeHs=False) 
    return [m for m in suppl if m is not None] 

def csv_generator(df, fname: str, index: Optional[bool]=False):
    filename = f"{fname}.csv" 
    csvdf = df.to_csv(filename, index=index) # include index positional argument for to_csv() 
    return csvdf

def check_reverse(ref1, ref2, sep='-') -> bool:
    reverse_string = sep.join(ref1.split(sep)[::-1])
    if reverse_string == ref2:
        return None
    
def rev_tor_label(label: str, sep='-') -> str:
    reverse_label = sep.join(label.split(sep)[::-1])
    return reverse_label

