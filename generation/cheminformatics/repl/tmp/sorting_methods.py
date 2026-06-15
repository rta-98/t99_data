import os 
os.chdir('/home/yang/projects/t99_calc/data/generation/cheminformatics')
from modules.app.services import *
from modules.app.display import * 
from modules.data.bridge import *
from typing import Optional, List 
from dataclasses import dataclass 
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import rdMolTransforms, rdMolDescriptors
from rdkit.Chem.rdchem import Mol, HybridizationType
from pathlib import Path 
from itertools import zip_longest, product 
from collections import Counter 
from collections.abc import Collection 
import inspect
import pandas as pd
import os
from pathlib import Path 
from rdkit import Chem 
from rdkit.Chem import MolFromSmiles 
from modules.data import bridge
from modules.data import sorting 
#|%%--%%| <C2CdsMCP6n|PgCmb6jmux>
os.chdir('/home/yang/projects/t99_calc/data/')
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

#|%%--%%| <PgCmb6jmux|rUh9h4qqeq>
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

#|%%--%%| <rUh9h4qqeq|gPpYHhvciX>
FORBIDDEN = frozenset({ "S", "N" }) # molecules to exclude  
TEMPLATES = {
        "Torsional Axes": '[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]',
        }

def has_atom(
        mol: Optional[Mol] = None, 
        smiles: Optional[str] = None, 
        forbidden: Collection[str] = FORBIDDEN) -> bool: 
    return any(atom.GetSymbol() in forbidden for atom in mol.GetAtoms()) 

def has_rot(
        mol: Optional[Mol] = None) -> bool: 
    try: 
        n_rot = rdMolDescriptors.CalcNumRotatableBonds(
            mol, rdMolDescriptors.NumRotatableBondsOptions.Strict
        ) 
    except AttributeError:
        n_rot = rdMolDescriptors.CalcNumRotatableBonds(mol, strict=True) 
    return n_rot > 0 

def calc_rot(
        mol: Optional[Mol] = None) -> int: 
    try: 
        n_rot = rdMolDescriptors.CalcNumRotatableBonds(
            mol, rdMolDescriptors.NumRotatableBondsOptions.Strict
        ) 
    except AttributeError:
        n_rot = rdMolDescriptors.CalcNumRotatableBonds(mol, strict=True) 
    return n_rot 

def analyzer(name, smiles, mol):
    analyzer_dict = {}
    analyzer_dict["Molecule"] = name
    analyzer_dict["SMILES"] = smiles
    stored_torsions = {
            "results" : [], 
            } 

    # parsing count_atoms() dict output; appending to analyzer dict
    num_atoms_dict = count_atoms(mol) 
    analyzer_dict["Global: Rot Bonds"] = calc_rot(mol) 
    for atom, atom_count in num_atoms_dict.items():
        analyzer_dict[f"Global {atom} Count"] = atom_count

    # parsing count_motif() dict output; appending to analyzer dict
    motif_dict = count_motif(mol) 
    for bond, bond_vals in motif_dict.items():
        analyzer_dict[bond] = bond_vals
   
    # passing the Mol object, and torsional templates key and value to count_dihedral()
    for key, val in TEMPLATES.items(): 
        torsions_dict = count_dihedral(template_key=key, template_val=val, mol=mol)
        stored_torsions["results"].append(torsions_dict)

    # parsing count_dihedral() dict output; appending to analyzer dict
    for results in stored_torsions["results"]:
        for tor, tor_vals in results.items():
            analyzer_dict[tor] = tor_vals

    return analyzer_dict

def count_atoms(mol): 
    cnt = Counter() 
    for atom in mol.GetAtoms():
        Z = atom.GetAtomicNum()
        hybrid = atom.GetHybridization()
        hybrid_str = f"{hybrid}".lower()
        if Z > 0: 
            E = Chem.GetPeriodicTable().GetElementSymbol(Z)
        else: 
            E = f"query({atom.GetSmarts()})" 
        cnt[f"{E}{hybrid_str}"] += 1

    return dict(cnt)  

def count_motif(mol): 
    motif_result = {} 
    bond_dict = {} 

    dash_bond_map = {
        Chem.BondType.SINGLE: "-", 
        Chem.BondType.DOUBLE: "=", 
        Chem.BondType.TRIPLE: "≡",
    }

    mol_h = Chem.AddHs(mol) # mol but with added H's
    for mol_bond in mol_h.GetBonds(): 
        bond_type = mol_bond.GetBondType()
        bond_type_str = f"Global: " + f"{bond_type}".capitalize() 
        dash_bond = dash_bond_map.get(mol_bond.GetBondType(), "Unk")

        atom1 = mol_bond.GetBeginAtom() 
        atom2 = mol_bond.GetEndAtom() 

        z1 = atom1.GetAtomicNum()
        z2 = atom1.GetAtomicNum()
        
        e1 = Chem.GetPeriodicTable().GetElementSymbol(z1) if z1 > 0 else atom1.GetSmarts()
        e2 = Chem.GetPeriodicTable().GetElementSymbol(z2) if z2 > 0 else atom1.GetSmarts()

        hyb1 = atom1.GetHybridization() 
        hyb2 = atom2.GetHybridization() 
        
        hyb1_str = f"{hyb1}".lower()
        hyb2_str = f"{hyb2}".lower()

        bond = f"{e1}{hyb1_str}{dash_bond}{e2}{hyb2_str}"
        bond_dict[bond] = bond_dict.get(bond, 0) + 1
        bond_dict[bond_type_str] = bond_dict.get(bond_type_str, 0) + 1
        
    return bond_dict

def canonical_torsion(tor_str: str):
    sep = "-"
    fwd = sep.join(tor_str.split(sep))
    rev = sep.join(tor_str.split(sep)[::-1])
    return min(fwd, rev) 

def count_dihedral(mol, template_key: str, template_val: str):
    # Template for rotatable bonds
    dash_bond_map = {
        Chem.BondType.SINGLE: "-", 
        Chem.BondType.DOUBLE: "=", 
        Chem.BondType.TRIPLE: "≡",
    }

    ROT_BONDS_SMARTS = Chem.MolFromSmarts(template_val)
    rot_matches = mol.GetSubstructMatches(ROT_BONDS_SMARTS)
    confs = mol.GetConformer() 
    traversed = set() 
    unique_rot_matches = []
    hybs = "central"
    canon_tor = True 
    for j, k in rot_matches:  
        bond = (min(j, k), max(j, k)) # e.g., min(7, 3), max(7, 3) -> (3, 7) 
        if bond not in traversed:
            traversed.add(bond) # object of type set naturally removes duplicates
            unique_rot_matches.append(bond) # now append to list  

    torsions = []
    bonds = []
    torsion_counts = {} 
    normal_dict = {} 
    num_torsions = 0
    for j, k in unique_rot_matches:
        atom_j = mol.GetAtomWithIdx(j)
        atom_k = mol.GetAtomWithIdx(k)
        hybrid_j = atom_j.GetHybridization()
        hybrid_k = atom_k.GetHybridization()
        hybrid_j_str = f"{hybrid_j}".lower()
        hybrid_k_str = f"{hybrid_k}".lower()

        # Neighbor atoms not including k (left side)
        i = [n.GetIdx() for n in atom_j.GetNeighbors() if n.GetIdx() != k and n.GetAtomicNum()]
        # Neighbor atoms not including j (right side)
        l = [n.GetIdx() for n in atom_k.GetNeighbors() if n.GetIdx() != j and n.GetAtomicNum()]
        num_tbond = 0 # counter for num. torsions around central bond

        for m, n in product(i, l): # All possible combinations via cartesian product 
            if m != n: # skips identical indices 
                num_tbond += 1 
        for m, n in product(i, l):
            if m == n:
                continue

            num_torsions += 1
            atom_idxs = [m, j, k, n]
            atoms = [mol.GetAtomWithIdx(idx) for idx in atom_idxs]
            atom_symbols = [a.GetSymbol() for a in atoms]
            atom_symbols_str = '-'.join(atom_symbols)
            hybrids = [atom.GetHybridization() for atom in atoms]
            hybrids_str = [f"{hybrid}".lower() for hybrid in hybrids]

            x = 0
            y = 2
            dash_bond_list = []
            for atom in atom_idxs: 
                pair = atom_idxs[x:y]
                if y < 4:
                    x += 1 
                    y += 1
                atom1 = int(pair[0])
                atom2 = int(pair[1])
                bond = mol.GetBondBetweenAtoms(atom1, atom2)
                dash_bond_rep = dash_bond_map.get(bond.GetBondType(), "Unk") 
                dash_bond_list.append(dash_bond_rep)

            # hybridization options 
            if hybs == "all":
                torsion_string = f"{atom_symbols[0]}{hybrids_str[0]}{dash_bond_list[0]}{atom_symbols[1]}{hybrids_str[1]}{dash_bond_list[1]}{atom_symbols[2]}{hybrids_str[2]}{dash_bond_list[2]}{atom_symbols[3]}{hybrids_str[3]}"
                print(torsion_string)
            elif hybs == "central":
                torsion_string = f"{atom_symbols[0]}{dash_bond_list[0]}{atom_symbols[1]}{hybrids_str[1]}{dash_bond_list[1]}{atom_symbols[2]}{hybrids_str[2]}{dash_bond_list[2]}{atom_symbols[3]}"
            elif hybs == "none":
                torsion_string = f"{atom_symbols[0]}{dash_bond_list[0]}{atom_symbols[1]}{dash_bond_list[1]}{atom_symbols[2]}{dash_bond_list[2]}{atom_symbols[3]}"

            # normalization options
            if canon_tor is True:
                canon_tor = canonical_torsion(tor_str=torsion_string) 
                torsions = canon_tor
            elif canon_tor is False:
                torsions = torsion_string 

            torsion_counts[torsion_string] = torsion_counts.get(torsion_string, 0) + 1 
    
    torsions_result = {}
    for label, val in torsion_counts.items():
        torsions_result[f"{label} {template_key}"] = val 
    if num_torsions != 0:
        torsions_result[f"Global: {template_key}"] = num_torsions

    return torsions_result 
        
#|%%--%%| <gPpYHhvciX|P3VTxeLmfm>
mol_10 = mol_202_list[0:2]
for mol in mol_10:
    smiles = Chem.MolToSmiles(mol)
    print(smiles)
"""---------------------------------"""
atoms_list = []
motif_list = [] 
"""---------------------------------"""
for mol in mol_10:
    atom_rd = count_atoms(mol)
    motif_rd = count_motif(mol)
    dihedral_rd = count_dihedral(mol=mol, template_key=key, template_val=val)

    atoms_list.append(atom_rd)
    motif_list.append(motif_rd)
"""---------------------------------"""
motif_list
dihedral_rd
