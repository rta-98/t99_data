from modules.app.services import *
from modules.app.display import * 
from modules.data.bridge import *
from typing import Optional, List 
from dataclasses import dataclass 
from rdkit import Chem
from rdkit.Chem import rdMolTransforms, rdMolDescriptors
from rdkit.Chem.rdchem import Mol
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
        analyzer_dict[f"{atom} Count"] = atom_count

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
        if Z > 0: 
            E = Chem.GetPeriodicTable().GetElementSymbol(Z)
        else: 
            E = f"query({atom.GetSmarts()})" 
        cnt[E] += 1
    return dict(cnt) 

def count_motif(mol): 
    # creates a map
    bo_map = {
        Chem.BondType.SINGLE: "bo 1.0", 
        Chem.BondType.DOUBLE: "bo 2.0", 
        Chem.BondType.TRIPLE: "bo 3.0",
        Chem.BondType.AROMATIC: "bo 1.5",
    }

    # creates first dict. (the inner-most, ultimately set to nest inside result)
    bonds_in_mol = {
        "bo 1.0": [],
        "bo 2.0": [],
        "bo 3.0" : [],
        "bo 1.5": [],
        "Unk": [],
    } 

    mol_h = Chem.AddHs(mol) # mol but with added H's
    for mol_bond in mol_h.GetBonds():

        # GetBondType() also produces: GetBeginAtom() and GetEndAtom() as options
        bond_type_mapped = bo_map.get(mol_bond.GetBondType(), "Unk") 

        # defining two points
        atom1 = mol_bond.GetBeginAtom()
        atom2 = mol_bond.GetEndAtom() 
        
        # convert two points to atomic numbers (Z)
        z1 = atom1.GetAtomicNum() 
        z2 = atom2.GetAtomicNum() 

        # convert atomic numbers (Z) to atomic symbols (E) 
        e1 = Chem.GetPeriodicTable().GetElementSymbol(z1) if z1 > 0 else atom1.GetSmarts() 
        e2 = Chem.GetPeriodicTable().GetElementSymbol(z2) if z2 > 0 else atom2.GetSmarts() 
        
        # append atomic symbols (E) to the bonds_in_mol dict.
        bonds_in_mol[bond_type_mapped].append({
            "Bond Pair ID": "-".join(sorted((e1, e2))),
        }) 
    
    # The nested dictionary result is initialized
    motif_result = {} 

    # Loop over the key (bond) and the value (bond_idx) in bonds_in_mol dict.  
    # first for loop simply takes the list associated with one of the 4 bondization 
    # values, and iterates over the total number of entries, i.e., 
    # for a given sp3 you get 
    # [Bond Pair ID: C-C, Bond Pair ID: C-C, Bond Pair ID: C-N], which has length 3.

    for bond_key, bond_list in bonds_in_mol.items():
        if bond_list:
            motif_result[f"Global: {bond_key}"] = len(bond_list) 
            pair_list = []
            pair_count_dict = {}
            for bonds_in_mol_dict in bond_list: 
                name = bonds_in_mol_dict["Bond Pair ID"] 
                pair_list.append(name)
                if name in pair_list:
                    pair_count_dict[name] = pair_count_dict.get(name, 0) + 1 # alternative to dict[key] += 1 
                else:
                    pair_count_dict[name] = 0

            for pair_name, pair_count in pair_count_dict.items():
                motif_result[f"Local: {bond_key} {pair_name}"] = pair_count

    return motif_result 

def count_dihedral(mol, template_key: str, template_val: str):
    # Template for rotatable bonds
    ROT_BONDS_SMARTS = Chem.MolFromSmarts(template_val)
    rot_matches = mol.GetSubstructMatches(ROT_BONDS_SMARTS)
    confs = mol.GetConformer() 
    traversed = set() 
    unique_rot_matches = []
    for j, k in rot_matches:  
        bond = (min(j, k), max(j, k)) # e.g., min(7, 3), max(7, 3) -> (3, 7) 
        if bond not in traversed:
            traversed.add(bond) # object of type set naturally removes duplicates
            unique_rot_matches.append(bond) # now append to list  

    torsions = []
    bonds = []
    torsion_counts = {} 
    num_torsions = 0
    for j, k in unique_rot_matches:
        atom_j = mol.GetAtomWithIdx(j)
        atom_k = mol.GetAtomWithIdx(k)
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
            atoms = [mol.GetAtomWithIdx(idx) for idx in (m, j, k, n)]
            atom_symbols = tuple(a.GetSymbol() for a in atoms) 

            # label being X-X-X-X, a key for torsion_counts dict. 
            label = "-".join(atom_symbols) 

            # if label/key exists, append increment
            torsion_counts[label] = torsion_counts.get(label, 0) + 1 

    # Combining superimposable molecular labels 
    stored_HCCC_labels = {} 
    for label, val in torsion_counts.items():
        if label == 'H-C-C-C' or label == 'C-C-C-H':
            stored_HCCC_labels['*H-C-C-C'] = stored_HCCC_labels.get('*H-C-C-C', 0) + val 

    stored_FCCC_labels = {} 
    for label, val in torsion_counts.items():
        if label == 'F-C-C-C' or label == 'C-C-C-F':
            stored_FCCC_labels['*F-C-C-C'] = stored_FCCC_labels.get('*F-C-C-C', 0) + val 

    torsion_counts.update(stored_FCCC_labels) 
    torsion_counts.update(stored_HCCC_labels) 

    torsions_result = {}
    for label, val in torsion_counts.items():
        torsions_result[f"Local: {label} {template_key}"] = val 

    if num_torsions != 0:
        torsions_result[f"Global: {template_key}"] = num_torsions

    return torsions_result

#|%%--%%| <gPpYHhvciX|C2CdsMCP6n>
mol_10 = mol_202_list[1:10]

atoms_list = []
motif_list = [] 

for mol in mol_10:
    atom_rd = count_atoms(mol)
    motif_rd = count_motif(mol)
    atoms_list.append(atom_rd)
    motif_list.append(motif_rd)

atoms_list
motif_list
