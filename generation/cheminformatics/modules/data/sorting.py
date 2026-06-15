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
from typing import Literal 

FORBIDDEN = frozenset({ "S", "N" }) # molecules to exclude  
TEMPLATES = {
        "Torsional Axes": '[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]',
        }

TEMPLATES = {
        "Torsional Axes": '[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]',
    }

class BytesPDB:
    """ INPUT: file name, mol objects, and SMILES as separate, equally sized lists. 
        FLAGS: 
            many_one: "many" or four different dictionaries sorted by 
                        conditions specified in analyze_all()
                      "one" or a single dictionary that isn't sorted. 
            hybs: "all" for all hybridization specified for atoms in torsion string. 
                  "central" for hybridization specified only for central atoms in torsion string. 
                  "none" for no hybridization specified for atoms in a torsion string. 
                        
            canon_tor: "True" for the normalization of torsion strings. 
                       "False" for the unnormalization of torsion strings.
        OUTPUT: flat dictionary of molecular info. """
    def __init__(self, 
                 name: list, # Usually derived from Path(<file>).stem and contained in "Molecules"  
                 smiles: list, # Accepts non-canonical SMILES
                 mol: list,
                 many_one: str,
                 hybs: str, 
                 canon_tor: Optional[bool] = True): # Mol objects stored in RAM  
                
        # --------------------------------- 
        self.name = name
        self.smiles = smiles
        self.mol = mol
        self.many_one = many_one
        self.hybs = hybs
        self.canon_tor = canon_tor
        # ---------------------------------
        self.smiles_names_mols_dict = {
                "Name": [],
                "SMILES": [],
                "Mol Object": []
        }

    def bookeeper(self) -> dict:
        for name, smiles, mol in zip(self.name, self.smiles, self.mol):
            self.smiles_names_mols_dict["Name"].append(name) 
            self.smiles_names_mols_dict["SMILES"].append(InternalValid.validator(smiles)) 
            self.smiles_names_mols_dict["Mol Object"].append(mol) 
        return self.smiles_names_mols_dict 

class MoleculeSorter: 
    def __init__(self, molecule_sorter: BytesPDB): 
        self.molecule_sorter: BytesPDB = molecule_sorter 
        self.imported_mol_data: dict = molecule_sorter.bookeeper() 
        self.many_one: str = molecule_sorter.many_one
        self.hybs: str = molecule_sorter.hybs
        self.canon_tor: bool = molecule_sorter.canon_tor 
        self.name = self.imported_mol_data["Name"]
        self.smiles = self.imported_mol_data["SMILES"] 
        self.mol = self.imported_mol_data["Mol Object"]
        self.mol_sorted_dict: dict = {} # rotatable bonds; S and N devoid
        self.forbidden_dict: dict = {} # rotatable bonds; S and N containing 
        self.non_rot_dict: dict = {} # non-rotatable bonds; S and N devoid
        self.custom_dict: dict = {} # 5/19: Viet wishes to combine non-rotatable bonds, with Sulfur & Nitrogen devoid dataset
        self.dud_list = [] # err 


    def analyze_all(self): 
        for i, (name, smiles, mol) in enumerate(
                zip(self.name, self.smiles, self.mol, strict=True)
        ):
            if not self.has_rot(mol) or not self.has_atom(mol):
                self.custom_dict[name] = self.analyzer(name, smiles, mol)
            if not self.has_rot(mol):
                self.non_rot_dict[name] = self.analyzer(name, smiles, mol) 
                continue 
            if not self.has_atom(mol):
                self.mol_sorted_dict[name] = self.analyzer(name, smiles, mol) 
            elif self.has_atom(mol):  
                self.forbidden_dict[name] = self.analyzer(name, smiles, mol)
        if self.many_one == "many":
            result = (self.mol_sorted_dict, self.forbidden_dict, self.non_rot_dict, self.custom_dict)
        if self.many_one == "one":
            result = self.custom_dict

        return result
    
    def has_atom(
            self,
            mol: Optional[Mol] = None, 
            smiles: Optional[str] = None, 
            forbidden: Collection[str] = FORBIDDEN) -> bool: 

        return any(atom.GetSymbol() in forbidden for atom in mol.GetAtoms()) 

    def has_rot(
            self, 
            mol: Optional[Mol] = None) -> bool: 
        try: 
            n_rot = rdMolDescriptors.CalcNumRotatableBonds(
                mol, rdMolDescriptors.NumRotatableBondsOptions.Strict
            ) 
        except AttributeError:
            n_rot = rdMolDescriptors.CalcNumRotatableBonds(mol, strict=True) 
        return n_rot > 0 

    def calc_rot(
            self, 
            mol: Optional[Mol] = None) -> int: 
        try: 
            n_rot = rdMolDescriptors.CalcNumRotatableBonds(
                mol, rdMolDescriptors.NumRotatableBondsOptions.Strict
            ) 
        except AttributeError:
            n_rot = rdMolDescriptors.CalcNumRotatableBonds(mol, strict=True) 
        return n_rot 

    def analyzer(self, name, smiles, mol):
        analyzer_dict = {}
        analyzer_dict["Molecule"] = name
        analyzer_dict["SMILES"] = smiles
        stored_torsions = {
                "results" : [], 
                } 

        # parsing count_atoms() dict output; appending to analyzer dict
        num_atoms_dict = self.count_atoms(mol) 
        analyzer_dict["Global: Rot Bonds"] = self.calc_rot(mol) 
        for atom, atom_count in num_atoms_dict.items():
            analyzer_dict[f"{atom} Count"] = atom_count

        # parsing count_motif() dict output; appending to analyzer dict
        motif_dict = self.count_motif(mol) 
        for bond, bond_vals in motif_dict.items():
            analyzer_dict[bond] = bond_vals
       
        # passing the Mol object, and torsional templates key and value to count_dihedral()
        for key, val in TEMPLATES.items(): 
            torsions_dict = self.count_dihedral(template_key=key, template_val=val, mol=mol)
            stored_torsions["results"].append(torsions_dict)

        # parsing count_dihedral() dict output; appending to analyzer dict
        for results in stored_torsions["results"]:
            for tor, tor_vals in results.items():
                analyzer_dict[tor] = tor_vals

        return analyzer_dict

    def count_atoms(self, mol): 
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

    def count_motif(self, mol): 
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
            
            # obtain hybridization of atoms in motif 
            hybrid1 = atom1.GetHybridization()
            hybrid2 = atom2.GetHybridization() 
            hyb_atom_dict = {e1: hybrid1, e2: hybrid2}
            a, b = sorted((e1, e2))
            hyb1_str = f"{hyb_atom_dict[a]}".lower()
            hyb2_str = f"{hyb_atom_dict[b]}".lower()

            # append atomic symbols (E) to the bonds_in_mol dict.
            bonds_in_mol[bond_type_mapped].append({
                "Bond Pair ID": f"{a}{hyb1_str}-{b}{hyb2_str}"
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
                    motif_result[f"{bond_key} {pair_name}"] = pair_count

        return motif_result 

    def canonical_torsion(self, tor_str: str):
        sep = "-"
        fwd = sep.join(tor_str.split(sep))
        rev = sep.join(tor_str.split(sep)[::-1])
        return min(fwd, rev) 

    def count_dihedral(self, mol, template_key: str, template_val: str):
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
                atoms = [mol.GetAtomWithIdx(idx) for idx in (m, j, k, n)]
                atom_symbols = [a.GetSymbol() for a in atoms]
                atom_symbols_str = '-'.join(atom_symbols)
                hybrids = [atom.GetHybridization() for atom in atoms]
                hybrids_str = [f"{hybrid}".lower() for hybrid in hybrids]

                # hybridization options 
                if self.hybs == "all":
                    torsion_string = f"{atom_symbols[0]}{hybrids_str[0]}-{atom_symbols[1]}{hybrids_str[1]}-{atom_symbols[2]}{hybrids_str[2]}-{atom_symbols[3]}{hybrids_str[3]}"
                elif self.hybs == "central":
                    torsion_string = f"{atom_symbols[0]}-{atom_symbols[1]}{hybrids_str[1]}-{atom_symbols[2]}{hybrids_str[2]}-{atom_symbols[3]}"
                elif self.hybs == "none":
                    torsion_string = f"{atom_symbols[0]}-{atom_symbols[1]}-{atom_symbols[2]}-{atom_symbols[3]}"

                # normalization options
                if self.canon_tor is True:
                    canon_tor = self.canonical_torsion(torsion_string) 
                    torsions = canon_tor
                elif self.canon_tor is False:
                    torsions = torsion_string 

                torsion_counts[torsions] = torsion_counts.get(torsions, 0) + 1 
        
        torsions_result = {}
        for label, val in torsion_counts.items():
            torsions_result[f"{label} {template_key}"] = val 
        if num_torsions != 0:
            torsions_result[f"Global: {template_key}"] = num_torsions

        return torsions_result 
