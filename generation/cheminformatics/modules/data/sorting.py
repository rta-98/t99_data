from modules.data.services import *
from modules.app.display import * 
from modules.data.bridge import *
from typing import Optional, List 
from dataclasses import dataclass 
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem import rdMolTransforms, rdMolDescriptors
from rdkit.Chem.rdchem import Mol
from pathlib import Path 
from itertools import zip_longest, product 
from collections import Counter 
from collections.abc import Collection 
from typing import Literal 
import pubchempy as pcp 

FORBIDDEN = frozenset({ "S", "N" }) # molecules to exclude  
TEMPLATES = {"Torsional Axes": '[!$(*#*)&!D1]-,=&!@[!$(*#*)&!D1]'}

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

            #dash_bo: "all" for all atoms
                      "central" for the middle two 
        OUTPUT: flat dictionary of molecular info. """

    def __init__(self, 
                 name: list, # Usually derived from Path(<file>).stem and contained in "Molecules"  
                 smiles: list, # Accepts non-canonical SMILES
                 mol: list,
                 many_one: str,
                 hybs: str, 
                 img_path: Optional[Path] = None,
                 pdb_path: Optional[Path] = None,
                 ysn: Optional[bool] = False, 
                 canon_tor: Optional[bool] = True): # Mol objects stored in RAM  
        # --------------------------------- 
        self.name = name
        self.smiles = smiles
        self.mol = mol
        self.many_one = many_one
        self.hybs = hybs
        self.ysn = ysn 
        self.canon_tor = canon_tor
        self.img_path = img_path
        self.pdb_path = pdb_path 
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
        self.ysn: bool = molecule_sorter.ysn 
        self.canon_tor: bool = molecule_sorter.canon_tor 
        self.img_path: Path = molecule_sorter.img_path
        self.pdb_path: Path = molecule_sorter.pdb_path
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
            if self.ysn is True: 
                self.custom_dict[name] = self.analyzer(name, smiles, mol)
            elif self.ysn is False: 
                if not self.has_atom(mol):
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

    def append_png(self, name: str, smiles: Optional[str] = None) -> Path:
        """ Saves to a specified dir; then prints the path for the final df """ 
        mol_canon = Chem.MolFromSmiles(InternalValid.validator(smiles))
        SIZE = (200,200)
        BG_COLOR = (.29, .31, .33)
        img_path = f"{name}.png"
        mol_h = mol_canon
        AllChem.Compute2DCoords(mol_h)
        drawer = rdMolDraw2D.MolDraw2DCairo(200, 200)
        live_ops = drawer.drawOptions()
        live_ops.setBackgroundColour(BG_COLOR) 
        live_ops.bracketsAroundAtomLists = False 
        drawer.DrawMolecule(mol_h)
        drawer.FinishDrawing() 
        drawer.WriteDrawingText(self.img_path / img_path)
        return img_path

#    def get_pubdata(self, smiles: str):
#        try: 
#            result = pcp.get_compounds(f"{smiles}", "smiles") 
#        except Exception as e:
#            result == "Unk." 
#        return result 

    def append_pdb(self, name: str): 
        for file in self.pdb_path.iterdir(): 
            if Path(file).stem == name: 
                pdb_name = Path(file).name 
                return pdb_name

    def get_subst(self, smiles: Optional[str] = None) -> bool:
        mol = Chem.MolFromSmiles(InternalValid.validator(smiles))
        mol_h = Chem.AddHs(mol)
        matcher = SubstructMatch() 
        cat = matcher.classifyMotif(mol=mol_h) 
        return str(cat)

    def count_cf(self, smiles):
        mol = Chem.MolFromSmiles(InternalValid.validator(smiles))
        mol_h = Chem.AddHs(mol)
        matcher = SubstructMatch() 
        cf_dict = matcher.classifyTail(mol=mol_h) 
        return cf_dict
        
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
        
        # appending substructure (motif) match column to analyzer_dict
        analyzer_dict["Motif"] = self.get_subst(smiles)

        # appending substructure (tail) match column to analyzer_dict 
        cf_dict = self.count_cf(smiles)
        for k, v in cf_dict.items():
            analyzer_dict[k] = v 

        # appending pdb path match column to analyzer_dict 
        analyzer_dict["pdb"] = self.append_pdb(name) 

        # appending img paths column to analyzer_dict 
        analyzer_dict["img"] = self.append_png(name, smiles)
         
        # parsing count_atoms() dict output; appending to analyzer dict
        num_atoms_dict = self.count_atoms(mol) 
        analyzer_dict["Global: Rot Bonds"] = self.calc_rot(mol) 
        for atom, atom_count in num_atoms_dict.items():
            analyzer_dict[f"{atom} Count"] = atom_count

        # parsing count_bond() dict output; appending to analyzer dict
        motif_dict = self.count_bond(mol) 
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

    def count_bond(self, mol): 
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
            z2 = atom2.GetAtomicNum()
            
            e1 = Chem.GetPeriodicTable().GetElementSymbol(z1) if z1 > 0 else atom1.GetSmarts()
            e2 = Chem.GetPeriodicTable().GetElementSymbol(z2) if z2 > 0 else atom2.GetSmarts()

            hyb1 = atom1.GetHybridization() 
            hyb2 = atom2.GetHybridization() 
            
            hyb1_pair = f"{e1}{hyb1}".capitalize()
            hyb2_pair = f"{e2}{hyb2}".capitalize()

            a, b = sorted((hyb1_pair, hyb2_pair))

            bond = f"{a}{dash_bond}{b}"
            bond_dict[bond] = bond_dict.get(bond, 0) + 1
            bond_dict[bond_type_str] = bond_dict.get(bond_type_str, 0) + 1
            
        return bond_dict



    def canonical_torsion(self, tor_str: str):
        parts = re.findall(r"[^-=]+|[-=]", tor_str)
        fwd = "".join(parts)
        rev = "".join(parts[::-1])
        return min(fwd, rev) 

    def count_dihedral(self, mol, template_key: str, template_val: str):
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
                if self.hybs == "all": # dash_bo = "all"
                    torsion_string = f"{atom_symbols[0]}{hybrids_str[0]}{dash_bond_list[0]}{atom_symbols[1]}{hybrids_str[1]}{dash_bond_list[1]}{atom_symbols[2]}{hybrids_str[2]}{dash_bond_list[2]}{atom_symbols[3]}{hybrids_str[3]}"
                elif self.hybs == "central": # dash_bo = "central"
                    torsion_string = f"{atom_symbols[0]}-{atom_symbols[1]}{hybrids_str[1]}{dash_bond_list[1]}{atom_symbols[2]}{hybrids_str[2]}-{atom_symbols[3]}"
                elif self.hybs == "none": # dash_bo = "central"
                    torsion_string = f"{atom_symbols[0]}-{atom_symbols[1]}{dash_bond_list[1]}{atom_symbols[2]}-{atom_symbols[3]}"

                if self.canon_tor is True:
                    torsion_string = self.canonical_torsion(torsion_string)
                else:
                    torsion_string 
                torsion_counts[torsion_string] = torsion_counts.get(torsion_string, 0) + 1 
        
        torsions_result = {}
        for label, val in torsion_counts.items():
            torsions_result[f"{label} {template_key}"] = val 
        if num_torsions != 0:
            torsions_result[f"Global: {template_key}"] = num_torsions

        return torsions_result 
