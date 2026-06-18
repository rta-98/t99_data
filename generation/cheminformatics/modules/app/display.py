import os
import re 
from rdkit import Chem
from rdkit.Chem import Draw 
from openbabel import openbabel as obab 
from pathlib import Path
from rdkit.Chem.Draw import rdMolDraw2D 

class FileParser:
    def __init__(self, dir_path):
        self.dir = dir_path
        self.smiles = []
        self.mols = [] 
        self.pdb_fname = [] 

    def parse(self):
        for fname in os.listdir(self.dir): 
            if fname.endswith('.pdb'):
                fpath = os.path.join(self.dir, fname) 
                try: 
                    smi, mol_obj = self.pdb_to_smi(fpath)
                except Exception as e:
                    raise RuntimeError(f'Error in parse(): {e}')  
                parsed_smiles = self.smiles.append(smi) 
                self.mols.append(Chem.MolFromSmiles(smi)) 
                self.pdb_fname.append(fname) 

    def pdb_to_smi(self, fpath): 
        try:
            try:
                mol_obj = obab.OBMol() # parse mol from pdb  
                conv = obab.OBConversion() # conversion object 
                conv.SetInAndOutFormats('pdb', 'smi')
            except Exception as e:
                raise ValueError(f'Something is not right: {e}') 
            mol_pdb = conv.ReadFile(mol_obj, fpath) 
            if mol_pdb:
                smi_pdb = conv.WriteString(mol_obj).strip() 
                return smi_pdb.split()[0], mol_obj  
            else: 
                raise RuntimeError(f'Failed to read PDB {mol_obj}')
        except Exception as e: 
            raise RuntimeError(f'Error in pdb_to_smi(): {e}') 
        
    def smi_to_png(self, molsPerRow=5, subImgSize=(200,200)): 
        try: 
            grid_png = Draw.MolsToGridImage(
                self.mols,
                molsPerRow=molsPerRow,
                subImgSize=subImgSize,
                legends=self.smiles[:len(self.mols)] 
            )
            return grid_png 
        except Exception as e:
            raise RuntimeError(f'Error in smi_to_png(): {e}') 
    
