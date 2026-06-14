from fastapi import FastAPI, Form, Request, HTTPException
from fastapi.responses import HTMLResponse, RedirectResponse
from fastapi.templating import Jinja2Templates
from fastapi.staticfiles import StaticFiles 
from pydantic import BaseModel, ConfigDict, Field, field_validator 
from typing import Optional, List  
import io 
import csv
import os
import rdkit
import requests
import glob
import pandas as pd
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdchem, Draw, PandasTools
import subprocess
import re
import pubchempy as pubpy
from pathlib import Path
from openpyxl import Workbook
from openpyxl.drawing.image import Image as pyxl_img
from openpyxl.utils import get_column_letter as col_char
from PIL import Image as pil_img
from rdkit import DataStructs
from rdkit.Chem.MolStandardize import rdMolStandardize 
import json

class InputValid(BaseModel):

    smi_in: str = Field(..., min_length=1)

    @field_validator('smi_in')
    @classmethod 
    def veri_smi(cls, v):
        v = v.strip()
        if not v: 
            raise TypeError('SMILES are of type str, not None') 
        smi_2_mol = Chem.MolFromSmiles(v)
        if smi_2_mol is not None:
            cannon_smi = Chem.MolToSmiles(smi_2_mol) 
            return cannon_smi
        else: 
            raise ValueError(f'Invalid SMILES string: {v}')
            
class MolData:

    def __init__(self): 
        self.mol_data = []

    def img_2_bytes(self, smi_in):
        mol_obj = Chem.MolFromSmarts(smi_in) 
        if mol_obj is None: 
            raise ValueError(f'Failed to parse SMILES: {smi_in}')
        try: 
            mol_img = Draw.MolToImage(mol_obj, size=(100,100)) 
            img_bytes = io.BytesIO()
            mol_img.save(img_bytes, format="PNG")
            img_bytes.seek(0) # moves the file pointer to the beginning of the buffer 
            return img_bytes.getvalue()
        except Exception as e:
            raise RuntimeError(f'General Error: {e}')

    def get_iupac(self, smi_in) -> Optional[str]:
        try:
            smi_iupac = pubpy.get_compounds(smi_in, 'smiles')
            return smi_iupac[0].iupac_name
        except Exception as e:
            raise RuntimeError(f'Failed to determine IUPAC from SMILES: {e}') 

    def append_mol_data(self, smi_in: str) -> dict:
        try: 
            col_iupac = self.get_iupac(smi_in)
            col_img_bytes = self.img_2_bytes(smi_in)
            
            dict_mol_data = {
                "smi_in": smi_in,
                "iupac": col_iupac or "Unk",
                "img_bytes": col_img_bytes or "Unk"
            }

            self.mol_data.append(dict_mol_data)
            return dict_mol_data 
        except Exception as e:
            raise RuntimeError(f'Error in append_mol_data: {e}') 

    def clear_mol_data(self):
        self.mol_data = []

class InternalValid: 
     @staticmethod
     def validator(non_canon): 
           custom_SO3 = (
               "SULFONIC_ACID\t"
               "[S:1]([O:2])([O:3])([O:4])[#6:5]>>"
               "[S:1](=[O:2])(=[O:3])([O:4])[#6:5]\n"
           )
           params = rdMolStandardize.CleanupParameters()
           norm_SO3 = rdMolStandardize.NormalizerFromData(custom_SO3,
           params)
           if not isinstance(non_canon, str):
               raise TypeError("SMILES must be a string")
           smi = non_canon.strip()
           mol = Chem.MolFromSmiles(smi) 
           if mol is None:
               raise ValueError(f"Bad SMILES: {smi}")
           mol = norm_SO3.normalize(mol)
           mol = rdMolStandardize.Cleanup(mol)
           return Chem.MolToSmiles(mol, canonical=True)


class SubstMatch(InternalValid):

    sub_mol_pfeca = Chem.MolFromSmarts('[*]-[#8]-[#6](-[#6](=[#8])-[#8]-[#1])(-[#9])-[*]') 
    sub_mol_pfsa = Chem.MolFromSmarts('[#8]=[#16](-[#8]-[#1])(=[#8])-[#6](-[#9])(-[#9])-[*]')
    sub_mol_ftoh = Chem.MolFromSmarts('[#8](-[#6](-[#6](-[#1])(-[#1])-[*])(-[#1])-[#1])-[#1]')
    sub_mol_pfoh = Chem.MolFromSmarts('[#8](-[#6](-[#6](-[*])(-[#9])-[#9])(-[#9])-[#9])-[#1]')
    sub_mol_mefasaa = Chem.MolFromSmarts('[#8](-[#1])-[#6](=[#8])-[#6](-[#1])(-[#1])-[#7](-[#1])-[#16](=[#8])(=[#8])-[#6](-[#9])(-[#9])-[*]')
    sub_mol_ftca = Chem.MolFromSmarts('[#6](-[#1])(-[#1])(-[#6](-[#6](-[#8]-[#1])=[#8])(-[#1])-[#1])-[*]')
    sub_mol_fts = Chem.MolFromSmarts('[#6](-[#1])(-[#1])(-[#6](-[#16](=[#8])(-[#8]-[#1])=[#8])(-[#1])-[#1])-[*]') 
    sub_mol_pfca = Chem.MolFromSmarts('[#6](=[#8])(-[#8]-[#1])-[#6](-[*])(-[#9])-[#9]')
    sub_mol_pfasa = Chem.MolFromSmarts('[#7](-[#16](=[#8])(=[#8])-[#6](-[#9])(-[#9])-[*])(-[#1])-[#1]')
    sub_mol_fasa = Chem.MolFromSmarts('[#7](-[#16](=[#8])(=[#8])-[#6](-[#6](-[*])(-[#1])-[#1])(-[#1])-[#1])(-[#1])-[#1]')
    sub_mol_pfal = Chem.MolFromSmarts('[#6](=[#8])(-[#9])-[#6](-[#9])(-[#9])-[*]') 
 
#     sub_mol_pfeca = Chem.MolFromSmarts('[O]C(C(O)=O)F')
#     sub_mol_pfsa = Chem.MolFromSmarts('O=[S](O)=O')
#     sub_mol_ftoh = Chem.MolFromSmarts('OC[CH2]')
#     sub_mol_mefasaa = Chem.MolFromSmarts('CN(CC(O)=O)[S](=O)=O')
#     sub_mol_ftca = Chem.MolFromSmarts('[CH2]CC(O)=O')
#     sub_mol_fts = Chem.MolFromSmarts('[CH2]CS(=O)(O)=O')
#     sub_mol_pfca = Chem.MolFromSmarts('O=[C]O')
#     sub_mol_fasa = Chem.MolFromSmarts('N[S](=O)=O')

    if sub_mol_pfeca is None: 
        print("ERROR") 

    def __init__(self) -> None:
        self.pfeca: List[tuple] = []
        self.fasa: List[tuple] = []
        self.pfasa: List[tuple] = []
        self.pfca: List[tuple] = []
        self.fts: List[tuple] = []
        self.ftca: List[tuple] = []
        self.mefasaa: List[tuple] = []
        self.ftoh: List[tuple] = []
        self.pfoh: List[tuple] = []
        self.pfsa: List[tuple] = []
        self.pfal: List[tuple] = []
#         self.CF_chain: List[tuple] = []

    def match_pfeca(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_pfeca):
            matched_smiles = Chem.MolToSmiles(mol_in)
            self.pfeca.append((matched_smiles, mol_in))
            return True
        return False

    def match_fasa(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_fasa):
            matched_smiles = Chem.MolToSmiles(mol_in)
            self.fasa.append((matched_smiles, mol_in))
            return True
        return False

    def match_pfasa(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_pfasa):
            matched_smiles = Chem.MolToSmiles(mol_in)
            self.pfasa.append((matched_smiles, mol_in))
            return True
        return False

    def match_pfca(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_pfca):
            matched_smiles = Chem.MolToSmiles(mol_in)
            self.pfca.append((matched_smiles, mol_in))
            return True
        return False

    def match_fts(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_fts):
            matched_smiles = Chem.MolToSmiles(mol_in)
            self.fts.append((matched_smiles, mol_in))
            return True
        return False

    def match_ftca(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_ftca):
            matched_smiles = Chem.MolToSmiles(mol_in)
            self.ftca.append((matched_smiles, mol_in))
            return True
        return False

    def match_mefasaa(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_mefasaa):
            matched_smiles = Chem.MolToSmiles(mol_in)
            self.mefasaa.append((matched_smiles, mol_in))
            return True
        return False

    def match_ftoh(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_ftoh):
            matched_smiles = Chem.MolToSmiles(mol_in)
            self.ftoh.append((matched_smiles, mol_in))
            return True
        return False

    def match_pfoh(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_pfoh):
            matched_smiles = Chem.MolToSmiles(mol_in)
            self.pfoh.append((matched_smiles, mol_in))
            return True
        return False

    def match_pfsa(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_pfsa):
            matched_smiles = Chem.MolToSmiles(mol_in)
            self.pfsa.append((matched_smiles, mol_in))
            return True
        return False

    def match_pfal(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_pfal):
            matched_smiles = Chem.MolToSmiles(mol_in)
            self.pfal.append((matched_smiles, mol_in))
            return True
        return False

#     def match_CF_chain(self, mol_in) -> bool:
#         if mol_in and mol_in.HasSubstructMatch(self.sub_mol_CF_chain):
#             matched_smiles = Chem.MolToSmiles(mol_in)
#             self.CF_chain.append((matched_smiles, mol_in))
#             return True
#         return False

    def classify(self, mol_in) -> List[str]:
        mol_cats: List[str] = []
        if self.match_pfeca(mol_in):
            mol_cats.append("PFECA")
        if self.match_pfasa(mol_in):
            mol_cats.append("PFASA") 
        if self.match_fasa(mol_in):
            mol_cats.append("FASA")
        if self.match_pfca(mol_in):
            mol_cats.append("PFCA")
        if self.match_fts(mol_in):
            mol_cats.append("FTS")
        if self.match_ftca(mol_in):
            mol_cats.append("FTCA")
        if self.match_mefasaa(mol_in):
            mol_cats.append("MeFASAA")
        if self.match_ftoh(mol_in):
            mol_cats.append("FTOH")
        if self.match_pfoh(mol_in):
            mol_cats.append("PFOH")
        if self.match_pfsa(mol_in):
            mol_cats.append("PFSA")
        if self.match_pfal(mol_in):
            mol_cats.append("PFAL")
#         if self.match_CF_chain(mol_in):
#             mol_cats.append("CF")
        return mol_cats
         
    def create_images(self, mol_array):
        try: 
            for smi, mol in self.mol_array: 
                fname_match = f"{smi}.png" 
                path = OUTPUT_DIR / fname_match 
                img_match = Draw.MolToImage(mol, size=SIZE)
                img_match.save(path)
        except Exception as e: 
            raise RuntimeError(f'Error in create_images(): {e}') 
