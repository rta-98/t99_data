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


class SubstructMatch(InternalValid):
    """
        Motif Smarts 
    """
    sub_mol_pfeca = Chem.MolFromSmarts('[*]-[#8]-[#6](-[#6](=[#8])-[#8]-[#1])(-[#9])-[*]') # perfluoroether carboxylic acids

    sub_mol_ftsm = Chem.MolFromSmarts('[#7](-[#16](=[#8])(=[#8])-[#6](-[#6](-[*])(-[#1])-[#1])(-[#1])-[#1])(-[#1])-[#1]') # fluorotelomer sulfonamides 
    sub_mol_fasa = Chem.MolFromSmarts('[#7](-[#16](=[#8])(=[#8])-[#6](-[#9])(-[#9])-[*])(-[#1])-[#1]') # perfluoroalkane sulfonamides 

    sub_mol_pfca = Chem.MolFromSmarts('[#6](=[#8])(-[#8]-[#1])-[#6](-[*])(-[#9])-[#9]') # perfluorocarboxylic acid 
    sub_mol_ftca = Chem.MolFromSmarts('[#6](-[#1])(-[#1])(-[#6](-[#6](-[#8]-[#1])=[#8])(-[#1])-[#1])-[*]') # fluorotelomer carboxylic acid

    sub_mol_pfsa = Chem.MolFromSmarts('[#8]=[#16](-[#8]-[#1])(=[#8])-[#6](-[#9])(-[#9])-[*]') # perfluoroalkane sulfonic acids 
    sub_mol_fts = Chem.MolFromSmarts('[#6](-[#1])(-[#1])(-[#6](-[#16](=[#8])(-[#8]-[#1])=[#8])(-[#1])-[#1])-[*]') # fLuorotelomer sulfonic acid 

    sub_mol_ftoh = Chem.MolFromSmarts('[#8](-[#6](-[#6](-[#1])(-[#1])-[*])(-[#1])-[#1])-[#1]') # fluorotelomer alcohols 
    sub_mol_pfoh = Chem.MolFromSmarts('[#8](-[#6](-[#6](-[*])(-[#9])-[#9])(-[#9])-[#9])-[#1]') # perfluoroalcohols 

    sub_mol_mefasaa = Chem.MolFromSmarts('[#8](-[#1])-[#6](=[#8])-[#6](-[#1])(-[#1])-[#7](-[#1])-[#16](=[#8])(=[#8])-[#6](-[#9])(-[#9])-[*]') # N-methyl perfluoroalkane sulfonamido acetic acid 

    sub_mol_pfal = Chem.MolFromSmarts('[#6](=[#8])(-[#9])-[#6](-[#9])(-[#9])-[*]') # perfluoroaldehydes  
 
    """
        Tail Smarts
    """

    sub_mol_CF3 = Chem.MolFromSmarts('[#6](-[*])(-[#9])(-[#9])-[#9]')

    sub_mol_CF2 = Chem.MolFromSmarts('[#6](-[*])(-[#9])(-[#9])-[*]')

    if sub_mol_pfeca is None: 
        print("ERROR") 

    def match_pfeca(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_pfeca):
            return True
        return False

    def match_ftsm(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_ftsm):
            return True
        return False

    def match_fasa(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_fasa):
            return True
        return False

    def match_pfca(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_pfca):
            return True
        return False

    def match_ftca(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_ftca):
            return True
        return False

    def match_pfsa(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_pfsa):
            return True
        return False

    def match_fts(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_fts):
            return True
        return False

    def match_ftoh(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_ftoh):
            return True
        return False

    def match_pfoh(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_pfoh):
            return True
        return False

    def match_mefasaa(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_mefasaa):
            return True
        return False

    def match_pfal(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_pfal):
            return True
        return False

    def classifyMotif(self, mol) -> str:
        cond_flag = False 
        if self.match_pfeca(mol):
            mol_cats = "PFECA"
            cond_flag = True 
        if self.match_ftsm(mol):
            mol_cats = "FTSm" 
            cond_flag = True 
        if self.match_fasa(mol):
            mol_cats = "FASA"
            cond_flag = True 
        if self.match_pfca(mol):
            mol_cats = "PFCA"
            cond_flag = True 
        if self.match_ftca(mol):
            mol_cats = "FTCA"
            cond_flag = True 
        if self.match_pfsa(mol):
            mol_cats = "PFSA"
            cond_flag = True 
        if self.match_fts(mol):
            mol_cats = "FTS"
            cond_flag = True 
        if self.match_ftoh(mol):
            mol_cats = "FTOH"
            cond_flag = True 
        if self.match_pfoh(mol):
            mol_cats = "PFOH"
            cond_flag = True 
        if self.match_mefasaa(mol):
            mol_cats = "MeFASAA"
            cond_flag = True 
        if self.match_pfal(mol):
            mol_cats = "PFAL"
            cond_flag = True 
        if cond_flag is False:
            mol_cats = "Unk" 

        return str(mol_cats)

    def classifyTail(self, mol) -> dict:
        cf2_matches = mol.GetSubstructMatches(self.sub_mol_CF2)
        cf3_matches = mol.GetSubstructMatches(self.sub_mol_CF3)
        matches = {
            "Global: CF2": f"{len(cf2_matches)}",
            "Global: CF3": f"{len(cf3_matches)}"
        }
        return matches


