from . import services 
import subprocess
from pathlib import Path 
from services import InternalValid
from rdkit import Chem 
from rdkit.Chem import rdchem
from rdkit.Chem.MolStandardize import rdMolStandardize 

class SmileFileParser(InternalValid): 

    def __init__(self, dir_path):  
        self.dir = dir_path
        self.smiles_list = []
        self.mols_list = []
        self.smi_fname = []
        self.smiles_dud_list = []
        self.pdb_files = []
        self.log_files = []
        self.log_mols = [] 
        self.log_files_mols_dict = {
                ".log path": [],
                "mol": []
        }

    def log_parse(self) -> dict: 
    #    results = {"orca": [], "g16": [], "unk": []} 
        for idx, i in enumerate(files): 
            head = i.read_text(errors="ignore", encoding="utf-8")[:20000]
            if "O R C A" in head or "ORCA" in head: 
        #            results["orca"].append(tuple([idx,i])) 
                fmt = "orca" 
            elif "Entering Gaussian" in head or "G16" in head or "Gaussian Inc." in head:
        #            results["g16"].append(tuple([idx,i])) 
                fmt = "g09" 
            else: 
        #            results["unk"].append(tuple([idx,i])) 
                fmt = "Unk" 
            for log_file in Path(self.dir).glob('*.log'): 
                sdf = log_file.with_suffix('.sdf')
                subprocess.run(
                    ["obabel", f"-i{fmt}", str(log_file), "-osdf", "-O", str(sdf)],
                    check=True
                )
                mol = Chem.SDMolSupplier(str(sdf), removeHs=False)[0]
                self.log_files_mols_dict[".log path"].append(log_file.stem) 
                self.log_files_mols_dict["mol"].append(mol) 

        return self.log_files_mols_dict 

    def smi_populate(self):
        for smi_file in Path(self.dir).rglob('*.smi'): 
            with open(smi_file) as f: 
                for line in f: 
                    tmp_smile = line.split()
                    if not tmp_smile: 
                        continue 
                    smiles = tmp_smile[0]
                    try:
                        canon_smi = InternalValid.validator(smiles)  
                        if canon_smi is None: 
                            self.smiles_dud_list.append(smiles)
                            continue 
                        mol = Chem.MolFromSmiles(canon_smi)
                        mol = Chem.AddHs(mol)
                        if mol is not None:
                            self.smiles_list.append(canon_smi)
                            self.mols_list.append(mol)
                            self.smi_fname.append(smi_file) 
                        else: 
                            self.smiles_dud_list.append(smiles)
                    except Exception as e: 
                        self.smiles_dud_list.append(smiles)
                        print(f"Error in smi_populate() {smiles}: {e}") 
                        continue 
