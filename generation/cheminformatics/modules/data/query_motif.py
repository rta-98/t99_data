from openbabel import openbabel as ob 
from openbabel import pybel
from pathlib import Path
#|%%--%%| <bBICAYgKS3|KzbsRcku5o>
pdb = Path("/home/yang/projects/t99_calc/data/generation/cheminformatics/script/pdb")
tfa_pdb = pdb / "tfa.pdb"
#|%%--%%| <KzbsRcku5o|vbbSdYfn3f>
mol = next(pybel.readfile("pdb", f"{tfa_pdb}"))

    def smi_populate(self):
        for smi_file in Path(self.dir).rglob('*.smi'): 
            pdb_file = smi_file.with_suffix('.pdb') 
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
                            self.pdb_files.append(pdb_file) 
                        else: 
                            self.smiles_dud_list.append(smiles)
                    except Exception as e: 
                        self.smiles_dud_list.append(smiles)
                        print(f"Error in smi_populate() {smiles}: {e}") 
                        continue 
