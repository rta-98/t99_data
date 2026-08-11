from rdkit import Chem 
from typing import Optional, List  

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
 
    if sub_mol_pfeca is None: 
        print("ERROR") 

    def match_pfeca(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_pfeca):
            return True
        return False

    def match_fasa(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_fasa):
            return True
        return False

    def match_pfasa(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_pfasa):
            return True
        return False

    def match_pfca(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_pfca):
            return True
        return False

    def match_fts(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_fts):
            return True
        return False

    def match_ftca(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_ftca):
            return True
        return False

    def match_mefasaa(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_mefasaa):
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

    def match_pfsa(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_pfsa):
            return True
        return False

    def match_pfal(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_pfal):
            return True
        return False

    def classify(self, mol) -> List[str]:
        cond_flag = False 
        mol_cats: List[str] = []
        if self.match_pfeca(mol):
            mol_cats.append("PFECA")
            cond_flag = True 
        if self.match_pfasa(mol):
            mol_cats.append("PFASA") 
            cond_flag = True 
        if self.match_fasa(mol):
            mol_cats.append("FASA")
            cond_flag = True 
        if self.match_pfca(mol):
            mol_cats.append("PFCA")
            cond_flag = True 
        if self.match_fts(mol):
            mol_cats.append("FTS")
            cond_flag = True 
        if self.match_ftca(mol):
            mol_cats.append("FTCA")
            cond_flag = True 
        if self.match_mefasaa(mol):
            mol_cats.append("MeFASAA")
            cond_flag = True 
        if self.match_ftoh(mol):
            mol_cats.append("FTOH")
            cond_flag = True 
        if self.match_pfoh(mol):
            mol_cats.append("PFOH")
            cond_flag = True 
        if self.match_pfsa(mol):
            mol_cats.append("PFSA")
            cond_flag = True 
        if self.match_pfal(mol):
            mol_cats.append("PFAL")
            cond_flag = True 
        if cond_flag is False:
            mol_cats.append("Unk") 

        return str(mol_cats)

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

