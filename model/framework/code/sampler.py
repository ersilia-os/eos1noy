import requests
from standardiser import standardise
from rdkit import Chem
import urllib
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem

def calculate_similarity(ref_mol, mol_list):
    ref_fp = AllChem.GetMorganFingerprint(ref_mol, 2)
    mol_fps = [AllChem.GetMorganFingerprint(mol, 2) for mol in mol_list]
    similarities = [
        DataStructs.TanimotoSimilarity(ref_fp, mol_fp) for mol_fp in mol_fps
    ]
    return similarities


def sort_molecules_by_similarity(ref_mol, mol_list):
    similarities = calculate_similarity(ref_mol, mol_list)
    paired = list(zip(mol_list, similarities))
    sorted_mols = sorted(paired, key=lambda x: x[1], reverse=True)
    print("Sorted mols", len(sorted_mols))
    return [mol for mol, _ in sorted_mols]


class ChemblSampler(object):
    def __init__(self):
        pass

    def _run_chembl_sampler(self, origin_smiles: str, similarity: float = 40):
        url_smiles = urllib.parse.quote(origin_smiles)
        url = f"https://www.ebi.ac.uk/chembl/api/data/similarity/{url_smiles}/{similarity}?format=json"
        try:
            reply = requests.get(url)
        except requests.exceptions.Timeout:
            print("ChEMBL seems to be down right now")
            return []
        try:
            data = reply.json()
        except:
            return []
        smiles = [d["molecule_structures"]["canonical_smiles"] for d in data["molecules"]]
        smiles = list(set(smiles))
        return smiles

    def _sample(self, origin_smiles):
        sampled_smiles = self._run_chembl_sampler(origin_smiles=origin_smiles)
        std_smiles = []
        for smi in sampled_smiles:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            try:
                mol = standardise.run(mol)
            except:
                continue
            std_smiles += [Chem.MolToSmiles(mol)]
        return std_smiles

    def _sort_smiles(self, origin_smiles, sampled_smiles):
        ref_mol = Chem.MolFromSmiles(origin_smiles)
        if len(sampled_smiles) == 0:
            return []
        mol_list = []
        for smi in sampled_smiles:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            mol_list += [mol]
        sorted_mols = sort_molecules_by_similarity(ref_mol, mol_list)
        sorted_smiles = []
        for m in sorted_mols:
            if m is None:
                continue
            sorted_smiles += [Chem.MolToSmiles(m)]
        return sorted_smiles

    def sample(self, origin_smiles):
        sampled_smiles = []
        sampled_smiles += self._sample(origin_smiles)
        sampled_smiles = list(set(sampled_smiles))
        sampled_smiles_ordered = self._sort_smiles(origin_smiles, sampled_smiles)
        sampled_smiles_selected = sampled_smiles_ordered[:100]
        if len(sampled_smiles_selected) < 100:
            sampled_smiles_selected.extend([None] * (100 - len(sampled_smiles_selected)))
        return sampled_smiles_selected