import requests
from standardiser import standardise
from rdkit import Chem
import urllib

def run_chembl_sampler(origin_smiles: str, similarity: float = 40):
    """Function adapted from Andrew White's Exmol"""

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


class ChemblSampler(object):
    def __init__(self):
        pass

    def _sample(self, smiles):
        smiles = run_chembl_sampler(origin_smiles=smiles)
        std_smiles = []
        for smi in smiles:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            try:
                mol = standardise.run(mol)
            except:
                continue
            std_smiles += [Chem.MolToSmiles(mol)]
        return std_smiles

    def sample(self, smiles):
        sampled_smiles = []
        sampled_smiles += self._sample(smiles)
        sampled_smiles = list(set(sampled_smiles))
        sampled_smiles = sampled_smiles[:100]
        return sampled_smiles