from rdkit import Chem
from rdkit.Chem import Draw

def substructure_search(mols, mol):
    substructure = Chem.MolFromSmiles(mol)
    matches = []
    try:
        for mol in mols:
            molecule = Chem.MolFromSmiles(mol)

            if molecule.HasSubstructMatch(substructure):
                matches.append(mol)
        return matches
    except:
        return f"Error: Invalid SMILES string"



# Perform a substructure search
final = substructure_search(["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "c1cccc1")
print(final)