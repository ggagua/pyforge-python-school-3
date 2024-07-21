from rdkit import Chem


def substructure_search(mols, mol):
    """
      Func returns a list of all molecules from argument 1 that contain substructure from argument 2.

      Arguments:
          mols (list): A list of SMILES strings containing the molecules to search.
          mol (str): A SMILES string which has the substructure to search for.

      Returns:
          list: A list of SMILES strings representing the molecules that contain the substructure.
          str: An error message if an exception occurs.
      """
    try:
        substructure = Chem.MolFromSmiles(mol)
        if substructure is None:
            raise ValueError("Unknown substructure SMILES string")

        matches = []
        for mol_smiles in mols:
            try:
                molecule = Chem.MolFromSmiles(mol_smiles)
                if molecule is None:
                    raise ValueError(f"Invalid molecule SMILES string: {mol_smiles}")

                if molecule.HasSubstructMatch(substructure):
                    matches.append(mol_smiles)
            except ValueError as ve:
                print(ve)

        return matches
    except Exception as e:
        return f"Error: {str(e)}"


final = substructure_search(["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "c1ccccc1")
print(final)
