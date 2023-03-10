from src.main import gen_smiles
from rdkit import Chem

def test_methane():
    molfile = "tests/mol/methane.mol"
    smiles_to_test = gen_smiles(molfile)
    expected_smiles = "C"
    assert smiles_to_test == expected_smiles

def test_ethane():
    molfile = "tests/mol/ethane.mol"
    smiles_to_test = gen_smiles(molfile)
    expected_smiles = "CC"
    assert smiles_to_test == expected_smiles

def test_propane():
    molfile = "tests/mol/propane.mol"
    smiles_to_test = gen_smiles(molfile)
    expected_smiles = "CCC"
    assert smiles_to_test == expected_smiles

def test_hexane():
    molfile = "tests/mol/hexane.mol"
    smiles_to_test = gen_smiles(molfile)
    expected_smiles = "CCCCCC"
    assert smiles_to_test == expected_smiles

def test_isobutane():
    molfile = "tests/mol/isobutane.mol"
    smiles_to_test = gen_smiles(molfile)
    expected_smiles = "CC(C)C"
    assert smiles_to_test == expected_smiles
    
def test_dimethylethylphenol():
    molfile = "tests/mol/dimethylethylphenol.mol"
    smiles_to_test = gen_smiles(molfile)
    mol = Chem.MolFromSmiles(smiles_to_test)
    assert mol is not None
    truth = Chem.MolFromMolFile(molfile)
    assert Chem.MolToSmiles(mol) == Chem.MolToSmiles(truth)

def test_ethylmethylketone():
    molfile = "tests/mol/ethylmethylketone.mol"
    smiles_to_test = gen_smiles(molfile)
    mol = Chem.MolFromSmiles(smiles_to_test)
    assert mol is not None
    truth = Chem.MolFromMolFile(molfile)
    assert Chem.MolToSmiles(mol) == Chem.MolToSmiles(truth)

def test_toluene():
    molfile = "tests/mol/toluene.mol"
    smiles_to_test = gen_smiles(molfile)
    mol = Chem.MolFromSmiles(smiles_to_test)
    assert mol is not None
    truth = Chem.MolFromMolFile(molfile)
    assert Chem.MolToSmiles(mol) == Chem.MolToSmiles(truth)

def test_acetaminophen():
    molfile = "tests/mol/acetaminophen.mol"
    smiles_to_test = gen_smiles(molfile)
    mol = Chem.MolFromSmiles(smiles_to_test)
    assert mol is not None
    truth = Chem.MolFromMolFile(molfile)
    assert Chem.MolToSmiles(mol) == Chem.MolToSmiles(truth)

def test_diphenyl():
    molfile = "tests/mol/diphenyl.mol"
    smiles_to_test = gen_smiles(molfile)
    mol = Chem.MolFromSmiles(smiles_to_test)
    assert mol is not None
    truth = Chem.MolFromMolFile(molfile)
    assert Chem.MolToSmiles(mol) == Chem.MolToSmiles(truth)

def test_3_amino_2_naphtoic_acid():
    molfile = "tests/mol/3-amino-2-naphthoic_acid.mol"
    smiles_to_test = gen_smiles(molfile)
    mol = Chem.MolFromSmiles(smiles_to_test)
    assert mol is not None
    truth = Chem.MolFromMolFile(molfile)
    assert Chem.MolToSmiles(mol) == Chem.MolToSmiles(truth)
