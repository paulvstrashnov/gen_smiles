import src.gen_smiles as gen_smiles
import os
import pytest
from rdkit import Chem

def test_methane():
    smiles = gen_smiles("tests/mol/methane.mol2")
    expected_smiles = "C"
    assert Chem.MolToSmiles(mol) == smiles

def test_ethane():
    smiles = gen_smiles("tests/mol/ethane.mol2")
    expected_smiles = "CC"
    mol = Chem.MolFromSmiles(expected_smiles)
    assert Chem.MolToSmiles(mol) == smiles

def test_propane():
    smiles = gen_smiles("tests/mol/propane.mol2")
    expected_smiles = "CCC"
    mol = Chem.MolFromSmiles(expected_smiles)
    assert Chem.MolToSmiles(mol) == smiles

def test_hexane():
    smiles = gen_smiles("tests/mol/hexane.mol2")
    expected_smiles = "CCCCCC"
    mol = Chem.MolFromSmiles(expected_smiles)
    assert Chem.MolToSmiles(mol) == smiles

def test_isobutane():
    smiles = gen_smiles("tests/mol/isobutane.mol2")
    expected_smiles = "CC(C)C"
    mol = Chem.MolFromSmiles(expected_smiles)
    assert Chem.MolToSmiles(mol) == smiles
    
 def test_dimethylethylphenol():
    smiles = gen_smiles("tests/mol/dimethylethylphenol.mol2")
    assert Chem.MolToSmiles(smiles) is not None

def test_benzene():
    smiles = gen_smiles("tests/mol/benzene.mol2")
    assert Chem.MolToSmiles(smiles) is not None

def test_toluene():
    smiles = gen_smiles("tests/mol/toluene.mol2")
    assert Chem.MolToSmiles(smiles) is not None

def test_acetaminophen():
    smiles = gen_smiles("tests/mol/acetaminophen.mol2")
    assert Chem.MolToSmiles(smiles) is not None

def test_diphenyl():
    smiles = gen_smiles("tests/mol/diphenyl.mol2")
    assert Chem.MolToSmiles(smiles) is not None
