#!/usr/bin/env python3
"""
Example usage of Osmordred descriptors from RDKit

This script demonstrates how to use the CalcOsmordred and Calculate functions
that are now available directly from rdkit.Chem.Osmordred
"""

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Osmordred

def example_single_molecule():
    """Example of calculating descriptors for a single molecule"""
    print("=== Single Molecule Example ===")
    
    # Example SMILES
    smiles = "CCO"  # Ethanol
    
    # Calculate descriptors with names
    descriptors, names = Osmordred.CalcOsmordred(smiles, version=2, names=True)
    
    print(f"SMILES: {smiles}")
    print(f"Number of descriptors: {len(descriptors)}")
    print(f"First 10 descriptors:")
    for i, (name, value) in enumerate(zip(names[:10], descriptors[:10])):
        print(f"  {name}: {value}")
    
    # Calculate without names (just values)
    descriptors_only = Osmordred.CalcOsmordred(smiles, version=2, names=False)
    print(f"\nDescriptors array shape: {descriptors_only.shape}")
    
    return descriptors, names

def example_from_mol_object():
    """Example of calculating descriptors from RDKit molecule object"""
    print("\n=== Molecule Object Example ===")
    
    # Create molecule from SMILES
    mol = Chem.MolFromSmiles("c1ccccc1")  # Benzene
    
    # Calculate descriptors directly from molecule
    descriptors = Osmordred.CalcOsmordredFromMol(mol, version=2, names=False)
    
    print(f"Molecule: {Chem.MolToSmiles(mol)}")
    print(f"Number of descriptors: {len(descriptors)}")
    print(f"First 5 descriptor values: {descriptors[:5]}")

def example_batch_processing():
    """Example of batch processing multiple molecules"""
    print("\n=== Batch Processing Example ===")
    
    # List of SMILES strings
    smiles_list = [
        "CCO",           # Ethanol
        "c1ccccc1",      # Benzene
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",  # Ibuprofen
        "CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F",  # Celecoxib
        "C1=CC=C(C=C1)CC2=CC=CC=C2"  # Biphenyl
    ]
    
    # Get descriptor names first
    _, descriptor_names = Osmordred.CalcOsmordred(smiles_list[0], version=2, names=True)
    
    # Process all molecules
    print(f"Processing {len(smiles_list)} molecules...")
    results_df = Osmordred.Calculate(
        smiles_list=smiles_list,
        ids=[f"mol_{i}" for i in range(len(smiles_list))],
        n_jobs=2,  # Use 2 parallel processes
        version=2,
        names=False,
        mynames=descriptor_names
    )
    
    print(f"Results DataFrame shape: {results_df.shape}")
    print(f"Columns: {list(results_df.columns[:5])}...")  # Show first 5 columns
    print(f"Index: {list(results_df.index)}")
    
    # Show first few rows
    print("\nFirst few rows:")
    print(results_df.head())
    
    return results_df

def example_version_comparison():
    """Example comparing different versions"""
    print("\n=== Version Comparison Example ===")
    
    smiles = "CCO"
    
    # Version 1
    desc_v1, names_v1 = Osmordred.CalcOsmordred(smiles, version=1, names=True)
    print(f"Version 1 descriptors: {len(desc_v1)}")
    
    # Version 2
    desc_v2, names_v2 = Osmordred.CalcOsmordred(smiles, version=2, names=True)
    print(f"Version 2 descriptors: {len(desc_v2)}")
    
    print(f"Version 2 has {len(desc_v2) - len(desc_v1)} additional descriptors")
    
    # Show some new descriptors in version 2
    new_descriptors = set(names_v2) - set(names_v1)
    print(f"New descriptors in version 2: {list(new_descriptors)[:5]}...")

def example_error_handling():
    """Example of error handling with invalid SMILES"""
    print("\n=== Error Handling Example ===")
    
    # Valid SMILES
    valid_smiles = "CCO"
    result_valid = Osmordred.CalcOsmordred(valid_smiles, version=2)
    print(f"Valid SMILES '{valid_smiles}': {len(result_valid)} descriptors")
    
    # Invalid SMILES
    invalid_smiles = "invalid_smiles"
    result_invalid = Osmordred.CalcOsmordred(invalid_smiles, version=2)
    print(f"Invalid SMILES '{invalid_smiles}': {result_invalid}")

def example_get_descriptor_names():
    """Example of getting descriptor names"""
    print("\n=== Descriptor Names Example ===")
    
    # Get names for version 2
    names_v2 = Osmordred.GetDescriptorNames(version=2)
    print(f"Version 2 has {len(names_v2)} descriptor names")
    print(f"First 10 names: {names_v2[:10]}")
    
    # Get names for version 1
    names_v1 = Osmordred.GetDescriptorNames(version=1)
    print(f"Version 1 has {len(names_v1)} descriptor names")

if __name__ == "__main__":
    print("Osmordred Descriptors Example")
    print("=" * 50)
    
    try:
        # Run all examples
        example_single_molecule()
        example_from_mol_object()
        example_batch_processing()
        example_version_comparison()
        example_error_handling()
        example_get_descriptor_names()
        
        print("\n" + "=" * 50)
        print("All examples completed successfully!")
        
    except Exception as e:
        print(f"Error running examples: {e}")
        import traceback
        traceback.print_exc() 