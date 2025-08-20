#!/usr/bin/env python3
"""
Usage example for Osmordred descriptors from RDKit

This script demonstrates how to use the new Osmordred functions
in a way similar to the original CalcOsmordred script.
"""

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Osmordred

def main():
    """Main function demonstrating Osmordred usage"""
    
    print("Osmordred Descriptors Usage Example")
    print("=" * 50)
    
    # Example 1: Single molecule calculation (similar to your original CalcOsmordred)
    print("\n1. Single Molecule Calculation:")
    smiles = "CCO"  # Ethanol
    
    # Get descriptors with names (similar to your original function)
    descriptors, names = Osmordred.CalcOsmordred(smiles, version=2, names=True)
    print(f"SMILES: {smiles}")
    print(f"Number of descriptors: {len(descriptors)}")
    print(f"First 5 descriptors:")
    for i in range(5):
        print(f"  {names[i]}: {descriptors[i]}")
    
    # Example 2: Batch processing (similar to your original Calculate function)
    print("\n2. Batch Processing:")
    
    # Sample SMILES list (you can replace with your actual data)
    smiles_list = [
        "CCO",           # Ethanol
        "c1ccccc1",      # Benzene
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",  # Ibuprofen
        "CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F",  # Celecoxib
        "C1=CC=C(C=C1)CC2=CC=CC=C2"  # Biphenyl
    ]
    
    # Get descriptor names for column headers
    _, descriptor_names = Osmordred.CalcOsmordred(smiles_list[0], version=2, names=True)
    
    # Process all molecules (similar to your original Calculate function)
    print(f"Processing {len(smiles_list)} molecules...")
    results_df = Osmordred.Calculate(
        smiles_list=smiles_list,
        ids=[f"mol_{i}" for i in range(len(smiles_list))],
        n_jobs=2,  # Number of parallel processes
        version=2,
        names=False,
        mynames=descriptor_names
    )
    
    print(f"Results shape: {results_df.shape}")
    print(f"Columns: {list(results_df.columns[:5])}...")
    
    # Show first few rows
    print("\nFirst few rows:")
    print(results_df.head())
    
    # Example 3: Working with molecule objects
    print("\n3. Working with RDKit Molecule Objects:")
    mol = Chem.MolFromSmiles("c1ccccc1")  # Benzene
    descriptors_from_mol = Osmordred.CalcOsmordredFromMol(mol, version=2, names=False)
    print(f"Molecule: {Chem.MolToSmiles(mol)}")
    print(f"Descriptors: {len(descriptors_from_mol)}")
    
    # Example 4: Version comparison
    print("\n4. Version Comparison:")
    desc_v1, names_v1 = Osmordred.CalcOsmordred(smiles, version=1, names=True)
    desc_v2, names_v2 = Osmordred.CalcOsmordred(smiles, version=2, names=True)
    
    print(f"Version 1: {len(desc_v1)} descriptors")
    print(f"Version 2: {len(desc_v2)} descriptors")
    print(f"Additional descriptors in v2: {len(desc_v2) - len(desc_v1)}")
    
    # Example 5: Error handling
    print("\n5. Error Handling:")
    invalid_smiles = "invalid_smiles"
    result = Osmordred.CalcOsmordred(invalid_smiles, version=2)
    print(f"Invalid SMILES result: {result}")
    
    # Example 6: Getting descriptor names
    print("\n6. Getting Descriptor Names:")
    names_v2 = Osmordred.GetDescriptorNames(version=2)
    print(f"Version 2 descriptor names: {len(names_v2)}")
    print(f"First 10 names: {names_v2[:10]}")
    
    # Example 7: Saving results to CSV (similar to your original script)
    print("\n7. Saving Results to CSV:")
    # Save the results DataFrame
    output_file = "osmordred_results.csv"
    results_df.to_csv(output_file)
    print(f"Results saved to: {output_file}")
    
    print("\n" + "=" * 50)
    print("Usage example completed successfully!")
    
    return results_df

if __name__ == "__main__":
    try:
        results = main()
        print(f"\nFinal results DataFrame shape: {results.shape}")
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc() 