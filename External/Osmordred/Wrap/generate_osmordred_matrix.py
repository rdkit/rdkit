#!/usr/bin/env python3
"""
Script to generate newosmordred.csv file by calculating Osmordred descriptors
for all molecules in testosmordred.csv
"""

import sys
import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Osmordred

# Add RDKit paths
rdkit_build_path = "/Users/guillaume-osmo/Github/rdkit/mac-build"
if rdkit_build_path not in sys.path:
    sys.path.insert(0, rdkit_build_path)

def calculate_descriptors(smiles):
    """Calculate all Osmordred descriptors for a given SMILES string."""
    try:
        # Calculate descriptors directly from SMILES
        descriptors = Osmordred.CalcOsmordred(smiles)
        return descriptors
    except Exception as e:
        print(f"Error processing molecule: {e}")
        return None

def main():
    # Read reference data
    ref_df = pd.read_csv("../testosmordred.csv", na_values=['-', 'nan', 'NaN', ''])
    
    # Get SMILES and descriptor columns
    smiles_list = ref_df['smiles'].tolist()
    descriptor_cols = [col for col in ref_df.columns if col not in ['smiles', 'RECORDID']]
    
    # Calculate descriptors for each molecule
    results = []
    for i, smiles in enumerate(smiles_list):
        print(f"Processing molecule {i+1}/{len(smiles_list)}", end='\r')
        descriptors = calculate_descriptors(smiles)
        if descriptors is not None:
            results.append(descriptors)
        else:
            # If calculation fails, add NaN values
            results.append([np.nan] * len(descriptor_cols))
    
    print("\nCalculations complete!")
    
    # Create DataFrame with results
    results_df = pd.DataFrame(results, columns=descriptor_cols)
    results_df.insert(0, 'smiles', smiles_list)
    
    # Save to CSV
    results_df.to_csv('newosmordred.csv', index=False)
    print("Results saved to newosmordred.csv")

if __name__ == "__main__":
    main()