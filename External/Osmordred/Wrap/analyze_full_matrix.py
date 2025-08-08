#!/usr/bin/env python3
"""
Comprehensive analysis of Osmordred implementation vs reference data.
Compares full 179x3585 matrix and reports percentage of wrong columns per descriptor type.
"""

import sys
import os
import pandas as pd
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
# import seaborn as sns  # Optional for styling

# Add RDKit paths
rdkit_build_path = "/Users/guillaume-osmo/Github/rdkit/mac-build"
rdkit_chem_path = "/Users/guillaume-osmo/Github/rdkit/mac-build/rdkit/Chem"
if rdkit_build_path not in sys.path:
    sys.path.insert(0, rdkit_build_path)
if rdkit_chem_path not in sys.path:
    sys.path.insert(0, rdkit_chem_path)

def load_reference_data():
    """Load the reference CSV data."""
    csv_path = "../testosmordred.csv"
    print(f"Loading reference data from: {csv_path}")
    
    # Read with proper na_values to handle '-' values
    df = pd.read_csv(csv_path, na_values=['-', 'nan', 'NaN', ''])
    print(f"Reference data shape: {df.shape}")
    print(f"Columns: {list(df.columns[:5])} ... {list(df.columns[-5:])}")
    
    return df

def calculate_our_descriptors():
    """Calculate descriptors using our Osmordred implementation."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Osmordred
        
        # Load reference data to get SMILES
        ref_df = load_reference_data()
        smiles_list = ref_df['smiles'].tolist()
        
        print(f"Calculating descriptors for {len(smiles_list)} molecules...")
        
        # Calculate descriptors for each molecule
        our_results = []
        successful_mols = 0
        
        for i, smiles in enumerate(smiles_list):
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    print(f" Failed to parse SMILES {i+1}: {smiles}")
                    our_results.append([np.nan] * 3585)  # Fill with NaN
                    continue
                
                # Calculate descriptors using SMILES string
                descriptors = Osmordred.CalcOsmordred(smiles)
                our_results.append(descriptors)
                successful_mols += 1
                
                if (i + 1) % 20 == 0:
                    print(f" Processed {i+1}/{len(smiles_list)} molecules...")
                    
            except Exception as e:
                print(f" Error processing molecule {i+1}: {e}")
                our_results.append([np.nan] * 3585)  # Fill with NaN
        
        print(f"Successfully calculated descriptors for {successful_mols}/{len(smiles_list)} molecules")
        
        # Convert to DataFrame
        our_df = pd.DataFrame(our_results)
        
        # Get column names from reference (excluding smiles and RECORDID)
        ref_columns = [col for col in ref_df.columns if col not in ['smiles', 'RECORDID']]
        our_df.columns = ref_columns
        
        # Add smiles and RECORDID columns to match reference format
        our_df.insert(0, 'RECORDID', ref_df['RECORDID'])
        our_df.insert(0, 'smiles', ref_df['smiles'])
        
        # Save our matrix to CSV
        output_csv = "newosmordred.csv"
        our_df.to_csv(output_csv, index=False)
        print(f" Our calculated matrix saved to: {output_csv}")
        print(f"   Shape: {our_df.shape}")
        print(f"   Columns: {list(our_df.columns[:5])} ... {list(our_df.columns[-5:])}")
        
        return our_df, ref_df
        
    except ImportError as e:
        print(f"Import error: {e}")
        return None, None

def analyze_descriptor_types():
    """Analyze which columns belong to which descriptor types."""
    ref_df = load_reference_data()
    
    # Get descriptor columns (exclude smiles and RECORDID)
    descriptor_cols = [col for col in ref_df.columns if col not in ['smiles', 'RECORDID']]
    
    # Group columns by descriptor type
    descriptor_groups = defaultdict(list)
    
    for col in descriptor_cols:
        # Extract descriptor type from column name
        # Most columns follow pattern: DescriptorType_Index
        parts = col.split('_')
        if len(parts) >= 2:
            desc_type = parts[0]
            descriptor_groups[desc_type].append(col)
        else:
            # Handle single word descriptors
            descriptor_groups[col].append(col)
    
    return descriptor_groups

def compare_matrices():
    """Compare our matrix with reference matrix and analyze differences per descriptor type."""
    print("=" * 80)
    print("FULL MATRIX COMPARISON ANALYSIS")
    print("=" * 80)
    
    # Load data
    our_df, ref_df = calculate_our_descriptors()
    if our_df is None:
        print(" Failed to calculate our descriptors")
        return
    
    # Get descriptor columns
    descriptor_cols = [col for col in ref_df.columns if col not in ['smiles', 'RECORDID']]
    
    # Extract reference descriptor data
    ref_descriptors = ref_df[descriptor_cols]
    our_descriptors = our_df[descriptor_cols]
    
    print(f"Reference matrix shape: {ref_descriptors.shape}")
    print(f"Our matrix shape: {our_descriptors.shape}")
    
    # Analyze descriptor types
    descriptor_groups = analyze_descriptor_types()
    
    print(f"\nFound {len(descriptor_groups)} descriptor types:")
    for desc_type, cols in sorted(descriptor_groups.items()):
        print(f"  {desc_type}: {len(cols)} descriptors")
    
    # Compare matrices with different precision levels
    precision_levels = [1e-3, 1e-4, 1e-5, 1e-6]
    
    for precision in precision_levels:
        print(f"\n{'='*60}")
        print(f"ANALYSIS WITH PRECISION: {precision}")
        print(f"{'='*60}")
        
        # Convert to numeric and handle non-numeric values
        ref_numeric = ref_descriptors.apply(pd.to_numeric, errors='coerce')
        our_numeric = our_descriptors.apply(pd.to_numeric, errors='coerce')
        
        # Calculate differences
        differences = np.abs(ref_numeric.values - our_numeric.values) > precision
        
        # Overall statistics
        total_elements = differences.size
        total_differences = np.sum(differences)
        overall_error_rate = (total_differences / total_elements) * 100
        
        print(f"Overall error rate: {overall_error_rate:.2f}% ({total_differences:,}/{total_elements:,} elements)")
        
        # Analyze per descriptor type
        print(f"\nError rates per descriptor type:")
        print(f"{'Descriptor Type':<25} {'Error %':<10} {'Wrong Cols':<12} {'Total Cols':<12}")
        print("-" * 65)
        
        desc_type_errors = {}
        
        for desc_type, cols in sorted(descriptor_groups.items()):
            if not cols:
                continue
                
            # Get columns for this descriptor type
            type_cols = [col for col in cols if col in descriptor_cols]
            if not type_cols:
                continue
            
            # Calculate error rate for this descriptor type
            type_ref = ref_descriptors[type_cols].apply(pd.to_numeric, errors='coerce').values
            type_our = our_descriptors[type_cols].apply(pd.to_numeric, errors='coerce').values
            type_differences = np.abs(type_ref - type_our) > precision
            
            total_type_elements = type_differences.size
            type_differences_count = np.sum(type_differences)
            type_error_rate = (type_differences_count / total_type_elements) * 100
            
            # Count columns with any errors
            cols_with_errors = np.sum(np.any(type_differences, axis=0))
            total_cols = len(type_cols)
            col_error_rate = (cols_with_errors / total_cols) * 100
            
            desc_type_errors[desc_type] = {
                'error_rate': type_error_rate,
                'col_error_rate': col_error_rate,
                'wrong_cols': cols_with_errors,
                'total_cols': total_cols,
                'total_elements': total_type_elements,
                'wrong_elements': type_differences_count
            }
            
            print(f"{desc_type:<25} {col_error_rate:<10.2f} {cols_with_errors:<12} {total_cols:<12}")
        
        # Find descriptor types with >60% error rate
        high_error_types = {k: v for k, v in desc_type_errors.items() 
                           if v['col_error_rate'] > 60}
        
        if high_error_types:
            print(f"\n DESCRIPTOR TYPES WITH >60% ERROR RATE:")
            print(f"{'Descriptor Type':<25} {'Error %':<10} {'Wrong Cols':<12} {'Total Cols':<12}")
            print("-" * 65)
            for desc_type, stats in sorted(high_error_types.items(), 
                                         key=lambda x: x[1]['col_error_rate'], reverse=True):
                print(f"{desc_type:<25} {stats['col_error_rate']:<10.2f} {stats['wrong_cols']:<12} {stats['total_cols']:<12}")
        
        # Save detailed results
        results_df = pd.DataFrame([
            {
                'Descriptor_Type': desc_type,
                'Error_Rate_Percent': stats['col_error_rate'],
                'Wrong_Columns': stats['wrong_cols'],
                'Total_Columns': stats['total_cols'],
                'Wrong_Elements': stats['wrong_elements'],
                'Total_Elements': stats['total_elements'],
                'Element_Error_Rate': (stats['wrong_elements'] / stats['total_elements']) * 100
            }
            for desc_type, stats in desc_type_errors.items()
        ])
        
        results_df = results_df.sort_values('Error_Rate_Percent', ascending=False)
        output_file = f"osmordred_analysis_precision_{precision}.csv"
        results_df.to_csv(output_file, index=False)
        print(f"\n Detailed results saved to: {output_file}")
        
        # Create visualization
        plt.figure(figsize=(15, 10))
        
        # Top 20 descriptor types by error rate
        top_20 = results_df.head(20)
        
        plt.subplot(2, 1, 1)
        bars = plt.barh(top_20['Descriptor_Type'], top_20['Error_Rate_Percent'])
        plt.xlabel('Error Rate (%)')
        plt.title(f'Top 20 Descriptor Types by Error Rate (Precision: {precision})')
        plt.gca().invert_yaxis()
        
        # Color bars based on error rate
        for bar, rate in zip(bars, top_20['Error_Rate_Percent']):
            if rate > 60:
                bar.set_color('red')
            elif rate > 30:
                bar.set_color('orange')
            else:
                bar.set_color('green')
        
        # Error rate distribution
        plt.subplot(2, 1, 2)
        plt.hist(results_df['Error_Rate_Percent'], bins=20, alpha=0.7, color='skyblue')
        plt.axvline(x=60, color='red', linestyle='--', label='60% threshold')
        plt.xlabel('Error Rate (%)')
        plt.ylabel('Number of Descriptor Types')
        plt.title('Distribution of Error Rates Across Descriptor Types')
        plt.legend()
        
        plt.tight_layout()
        plot_file = f"osmordred_analysis_precision_{precision}.png"
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        print(f" Visualization saved to: {plot_file}")
        plt.close()

def main():
    """Main analysis function."""
    print(" OSMORDRED FULL MATRIX ANALYSIS")
    print("=" * 80)
    
    # Set environment
    os.environ['PYTHONPATH'] = '/Users/guillaume-osmo/Github/rdkit/mac-build'
    
    # Run comparison
    compare_matrices()
    
    print("\n Analysis complete!")

if __name__ == "__main__":
    main() 