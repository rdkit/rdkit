#!/usr/bin/env python3
"""
Comprehensive analysis of all 60 descriptor types using corrected relative difference.
Uses 1e-6 precision tolerance for zero cases.
"""

import sys
import os
import pandas as pd
import numpy as np
from collections import defaultdict

# Add RDKit paths
rdkit_build_path = "/Users/guillaume-osmo/Github/rdkit/mac-build"
rdkit_chem_path = "/Users/guillaume-osmo/Github/rdkit/mac-build/rdkit/Chem"
if rdkit_build_path not in sys.path:
    sys.path.insert(0, rdkit_build_path)
if rdkit_chem_path not in sys.path:
    sys.path.insert(0, rdkit_chem_path)

def calculate_relative_difference_correct(ref_val, our_val, tolerance=1e-6):
    """Calculate relative difference with 1e-6 tolerance for zero cases."""
    if pd.isna(ref_val) or pd.isna(our_val):
        return np.nan
    
    abs_diff = abs(our_val - ref_val)
    
    # If absolute difference is within tolerance of zero, relative difference is 0%
    if abs_diff <= tolerance:
        return 0.0
    
    # If reference is zero, handle specially
    if abs(ref_val) <= tolerance:
        if abs(our_val) <= tolerance:
            return 0.0  # Both effectively zero = no difference
        else:
            return np.inf  # Infinite relative difference
    
    return (abs_diff / abs(ref_val)) * 100

def extract_descriptor_type(column_name):
    """Extract descriptor type from column name."""
    # Remove _1, _2, etc. suffix
    base_name = column_name.split('_')[0]
    return base_name

def analyze_all_60_types():
    """Analyze all 60 descriptor types with corrected relative difference."""
    print("🔬 COMPREHENSIVE ANALYSIS OF ALL 60 DESCRIPTOR TYPES")
    print("=" * 80)
    
    # Load data
    ref_df = pd.read_csv("../testosmordred.csv", na_values=['-', 'nan', 'NaN', ''])
    our_df = pd.read_csv("newosmordred.csv")
    
    print(f"Reference data shape: {ref_df.shape}")
    print(f"Our data shape: {our_df.shape}")
    
    # Get descriptor columns
    descriptor_cols = [col for col in ref_df.columns if col not in ['smiles', 'RECORDID']]
    
    # Convert to numeric
    ref_numeric = ref_df[descriptor_cols].apply(pd.to_numeric, errors='coerce')
    our_numeric = our_df[descriptor_cols].apply(pd.to_numeric, errors='coerce')
    
    # Group descriptors by type
    descriptor_types = defaultdict(list)
    for col in descriptor_cols:
        desc_type = extract_descriptor_type(col)
        descriptor_types[desc_type].append(col)
    
    print(f"\nFound {len(descriptor_types)} descriptor types")
    
    # Test different relative difference thresholds
    relative_thresholds = [0.1, 0.5, 1.0, 2.0, 5.0]
    
    for threshold in relative_thresholds:
        print(f"\n{'='*60}")
        print(f"RELATIVE DIFFERENCE THRESHOLD: {threshold}% (1e-6 tolerance)")
        print(f"{'='*60}")
        
        # Analyze each descriptor type
        type_results = []
        
        for desc_type, cols in descriptor_types.items():
            total_matches = 0
            total_valid = 0
            zero_tolerance_matches = 0
            all_relative_diffs = []
            
            for col in cols:
                ref_vals = ref_numeric[col]
                our_vals = our_numeric[col]
                
                # Count valid comparisons (both not NaN and ref_val != 0)
                valid_mask = ref_vals.notna() & our_vals.notna() & (ref_vals != 0)
                
                if valid_mask.sum() > 0:
                    for idx in valid_mask[valid_mask].index:
                        ref_val = ref_vals[idx]
                        our_val = our_vals[idx]
                        rel_diff = calculate_relative_difference_correct(ref_val, our_val)
                        
                        if rel_diff <= threshold:
                            total_matches += 1
                        
                        total_valid += 1
                        all_relative_diffs.append(rel_diff)
                        
                        # Count zero tolerance matches
                        if rel_diff == 0.0:
                            zero_tolerance_matches += 1
            
            if total_valid > 0:
                match_rate = (total_matches / total_valid) * 100
                avg_relative_diff = np.mean(all_relative_diffs) if all_relative_diffs else 0
                max_relative_diff = np.max(all_relative_diffs) if all_relative_diffs else 0
                min_relative_diff = np.min(all_relative_diffs) if all_relative_diffs else 0
                
                type_results.append({
                    'type': desc_type,
                    'match_rate': match_rate,
                    'matches': total_matches,
                    'total_valid': total_valid,
                    'num_columns': len(cols),
                    'zero_tolerance_matches': zero_tolerance_matches,
                    'avg_relative_diff': avg_relative_diff,
                    'max_relative_diff': max_relative_diff,
                    'min_relative_diff': min_relative_diff
                })
        
        # Sort by match rate (worst first)
        type_results.sort(key=lambda x: x['match_rate'])
        
        print(f"\nDescriptor Types Performance (threshold: {threshold}%)")
        print("-" * 100)
        print(f"{'Type':<25} {'Match%':<8} {'Matches':<8} {'Total':<8} {'Columns':<8} {'ZeroTol':<8} {'AvgRel%':<8} {'MaxRel%':<8}")
        print("-" * 100)
        
        for result in type_results:
            print(f"{result['type']:<25} {result['match_rate']:<8.2f} {result['matches']:<8d} {result['total_valid']:<8d} "
                  f"{result['num_columns']:<8d} {result['zero_tolerance_matches']:<8d} "
                  f"{result['avg_relative_diff']:<8.3f} {result['max_relative_diff']:<8.3f}")
        
        # Summary statistics
        total_types = len(type_results)
        perfect_types = sum(1 for r in type_results if r['match_rate'] == 100.0)
        excellent_types = sum(1 for r in type_results if r['match_rate'] >= 95.0)
        good_types = sum(1 for r in type_results if r['match_rate'] >= 90.0)
        poor_types = sum(1 for r in type_results if r['match_rate'] < 90.0)
        
        print(f"\nSummary:")
        print(f"  Total descriptor types: {total_types}")
        print(f"  Perfect (100%): {perfect_types}")
        print(f"  Excellent (≥95%): {excellent_types}")
        print(f"  Good (≥90%): {good_types}")
        print(f"  Poor (<90%): {poor_types}")
        print(f"  Average match rate: {np.mean([r['match_rate'] for r in type_results]):.2f}%")

def show_detailed_analysis():
    """Show detailed analysis of performance categories."""
    print(f"\n🔍 DETAILED ANALYSIS BY PERFORMANCE CATEGORY")
    print("=" * 80)
    
    # Load data
    ref_df = pd.read_csv("../testosmordred.csv", na_values=['-', 'nan', 'NaN', ''])
    our_df = pd.read_csv("newosmordred.csv")
    
    # Get descriptor columns
    descriptor_cols = [col for col in ref_df.columns if col not in ['smiles', 'RECORDID']]
    
    # Convert to numeric
    ref_numeric = ref_df[descriptor_cols].apply(pd.to_numeric, errors='coerce')
    our_numeric = our_df[descriptor_cols].apply(pd.to_numeric, errors='coerce')
    
    # Group descriptors by type
    descriptor_types = defaultdict(list)
    for col in descriptor_cols:
        desc_type = extract_descriptor_type(col)
        descriptor_types[desc_type].append(col)
    
    threshold = 0.1  # Use 0.1% threshold
    
    # Analyze each descriptor type
    type_results = []
    
    for desc_type, cols in descriptor_types.items():
        total_matches = 0
        total_valid = 0
        zero_tolerance_matches = 0
        all_relative_diffs = []
        
        for col in cols:
            ref_vals = ref_numeric[col]
            our_vals = our_numeric[col]
            
            # Count valid comparisons
            valid_mask = ref_vals.notna() & our_vals.notna() & (ref_vals != 0)
            
            if valid_mask.sum() > 0:
                for idx in valid_mask[valid_mask].index:
                    ref_val = ref_vals[idx]
                    our_val = our_vals[idx]
                    rel_diff = calculate_relative_difference_correct(ref_val, our_val)
                    
                    if rel_diff <= threshold:
                        total_matches += 1
                    
                    total_valid += 1
                    all_relative_diffs.append(rel_diff)
                    
                    # Count zero tolerance matches
                    if rel_diff == 0.0:
                        zero_tolerance_matches += 1
        
        if total_valid > 0:
            match_rate = (total_matches / total_valid) * 100
            avg_relative_diff = np.mean(all_relative_diffs) if all_relative_diffs else 0
            max_relative_diff = np.max(all_relative_diffs) if all_relative_diffs else 0
            
            type_results.append({
                'type': desc_type,
                'match_rate': match_rate,
                'matches': total_matches,
                'total_valid': total_valid,
                'num_columns': len(cols),
                'zero_tolerance_matches': zero_tolerance_matches,
                'avg_relative_diff': avg_relative_diff,
                'max_relative_diff': max_relative_diff,
                'columns': cols
            })
    
    # Sort by match rate
    type_results.sort(key=lambda x: x['match_rate'])
    
    # Categorize results
    perfect_types = [r for r in type_results if r['match_rate'] == 100.0]
    excellent_types = [r for r in type_results if 95.0 <= r['match_rate'] < 100.0]
    good_types = [r for r in type_results if 90.0 <= r['match_rate'] < 95.0]
    poor_types = [r for r in type_results if r['match_rate'] < 90.0]
    
    print(f"Perfect Types (100%): {len(perfect_types)}")
    print("-" * 50)
    for result in perfect_types:
        print(f"  {result['type']:<25} {result['match_rate']:<8.2f}% ({result['matches']}/{result['total_valid']}) - {result['num_columns']} columns")
    
    print(f"\nExcellent Types (95-99.99%): {len(excellent_types)}")
    print("-" * 50)
    for result in excellent_types:
        print(f"  {result['type']:<25} {result['match_rate']:<8.2f}% ({result['matches']}/{result['total_valid']}) - {result['num_columns']} columns")
    
    print(f"\nGood Types (90-94.99%): {len(good_types)}")
    print("-" * 50)
    for result in good_types:
        print(f"  {result['type']:<25} {result['match_rate']:<8.2f}% ({result['matches']}/{result['total_valid']}) - {result['num_columns']} columns")
    
    print(f"\nPoor Types (<90%): {len(poor_types)}")
    print("-" * 50)
    for result in poor_types:
        print(f"  {result['type']:<25} {result['match_rate']:<8.2f}% ({result['matches']}/{result['total_valid']}) - {result['num_columns']} columns")
        print(f"    Zero tolerance matches: {result['zero_tolerance_matches']}")
        print(f"    Avg rel diff: {result['avg_relative_diff']:.3f}%, Max: {result['max_relative_diff']:.3f}%")

def main():
    """Main analysis function."""
    print("🔬 COMPREHENSIVE ANALYSIS OF ALL 60 DESCRIPTOR TYPES")
    print("=" * 80)
    
    # Set environment
    os.environ['PYTHONPATH'] = '/Users/guillaume-osmo/Github/rdkit/mac-build'
    
    # Run analyses
    analyze_all_60_types()
    show_detailed_analysis()
    
    print("\n✅ Comprehensive analysis complete!")

if __name__ == "__main__":
    main() 