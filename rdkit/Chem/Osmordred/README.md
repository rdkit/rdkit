# Osmordred Descriptors for RDKit

This module provides easy access to Osmordred molecular descriptors directly from RDKit. It includes both single-molecule and batch processing capabilities with parallel execution support.

## Installation

The Osmordred module is included with RDKit when built with Osmordred support enabled. Make sure you have the following dependencies:

```bash
conda install -c conda-forge tqdm numpy pandas
```

## Quick Start

```python
from rdkit import Chem
from rdkit.Chem import Osmordred

# Calculate descriptors for a single molecule
smiles = "CCO"
descriptors, names = Osmordred.CalcOsmordred(smiles, version=2, names=True)
print(f"Got {len(descriptors)} descriptors")

# Batch processing
smiles_list = ["CCO", "c1ccccc1", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"]
results_df = Osmordred.Calculate(smiles_list, n_jobs=4, version=2)
print(f"Results shape: {results_df.shape}")
```

## Functions

### `CalcOsmordred(smiles, version=2, names=False, mynames=None)`

Calculate Osmordred descriptors for a single SMILES string.

**Parameters:**
- `smiles` (str): SMILES string of the molecule
- `version` (int): Version of descriptors (1 or 2, default: 2)
- `names` (bool): Whether to return descriptor names along with values
- `mynames` (List[str], optional): Custom descriptor names list

**Returns:**
- If `names=False`: numpy array of descriptor values
- If `names=True`: tuple of (descriptor_values, descriptor_names)

**Example:**
```python
# Get descriptors with names
descriptors, names = Osmordred.CalcOsmordred("CCO", version=2, names=True)
print(f"First descriptor: {names[0]} = {descriptors[0]}")

# Get descriptors without names
descriptors = Osmordred.CalcOsmordred("CCO", version=2, names=False)
print(f"Number of descriptors: {len(descriptors)}")
```

### `Calculate(smiles_list, ids=None, n_jobs=4, version=2, names=False, mynames=None)`

Compute molecular descriptors for a list of SMILES with parallel processing.

**Parameters:**
- `smiles_list` (List[str]): List of SMILES strings
- `ids` (List, optional): List of unique identifiers (same length as smiles_list)
- `n_jobs` (int): Number of parallel processes (default: 4)
- `version` (int): Version of descriptors (default: 2)
- `names` (bool): Whether to include names
- `mynames` (List[str], optional): Custom names list

**Returns:**
- `pd.DataFrame`: DataFrame with results, indexed by ID with SMILES as first column

**Example:**
```python
smiles_list = ["CCO", "c1ccccc1", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"]
ids = ["mol_1", "mol_2", "mol_3"]

results_df = Osmordred.Calculate(
    smiles_list=smiles_list,
    ids=ids,
    n_jobs=4,
    version=2
)
print(results_df.head())
```

### `CalcOsmordredFromMol(mol, version=2, names=False, mynames=None)`

Calculate Osmordred descriptors from an RDKit molecule object.

**Parameters:**
- `mol`: RDKit molecule object
- `version` (int): Version of descriptors (1 or 2, default: 2)
- `names` (bool): Whether to return descriptor names along with values
- `mynames` (List[str], optional): Custom descriptor names list

**Returns:**
- If `names=False`: numpy array of descriptor values
- If `names=True`: tuple of (descriptor_values, descriptor_names)

**Example:**
```python
mol = Chem.MolFromSmiles("c1ccccc1")
descriptors = Osmordred.CalcOsmordredFromMol(mol, version=2)
print(f"Descriptors: {len(descriptors)}")
```

### `GetDescriptorNames(version=2)`

Get the list of descriptor names for a given version.

**Parameters:**
- `version` (int): Version of descriptors (1 or 2, default: 2)

**Returns:**
- `List[str]`: List of descriptor names

**Example:**
```python
names = Osmordred.GetDescriptorNames(version=2)
print(f"Version 2 has {len(names)} descriptors")
print(f"First 5: {names[:5]}")
```

## Version Differences

- **Version 1**: Original Mordred descriptors (1,599 descriptors)
- **Version 2**: Extended descriptors with additional features (3,585 descriptors)

Version 2 includes additional descriptors such as:
- Extended EState descriptors
- Additional fragment descriptors
- Enhanced information content calculations
- New topological indices

## Error Handling

The functions handle errors gracefully:

```python
# Invalid SMILES returns None
result = Osmordred.CalcOsmordred("invalid_smiles", version=2)
print(result)  # None

# Functions that fail return NaN arrays
# (handled automatically in batch processing)
```

## Performance Tips

1. **Use batch processing** for multiple molecules:
   ```python
   # Good: Use Calculate for multiple molecules
   results = Osmordred.Calculate(smiles_list, n_jobs=4)
   
   # Less efficient: Loop through individual calls
   for smiles in smiles_list:
       result = Osmordred.CalcOsmordred(smiles)
   ```

2. **Adjust n_jobs** based on your system:
   ```python
   # For CPU-intensive tasks, use number of CPU cores
   import multiprocessing
   n_jobs = multiprocessing.cpu_count()
   results = Osmordred.Calculate(smiles_list, n_jobs=n_jobs)
   ```

3. **Reuse descriptor names** for batch processing:
   ```python
   # Get names once
   _, names = Osmordred.CalcOsmordred(smiles_list[0], names=True)
   
   # Use in batch processing
   results = Osmordred.Calculate(smiles_list, mynames=names)
   ```

## Examples

### Basic Usage
```python
from rdkit import Chem
from rdkit.Chem import Osmordred

# Single molecule
descriptors = Osmordred.CalcOsmordred("CCO", version=2)
print(f"Ethanol has {len(descriptors)} descriptors")

# Batch processing
smiles_list = ["CCO", "c1ccccc1", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"]
results_df = Osmordred.Calculate(smiles_list, n_jobs=2)
print(f"Processed {len(results_df)} molecules")
```

### Working with DataFrames
```python
import pandas as pd

# Load your data
df = pd.read_csv("molecules.csv")
smiles_list = df['smiles'].tolist()

# Calculate descriptors
results = Osmordred.Calculate(smiles_list, ids=df['id'].tolist())

# Merge with original data
final_df = pd.merge(df, results, left_on='id', right_index=True)
final_df.to_csv("molecules_with_descriptors.csv")
```

### Custom Processing
```python
# Get descriptor names for custom processing
_, names = Osmordred.CalcOsmordred("CCO", names=True)

# Process with custom names
results = Osmordred.Calculate(
    smiles_list=smiles_list,
    mynames=names,
    n_jobs=4
)

# Filter specific descriptors
abc_descriptors = [col for col in results.columns if col.startswith('ABCIndex')]
abc_results = results[abc_descriptors]
```

## File Structure

```
rdkit/Chem/Osmordred/
├── __init__.py          # Main module with functions
├── example.py           # Comprehensive examples
├── usage_example.py     # Usage demonstration
└── README.md           # This file
```

## Troubleshooting

### Import Errors
If you get import errors, ensure:
1. RDKit is built with Osmordred support
2. The `rdOsmordred.so` module is available
3. All dependencies are installed (`tqdm`, `numpy`, `pandas`)

### Performance Issues
- Use `n_jobs=1` for debugging
- Reduce batch size for memory-constrained systems
- Monitor memory usage with large datasets

### Version Compatibility
- Version 1: 1,599 descriptors
- Version 2: 3,585 descriptors
- Use version 2 for maximum descriptor coverage

## Contributing

The Osmordred module is part of RDKit. For issues or contributions:
1. Check the RDKit issue tracker
2. Ensure Osmordred support is enabled in your build
3. Test with both versions (1 and 2) of descriptors 