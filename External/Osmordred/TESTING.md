# Osmordred Testing Guide

This document describes the comprehensive test suite for Osmordred functionality in RDKit.

## Overview

Osmordred provides molecular descriptors for cheminformatics applications. The test suite covers both C++ and Python implementations to ensure functionality and correctness across all platforms.

## Test Structure

### C++ Tests (`test_osmordred.cpp`)

The C++ test suite is located in `External/Osmordred/test_osmordred.cpp` and uses the Catch2 testing framework.

#### Test Categories

1. **Basic Functionality Tests**
   - Simple molecule tests (ethanol, benzene)
   - Atom and bond counting
   - Basic molecular properties

2. **Matrix Descriptor Tests**
   - Distance matrix descriptors
   - Adjacency matrix descriptors
   - Detour matrix descriptors
   - Barysz matrix descriptors
   - Version parameter testing

3. **Carbon Type Tests**
   - Carbon type descriptors
   - Version parameter testing

4. **EState Descriptor Tests**
   - EState descriptors
   - Extended EState descriptors

5. **Chi Descriptor Tests**
   - All Chi descriptors
   - Individual Chi descriptor types

6. **Count Function Tests**
   - Acidic/basic group counts
   - Aromatic atom/bond counts

7. **Information Content Tests**
   - Information content descriptors
   - Different radius parameters

8. **Edge Case Tests**
   - Empty molecules
   - Single atom molecules
   - Large complex molecules

9. **Validation Tests**
   - Consistency checks
   - Physical property validation

### Python Tests (`test_osmordred.py`)

The Python test suite is located in `External/Osmordred/Wrap/test_osmordred.py` and uses the unittest framework.

#### Test Categories

1. **Basic Functionality Tests**
   - Simple molecule tests
   - Atom and bond counting
   - Molecular weight calculations

2. **Aromatic Functionality Tests**
   - Aromatic counts
   - Ring counts
   - Aromatic atom/bond counts

3. **Matrix Descriptor Tests**
   - Distance matrix descriptors
   - Adjacency matrix descriptors
   - Detour matrix descriptors
   - Barysz matrix descriptors

4. **Physical Property Tests**
   - Molecular weight
   - Van der Waals volume
   - McGowan volume
   - Polarizability
   - SLogP and LogS

5. **Topological Descriptor Tests**
   - Wiener index
   - Balaban J index
   - Bertz CT index
   - Zagreb index
   - Eccentric connectivity index

6. **Chi Descriptor Tests**
   - All Chi descriptors
   - Individual Chi descriptor types

7. **EState Descriptor Tests**
   - EState descriptors
   - BEState descriptors
   - HEState descriptors

8. **Count Function Tests**
   - Acidic/basic group counts
   - Aromatic atom/bond counts

9. **Complex Molecule Tests**
   - Large molecule functionality
   - Error handling

10. **Edge Case Tests**
    - None molecules
    - Empty SMILES
    - Single atoms

11. **Consistency Tests**
    - Multiple call consistency
    - Result validation

12. **Information Content Tests**
    - Information content descriptors

13. **Function Availability Tests**
    - All expected functions exist

## Running Tests

### Using the Test Runner Script

The main test runner script is located at `Scripts/test_osmordred.py`:

```bash
# Run all tests
python Scripts/test_osmordred.py

# Run only C++ tests
python Scripts/test_osmordred.py --cpp-only

# Run only Python tests
python Scripts/test_osmordred.py --python-only

# Quick validation
python Scripts/test_osmordred.py --quick

# Specify build directory
python Scripts/test_osmordred.py --build-dir /path/to/build
```

### Manual C++ Testing

```bash
# Build the tests
cd build
make testOsmordred

# Run the tests
./External/Osmordred/testOsmordred
```

### Manual Python Testing

```bash
# Run Python tests
cd External/Osmordred/Wrap
python -m unittest test_osmordred

# Or run specific test
python -c "
import sys
sys.path.insert(0, 'External/Osmordred/Wrap')
from test_osmordred import TestOsmordred
import unittest
suite = unittest.TestLoader().loadTestsFromTestCase(TestOsmordred)
unittest.TextTestRunner(verbosity=2).run(suite)
"
```

## Test Molecules

The test suite uses several standard test molecules:

1. **Ethanol (CCO)**: Simple aliphatic molecule
2. **Benzene (c1ccccc1)**: Simple aromatic molecule
3. **Complex Molecule**: Large pharmaceutical-like molecule
4. **t-Butyl Benzene**: Medium-sized molecule with branching

## Expected Results

### C++ Tests

All C++ tests should:
- Compile without errors
- Run without crashes
- Return non-empty results for all descriptor functions
- Maintain consistency across multiple calls
- Handle edge cases gracefully

### Python Tests

All Python tests should:
- Import without errors
- Return expected data types (list for descriptors, int for counts)
- Handle None and invalid inputs appropriately
- Maintain consistency across multiple calls
- Provide meaningful error messages

## Coverage

The test suite covers:

- **All 100+ Osmordred functions**
- **Multiple molecule types** (aliphatic, aromatic, complex)
- **Edge cases** (empty, single atom, large molecules)
- **Parameter variations** (different versions, extended modes)
- **Error handling** (invalid inputs, None values)
- **Consistency validation** (reproducible results)

## Troubleshooting

### Common Issues

1. **Build Errors**
   - Ensure LAPACK and LAPACKE are installed
   - Check Eigen3 installation
   - Verify CMake configuration

2. **Import Errors**
   - Ensure RDKit is properly installed
   - Check Python path configuration
   - Verify Osmordred module compilation

3. **Test Failures**
   - Check molecule validity
   - Verify descriptor function availability
   - Review error messages for specific issues

### Debug Mode

For debugging, run tests with verbose output:

```bash
# C++ tests with verbose output
./External/Osmordred/testOsmordred --verbose

# Python tests with verbose output
python Scripts/test_osmordred.py --verbose
```

## Contributing

When adding new Osmordred functionality:

1. Add corresponding C++ tests in `test_osmordred.cpp`
2. Add corresponding Python tests in `test_osmordred.py`
3. Update this documentation
4. Ensure tests pass on all platforms

## Performance Considerations

- C++ tests should complete within 5 minutes
- Python tests should complete within 2 minutes
- Large molecule tests may take longer but should not timeout
- Memory usage should remain reasonable

## Platform Support

Tests are designed to work on:
- Linux (GCC/Clang)
- macOS (Clang)
- Windows (MSVC)
- Python 3.7+

## Continuous Integration

The test suite is integrated into RDKit's CI pipeline and runs:
- On every pull request
- On every merge to main
- On all supported platforms
- With both Debug and Release builds 