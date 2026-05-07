# RDKit 2D Coordinate Generation Timing Script

## Overview

`time_coord_gen.sh` benchmarks 2D coordinate generation using the locally built RDKit (not system RDKit).

## Usage

```bash
# Time coordinate generation from a file
./time_coord_gen.sh test_smiles.txt

# Or pipe from stdin
cat smiles.txt | ./time_coord_gen.sh -
```

## Input Format

The input should be a text file with one SMILES per line:
```
c1ccccc1
c1ccc2c(c1)ccc1ccccc12
CC(C)Cc1ccc(cc1)C(C)C(=O)O
```

## Output

The script outputs:
- Per-molecule timing results
- Summary statistics (mean, median, std dev, min, max)
- Total throughput in molecules/sec
- Success/failure counts

## Example Output

```
[1/10] ✓ 0.85ms - c1ccccc1
[2/10] ✓ 48.64ms - c1ccc2c(c1)ccc1ccccc12
...

======================================================================
SUMMARY
======================================================================
Total molecules:     10
Successful:          9
Failed:              1

Total time:          0.052s
Mean time:           5.76ms
Median time:         0.48ms
Std dev:             16.08ms
Min time:            0.11ms
Max time:            48.64ms
Throughput:          173.5 molecules/sec
```

## Implementation Details

- Uses local RDKit build from `$SCHRODINGER/rdkit`
- Runs with templates enabled (`useRingTemplates=True`)
- Python script is in `~/bin/time_coord_gen.py` to avoid circular imports
- Uses `$SCHRODINGER/run python3` for proper environment setup
- Measures only coordinate generation time, not parsing or output
