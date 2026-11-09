# RDKit Benchmarks

To run:

```bash
mkdir build
cd build
cmake ..
cmake --build . --target bench -j "$(nproc)"
# see `./Code/Bench/bench --help` for options
export RDBASE=".."
./Code/Bench/bench
```

To run the benchmarks that use the sample molecules on your own molecules
instead, pass a SMILES file (one molecule per line, SMILES in the first
column; lines that fail to parse are skipped):

```bash
./Code/Bench/bench --smiles my_molecules.smi "[forcefield]"
```
