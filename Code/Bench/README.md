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
