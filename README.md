# AtomicFingerprint
Numerical fingerprinting tools for RNA

## Quick Start
### Install
```
make clean
make
```

### Example
- Cavity Distances
```
label=native
./bin/test -mol2 tests/1aju/cavity_${label}.mol2 tests/1aju/complex_cavity_${label}.pdb > tests/output_${label}.txt
label=decoy
./bin/test -mol2 tests/1aju/cavity_${label}.mol2 tests/1aju/complex_cavity_${label}.pdb > tests/output_${label}.txt
python py/process_dists.py
```
