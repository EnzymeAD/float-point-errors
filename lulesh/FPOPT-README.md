1. `make ser-single-forward.exe` to get an original LULESH binary.
2. `make ser-single-gradient.exe` to get an instrumented LULESH binary.
3. `script -c "time ./ser-single-gradient.exe" lulesh.txt` to get the profile (`lulesh.txt`).
4. `script -c "time make ser-single-forward-fpopt.exe" out` to get an optimized LULESH with 0 computation cost (DP table and critical cost budgets saved in `out`).
5. WIP: Evaluation of optimized binaries

## WIP Evaluation
1. First generate a ground truth LULESH binary and use it to get reference results (final original energy? also there are some symmetry checks inside LULESH)
2. Generate N optimized binaries using a subset of critical cost budgets.
3. Compare results from optimized binaries with reference results.
4. Measure runtimes of optimized binaries and the original binary.
5. Plot results (see FPBench scripts)