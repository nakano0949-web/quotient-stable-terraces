# Quotient-stable directed terraces of Z_{2m} — code

Python scripts accompanying the paper *Quotient-stable directed terraces of Z_{2m} and complete Latin squares*
(Tatsuya Nakano). Only the Python 3 standard library is used. Keep the four files in one folder.

| File | What it does | Typical call |
|---|---|---|
| `quotient_stable_terraces.py` | Enumerates quotient-stable directed terraces of Z_n (n = 2m, 4 \| n) with t_0 = 0: pruned depth-first search (cross-checked against full enumeration for n <= 12), counts of stable and translation-type terraces, orbits under units and reversal | `python3 quotient_stable_terraces.py 4 8 12 16` |
| `stable_terraces_fibers.py` | Structure theorem: sums the solution counts of the GF(2) systems over all pairs (P', R') and compares with the direct counts; prints the distribution of fibre dimensions; compares fibres with brute force for m = 4, 6 | `python3 stable_terraces_fibers.py 2 4 6 8 10` |
| `stable_terraces_lift.py` | Lifting algorithm Z_m -> Z_{2m} (Gaussian elimination over GF(2)); towers Z_4 -> ... -> Z_256 with checks at every step; mixed lifts at m = 8 | `python3 stable_terraces_lift.py` |
| `bop_theorem3_check.py` | Checks the Bailey–Ollis–Preece half-projection/lift counts for symmetric directed terraces of Z_n (n = 4..16) and counts how many projections are directed | `python3 bop_theorem3_check.py` |

## Expected results (as reported in the paper)

* Quotient-stable terraces with t_0 = 0 for 2m = 4, 8, 12, 16, 20: 2, 8, 56, 896, 106160
  (the case 2m = 20 takes a few minutes in pure Python; the fibre computation reproduces the same total independently).
* Non-symmetric ones: 0, 0, 24, 512, 96944.
* `bop_theorem3_check.py`: numbers of symmetric directed terraces 2, 4, 8, 48, 176, 576, 3712 for n = 4, 6, ..., 16,
  of which 2, 0, 8, 0, 32, 0, 384 have a directed half-projection.

## Reproducibility notes

* The searches for n = 16 and n = 20 are pruned depth-first searches; for n <= 12 they were compared with full enumeration.
* In Jupyter, call the functions directly, e.g. `main([4, 8, 12])`, instead of using command-line arguments.

## Citation

Please cite the paper and the archived version of this code (see `CITATION.cff`).

## License

See `LICENSE` (to be chosen by the author).
