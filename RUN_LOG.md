# RUN_LOG: vec-t2-zero-data-starter

Recorded 2026-11-18 on Phoenix Forge (Linux), CPU only, pinned to 3 cores (`taskset -c 0-2`, `OMP_NUM_THREADS=3`). Fresh `uv sync` env: Python 3.14.7, veckit 0.1.2, anndata 0.13.4, numpy 2.5.5.

## Default run (synthetic gene names)
```
$ OMP_NUM_THREADS=3 taskset -c 0-2 uv run python -W ignore t2_starter.py --out /tmp/ve-b2-t2run
SYNTHETIC data: 500 genes, stages E8.25=1200 cells, E8.75=2018 cells, E9.5=3070 cells, E8.5=1626 cells, E10.5=4309 cells -> /tmp/ve-b2-t2run

SYNTHETIC scores (copy_last floor = 50, half-vs-half ceiling = 100; NOT predictive of the real board)

  board extrap (target E10.5, anchor stage E9.5)
    method           valid headline |   de_score de_directi      mmd_u  variogram   d2_shape occupancy_ scale_log_ neighborho
    copy_last         True    50.00 |      50.00      50.00      50.00      50.00      50.00      50.00      50.00      50.00
    grow_and_shift    True    66.79 |      65.22      81.13      61.61      68.87      57.00      81.37      60.37      63.21
    same prediction, randomly rotated + translated: headline 66.79
    same prediction, all-zero coordinates: valid=True headline 41.56 warnings=1

  board interp (target E8.5, anchor stage E8.25)
    method           valid headline |   de_score de_directi      mmd_u  variogram   d2_shape occupancy_ scale_log_ neighborho
    copy_last         True    50.00 |      50.00      50.00      50.00      50.00      50.00      50.00      50.00      50.00
    grow_and_shift    True    71.91 |      73.33      67.96      58.85      58.79      86.29      98.31      86.79      67.72
    same prediction, randomly rotated + translated: headline 71.91
    same prediction, all-zero coordinates: valid=True headline 34.60 warnings=1

Silent killers:
  no_coords      -> valid=False errors=['obsm["spatial_3D"] missing: expression without coordinates is rejected']
  centred_X      -> valid=False errors=['.X has negative values (centred or scaled upstream?)']
  shuffled_genes -> valid=False errors=['var_names: right genes, WRONG ORDER (rejected just the same)']
elapsed 22.43 s, max RSS 1306700 KB
```

## Public panel gene names (T2__heart__val_extrap.genes.txt, 500 names, public metadata)
Same seed, so the same numbers: only the gene names change. This is a format test of the full-width panel.
```
$ taskset -c 0-2 uv run python -W ignore t2_starter.py --genes-file T2__heart__val_extrap.genes.txt --out /tmp/ve-b2-t2pub
SYNTHETIC data: 500 genes, stages E8.25=1200 cells, E8.75=2018 cells, E9.5=3070 cells, E8.5=1626 cells, E10.5=4309 cells -> /tmp/ve-b2-t2pub
SYNTHETIC scores (copy_last floor = 50, half-vs-half ceiling = 100; NOT predictive of the real board)
  board extrap (target E10.5, anchor stage E9.5)
    copy_last         True    50.00 |      50.00      50.00      50.00      50.00      50.00      50.00      50.00      50.00
    grow_and_shift    True    66.79 |      65.22      81.13      61.61      68.87      57.00      81.37      60.37      63.21
  board interp (target E8.5, anchor stage E8.25)
    copy_last         True    50.00 |      50.00      50.00      50.00      50.00      50.00      50.00      50.00      50.00
    grow_and_shift    True    71.91 |      73.33      67.96      58.85      58.79      86.29      98.31      86.79      67.72
```

## Tests
```
$ taskset -c 0-2 uv run pytest -q
9 passed, 107 warnings in 20.08s
```
(The 107 warnings are anndata/pandas/sklearn deprecation notices from dependencies, not test problems.)
