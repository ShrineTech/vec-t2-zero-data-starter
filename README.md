# T2 zero-data starter: a first valid Task 2 (3D spatial) submission, built on synthetic data

A community tutorial for the [Virtual Embryo Challenge](https://virtualembryo.ai/challenge), Task 2. It is not affiliated with or endorsed by the organisers.

## What it is
A single script, also shipped as a Colab or Jupyter notebook, that takes you from nothing to scored, validated Task 2 predictions on **both kinds of Task 2 board, extrapolation and interpolation**, without registering or downloading the release:

1. It generates a **synthetic** series of 3D MERFISH-shaped embryos in the documented Task 2 layout:
   - `.h5ad` files with `.X` float32 log-normalised on a 500-gene panel, `obs["celltype"]` and `obsm["spatial_3D"]` float32 [n, 3];
   - the heart-setting stage plan: train E8.25, E8.75 and E9.5, interpolation target E8.5, extrapolation target E10.5;
   - a tissue that grows, gains cells and elongates, with spatial cell-type bands that shift over time. This mimics the shape of the problem, not the biology;
   - every stage in its **own** random rotation and translation, because the documented frame is "per-embryo local and is not registered across time points";
   - optionally the real public panel gene names (for example `T2__heart__val_extrap.genes.txt`, public metadata).
2. It runs two baselines on each board:
   - `copy_last`, the official floor;
   - `grow_and_shift`. It fits the cell count and the spread along each principal axis, each measured in the stage's **own** PCA frame, so no registration is needed. It then stretches the anchor stage's coordinates along its own axes, resamples its cells to the target count, and moves each cell type's mean expression along its trend.
3. It validates each prediction against the documented contract and demonstrates the **spatial silent killers**:
   - **Dropping the coordinates.** "A file with expression but no coordinates is rejected before scoring."
   - **A centred matrix.** "Must be finite and non-negative; a negative entry usually means the matrix was centred somewhere upstream."
   - **Placeholder coordinates** (all zeros). These pass every format check, then lose most of the spatial score. The validator warns.
   - **Wrong gene order**, which is rejected even when the genes themselves are right.
4. It shows that **the frame is yours**. The same prediction, randomly rotated and translated, gets the same score, because "every spatial term is invariant to translation and rotation" (Evaluation page). You never need to register your output to the target.
5. It scores with the **official `veckit` scorer** (`veckit.score(task="T2", ...)`), the documented skill map and the T2 weights:
   - expression change 25%;
   - cell-state distribution 25%;
   - tissue shape and growth 25%;
   - local spatial organisation 25%.

   Anchors are computed on the synthetic target exactly as documented: floor = `copy_last`, ceiling = one half of the target against the other half. `scale_log_ratio` is two-sided.

## Why it helps participants
- Among the ~716 contributions on the Community resources page (2026-11-18 register), no starter takes you through a first **Task 2** submission. The Task 2 material that exists covers the scorer (veckit PR #4, synthetic conformance checks), viewers and diagnostics, or needs the real ~1 GB of data.
- Task 2 is the task with coordinates, and that is where first submissions go wrong:
  - a missing `obsm` key;
  - a centred matrix;
  - placeholder coordinates;
  - time spent registering frames that the metrics ignore.
- You can build and debug your loader, writer, validator and scoring loop for both board types before registering. It is CPU only, takes about 25 s and fits Colab.
- The Challenge is annual, so the starter stays useful for the next round after this one closes.

## How to run
```bash
uv sync && uv run python t2_starter.py --out synthetic_t2          # or: pip install anndata scipy pandas veckit==0.1.2
curl -O https://virtualembryo.ai/challenge/panels/T2__heart__val_extrap.genes.txt
uv run python t2_starter.py --genes-file T2__heart__val_extrap.genes.txt   # the real public panel's gene names
uv run pytest -q
```
For Colab, open `t2_starter.ipynb`. Its first code cell installs the dependencies. The notebook is generated from the script by `python make_notebook.py`, and a test keeps the two in sync.

Observed output (2026-11-18, CPU pinned to 3 cores, 500 synthetic genes; about 22 s and 1.3 GB RSS; full log in `RUN_LOG.md`):
```
  board extrap (target E10.5, anchor stage E9.5)
    method           valid headline |   de_score de_directi      mmd_u  variogram   d2_shape occupancy_ scale_log_ neighborho
    copy_last         True    50.00 |      50.00      50.00      50.00      50.00      50.00      50.00      50.00      50.00
    grow_and_shift    True    66.79 |      65.22      81.13      61.61      68.87      57.00      81.37      60.37      63.21
    same prediction, randomly rotated + translated: headline 66.79
    same prediction, all-zero coordinates: valid=True headline 41.56 warnings=1

  board interp (target E8.5, anchor stage E8.25)
    copy_last         True    50.00 |      50.00      50.00      50.00      50.00      50.00      50.00      50.00      50.00
    grow_and_shift    True    71.91 |      73.33      67.96      58.85      58.79      86.29      98.31      86.79      67.72
    same prediction, randomly rotated + translated: headline 71.91
    same prediction, all-zero coordinates: valid=True headline 34.60 warnings=1

Silent killers:
  no_coords      -> valid=False errors=['obsm["spatial_3D"] missing: expression without coordinates is rejected']
  centred_X      -> valid=False errors=['.X has negative values (centred or scaled upstream?)']
  shuffled_genes -> valid=False errors=['var_names: right genes, WRONG ORDER (rejected just the same)']
```

**Synthetic scores say nothing about the real board.** We wrote the generator, so a baseline that models growth is expected to help here. One lesson does carry over: resampling cells **by cell type** to fix the composition sounds natural, but it moves mass around the tissue and costs shape points. The baseline therefore resamples uniformly. Try it both ways.

## Tests (`uv run pytest -q`)
The tests check that:

- `copy_last` scores exactly 50 on both boards;
- `grow_and_shift` scores above 55 on both boards on the synthetic data;
- a random rotation and translation of a prediction leaves its score unchanged;
- all-zero coordinates pass the format check, trigger the placeholder warning and lose at least 10 points;
- both predictions validate and are float32 `.h5ad` files with float32 `obsm["spatial_3D"]` of shape [n, 3];
- the missing-coordinates, centred-matrix and wrong-gene-order killers are caught;
- the data is marked SYNTHETIC and the stages are not registered to each other;
- the skill map's anchor points are correct, including the two-sided `scale_log_ratio`, and the T2 weights sum to 1;
- the notebook is in sync with the script.

## Known limits (read these)
- **Synthetic only, toy biology.** The tissue is an ellipsoid with banded cell types; real hearts and embryos are not.
- **Structure gate.** The Baselines page mentions a Task 2 "structure gate" (shape credit multiplied by neighbourhood-structure skill) that the scoring page does not specify. The headline here is the plain weighted sum, without it.
- **Size metric.** The ranked size term follows the anchor list in `panels/index.json` (`scale_log_ratio`, two-sided). veckit 0.1.2's T2 panel comment says `size_fidelity` "is what ranks". If the board uses `size_fidelity`, change one key in `T2_WEIGHTS`.
- **Settings.** Only the heart-setting stage plan is simulated. The embryo setting (train E6.75, E7.25, E8.0; interpolation targets E7.5 and E7.75) "share[s] this schema exactly", so the same code applies with different stage names. Note that `T2:embryo:val_interp` has **498** genes, not 500.
- **Cell limits.** The synthetic run scales the minimum cell count down. The real boards allow, for example, 1,000-25,179 cells (`T2:heart:val_extrap`) and 583-5,000 (`T2:embryo:val_interp`); see `panels/index.json`.
- **Floor stage.** For interpolation the starter's floor copies the preceding stage. On the real boards, use the published anchors in `panels/index.json`, which are the authority.
- The portal's own format check is the authority (Rules s.11). Run it too.

## Related work
- [ShrineTech/vec-t1-zero-data-starter](https://github.com/ShrineTech/vec-t1-zero-data-starter): the Task 1 counterpart (same layout, same skill map).
- `vec-t3-zero-data-starter`: the Task 3 sibling (knock-outs with two replicates), published alongside this one.
- veckit PR #4 ("synthetic T2 conformance checks") tests the scorer itself. Shashwat Srivastav's Submission Viewer and the random-guy-05 SpatialQC tool inspect predictions you already have. This repo is the participant-side path to a first Task 2 file.

## Data and AI disclosure
- It contains **no challenge data**. Every cell is generated. The only real input is the optional public gene-name list.
- It was written with AI assistance (Claude Code, by Anthropic) and reviewed and filed by a human.
- `veckit` is the organisers' MIT-licensed scorer.
- Licence: MIT.
