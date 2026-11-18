# FILING: vec-t2-zero-data-starter (Community Contribution Award, up to $200)

Drafted 2026-11-18 by the fleet. **Nothing has been published, pushed or filed.** This folder is a local git repo on branch `fleet/ve-community-batch2-20261118`, with no remote.

| Item | Value |
|---|---|
| Suggested repo | `ShrineTech/vec-t2-zero-data-starter` (public) |
| Kind (form) | teach |
| Tests (2026-11-18) | 9 passed (fresh `uv sync` env, 3 CPU cores); see `RUN_LOG.md` |
| Overlap risk | **Low to medium.** No Task 2 first-submission starter among the ~716 filed contributions (2026-11-18 register). Nearby: veckit PR #4 (synthetic T2 conformance of the *scorer*), Shashwat Srivastav's Submission Viewer and random-guy-05's SpatialQC (inspect an existing prediction). None generates data, runs baselines and scores a first T2 file end to end without the download. |

## 1. Publish (about 5 min)
1. [ ] Read `README.md` yourself, so that "reviewed" in the AI disclosure stays true.
2. [ ] On GitHub, create an **empty public** repo `ShrineTech/vec-t2-zero-data-starter`: no README, licence or .gitignore, since they are already here.
3. [ ] From this folder, run:
   ```bash
   git branch -m main
   git remote add origin https://github.com/ShrineTech/vec-t2-zero-data-starter.git
   git push -u origin main
   ```
4. [ ] Open the repo page. Check that the README renders and that `t2_starter.ipynb` opens. Optional: open it in Colab and run all (about 1 min on a free CPU runtime; untested by the fleet).
5. [ ] Once the T3 repo is public too, optionally turn the plain-text `vec-t3-zero-data-starter` mention under "Related work" into a link (and the same in the T3 README).

## 2. File (about 5 min)
Go to https://virtualembryo.ai/challenge/community, sign in, and open the Community form. Fill in:

- **kind:** teach
- **title:** Task 2 zero-data starter: a first valid 3D spatial submission on synthetic data, extrapolation and interpolation (script + Colab notebook)
- **url:** https://github.com/ShrineTech/vec-t2-zero-data-starter
- **summary:** This lets you build and test a whole Task 2 pipeline before you register or download anything. It generates synthetic 3D MERFISH-shaped embryos in the documented T2 layout: 500-gene panel (optionally the public T2 gene list), obs celltype, obsm spatial_3D, and every stage in its own unregistered frame, on the heart-setting stage plan. It runs copy_last (the official floor, exactly 50) and a frame-free grow_and_shift baseline on both an extrapolation and an interpolation board. The baseline fits cell count and principal-axis growth in each stage's own PCA frame, so no registration is needed. It demonstrates the spatial silent killers: missing coordinates, a centred matrix, placeholder coordinates (which pass the format check but lose most of the spatial score) and wrong gene order. It shows that a rotated copy of a prediction scores the same. It scores with the official veckit scorer, the documented skill map and the T2 weights, and flags two open points: the unspecified structure gate, and scale_log_ratio vs size_fidelity. CPU only, about 25 s. Companion to our T1 starter and the T3 sibling. Synthetic data only, no challenge data. Written with AI assistance (Claude Code, by Anthropic) under human direction; reviewed, tested and filed by Chuck Gaffney for ShrineTech.
- **reproduce:** `uv sync && uv run python t2_starter.py` (or open `t2_starter.ipynb` in Colab) ; tests: `uv run pytest -q`
- **display_name (Credited as):** ShrineTech

## 3. Known gaps (fleet, honest)
- Only the heart-setting stage plan is simulated. The embryo setting uses the same schema, but its interpolation board has 498 genes.
- The T2 "structure gate" is unspecified by the organisers. The headline is computed without it, as the README states.
- The Colab run is untested. The notebook is generated from the script, and a test checks that the two match.
