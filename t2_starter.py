# %% [markdown]
# # Virtual Embryo Challenge, Task 2: zero-data starter for 3D spatial prediction (SYNTHETIC data only)
#
# Build and test a whole Task 2 pipeline **before** you register or download anything:
#
# 1. generate synthetic 3D MERFISH-shaped embryos in the documented T2 file layout (`.X` on a 500-gene panel
#    plus `obsm["spatial_3D"]`), each stage in its **own local frame**, as the real release is;
# 2. run two baselines on both kinds of Task 2 board: `copy_last` (the official floor) and `grow_and_shift`
#    (a simple size, cell-count and composition-aware baseline), for **extrapolation** and **interpolation**;
# 3. validate each prediction against the documented contract, and demonstrate the spatial silent killers:
#    missing coordinates, a centred matrix, placeholder coordinates, wrong gene order;
# 4. show that the spatial metrics do not care about your coordinate frame (no registration needed);
# 5. score with the official `veckit` scorer and the documented 0-100 skill map and T2 weights.
#
# **Nothing here is challenge data.** Synthetic numbers say nothing about the real board.
# Runs as a script (`python t2_starter.py`) or cell by cell in Jupyter / Colab / VS Code (`# %%` cells).

# %%
from __future__ import annotations

import argparse
import json
import math
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial.transform import Rotation

# %% [markdown]
# ## 1. Synthetic T2-shaped data
# Documented layout (Data page): AnnData `.h5ad`, one file per stage; `.X` float32 log-normalised on the
# 500-gene MERFISH panel (finite, non-negative); `obs["celltype"]` on training stages; var index = the panel in
# panel order; `obsm["spatial_3D"]` float32 `[n_obs, 3]`, "per-embryo local and is not registered across time
# points". The stages below follow the **heart setting** (train E8.25, E8.75, E9.5; interpolation target E8.5;
# extrapolation target E10.5). The generative story is a toy: an ellipsoid tissue that grows, gains cells, and
# whose spatial cell-type bands shift and drift over time. It mimics the shape of the problem, not the biology.

# %%
TRAIN = {"E8.25": 8.25, "E8.75": 8.75, "E9.5": 9.5}
TARGETS = {"E8.5": 8.5, "E10.5": 10.5}
# board -> (stages the baseline may use, target stage, the stage copy_last resubmits and DE is measured from)
BOARDS = {
    "extrap": (("E8.75", "E9.5"), "E10.5", "E9.5"),
    "interp": (("E8.25", "E8.75"), "E8.5", "E8.25"),
}


def _type_table(rng, G, K):
    base = rng.normal(-1.0, 1.2, size=G)
    prof, drift = [], []
    for _ in range(K):
        p = base.copy()
        prog = rng.choice(G, size=max(5, G // 12), replace=False)
        p[prog] += rng.normal(2.0, 0.5, size=prog.size)
        d = np.zeros(G)
        dg = rng.choice(G, size=max(5, G // 10), replace=False)
        d[dg] = rng.normal(0, 1.2, size=dg.size)
        prof.append(p); drift.append(d)
    grad = np.zeros(G)                                   # a few genes follow a spatial gradient
    gg = rng.choice(G, size=max(3, G // 25), replace=False)
    grad[gg] = rng.normal(0, 1.0, size=gg.size)
    return np.array(prof), np.array(drift), grad, [f"SynType{k:02d}" for k in range(K)]


def make_stage(rng, t, tbl, genes, n0=1200):
    prof, drift, grad, names = tbl
    K = len(names)
    dt = 2.5 * math.log1p((t - 8.25) / 0.5)              # a developmental clock: change is fastest early
    n = int(n0 * math.exp(0.3 * dt))                     # proliferation: more cells later
    size = 40.0 * math.exp(0.2 * dt)                     # growth: bigger tissue later
    axes = size * np.array([1.4 * math.exp(0.15 * dt), 1.0, 0.8 * math.exp(-0.1 * dt)])   # and it elongates
    u = rng.normal(size=(n, 3)); u /= np.linalg.norm(u, axis=1, keepdims=True)
    C = u * rng.uniform(0, 1, size=(n, 1)) ** (1 / 3) * axes          # uniform inside the ellipsoid
    s = C[:, 0] / axes[0]                                              # position along the long axis, -1..1
    edges = np.linspace(-1, 1, K + 1)[1:-1] + 0.4 * np.tanh(dt / 2)    # bands move with time: composition turns over
    ct = np.clip(np.searchsorted(edges, s + rng.normal(0, 0.08, size=n)), 0, K - 1)
    lam = prof[ct] + drift[ct] * dt + np.outer(s, grad) + rng.normal(0, 0.3, size=(n, len(genes)))
    lam = np.exp(lam - lam.max(1, keepdims=True)); lam /= lam.sum(1, keepdims=True)
    lib = np.exp(rng.normal(np.log(600), 0.4, size=(n, 1)))           # MERFISH-like low counts per cell
    counts = rng.poisson(lam * lib).astype(np.float32)
    tot = counts.sum(1, keepdims=True); tot[tot == 0] = 1
    X = np.log1p(1e3 * counts / tot).astype(np.float32)
    frame = Rotation.random(random_state=int(rng.integers(1 << 31)))   # per-embryo local frame, not registered
    C = (frame.apply(C) + rng.normal(0, 500, size=3)).astype(np.float32)
    obs = pd.DataFrame({"celltype": pd.Categorical([names[i] for i in ct])},
                       index=[f"SYN_{t}_{i}" for i in range(n)])
    A = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=list(genes)))
    A.obsm["spatial_3D"] = C
    A.uns["SYNTHETIC"] = "generated by t2_starter.py; contains no challenge data"
    return A


def generate(n_genes=500, n_types=8, n0=1200, seed=7, genes=None) -> dict:
    rng = np.random.default_rng(seed)
    genes = list(genes) if genes is not None else [f"Gsyn{i:03d}" for i in range(n_genes)]
    tbl = _type_table(rng, len(genes), n_types)
    return {s: make_stage(rng, t, tbl, genes, n0) for s, t in {**TRAIN, **TARGETS}.items()}


# %% [markdown]
# ## 2. Baselines
# `copy_last` is the official floor ("copy_last for Tasks 1 and 2"): resubmit an observed stage verbatim, which
# scores exactly 50 by construction. `grow_and_shift` uses only frame-free quantities, so it never needs to
# register stages: it fits the log cell count and the log of each **principal-axis spread** (measured in each
# stage's own PCA frame) over time, resamples the anchor stage's cells uniformly to the target count, stretches
# the anchor's own coordinates along its own principal axes, and moves each cell type's mean expression part of
# the way along its trend. (Resampling cells *by cell type* to fix composition looks natural but moves mass
# around the tissue and costs shape points: try it.)
# Extrapolation fits the last two stages forward; interpolation goes between the bracketing stages.

# %%
def _dense(X):
    return (X.toarray() if sparse.issparse(X) else np.asarray(X)).astype(np.float32, copy=False)


def rms_radius(C):
    C = np.asarray(C, float)
    return float(np.sqrt(((C - C.mean(0)) ** 2).sum(1).mean()))


def copy_last(a, b, *, anchor, **_):
    src = {"a": a, "b": b}[anchor]
    out = ad.AnnData(X=_dense(src.X).copy(), var=src.var.copy())
    out.obsm["spatial_3D"] = np.asarray(src.obsm["spatial_3D"], np.float32).copy()
    return out


def principal_sd(C):
    """Sorted principal standard deviations and axes of a point cloud: frame-free, so comparable across stages."""
    C = np.asarray(C, float); C = C - C.mean(0)
    w, V = np.linalg.eigh(np.cov(C.T))
    order = np.argsort(w)[::-1]
    return np.sqrt(np.maximum(w[order], 1e-12)), V[:, order]


def grow_and_shift(a, b, *, t_a, t_b, t_target, anchor, max_cells=None, shift_damp=0.6, seed=0):
    """a, b: the two stages used (time t_a < t_b); anchor: 'a' or 'b', whose cells and frame are reused."""
    from scipy.spatial import cKDTree
    rng = np.random.default_rng(seed)
    w = (t_target - t_a) / (t_b - t_a)                    # 0..1 interpolates, >1 extrapolates
    src = {"a": a, "b": b}[anchor]
    t_src = t_a if anchor == "a" else t_b
    # cell count: log-linear in time
    n = int(round(math.exp(math.log(a.n_obs) + w * (math.log(b.n_obs) - math.log(a.n_obs)))))
    if max_cells:
        n = min(n, max_cells)
    # shape and size: each principal axis grows log-linearly, measured in each stage's OWN frame
    sd_a, _ = principal_sd(a.obsm["spatial_3D"])
    sd_b, _ = principal_sd(b.obsm["spatial_3D"])
    sd_t = np.exp(np.log(sd_a) + w * (np.log(sd_b) - np.log(sd_a)))
    Cs = np.asarray(src.obsm["spatial_3D"], float)
    sd_s, V = principal_sd(Cs)
    mu = Cs.mean(0)
    # resample cells uniformly (keeps the anchor's spatial density); duplicates get a sub-neighbour offset
    idx = rng.choice(src.n_obs, size=n, replace=n > src.n_obs)
    nn = float(np.median(cKDTree(Cs).query(Cs, k=2)[0][:, 1]))
    C = Cs[idx] + rng.normal(0, 0.5 * nn, size=(n, 3))
    C = mu + ((C - mu) @ V) * (sd_t / sd_s) @ V.T
    # expression: each cell type's mean moves part of the way along its own trend between the two stages
    Xs = _dense(src.X)
    lab = np.asarray(src.obs["celltype"]).astype(str)
    Xa, Xb = _dense(a.X), _dense(b.X)
    la, lb = np.asarray(a.obs["celltype"]).astype(str), np.asarray(b.obs["celltype"]).astype(str)
    X, lab_out = Xs[idx].copy(), lab[idx]
    for k in set(lab_out):
        if (la == k).any() and (lb == k).any():
            slope = (Xb[lb == k].mean(0) - Xa[la == k].mean(0)) / (t_b - t_a)
            X[lab_out == k] += shift_damp * slope * (t_target - t_src)
    out = ad.AnnData(X=np.clip(X, 0, None).astype(np.float32), var=src.var.copy())
    out.obsm["spatial_3D"] = C.astype(np.float32)
    return out


# %% [markdown]
# ## 3. Validate against the documented contract
# What the Data page and the scorer require: var index = the board's panel **in order** (the right genes in the
# wrong order are "rejected just the same"); `.X` float32, finite, **non-negative** ("a negative entry usually
# means the matrix was centred somewhere upstream"); `obsm["spatial_3D"]` present, `[n_obs, >=3]`, finite ("A file
# with expression but no coordinates is rejected before scoring"); the cell count inside the board's limits
# (from `panels/index.json`, e.g. T2:heart:val_extrap allows 1,000-25,179). Two warnings catch files that pass
# but score badly: coordinates that are a placeholder (all equal / collapsed), and `.X` that looks like raw counts.

# %%
@dataclass
class Check:
    ok: bool = True
    errors: list = field(default_factory=list)
    warnings: list = field(default_factory=list)

    def err(self, m):
        self.ok = False
        self.errors.append(m)


def validate(a, genes, *, min_cells, max_cells) -> Check:
    c = Check()
    got = [str(g) for g in a.var_names]
    if got != list(genes):
        same_set = set(got) == set(genes)
        c.err("var_names: " + ("right genes, WRONG ORDER (rejected just the same)" if same_set else
                               f"expected {len(genes)} panel genes, got {len(got)}"))
    vals = a.X.data if sparse.issparse(a.X) else np.asarray(a.X)
    if np.asarray(vals).dtype != np.float32:
        c.err(f".X dtype {np.asarray(vals).dtype}, expected float32")
    if not np.isfinite(vals).all():
        c.err(".X has NaN/inf")
    if (vals < 0).any():
        c.err(".X has negative values (centred or scaled upstream?)")
    if not (min_cells <= a.n_obs <= max_cells):
        c.err(f"n_obs={a.n_obs} outside [{min_cells}, {max_cells}]")
    if "spatial_3D" not in a.obsm:
        c.err('obsm["spatial_3D"] missing: expression without coordinates is rejected')
    else:
        C = np.asarray(a.obsm["spatial_3D"])
        if C.ndim != 2 or C.shape[0] != a.n_obs or C.shape[1] < 3:
            c.err(f'obsm["spatial_3D"] shape {C.shape}, expected [{a.n_obs}, 3]')
        elif not np.isfinite(C[:, :3]).all():
            c.err('obsm["spatial_3D"] has NaN/inf')
        else:
            spread = np.ptp(C[:, :3], axis=0)
            if (spread < 1e-6).sum() >= 2 or len(np.unique(np.round(C[:, :3], 6), axis=0)) < 0.5 * a.n_obs:
                c.warnings.append('obsm["spatial_3D"] looks like a placeholder (collapsed or mostly duplicate '
                                  "points). This PASSES format checks but wrecks every shape and local metric.")
    d = _dense(a.X[: min(500, a.n_obs)])
    if d.size and (d.max() > 30 or np.allclose(d, np.round(d))):
        c.warnings.append(".X looks like raw counts, not log-normalised. This passes format checks but scores wrongly.")
    return c


def write_submission(a, path):
    out = ad.AnnData(X=np.asarray(_dense(a.X), dtype=np.float32), var=pd.DataFrame(index=list(a.var_names)))
    out.obs_names = [f"pred_{i}" for i in range(out.n_obs)]
    if "spatial_3D" in a.obsm:
        out.obsm["spatial_3D"] = np.asarray(a.obsm["spatial_3D"], dtype=np.float32)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    out.write_h5ad(path, compression="gzip")


# %% [markdown]
# ## 4. Score: the official veckit scorer + the documented skill map
# Raw metrics come from `veckit.score(task="T2", ...)` itself (`pip install veckit==0.1.2`, MIT): the same entry
# point as the `veckit` CLI. Anchors are computed on the synthetic target exactly as documented: floor =
# copy_last, ceiling = one half of the target scored against the other, disjoint half. skill =
# min(d_floor / (d_floor + d), 1); `scale_log_ratio` is two-sided (absolute values, anchors included).
# T2 weights from the scoring page: expression change 25% (DES .5, DCS .5), cell-state distribution 25%
# (MMD .6, variogram .4), tissue shape and growth 25% (d2_shape, occupancy_dice, scale_log_ratio, 1/3 each),
# local spatial organisation 25% (neighborhood_mmd).
# For real boards use the published anchors in `panels/index.json` instead.

# %%
T2_WEIGHTS = {"de_score": .125, "de_direction": .125, "mmd_u": .15, "variogram": .10,
              "d2_shape": .25 / 3, "occupancy_dice": .25 / 3, "scale_log_ratio": .25 / 3, "neighborhood_mmd": .25}
LOWER_IS_BETTER = {"de_score": False, "de_direction": False, "mmd_u": True, "variogram": True,
                   "d2_shape": True, "occupancy_dice": False, "scale_log_ratio": True, "neighborhood_mmd": True}
TWO_SIDED = {"scale_log_ratio"}


def skill(v, f, c, lower, two_sided=False):
    if v is None or not all(map(math.isfinite, (v, f, c))):
        return 0.0
    if two_sided:
        v, f, c = abs(v), abs(f), abs(c)
    if abs(c - f) < 1e-12:
        return 0.0
    d_floor, d = abs(c - f), (v - c) if lower else (c - v)
    return 1.0 if d_floor + d <= 1e-9 else min(d_floor / (d_floor + d), 1.0)


def headline(raw, anchors):
    sk = {m: skill(raw.get(m), anchors[m]["floor"], anchors[m]["ceiling"], LOWER_IS_BETTER[m], m in TWO_SIDED)
          for m in T2_WEIGHTS}
    return 100 * sum(T2_WEIGHTS[m] * sk[m] for m in T2_WEIGHTS), {m: 100 * s for m, s in sk.items()}


def veckit_raw(pred_path, truth_path, ref_path, seed=0):
    from veckit import score
    return score(task="T2", input=pred_path, target=truth_path, reference=ref_path, seed=seed)["metrics"]


def split_target(target, out_dir, seed=0):
    """Write the two disjoint halves of the target; half A is the truth, half B is the ceiling 'prediction'."""
    from common.core_metrics import split_half
    a, b = split_half(target.n_obs, seed=seed)
    pa, pb = Path(out_dir) / "truth_halfA.h5ad", Path(out_dir) / "ceiling_halfB.h5ad"
    target[a].copy().write_h5ad(pa)
    write_submission(target[b].copy(), pb)
    return pa, pb


# %% [markdown]
# ## 5. Run it

# %%
def run(out_dir, n_genes=500, n0=1200, genes_file=None, seed=7, quiet=False):
    log = (lambda *a: None) if quiet else print
    out = Path(out_dir)
    genes = [g.strip() for g in Path(genes_file).read_text().splitlines() if g.strip()] if genes_file else None
    data = generate(n_genes=n_genes, n0=n0, seed=seed, genes=genes)
    genes = list(data["E8.25"].var_names)
    for s, A in data.items():
        sub = "released" if s in TRAIN else "hidden"
        (out / sub).mkdir(parents=True, exist_ok=True)
        A.write_h5ad(out / sub / f"{s}.h5ad", compression="gzip")
    log(f"SYNTHETIC data: {len(genes)} genes, stages " +
        ", ".join(f"{s}={A.n_obs} cells" for s, A in data.items()) + f" -> {out}")
    results = {}
    for board, ((sa, sb), st, sref) in BOARDS.items():
        a, b, target = data[sa], data[sb], data[st]
        anchor = "b" if board == "extrap" else "a"
        bdir = out / board
        bdir.mkdir(parents=True, exist_ok=True)
        truth_p, ceil_p = split_target(target, bdir, seed=seed)
        ref_p = out / "released" / f"{sref}.h5ad"
        lo, hi = min(1000, target.n_obs // 4), 3 * max(target.n_obs, data[sref].n_obs)
        kw = dict(t_a=TRAIN[sa], t_b=TRAIN[sb], t_target=TARGETS[st], anchor=anchor)
        preds = {"copy_last": copy_last(a, b, **kw), "grow_and_shift": grow_and_shift(a, b, max_cells=hi, **kw)}
        raws = {}
        for name, pred in preds.items():
            p = bdir / f"{name}.h5ad"
            write_submission(pred, p)
            raws[name] = (veckit_raw(p, truth_p, ref_p), validate(pred, genes, min_cells=lo, max_cells=hi), p)
        floor_raw = raws["copy_last"][0]
        ceil_raw = veckit_raw(ceil_p, truth_p, ref_p)
        anchors = {m: {"floor": floor_raw[m], "ceiling": ceil_raw[m]} for m in T2_WEIGHTS}
        br = {"target": st, "anchor_stage": {"a": sa, "b": sb}[anchor], "anchors": anchors}
        for name, (raw, chk, p) in raws.items():
            h, sk = headline(raw, anchors)
            br[name] = {"valid": chk.ok, "errors": chk.errors, "warnings": chk.warnings, "headline": h,
                        "skill_x100": sk, "raw": {m: raw[m] for m in T2_WEIGHTS}, "file": str(p)}
        # Frame freedom: a random rotation + translation of the same prediction scores the same.
        g = preds["grow_and_shift"]
        moved = ad.AnnData(X=_dense(g.X).copy(), var=g.var.copy())
        rot = Rotation.random(random_state=seed + 1)
        moved.obsm["spatial_3D"] = (rot.apply(np.asarray(g.obsm["spatial_3D"], float)) + 1234.5).astype(np.float32)
        pm = bdir / "grow_and_shift_rotated.h5ad"
        write_submission(moved, pm)
        br["rotated_headline"] = headline(veckit_raw(pm, truth_p, ref_p), anchors)[0]
        # Placeholder coordinates: valid file, collapsed geometry.
        flat = ad.AnnData(X=_dense(g.X).copy(), var=g.var.copy())
        flat.obsm["spatial_3D"] = np.zeros((g.n_obs, 3), np.float32)
        pf = bdir / "grow_and_shift_zero_coords.h5ad"
        write_submission(flat, pf)
        br["zero_coords"] = {"check": vars(validate(flat, genes, min_cells=lo, max_cells=hi)),
                             "headline": headline(veckit_raw(pf, truth_p, ref_p), anchors)[0]}
        results[board] = br
    # The other silent killers, on the extrapolation grow_and_shift prediction:
    a, b = data["E8.75"], data["E9.5"]
    g = grow_and_shift(a, b, t_a=8.75, t_b=9.5, t_target=10.5, anchor="b")
    lo, hi = 1000, 10 ** 6
    no_coords = ad.AnnData(X=_dense(g.X).copy(), var=g.var.copy())
    centred = ad.AnnData(X=(_dense(g.X) - _dense(g.X).mean(0)).astype(np.float32), var=g.var.copy())
    centred.obsm["spatial_3D"] = g.obsm["spatial_3D"]
    shuffled = g[:, list(reversed(genes))].copy()
    results["_killers"] = {
        "no_coords": vars(validate(no_coords, genes, min_cells=lo, max_cells=hi)),
        "centred_X": vars(validate(centred, genes, min_cells=lo, max_cells=hi)),
        "shuffled_genes": vars(validate(shuffled, genes, min_cells=lo, max_cells=hi)),
    }
    log("\nSYNTHETIC scores (copy_last floor = 50, half-vs-half ceiling = 100; NOT predictive of the real board)")
    for board, br in results.items():
        if board.startswith("_"):
            continue
        log(f"\n  board {board} (target {br['target']}, anchor stage {br['anchor_stage']})")
        log(f"    {'method':<16}{'valid':>6}{'headline':>9} | " + " ".join(f"{m[:10]:>10}" for m in T2_WEIGHTS))
        for n in ("copy_last", "grow_and_shift"):
            r = br[n]
            log(f"    {n:<16}{str(r['valid']):>6}{r['headline']:>9.2f} | "
                + " ".join(f"{r['skill_x100'][m]:>10.2f}" for m in T2_WEIGHTS))
        log(f"    same prediction, randomly rotated + translated: headline {br['rotated_headline']:.2f}")
        log(f"    same prediction, all-zero coordinates: valid={br['zero_coords']['check']['ok']} "
            f"headline {br['zero_coords']['headline']:.2f} warnings={len(br['zero_coords']['check']['warnings'])}")
    k = results["_killers"]
    log("\nSilent killers:")
    for n in ("no_coords", "centred_X", "shuffled_genes"):
        log(f"  {n:<15}-> valid={k[n]['ok']} errors={k[n]['errors']}")
    (out / "results.json").write_text(json.dumps(results, indent=2, default=float))
    return results


# %%
if __name__ == "__main__" and "ipykernel" not in sys.modules:
    ap = argparse.ArgumentParser(description="VEC T2 zero-data starter (synthetic only)")
    ap.add_argument("--out", default="synthetic_t2")
    ap.add_argument("--genes", type=int, default=500)
    ap.add_argument("--cells", type=int, default=1200, help="cells in the first stage; later stages grow")
    ap.add_argument("--genes-file", help="public T2 panel (e.g. T2__heart__val_extrap.genes.txt) for a format test")
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    run(a.out, n_genes=a.genes, n0=a.cells, genes_file=a.genes_file, seed=a.seed)
elif "ipykernel" in sys.modules:
    _results = run(tempfile.mkdtemp(prefix="vec_t2_"))
