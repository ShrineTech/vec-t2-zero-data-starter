"""Tests for t2_starter (SYNTHETIC data only). Needs numpy, scipy, anndata, veckit==0.1.2.

    python -m pytest -q test_t2_starter.py
"""
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import anndata as ad  # noqa: E402
import t2_starter as st  # noqa: E402


class TestStarter(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.res = st.run(cls.tmp.name, n0=1000, quiet=True)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_copy_last_is_exactly_the_floor(self):
        for board in ("extrap", "interp"):
            self.assertAlmostEqual(self.res[board]["copy_last"]["headline"], 50.0, places=6, msg=board)

    def test_grow_and_shift_beats_floor_on_synthetic(self):
        for board in ("extrap", "interp"):
            self.assertGreater(self.res[board]["grow_and_shift"]["headline"], 55.0, msg=board)

    def test_frame_does_not_matter(self):
        for board in ("extrap", "interp"):
            r = self.res[board]
            self.assertAlmostEqual(r["rotated_headline"], r["grow_and_shift"]["headline"], delta=2.0, msg=board)

    def test_placeholder_coords_pass_format_but_warn_and_lose(self):
        for board in ("extrap", "interp"):
            z = self.res[board]["zero_coords"]
            self.assertTrue(z["check"]["ok"])
            self.assertTrue(any("placeholder" in w for w in z["check"]["warnings"]))
            self.assertLess(z["headline"], self.res[board]["grow_and_shift"]["headline"] - 10)

    def test_predictions_valid_and_readable(self):
        for board in ("extrap", "interp"):
            for n in ("copy_last", "grow_and_shift"):
                r = self.res[board][n]
                self.assertTrue(r["valid"], r["errors"])
                self.assertEqual(r["warnings"], [])
                a = ad.read_h5ad(r["file"])
                self.assertEqual(a.X.dtype.name, "float32")
                self.assertEqual(a.n_vars, 500)
                self.assertEqual(a.obsm["spatial_3D"].shape, (a.n_obs, 3))
                self.assertEqual(a.obsm["spatial_3D"].dtype.name, "float32")

    def test_silent_killers_caught(self):
        k = self.res["_killers"]
        self.assertFalse(k["no_coords"]["ok"])
        self.assertIn("spatial_3D", k["no_coords"]["errors"][0])
        self.assertFalse(k["centred_X"]["ok"])
        self.assertIn("negative", k["centred_X"]["errors"][0])
        self.assertFalse(k["shuffled_genes"]["ok"])
        self.assertIn("WRONG ORDER", k["shuffled_genes"]["errors"][0])

    def test_data_marked_synthetic_and_frames_unregistered(self):
        a = ad.read_h5ad(Path(self.tmp.name) / "released" / "E8.75.h5ad")
        b = ad.read_h5ad(Path(self.tmp.name) / "released" / "E9.5.h5ad")
        self.assertIn("SYNTHETIC", a.uns)
        self.assertTrue(all(n.startswith("SYN_") for n in a.obs_names))
        # each stage sits in its own frame: centroids far apart, principal axes not aligned
        self.assertGreater(np.linalg.norm(a.obsm["spatial_3D"].mean(0) - b.obsm["spatial_3D"].mean(0)), 50)

    def test_skill_map(self):
        self.assertEqual(st.skill(0.02455, 0.02455, 0.00011, True), 0.5)
        self.assertEqual(st.skill(0.00011, 0.02455, 0.00011, True), 1.0)
        self.assertEqual(st.skill(None, 0, 1, False), 0.0)
        self.assertEqual(st.skill(float("nan"), 0, 1, False), 0.0)
        # two-sided: -0.268 floor and +0.268 prediction are the same distance from ideal
        self.assertEqual(st.skill(0.268, -0.268, 0.0074, True, two_sided=True), 0.5)
        self.assertAlmostEqual(sum(st.T2_WEIGHTS.values()), 1.0)


class TestNotebookInSync(unittest.TestCase):
    def test_ipynb_matches_script(self):
        import json
        import make_notebook
        self.assertEqual(json.loads((HERE / "t2_starter.ipynb").read_text()), make_notebook.build(),
                         "run `python make_notebook.py` after editing t2_starter.py")


if __name__ == "__main__":
    unittest.main()
