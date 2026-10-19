"""Fast offline tests for the heatmap helper functions (no Cooler file required)."""
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from plot_heatmap import draw_heatmap, prepare_log_heatmap, resolve_region


class TestHeatmap(unittest.TestCase):
    def test_resolve_chromosome_aliases(self):
        self.assertEqual(resolve_region("3R:100-900", ["chr3R"]), "chr3R:100-900")
        self.assertEqual(resolve_region("chr3R:100-900", ["3R"]), "3R:100-900")
        self.assertEqual(resolve_region("3R:100-900", ["3R"]), "3R:100-900")

    def test_masks_zero_and_nan(self):
        mat, low, high = prepare_log_heatmap(np.array([[0, 1], [float('inf'), 9]]))
        self.assertTrue(np.isnan(mat[0, 0]))
        self.assertTrue(np.isnan(mat[1, 0]))
        self.assertEqual(low, 1.0)
        self.assertGreater(high, low)

    def test_empty_contacts(self):
        with self.assertRaisesRegex(ValueError, "no positive"):
            prepare_log_heatmap(np.zeros((3, 3)))

    def test_fixed_bounds(self):
        _mat, low, high = prepare_log_heatmap(np.ones((3, 3)), vmin=0.1, vmax=50)
        self.assertEqual((low, high), (0.1, 50.0))

    def test_invalid_bounds(self):
        with self.assertRaisesRegex(ValueError, "Invalid log color bounds"):
            prepare_log_heatmap(np.ones((3, 3)), vmin=10, vmax=1)

    def test_png_can_be_rendered_with_synthetic_matrix(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "synthetic_example.png"
            bounds = draw_heatmap(np.array([[1, 3, 0], [3, 5, 2], [0, 2, 1]]),
                                  str(out), region="chr3R:0-1200")
            self.assertTrue(out.exists())
            self.assertGreater(out.stat().st_size, 1000)
            self.assertGreater(bounds[1], bounds[0])


if __name__ == "__main__":
    unittest.main()
