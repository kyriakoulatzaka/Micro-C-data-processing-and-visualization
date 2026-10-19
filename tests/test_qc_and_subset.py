"""Small synthetic tests of NaN summarization and bin-ID remapping."""
import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from qc_and_subset import weight_stats, remap_pixel_chunk


class QcSubsetTests(unittest.TestCase):
    def test_weight_statistics(self):
        result = weight_stats([1., float("nan"), 0., float("inf"), 3.])
        self.assertEqual(result["bins"], 5)
        self.assertEqual(result["nan"], 1)
        self.assertEqual(result["inf"], 1)
        self.assertEqual(result["finite"], 3)
        self.assertEqual(result["zero"], 1)
        self.assertEqual(result["finite_mean"], 4 / 3)

    def test_reindex_noncontiguous_bins(self):
        mapping = np.array([0, -1, 1, -1, 2], dtype=np.int64)
        pixels = pd.DataFrame({
            "bin1_id": [0, 0, 1, 2, 2, 4],
            "bin2_id": [0, 2, 4, 2, 4, 4],
            "count": [1, 3, 5, 7, 11, 13],
        })
        result = remap_pixel_chunk(pixels, mapping)
        self.assertEqual(result["bin1_id"].tolist(), [0, 0, 1, 1, 2])
        self.assertEqual(result["bin2_id"].tolist(), [0, 1, 1, 2, 2])
        self.assertEqual(result["count"].tolist(), [1, 3, 7, 11, 13])

    def test_out_of_range(self):
        data = pd.DataFrame({"bin1_id": [0], "bin2_id": [5], "count": [1]})
        with self.assertRaises(ValueError):
            remap_pixel_chunk(data, np.array([0, 1, 2]))


if __name__ == "__main__":
    unittest.main()
