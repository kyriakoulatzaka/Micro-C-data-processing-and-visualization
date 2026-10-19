import unittest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from convert_contacts import convert_lines, match_chrom


class TestCoolToJuicer(unittest.TestCase):
    def test_add_prefix_only_when_needed(self):
        valid = {"chr3R", "chr2L"}
        self.assertEqual(match_chrom("3R", valid), "chr3R")
        self.assertEqual(match_chrom("chr3R", valid), "chr3R")

    def test_remove_prefix_only_when_needed(self):
        self.assertEqual(match_chrom("chr3R", {"3R"}), "3R")

    def test_zero_based_and_cis_only(self):
        rows = [
            "3R\t0\t400\t3R\t400\t800\t5\n",
            "3R\t0\t400\t2L\t0\t400\t3\n",
            "3R\t800\t1200\t3R\t1200\t1600\t0\n",
        ]
        self.assertEqual(
            list(convert_lines(rows, {"chr3R", "chr2L"}, True)),
            ["0\tchr3R\t1\t0\t16\tchr3R\t401\t1\t5"],
        )
        self.assertEqual(len(list(convert_lines(rows, {"chr3R", "chr2L"}, False))), 2)

    def test_reject_fractional_count(self):
        with self.assertRaises(ValueError):
            list(convert_lines(["3R\t0\t400\t3R\t400\t800\t0.5"], {"chr3R"}))

    def test_optional_header(self):
        rows = ["chrom1 start1 end1 chrom2 start2 end2 count\n", "chr3R 0 400 chr3R 400 800 2\n"]
        self.assertEqual(len(list(convert_lines(rows, {"chr3R"}))), 1)


if __name__ == "__main__":
    unittest.main()
