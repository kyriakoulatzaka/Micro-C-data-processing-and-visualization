#!/usr/bin/env python3
"""Print a few metadata fields before preparing a .hic file."""
import argparse
import cooler
import h5py


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cool", help="Input single-resolution .cool file")
    parser.add_argument("--expected-resolution", type=int, help="Stop if bin size does not match")
    args = parser.parse_args()

    c = cooler.Cooler(args.cool)
    print("Input:", args.cool)
    print("Assembly:", c.info.get("genome-assembly", "not recorded"))
    print("Bin size (bp):", c.binsize)
    print("Number of pixels:", c.info.get("nnz", "not recorded"))
    print("Number of bins:", c.info.get("nbins", "not recorded"))
    with h5py.File(args.cool, "r") as handle:
        print("Existing /bins/weight:", "weight" in handle["bins"])
    if args.expected_resolution is not None and c.binsize != args.expected_resolution:
        parser.error(f"Input bin size is {c.binsize}; expected {args.expected_resolution}")
    if c.binsize is None:
        parser.error("Variable-size bins are not supported by this pipeline")


if __name__ == "__main__":
    main()
