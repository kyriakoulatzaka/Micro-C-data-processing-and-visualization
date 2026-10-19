#!/usr/bin/env python3
"""Inspect Cooler balancing weights or subset a .cool to canonical chromosomes."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

MAIN_CHROMS = ("2L", "2R", "3L", "3R", "X", "4")


def weight_stats(weights):
    """Descriptive weights only; finite values are not a quality guarantee."""
    w = np.asarray(weights, dtype=float)
    finite = w[np.isfinite(w)]
    return {
        "bins": int(w.size),
        "nan": int(np.isnan(w).sum()),
        "inf": int(np.isinf(w).sum()),
        "zero": int((w == 0).sum()),
        "finite": int(np.isfinite(w).sum()),
        "nan_fraction": float(np.isnan(w).mean()) if w.size else float("nan"),
        "finite_min": float(finite.min()) if finite.size else None,
        "finite_max": float(finite.max()) if finite.size else None,
        "finite_mean": float(finite.mean()) if finite.size else None,
    }


def remap_pixel_chunk(pixels, bin_map):
    """Drop contacts to removed bins and renumber both bin IDs."""
    old_a = pixels["bin1_id"].to_numpy(dtype=np.int64)
    old_b = pixels["bin2_id"].to_numpy(dtype=np.int64)
    if ((old_a < 0) | (old_a >= bin_map.size) |
            (old_b < 0) | (old_b >= bin_map.size)).any():
        raise ValueError("Pixel contains an out-of-range bin ID")
    a, b = bin_map[old_a], bin_map[old_b]
    keep = (a >= 0) & (b >= 0)
    result = pixels.loc[keep, ["bin1_id", "bin2_id", "count"]].copy()
    result["bin1_id"] = a[keep]
    result["bin2_id"] = b[keep]
    return result.reset_index(drop=True)


def inspect(uri, window_bins):
    import cooler
    c = cooler.Cooler(uri)
    bins = c.bins()[:]
    print("Input:", uri)
    print("Assembly:", c.info.get("genome-assembly", "not recorded"))
    print("Resolution:", c.binsize)
    print("Bins:", len(bins), "Stored pixels:", c.info.get("nnz", "unknown"))
    if "weight" not in bins:
        print("NO_WEIGHT_COLUMN: calculate weights before balance=True")
        return
    for key, val in weight_stats(bins["weight"]).items():
        print(f"{key}: {val}")
    summary = (bins.assign(nan=bins["weight"].isna())
                  .groupby("chrom", observed=True)["nan"]
                  .agg(nan="sum", bins="count"))
    summary["nan_fraction"] = summary["nan"] / summary["bins"]
    print("\nChromosomes by absolute NaN count:")
    print(summary.sort_values(["nan", "nan_fraction"], ascending=False)
          .head(25).to_string())
    main = [x for x in summary.index if str(x).removeprefix("chr") in MAIN_CHROMS]
    print("\nMain chromosomes:")
    print(summary.loc[main].to_string() if main else "No matching chromosomes")
    if window_bins:
        stop = min(window_bins, int(c.info["nbins"]))
        raw = c.matrix(balance=False)[0:stop, 0:stop]
        bal = c.matrix(balance=True)[0:stop, 0:stop]
        print(f"\nFirst {stop} x {stop} bins:")
        print("Raw finite:", int(np.isfinite(raw).sum()), "sum:", float(np.nansum(raw)))
        print("Balanced finite:", int(np.isfinite(bal).sum()),
              "sum:", float(np.nansum(bal)))
        print("Raw and balanced sums are on different scales.")
    print("A successful balance=True read does not establish normalization quality.")


def subset(uri, output, keep_chroms, chunk_size):
    import cooler
    src = cooler.Cooler(uri)
    if "::" in output:
        raise ValueError("Output should be a plain .cool file path")
    dest = Path(output).expanduser()
    if dest.exists():
        raise FileExistsError(f"Refusing to overwrite {dest}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    bins = src.bins()[:]
    selected = (bins["chrom"].astype(str).str.removeprefix("chr")
                .isin(set(keep_chroms)).to_numpy())
    if not selected.any():
        raise ValueError("No bins matched requested chromosomes")
    bin_map = np.full(len(bins), -1, dtype=np.int64)
    bin_map[selected] = np.arange(selected.sum(), dtype=np.int64)
    # Drop old balancing weights. Rebalance new file separately, if desired.
    out_bins = bins.loc[selected, ["chrom", "start", "end"]].reset_index(drop=True)
    print(f"Keeping {selected.sum()}/{len(selected)} bins")
    nnz = int(src.info["nnz"])

    def pixel_chunks():
        kept = 0
        for start in range(0, nnz, chunk_size):
            original = src.pixels()[start:min(start + chunk_size, nnz)]
            filtered = remap_pixel_chunk(original, bin_map)
            kept += len(filtered)
            if len(filtered):
                yield filtered
        print("Kept pixels:", kept, "out of", nnz)

    cooler.create_cooler(str(dest), bins=out_bins,
                         pixels=pixel_chunks(), ordered=True)
    print("Wrote:", dest)
    print("Optional: cooler balance", dest)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="task", required=True)
    q = sub.add_parser("inspect", help="Summarize NaN weights and matrices")
    q.add_argument("cool")
    q.add_argument("--window-bins", type=int, default=0)
    s = sub.add_parser("subset", help="Subset chromosomes with renumbered bin IDs")
    s.add_argument("cool")
    s.add_argument("output")
    s.add_argument("--keep", nargs="+", default=list(MAIN_CHROMS))
    s.add_argument("--chunk-size", type=int, default=1_000_000)
    a = p.parse_args()
    if a.task == "inspect":
        if a.window_bins < 0:
            p.error("--window-bins must be nonnegative")
        inspect(a.cool, a.window_bins)
    else:
        if a.chunk_size < 1:
            p.error("--chunk-size must be positive")
        subset(a.cool, a.output, a.keep, a.chunk_size)


if __name__ == "__main__":
    main()
