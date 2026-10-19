#!/usr/bin/env python3
"""Generate a locus-specific Micro-C contact heatmap from a .cool/.mcool URI.

The default visual style follows the user's original ems heatmap example:
YlOrBr colors, white zero/missing pixels, logarithmic scale, no axis ticks.
Raw counts are used by default; balancing is explicitly opt-in.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np


def resolve_region(region: str, chromnames) -> str:
    """Accept 3R:... or chr3R:... when the matrix uses the other naming style."""
    chroms = set(chromnames)
    if ":" not in region:
        chromosome, suffix = region, ""
    else:
        chromosome, rest = region.split(":", 1)
        suffix = ":" + rest
    if chromosome in chroms:
        return region
    alternative = chromosome[3:] if chromosome.startswith("chr") else "chr" + chromosome
    if alternative in chroms:
        return alternative + suffix
    raise ValueError(f"Chromosome {chromosome!r} not found in input. Available: {', '.join(sorted(chroms)[:15])}")


def prepare_log_heatmap(matrix: np.ndarray, *, vmin=None, vmax=None):
    """Mask zeros/invalid values and choose positive log color limits.

    Without explicit limits, the lower bound is the smallest positive value,
    and the upper bound is the 99th percentile (or higher if required).
    """
    cleaned = np.asarray(matrix, dtype=float).copy()
    if cleaned.ndim != 2 or cleaned.size == 0:
        raise ValueError("Expected a non-empty 2D contact matrix")
    cleaned[~np.isfinite(cleaned) | (cleaned <= 0)] = np.nan
    positives = cleaned[np.isfinite(cleaned)]
    if not positives.size:
        raise ValueError("Region contains no positive contacts to plot")
    low = float(positives.min()) if vmin is None else float(vmin)
    if vmax is None:
        high = float(np.percentile(positives, 99))
        # A uniform matrix still needs a strictly positive log range.
        high = max(high, low * 1.1)
    else:
        high = float(vmax)
    if not np.isfinite(low) or not np.isfinite(high) or low <= 0 or high <= low:
        raise ValueError(f"Invalid log color bounds: vmin={low}, vmax={high}; require 0 < vmin < vmax")
    return cleaned, low, high


def draw_heatmap(matrix, output: str, *, region: str, balanced: bool = False,
                 vmin=None, vmax=None, cmap: str = "YlOrBr", dpi: int = 300,
                 title: str | None = None):
    """Render a static heatmap; returns the color scale bounds."""
    import matplotlib
    matplotlib.use("Agg")  # works in headless WSL, Linux and CI
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.ticker import LogLocator, LogFormatterMathtext

    values, low, high = prepare_log_heatmap(matrix, vmin=vmin, vmax=vmax)
    colormap = plt.get_cmap(cmap).copy()
    colormap.set_bad("white")
    colormap.set_under("white")
    fig, ax = plt.subplots(figsize=(6, 5), dpi=200)
    try:
        im = ax.imshow(values, origin="lower", interpolation="nearest",
                       aspect="equal", cmap=colormap,
                       norm=LogNorm(vmin=low, vmax=high))
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        if title:
            ax.set_title(title, fontsize=11)
        colorbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        colorbar.set_label("Balanced contact signal" if balanced else "Raw contact count")
        colorbar.locator = LogLocator(base=10)
        colorbar.formatter = LogFormatterMathtext(base=10)
        colorbar.update_ticks()
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.06)
    finally:
        plt.close(fig)
    return low, high


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cool", required=True,
                        help="Input single-resolution .cool or .mcool::/resolutions/400 URI")
    parser.add_argument("--region", required=True,
                        help="Genomic region, e.g. 3R:13778677-14027796")
    parser.add_argument("--out", required=True,
                        help="Output image (.png, .pdf or .svg)")
    parser.add_argument("--balanced", action="store_true",
                        help="Use existing Cooler balancing weights (default: raw counts)")
    parser.add_argument("--vmin", type=float,
                        help="Fixed positive logarithmic minimum for comparable figures")
    parser.add_argument("--vmax", type=float,
                        help="Fixed logarithmic maximum for comparable figures")
    parser.add_argument("--cmap", default="YlOrBr", help="Matplotlib colormap (default: YlOrBr)")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--title", default=None, help="Optional figure title (default: none)")
    parser.add_argument("--max-bins", type=int, default=2500,
                        help="Safety limit for each axis of the dense matrix (default: 2500)")
    args = parser.parse_args(argv)
    if args.dpi <= 0 or args.max_bins <= 0:
        parser.error("--dpi and --max-bins must be positive integers")
    if Path(args.out).suffix.lower() not in {".png", ".pdf", ".svg"}:
        parser.error("--out must end in .png, .pdf or .svg")

    try:
        import cooler
        c = cooler.Cooler(args.cool)
        region = resolve_region(args.region, c.chromnames)
        lo, hi = c.extent(region)
        if hi - lo > args.max_bins:
            raise ValueError(
                f"Selected {hi-lo} bins per axis (> --max-bins {args.max_bins}); "
                "reduce the region or explicitly raise the limit if enough RAM is available"
            )
        try:
            matrix = c.matrix(balance=args.balanced, sparse=False).fetch(region)
        except (ValueError, KeyError) as exc:
            if args.balanced:
                raise ValueError("Could not read a balanced matrix; check that /bins/weight exists and is valid") from exc
            raise
        low, high = draw_heatmap(matrix, args.out, region=region,
                                 balanced=args.balanced, vmin=args.vmin,
                                 vmax=args.vmax, cmap=args.cmap,
                                 dpi=args.dpi, title=args.title)
    except (OSError, ValueError, KeyError, RuntimeError) as exc:
        parser.exit(1, f"ERROR: {exc}\n")
    print(f"Saved: {args.out}")
    print(f"Region: {region} | Bins: {hi-lo} | Mode: {'balanced' if args.balanced else 'raw counts'}")
    print(f"Log color scale: vmin={low:g}, vmax={high:g}")


if __name__ == "__main__":
    main()
