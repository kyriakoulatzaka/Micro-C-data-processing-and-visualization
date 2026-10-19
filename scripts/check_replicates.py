#!/usr/bin/env python3
"""Validate compatible Cooler genomic axes before combining count matrices.

Accepts plain .cool paths or cooler URI syntax: sample.mcool::/resolutions/400
"""
import argparse
import sys
import cooler


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected-resolution', type=int, required=True)
    parser.add_argument('inputs', nargs='+', help='Two or more Cooler paths/URIs')
    args = parser.parse_args()
    if args.expected_resolution <= 0 or len(args.inputs) < 2:
        parser.error('Supply a positive resolution and at least two input Coolers')

    reference_bins = None
    reference_assembly = None
    metadata_counts = []
    for i, path in enumerate(args.inputs, 1):
        try:
            c = cooler.Cooler(path)
            if c.binsize != args.expected_resolution:
                raise ValueError(f'Input {i}: expected {args.expected_resolution} bp but got {c.binsize}')
            bins = c.bins()[:][['chrom', 'start', 'end']]
            assembly = c.info.get('genome-assembly')
            if reference_bins is None:
                reference_bins = bins
                reference_assembly = assembly
            elif not reference_bins.equals(bins):
                raise ValueError(f'Input {i}: genomic bins or ordering differ; refusing to merge')
            if reference_assembly and assembly and reference_assembly != assembly:
                raise ValueError(f'Input {i}: assembly metadata mismatch: {reference_assembly} vs {assembly}')
            metadata_counts.append(c.info.get('sum'))
            print(f'Input {i}: bins={c.info.get("nbins")}, nnz={c.info.get("nnz")}, '
                  f'sum={c.info.get("sum", "n/a")}, assembly={assembly or "unrecorded"}')
        except (OSError, ValueError, KeyError) as err:
            parser.exit(1, f'ERROR: {path}: {err}\n')
    print('PASS: all input matrices share identical bins at the requested resolution.')
    if all(v is not None for v in metadata_counts):
        print(f'Expected sum of input contact counts: {sum(metadata_counts)}')
    print('Note: merging does not test biological replicate concordance.')


if __name__ == '__main__':
    main()
