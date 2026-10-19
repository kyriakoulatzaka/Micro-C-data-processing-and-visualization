#!/usr/bin/env python3
"""Stream cooler --join pixel rows into Juicer shortWithScore rows.

Input columns: chrom1 start1 end1 chrom2 start2 end2 count
Output columns: str1 chr1 pos1 frag1 str2 chr2 pos2 frag2 score
"""
import argparse
from decimal import Decimal, InvalidOperation
import sys


def read_chrom_names(path):
    names = set()
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 2 and not line.startswith("#"):
                names.add(fields[0])
    if not names:
        raise ValueError("No chromosome names found in chromosome sizes file")
    return names


def match_chrom(name, valid_names):
    # Use exact chromosome names from the supplied genome assembly.
    if name in valid_names:
        return name
    candidates = ((name[3:],) if name.startswith("chr") else ("chr" + name,))
    for candidate in candidates:
        if candidate in valid_names:
            return candidate
    raise ValueError(f"Chromosome {name!r} is not in chromosome sizes (also checked optional 'chr' prefix)")


def convert_lines(lines, chrom_names, cis_only=False):
    for row_num, line in enumerate(lines, 1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.rstrip("\n").split()
        if fields[0] == "chrom1":  # optional cooler --header
            continue
        if len(fields) != 7:
            raise ValueError(f"Line {row_num}: expected 7 cooler --join columns, got {len(fields)}")
        source_chrom1, start1, _end1, source_chrom2, start2, _end2, score_text = fields
        chrom1 = match_chrom(source_chrom1, chrom_names)
        chrom2 = match_chrom(source_chrom2, chrom_names)
        if cis_only and chrom1 != chrom2:
            continue
        try:
            pos1 = int(start1) + 1  # zero-based cooler -> one-based Juicer
            pos2 = int(start2) + 1
            score = Decimal(score_text)
        except (ValueError, InvalidOperation) as exc:
            raise ValueError(f"Line {row_num}: invalid coordinate or count") from exc
        if pos1 < 1 or pos2 < 1:
            raise ValueError(f"Line {row_num}: negative genomic start coordinate")
        if not score.is_finite() or score < 0 or score != score.to_integral_value():
            raise ValueError(f"Line {row_num}: count must be a finite nonnegative integer; received {score_text!r}")
        if score == 0:
            continue
        # The two dummy fragment indices must differ (0 and 1).
        yield f"0\t{chrom1}\t{pos1}\t0\t16\t{chrom2}\t{pos2}\t1\t{int(score)}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chrom-sizes", required=True, help="Genome chrom.sizes matching the intended .hic output")
    parser.add_argument("--cis-only", action="store_true", help="Exclude inter-chromosomal contacts")
    args = parser.parse_args()
    try:
        valid = read_chrom_names(args.chrom_sizes)
        for line in convert_lines(sys.stdin, valid, args.cis_only):
            print(line)
    except (ValueError, OSError) as exc:
        parser.exit(1, f"ERROR: {exc}\n")


if __name__ == "__main__":
    main()
