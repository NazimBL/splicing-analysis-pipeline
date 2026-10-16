#!/usr/bin/env python3
"""Infer library strandedness from STAR GeneCounts (ReadsPerGene.out.tab).

Column 3 counts read pairs whose read 1 is on the gene's strand, column 4 those
on the opposite strand. f = col3 / (col3 + col4):
    f >= 0.8         fr-secondstrand  (read 1 = RNA strand)
    f <= 0.2         fr-firststrand   (read 1 = reverse of RNA, e.g. dUTP)
    0.4 <= f <= 0.6  fr-unstranded
Anything in between, or samples that disagree, is an error: it points to mixed
library preps or a sample mix-up, and rMATS must get one --libType per run.

Usage: infer_strandedness.py <out.tsv> <ReadsPerGene.out.tab>...
Writes a per-sample table and prints the shared rMATS --libType on stdout.
"""
import os
import sys


def classify(f):
    if f >= 0.8:
        return "fr-secondstrand"
    if f <= 0.2:
        return "fr-firststrand"
    if 0.4 <= f <= 0.6:
        return "fr-unstranded"
    return "ambiguous"


def main(out_tsv, tabs):
    rows = []
    for path in sorted(tabs):
        sample = os.path.basename(path).split(".ReadsPerGene")[0]
        same = opposite = 0
        with open(path) as fh:
            for line in fh:
                if line.startswith("N_"):
                    continue
                f = line.split("\t")
                same += int(f[2])
                opposite += int(f[3])
        if same + opposite == 0:
            sys.exit(f"{sample}: no reads assigned to genes in {path}")
        frac = same / (same + opposite)
        rows.append((sample, same, opposite, frac, classify(frac)))

    with open(out_tsv, "w") as out:
        out.write("sample\tsame_strand\topposite_strand\tfrac_same\tlib_type\n")
        for s, a, b, f, t in rows:
            out.write(f"{s}\t{a}\t{b}\t{f:.4f}\t{t}\n")

    types = {t for *_, t in rows}
    if len(types) != 1 or "ambiguous" in types:
        bad = [f"{s} {f:.3f} {t}" for s, _, _, f, t in rows]
        sys.exit("Strandedness is not consistent across samples:\n  " + "\n  ".join(bad))
    print(types.pop())


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2:])
