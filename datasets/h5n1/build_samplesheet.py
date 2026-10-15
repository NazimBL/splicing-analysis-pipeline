#!/usr/bin/env python3
"""Build datasets/h5n1/samplesheet.csv from the raw FASTQ file names.

File names look like  p24240-s001_1E9-7442-D0_S50_L004_R1_001.fastq.gz
                      <project>-<sample>_<group>-<animal>-<timepoint>_<S#>_<lane>_<read>_001

Rules (see docs/DECISIONS.md, 2026-10-01 and 2026-10-15):
  - keep timepoints D0, D1, D2, D4;
  - groups 1E7/1E9 are rhesus (dose in the name), V/VII are cynomolgus;
  - DVS31 and DV531 are one animal (typo), stored as DV531;
  - a sample without R2 is excluded (s036), so no contrast mixes layouts.

FASTQ paths are written relative to `raw_dir` in params.yaml.

Usage: build_samplesheet.py <raw_dir> > samplesheet.csv
"""
import csv
import re
import sys
from pathlib import Path

TIMEPOINTS = {"D0": 0, "D1": 1, "D2": 2, "D4": 4}
SPECIES = {"1E7": "rhesus", "1E9": "rhesus", "V": "cynomolgus", "VII": "cynomolgus"}
DOSE = {"1E7": "1e7", "1E9": "1e9", "V": "NA", "VII": "NA"}
ANIMAL_ALIAS = {"DVS31": "DV531"}
NAME = re.compile(r"^(?P<sample>p24240-s\d{3})_(?P<group>[^-]+)-(?P<animal>[^-]+)-"
                  r"(?P<tp>[^_]+)_S\d+_L\d+_R(?P<read>[12])_001\.fastq\.gz$")


def main(raw_dir):
    samples = {}
    for f in sorted(Path(raw_dir).glob("*.fastq.gz")):
        m = NAME.match(f.name)
        if not m:
            sys.exit(f"unexpected file name: {f.name}")
        s = samples.setdefault(m["sample"], {**m.groupdict(), "R1": "", "R2": ""})
        s[f"R{m['read']}"] = f.name

    w = csv.writer(sys.stdout, lineterminator="\n")
    w.writerow(["sample", "animal", "animal_raw", "species", "group", "dose", "tissue",
                "timepoint", "day", "fastq_1", "fastq_2"])
    for sid, s in sorted(samples.items()):
        if s["tp"] not in TIMEPOINTS:
            continue
        if not s["R2"]:
            print(f"excluded {sid}: no R2 file", file=sys.stderr)
            continue
        w.writerow([sid, ANIMAL_ALIAS.get(s["animal"], s["animal"]), s["animal"],
                    SPECIES[s["group"]], s["group"], DOSE[s["group"]], "whole_blood",
                    s["tp"], TIMEPOINTS[s["tp"]], s["R1"], s["R2"]])


if __name__ == "__main__":
    main(sys.argv[1])
