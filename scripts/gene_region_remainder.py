#!/usr/bin/env python3
"""Separate mapped but unclassified UCEs from UCEs absent from the reference."""

import shutil
import sys
from pathlib import Path


def main(copy_dir, output_dir, locations_file):
    copy_dir, output_dir, locations_file = map(Path, (copy_dir, output_dir, locations_file))
    mapped = {line.split("\t", 1)[0] + ".nexus"
              for line in locations_file.read_text().splitlines() if line}
    mapped_dir = output_dir / "unclassified-mapped"
    unmapped_dir = output_dir / "unassigned-regions"
    mapped_dir.mkdir(exist_ok=True)
    unmapped_dir.mkdir(exist_ok=True)
    for alignment in copy_dir.glob("*.nexus"):
        destination = mapped_dir if alignment.name in mapped else unmapped_dir
        shutil.move(str(alignment), destination / alignment.name)


if __name__ == "__main__":
    main(*sys.argv[1:])
