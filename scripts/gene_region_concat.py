#!/usr/bin/env python3
"""Run GeneRegion PHYLUCE concatenations without shell interpretation of names."""

import argparse
import concurrent.futures
import subprocess
import sys
from pathlib import Path


def jobs(mode, source, destination):
    if mode == "exon":
        directories = (p for p in source.glob("*/*") if p.is_dir())
    else:
        directories = (p for p in source.iterdir() if p.is_dir()) if source.exists() else ()
    for directory in sorted(directories):
        if mode == "exon":
            name = (directory.parent.name + "__" + directory.name).replace("-", "_")
        else:
            name = directory.name.replace("-", "_")
        yield directory, destination / name, name


def run_job(job, mode, phyluce, log_root):
    source, output, name = job
    job_log = log_root / name
    job_log.mkdir(parents=True, exist_ok=True)
    command = [str(phyluce), "--alignments", str(source), "--log", str(job_log),
               "--output", str(output)]
    if mode == "gene":
        command.extend(("--phylip", "--input-format", "nexus"))
        expected = (output / (name + ".phylip"), output / (name + ".charsets"))
    else:
        command.append("--nexus")
        expected = (output / (name + ".nexus"),)
    with (job_log / "console.log").open("w") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
    if result.returncode != 0 or any(not path.is_file() or path.stat().st_size == 0 for path in expected):
        raise RuntimeError(f"{mode} concatenation failed for {source}; see {job_log}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("exon", "intron", "gene"))
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("log_root", type=Path)
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be positive")
    args.destination.mkdir(parents=True, exist_ok=True)
    args.log_root.mkdir(parents=True, exist_ok=True)
    work = list(jobs(args.mode, args.source, args.destination))
    phyluce = Path(sys.executable).parent / "phyluce_align_concatenate_alignments"
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = [executor.submit(run_job, job, args.mode, phyluce, args.log_root) for job in work]
        for future in concurrent.futures.as_completed(futures):
            future.result()
    print(f"Completed {len(work)} {args.mode} concatenations", flush=True)


if __name__ == "__main__":
    main()
