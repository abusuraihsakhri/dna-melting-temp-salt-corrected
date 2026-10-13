#!/usr/bin/env python3
"""Sequence-based DNA melting-temperature CLI and CSV batch processor."""
import argparse
import csv
import json
import sys
from typing import Any, Dict

from mg2_correction import calculate_tm_with_mg2


def calculate_metrics(sequence: str, na_mM: float = 50.0, mg_mM: float = 1.5,
                      oligo_conc_nM: float = 250.0) -> Dict[str, Any]:
    """Calculate DNA/DNA duplex thermodynamics (never clinical risk scores)."""
    return calculate_tm_with_mg2(sequence, na_mM, mg_mM, oligo_conc_nM)


def process_batch(input_csv: str, output_csv: str, na_mM: float = 50.0,
                  mg_mM: float = 1.5, oligo_conc_nM: float = 250.0) -> int:
    """Evaluate a CSV containing a 'sequence' column; preserve invalid rows with errors."""
    with open(input_csv, newline="", encoding="utf-8-sig") as source:
        reader = csv.DictReader(source)
        fields = list(reader.fieldnames or [])
        if "sequence" not in fields:
            raise ValueError("CSV must contain a 'sequence' header")
        extra = ["tm_nn_C", "tm_na_corrected_C", "tm_mg_corrected_C", "gc_fraction", "error"]
        output_fields = fields + [f for f in extra if f not in fields]
        with open(output_csv, "w", newline="", encoding="utf-8") as target:
            writer = csv.DictWriter(target, fieldnames=output_fields, extrasaction="ignore")
            writer.writeheader()
            total = 0
            for row in reader:
                result = calculate_metrics(row.get("sequence") or "", na_mM, mg_mM, oligo_conc_nM)
                writer.writerow({**row, **{f: result.get(f, "") for f in extra}})
                total += 1
    return total


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Salt-corrected nearest-neighbor DNA melting temperature")
    subparsers = parser.add_subparsers(dest="command", required=True)
    single = subparsers.add_parser("single", help="Estimate melting temperature for one sequence")
    single.add_argument("--sequence", "-s", required=True, help="DNA sequence (5' to 3')")
    batch = subparsers.add_parser("batch", help="Evaluate a CSV with a sequence column")
    batch.add_argument("--input", "-i", required=True)
    batch.add_argument("--output", "-o", default="results.csv")
    for p in (single, batch):
        p.add_argument("--na", type=float, default=50.0, help="Na+ concentration, mM")
        p.add_argument("--mg", type=float, default=1.5, help="Free Mg2+ concentration, mM")
        p.add_argument("--oligo", type=float, default=250.0, help="Total oligo concentration, nM")
    args = parser.parse_args(argv)
    if args.command == "single":
        result = calculate_metrics(args.sequence, args.na, args.mg, args.oligo)
        print(json.dumps(result, indent=2))
        return 2 if "error" in result else 0
    try:
        count = process_batch(args.input, args.output, args.na, args.mg, args.oligo)
    except (OSError, ValueError, csv.Error) as exc:
        parser.error(str(exc))
    print(f"Processed {count} sequences -> {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
