#!/usr/bin/env python3
"""
Multi-Primer Pair Comparison Engine
Compare and rank multiple primer pairs by Tm matching, GC balance,
and secondary structure absence.
"""
import csv
import json
from typing import Dict, Any, List, Optional

from mg2_correction import calculate_tm_with_mg2
from secondary_structure import analyze_secondary_structures


def evaluate_primer_pair(
    forward: str,
    reverse: str,
    na_mM: float = 50.0,
    mg_mM: float = 1.5,
    oligo_conc_nM: float = 250.0,
) -> Dict[str, Any]:
    """Evaluate a single primer pair with all metrics.

    Args:
        forward: Forward primer sequence.
        reverse: Reverse primer sequence.
        na_mM: Na+ concentration in mM.
        mg_mM: Mg2+ concentration in mM.
        oligo_conc_nM: Oligo concentration in nM.

    Returns:
        Dictionary with all evaluation metrics.
    """
    fwd_result = calculate_tm_with_mg2(forward, na_mM, mg_mM, oligo_conc_nM)
    rev_result = calculate_tm_with_mg2(reverse, na_mM, mg_mM, oligo_conc_nM)

    fwd_struct = analyze_secondary_structures(forward)
    rev_struct = analyze_secondary_structures(reverse)

    tm_diff = abs(fwd_result["tm_mg_corrected_C"] - rev_result["tm_mg_corrected_C"])
    gc_diff = abs(fwd_result["gc_fraction"] - rev_result["gc_fraction"])

    # Scoring: lower is better
    tm_score = tm_diff * 10  # Heavy penalty for Tm mismatch
    gc_score = gc_diff * 50  # Penalty for GC imbalance
    struct_score = abs(fwd_struct["worst_dg_kcal_mol"]) + abs(rev_struct["worst_dg_kcal_mol"])

    total_score = tm_score + gc_score + struct_score

    # Quality flags
    flags = []
    if tm_diff > 2.0:
        flags.append(f"Tm difference ({tm_diff:.1f}°C) exceeds 2°C threshold")
    if fwd_result["gc_fraction"] < 0.4 or fwd_result["gc_fraction"] > 0.6:
        flags.append(f"Forward GC% ({fwd_result['gc_fraction']*100:.0f}%) outside 40-60% range")
    if rev_result["gc_fraction"] < 0.4 or rev_result["gc_fraction"] > 0.6:
        flags.append(f"Reverse GC% ({rev_result['gc_fraction']*100:.0f}%) outside 40-60% range")
    if fwd_struct["worst_dg_kcal_mol"] < -3.0:
        flags.append("Forward primer has stable secondary structure")
    if rev_struct["worst_dg_kcal_mol"] < -3.0:
        flags.append("Reverse primer has stable secondary structure")

    return {
        "forward": forward,
        "reverse": reverse,
        "fwd_tm_C": fwd_result["tm_mg_corrected_C"],
        "rev_tm_C": rev_result["tm_mg_corrected_C"],
        "tm_difference_C": round(tm_diff, 1),
        "fwd_gc_pct": round(fwd_result["gc_fraction"] * 100, 1),
        "rev_gc_pct": round(rev_result["gc_fraction"] * 100, 1),
        "fwd_hairpin_dg": fwd_struct["hairpin"].get("dg_total_kcal_mol", 0.0),
        "rev_hairpin_dg": rev_struct["hairpin"].get("dg_total_kcal_mol", 0.0),
        "fwd_self_dimer_dg": fwd_struct["self_dimer"].get("dg_total_kcal_mol", 0.0),
        "rev_self_dimer_dg": rev_struct["self_dimer"].get("dg_total_kcal_mol", 0.0),
        "total_score": round(total_score, 2),
        "flags": flags,
    }


def rank_primer_pairs(
    pairs: List[Dict[str, str]],
    na_mM: float = 50.0,
    mg_mM: float = 1.5,
    oligo_conc_nM: float = 250.0,
    top_n: int = 10,
) -> List[Dict[str, Any]]:
    """Rank multiple primer pairs by quality metrics.

    Args:
        pairs: List of dicts with 'forward' and 'reverse' keys.
        na_mM: Na+ concentration in mM.
        mg_mM: Mg2+ concentration in mM.
        oligo_conc_nM: Oligo concentration in nM.
        top_n: Number of top pairs to return.

    Returns:
        Ranked list of primer pair evaluations.
    """
    results = []
    for i, pair in enumerate(pairs):
        fwd = pair.get("forward", pair.get("fwd", ""))
        rev = pair.get("reverse", pair.get("rev", ""))
        if not fwd or not rev:
            continue
        eval_result = evaluate_primer_pair(fwd, rev, na_mM, mg_mM, oligo_conc_nM)
        eval_result["pair_id"] = pair.get("id", f"pair_{i+1}")
        results.append(eval_result)

    # Sort by total score (lower is better)
    results.sort(key=lambda x: x["total_score"])

    for rank, result in enumerate(results[:top_n], 1):
        result["rank"] = rank

    return results[:top_n]


def compare_from_csv(
    input_csv: str,
    output_csv: Optional[str] = None,
    na_mM: float = 50.0,
    mg_mM: float = 1.5,
) -> List[Dict[str, Any]]:
    """Compare primer pairs from a CSV file.

    Args:
        input_csv: Path to CSV with columns: id, forward, reverse.
        output_csv: Optional path for output CSV.
        na_mM: Na+ concentration in mM.
        mg_mM: Mg2+ concentration in mM.

    Returns:
        Ranked list of primer pair evaluations.
    """
    with open(input_csv, encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        pairs = list(reader)

    results = rank_primer_pairs(pairs, na_mM=na_mM, mg_mM=mg_mM)

    if output_csv:
        fieldnames = [
            "rank", "pair_id", "forward", "reverse",
            "fwd_tm_C", "rev_tm_C", "tm_difference_C",
            "fwd_gc_pct", "rev_gc_pct",
            "fwd_hairpin_dg", "rev_hairpin_dg",
            "fwd_self_dimer_dg", "rev_self_dimer_dg",
            "total_score", "flags",
        ]
        with open(output_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            for r in results:
                row = dict(r)
                row["flags"] = "; ".join(row.get("flags", []))
                writer.writerow(row)

    return results


if __name__ == "__main__":
    import json
    test_pairs = [
        {"id": "pair1", "forward": "ATCGATCGATCGATCG", "reverse": "GCTAGCTAGCTAGCTA"},
        {"id": "pair2", "forward": "GCGCGCGCGCGCGCGC", "reverse": "ATATATATATATATAT"},
    ]
    results = rank_primer_pairs(test_pairs)
    print(json.dumps(results, indent=2))
