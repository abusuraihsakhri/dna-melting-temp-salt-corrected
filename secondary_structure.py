#!/usr/bin/env python3
"""
Secondary Structure Predictor for Primer Hairpin and Self-Dimer Stability.
Implements Nussinov-style algorithm for hairpin prediction and
sliding-window complementarity for self-dimer detection.

Uses SantaLucia 1998 nearest-neighbor parameters for ΔG estimation.
"""
import math
from typing import Dict, Any, List, Tuple, Optional

# SantaLucia 1998 nearest-neighbor ΔG°37 (kcal/mol) in 1M NaCl
NN_DG37 = {
    "AA": -1.0, "TT": -1.0,
    "AT": -0.88,
    "TA": -0.58,
    "CA": -1.45, "TG": -1.45,
    "GT": -1.44, "AC": -1.44,
    "CT": -1.28, "AG": -1.28,
    "GA": -1.30, "TC": -1.30,
    "CG": -2.17,
    "GC": -2.24,
    "GG": -1.84, "CC": -1.84,
}

# Initiation parameters
INIT_DG37 = {"G": 0.98, "C": 0.98, "A": 1.03, "T": 1.03}

COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return "".join(COMPLEMENT.get(b, "N") for b in reversed(seq))


def complement(seq: str) -> str:
    """Return the complement of a DNA sequence."""
    return "".join(COMPLEMENT.get(b, "N") for b in seq)


def compute_duplex_dg(seq1: str, seq2: str) -> float:
    """Compute ΔG°37 for a perfect duplex alignment using nearest-neighbor model.

    Args:
        seq1: First strand (5' to 3').
        seq2: Second strand (3' to 5', i.e., complement of the pairing strand).

    Returns:
        ΔG°37 in kcal/mol (negative = stable).
    """
    if len(seq1) != len(seq2):
        return 0.0

    dg = 0.0
    for i in range(len(seq1) - 1):
        pair = seq1[i:i+2]
        dg += NN_DG37.get(pair, 0.0)

    # Initiation
    dg += INIT_DG37.get(seq1[0], 1.0)
    dg += INIT_DG37.get(seq1[-1], 1.0)

    return round(dg, 2)


def predict_hairpin(
    sequence: str,
    min_loop_size: int = 3,
    min_stem_length: int = 3,
) -> Dict[str, Any]:
    """Predict the most stable hairpin structure using a Nussinov-style approach.

    Scans for self-complementary regions that can form a stem with a loop.

    Args:
        sequence: DNA sequence (5' to 3').
        min_loop_size: Minimum number of bases in the loop.
        min_stem_length: Minimum number of base pairs in the stem.

    Returns:
        Dictionary with hairpin prediction results.
    """
    seq = sequence.upper().replace(" ", "")
    n = len(seq)
    best_dg = 0.0
    best_structure = None

    for loop_start in range(min_loop_size, n - min_stem_length):
        for stem_len in range(min_stem_length, min(loop_start, n - loop_start) + 1):
            # Check if the regions flanking the loop can form a stem
            left_start = loop_start - stem_len
            right_start = loop_start
            right_end = right_start + stem_len

            if right_end > n:
                continue

            # Check complementarity
            matches = 0
            for i in range(stem_len):
                left_base = seq[left_start + i]
                right_base = seq[right_end - 1 - i]
                if COMPLEMENT.get(left_base) == right_base:
                    matches += 1

            if matches < min_stem_length:
                continue

            # Compute ΔG for the stem
            stem_seq = seq[left_start:left_start + stem_len]
            stem_comp = reverse_complement(seq[right_start:right_start + stem_len])
            dg_stem = compute_duplex_dg(stem_seq, stem_comp)

            # Loop penalty: ~1.5 kcal/mol per loop base (simplified)
            loop_penalty = 1.5 * (loop_start - left_start - stem_len + 1) if loop_start - left_start > stem_len else 3.0
            dg_total = dg_stem + loop_penalty

            if dg_total < best_dg:
                best_dg = dg_total
                best_structure = {
                    "stem_start": left_start,
                    "stem_length": stem_len,
                    "loop_start": loop_start,
                    "loop_size": loop_start - left_start - stem_len if loop_start - left_start > stem_len else 0,
                    "matches": matches,
                    "dg_stem": round(dg_stem, 2),
                    "dg_total": round(dg_total, 2),
                }

    if best_structure is None:
        return {
            "sequence": seq,
            "hairpin_found": False,
            "dg_total": 0.0,
            "message": "No stable hairpin structure predicted.",
        }

    return {
        "sequence": seq,
        "hairpin_found": True,
        "stem_start": best_structure["stem_start"],
        "stem_length": best_structure["stem_length"],
        "loop_start": best_structure["loop_start"],
        "loop_size": best_structure["loop_size"],
        "matches": best_structure["matches"],
        "dg_stem_kcal_mol": best_structure["dg_stem"],
        "dg_total_kcal_mol": best_structure["dg_total"],
    }


def predict_self_dimer(
    sequence: str,
    min_overlap: int = 2,
) -> Dict[str, Any]:
    """Predict the most stable self-dimer structure.

    Slides the primer against its reverse complement and finds the best alignment.

    Args:
        sequence: DNA sequence (5' to 3').
        min_overlap: Minimum overlap length to consider.

    Returns:
        Dictionary with self-dimer prediction results.
    """
    seq = sequence.upper().replace(" ", "")
    n = len(seq)
    rc = reverse_complement(seq)

    best_dg = 0.0
    best_offset = 0
    best_matches = 0

    for offset in range(-(n - min_overlap), n - min_overlap + 1):
        matches = 0
        aligned_pairs = []

        for i in range(n):
            j = i + offset
            if 0 <= j < n:
                if COMPLEMENT.get(seq[i]) == seq[j]:
                    matches += 1
                    aligned_pairs.append((i, j))

        if matches < min_overlap:
            continue

        # Compute ΔG for the overlap region
        overlap_seq = ""
        for i, j in aligned_pairs:
            overlap_seq += seq[i]

        if len(overlap_seq) < 2:
            continue

        dg = compute_duplex_dg(overlap_seq, reverse_complement(overlap_seq))

        if dg < best_dg:
            best_dg = dg
            best_offset = offset
            best_matches = matches

    if best_dg >= 0:
        return {
            "sequence": seq,
            "self_dimer_found": False,
            "dg_total": 0.0,
            "message": "No stable self-dimer structure predicted.",
        }

    return {
        "sequence": seq,
        "self_dimer_found": True,
        "offset": best_offset,
        "matches": best_matches,
        "dg_total_kcal_mol": round(best_dg, 2),
    }


def analyze_secondary_structures(sequence: str) -> Dict[str, Any]:
    """Complete secondary structure analysis for a primer.

    Args:
        sequence: DNA primer sequence.

    Returns:
        Dictionary with hairpin and self-dimer analysis.
    """
    seq = sequence.upper().replace(" ", "")
    hairpin = predict_hairpin(seq)
    self_dimer = predict_self_dimer(seq)

    worst_dg = min(
        hairpin.get("dg_total_kcal_mol", 0.0),
        self_dimer.get("dg_total_kcal_mol", 0.0),
    )

    stability = "Stable" if worst_dg < -3.0 else "Marginal" if worst_dg < -1.0 else "Unstable"

    return {
        "sequence": seq,
        "length": len(seq),
        "hairpin": hairpin,
        "self_dimer": self_dimer,
        "worst_dg_kcal_mol": worst_dg,
        "stability_assessment": stability,
    }


if __name__ == "__main__":
    import json
    test_seq = "ATCGATCGATCGATCGATCG"
    result = analyze_secondary_structures(test_seq)
    print(json.dumps(result, indent=2))
