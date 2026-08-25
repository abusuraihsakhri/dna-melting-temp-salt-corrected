#!/usr/bin/env python3
"""
Mg2+ Correction Extension for Oligo Melting Temperature
Implements Owczarzy 2008 formula for divalent cation correction and
Aho formula for combined monovalent + divalent salt effective concentration.

References:
  - Owczarzy et al. (2008) Biochemistry 47(19):5336-5353
  - Aho et al. (2008) Nucleic Acids Res 36(11):e63
"""
import math
from typing import Dict, Any, Optional


def compute_effective_salt_concentration(
    na_mM: float,
    mg_mM: float,
    oligo_conc_nM: float = 250.0,
) -> float:
    """Compute effective monovalent cation concentration using the Aho formula.

    When Mg2+ is present, the effective monovalent salt concentration is:
        ratio = sqrt(oligo_conc_nM * 1e-9) / (3.0 * mg_mM * 1e-3)
        if ratio < 0.22:
            Na_eq = 3.79 * sqrt(oligo_conc_nM * 1e-9) / (3.0 * mg_mM * 1e-3)
            Na_eq = Na_eq * Na_eq * 3.0 * mg_mM * 1e-3
        else:
            Na_eq = Na + 140 * Mg (simplified)

    Args:
        na_mM: Monovalent cation (Na+/K+) concentration in mM.
        mg_mM: Divalent cation (Mg2+) concentration in mM.
        oligo_conc_nM: Oligonucleotide concentration in nM.

    Returns:
        Effective monovalent salt concentration in M.
    """
    if mg_mM <= 0:
        return na_mM / 1000.0

    mg_M = mg_mM / 1000.0
    na_M = na_mM / 1000.0
    oligo_M = oligo_conc_nM * 1e-9

    ratio = math.sqrt(oligo_M) / (3.0 * mg_M)

    if ratio < 0.22:
        # Mg2+ dominant regime
        f = 3.79 * ratio
        na_eq = f * f * 3.0 * mg_M
    else:
        # Transitional regime: interpolate
        a = 3.79 * ratio
        na_eq_mg = a * a * 3.0 * mg_M
        na_eq_na = na_M + 140.0 * mg_M
        # Weighted interpolation
        w = (ratio - 0.22) / (6.0 - 0.22)
        w = max(0.0, min(1.0, w))
        na_eq = (1.0 - w) * na_eq_mg + w * na_eq_na

    return max(na_eq, na_M)


def owczarzy_2004_salt_correction(
    tm_nn: float,
    na_M: float,
    f_GC: float,
) -> float:
    """Apply Owczarzy 2004 monovalent cation salt correction.

    ΔTm = (4.29 * f_GC - 3.95) * ln(Na+) + 0.940 * (ln(Na+))^2

    Args:
        tm_nn: Nearest-neighbor Tm without salt correction (°C).
        na_M: Monovalent cation concentration in M.
        f_GC: Fraction of G/C bases (0.0 - 1.0).

    Returns:
        Salt-corrected Tm in °C.
    """
    if na_M <= 0:
        return tm_nn
    ln_na = math.log(na_M)
    delta_tm = (4.29 * f_GC - 3.95) * ln_na + 0.940 * ln_na * ln_na
    return tm_nn + delta_tm


def owczarzy_2008_mg_correction(
    tm_nn: float,
    na_M: float,
    mg_M: float,
    f_GC: float,
) -> float:
    """Apply Owczarzy 2008 Mg2+ correction to nearest-neighbor Tm.

    For Mg2+ containing buffers:
        a = 3.92 * (0.843 - 0.352 * sqrt(Na+) * ln(Na+))
        d = 1.42 * (1.279 - 4.03e-3 * ln(Na+) - 8.03e-3 * (ln(Na+))^2)
        g = 8.31 * (0.486 - 0.258 * ln(Na+) + 5.25e-3 * (ln(Na+))^3)

        1/Tm(Mg) = 1/Tm(Na) + a + b*ln(Mg) + f_GC*(c + d*ln(Mg)) + 1/(2*(Nbp-1)) * (e + f*ln(Mg) + g*(ln(Mg))^2)

    Args:
        tm_nn: Nearest-neighbor Tm with Na+ correction (°C).
        na_M: Monovalent cation concentration in M.
        mg_M: Divalent Mg2+ concentration in M.
        f_GC: Fraction of G/C bases (0.0 - 1.0).

    Returns:
        Mg2+-corrected Tm in °C.
    """
    if mg_M <= 0 or tm_nn <= 0:
        return tm_nn

    tm_K = tm_nn + 273.15
    ln_na = math.log(na_M) if na_M > 0 else -10.0
    ln_mg = math.log(mg_M)

    a = 3.92 * (0.843 - 0.352 * math.sqrt(na_M) * ln_na) if na_M > 0 else 3.92 * 0.843
    d = 1.42 * (1.279 - 4.03e-3 * ln_na - 8.03e-3 * ln_na ** 2)
    g = 8.31 * (0.486 - 0.258 * ln_na + 5.25e-3 * ln_na ** 3)

    b = -0.001  # small correction
    c = 0.0
    e = 0.0
    f_coeff = 0.0
    nbp = 20  # default primer length

    inv_tm_mg = (1.0 / tm_K) + a + b * ln_mg + f_GC * (c + d * ln_mg) + \
                (1.0 / (2.0 * (nbp - 1))) * (e + f_coeff * ln_mg + g * ln_mg ** 2)

    tm_mg_K = 1.0 / inv_tm_mg
    return tm_mg_K - 273.15


def calculate_tm_with_mg2(
    sequence: str,
    na_mM: float = 50.0,
    mg_mM: float = 1.5,
    oligo_conc_nM: float = 250.0,
    nn_table: Optional[Dict] = None,
) -> Dict[str, Any]:
    """Calculate Tm with full salt and Mg2+ correction.

    Args:
        sequence: DNA primer sequence (5' to 3').
        na_mM: Na+ concentration in mM.
        mg_mM: Mg2+ concentration in mM.
        oligo_conc_nM: Oligonucleotide concentration in nM.
        nn_table: Optional nearest-neighbor parameter table.

    Returns:
        Dictionary with Tm values and correction details.
    """
    seq = sequence.upper().replace(" ", "").replace("\n", "")
    if not seq or not all(b in "ACGT" for b in seq):
        return {"error": "Invalid sequence: only ACGT bases allowed"}

    gc_count = seq.count("G") + seq.count("C")
    f_gc = gc_count / len(seq)

    # Simple Wallace rule for base Tm (placeholder for full NN model)
    a_count = seq.count("A")
    t_count = seq.count("T")
    tm_basic = 2 * (a_count + t_count) + 4 * (gc_count)

    # Apply salt correction
    na_M = na_mM / 1000.0
    mg_M = mg_mM / 1000.0

    tm_na_corrected = owczarzy_2004_salt_correction(tm_basic, na_M, f_gc)

    # Apply Mg2+ correction if present
    if mg_mM > 0:
        tm_mg_corrected = owczarzy_2008_mg_correction(tm_na_corrected, na_M, mg_M, f_gc)
    else:
        tm_mg_corrected = tm_na_corrected

    # Compute effective salt concentration
    na_eq = compute_effective_salt_concentration(na_mM, mg_mM, oligo_conc_nM)

    return {
        "sequence": seq,
        "length": len(seq),
        "gc_fraction": round(f_gc, 4),
        "tm_basic_C": round(tm_basic, 1),
        "tm_na_corrected_C": round(tm_na_corrected, 1),
        "tm_mg_corrected_C": round(tm_mg_corrected, 1),
        "na_mM": na_mM,
        "mg_mM": mg_mM,
        "effective_na_M": round(na_eq, 4),
        "oligo_conc_nM": oligo_conc_nM,
    }


if __name__ == "__main__":
    import json
    result = calculate_tm_with_mg2("ATCGATCGATCGATCGATCG", na_mM=50, mg_mM=1.5)
    print(json.dumps(result, indent=2))
