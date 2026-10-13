#!/usr/bin/env python3
"""DNA duplex melting temperatures from SantaLucia nearest neighbors and Owczarzy salt corrections.

Assumptions: perfect-match DNA/DNA duplex, non-self-complementary strand
concentration convention (Ct/4), and free Mg2+ supplied by the caller.
Values are estimates and are not a substitute for experimental calibration.
"""
import math
from typing import Any, Dict, Optional

# SantaLucia/Allawi nearest-neighbor DNA/DNA parameters:
# enthalpy in kcal/mol, entropy in cal/(mol K).
NN_DNA = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2), "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}
R_GAS = 1.987  # cal/(mol K)


def nearest_neighbor_tm(sequence: str, oligo_conc_nM: float = 250.0,
                        nn_table: Optional[Dict] = None) -> float:
    """Return uncorrected DNA/DNA Tm (°C), with Ct/4 strand convention."""
    seq = "".join(sequence.split()).upper()
    if len(seq) < 2 or any(base not in "ACGT" for base in seq):
        raise ValueError("sequence must contain at least two A/C/G/T bases")
    if not math.isfinite(oligo_conc_nM) or oligo_conc_nM <= 0:
        raise ValueError("oligo_conc_nM must be finite and positive")
    table = NN_DNA if nn_table is None else nn_table
    # Terminal initiation terms for A/T and G/C ends.
    delta_h = 0.0
    delta_s = 0.0
    for base in (seq[0], seq[-1]):
        dh, ds = (2.3, 4.1) if base in "AT" else (0.1, -2.8)
        delta_h += dh
        delta_s += ds
    for i in range(len(seq) - 1):
        dh, ds = table[seq[i:i + 2]]
        delta_h += dh
        delta_s += ds
    strand_factor = oligo_conc_nM * 1e-9 / 4.0
    denominator = delta_s + R_GAS * math.log(strand_factor)
    if denominator == 0:
        raise ValueError("invalid thermodynamic parameters")
    tm_kelvin = 1000.0 * delta_h / denominator
    if not math.isfinite(tm_kelvin) or tm_kelvin <= 0:
        raise ValueError("calculated temperature is outside the physical range")
    return tm_kelvin - 273.15


def compute_effective_salt_concentration(
    na_mM: float, mg_mM: float, oligo_conc_nM: float = 250.0
) -> float:
    """Diagnostic sodium-equivalent M by the von Ahsen approximation.

    This *is not* the concentration used by the Owczarzy 2008 correction.
    The oligo concentration argument is retained for API compatibility.
    """
    if not all(math.isfinite(v) for v in (na_mM, mg_mM, oligo_conc_nM)):
        raise ValueError("concentrations must be finite")
    if na_mM < 0 or oligo_conc_nM <= 0:
        raise ValueError("concentrations must not be negative")
    return (na_mM + 120.0 * math.sqrt(max(0.0, mg_mM))) / 1000.0


def _apply_inverse_kelvin_correction(tm_celsius: float, correction: float) -> float:
    tm_kelvin = tm_celsius + 273.15
    if not math.isfinite(tm_kelvin) or tm_kelvin <= 0:
        raise ValueError("uncorrected Tm must exceed absolute zero")
    denominator = 1.0 / tm_kelvin + correction
    if denominator <= 0 or not math.isfinite(denominator):
        raise ValueError("salt correction gives an invalid temperature")
    return 1.0 / denominator - 273.15


def owczarzy_2004_salt_correction(tm_nn: float, na_M: float,
                                   f_GC: float) -> float:
    """Owczarzy 2004 monovalent correction applied to inverse Tm in Kelvin."""
    if not all(map(math.isfinite, (tm_nn, na_M, f_GC))):
        raise ValueError("inputs must be finite")
    if not 0 <= f_GC <= 1 or na_M < 0:
        raise ValueError("invalid GC fraction or salt concentration")
    if na_M == 0:
        return tm_nn  # No monovalent correction can be evaluated.
    log_na = math.log(na_M)
    corr = (4.29 * f_GC - 3.95) * 1e-5 * log_na + 9.40e-6 * log_na ** 2
    return _apply_inverse_kelvin_correction(tm_nn, corr)


def owczarzy_2008_mg_correction(tm_nn: float, na_M: float, mg_M: float,
                                f_GC: float, nbp: int = 20) -> float:
    """Owczarzy 2008 salt model; tm_nn is *uncorrected* nearest-neighbor Tm.

    When Mg2+ is absent, returns tm_nn for backwards-compatible standalone use.
    Use owczarzy_2004_salt_correction explicitly for Na-only buffers.
    """
    if not all(map(math.isfinite, (tm_nn, na_M, mg_M, f_GC))):
        raise ValueError("inputs must be finite")
    if na_M < 0 or mg_M < 0 or not 0 <= f_GC <= 1:
        raise ValueError("invalid salt concentration or GC fraction")
    if mg_M == 0:
        return tm_nn
    if nbp < 2:
        raise ValueError("nbp must be >= 2")
    ratio = math.sqrt(mg_M) / na_M if na_M else math.inf
    if ratio < 0.22:
        return owczarzy_2004_salt_correction(tm_nn, na_M, f_GC)

    a, b, c, d, e, f, g = 3.92, -0.911, 6.26, 1.42, -48.2, 52.5, 8.31
    if ratio < 6.0 and na_M > 0:
        ln_na = math.log(na_M)
        a *= 0.843 - 0.352 * math.sqrt(na_M) * ln_na
        d *= 1.279 - 4.03e-3 * ln_na - 8.03e-3 * ln_na ** 2
        g *= 0.486 - 0.258 * ln_na + 5.25e-3 * ln_na ** 3
    ln_mg = math.log(mg_M)
    correction = (
        a + b * ln_mg + f_GC * (c + d * ln_mg)
        + (e + f * ln_mg + g * ln_mg ** 2) / (2.0 * (nbp - 1))
    ) * 1e-5
    return _apply_inverse_kelvin_correction(tm_nn, correction)


def calculate_tm_with_mg2(
    sequence: str, na_mM: float = 50.0, mg_mM: float = 1.5,
    oligo_conc_nM: float = 250.0, nn_table: Optional[Dict] = None,
) -> Dict[str, Any]:
    """Calculate a perfect-match DNA/DNA nearest-neighbor Tm.

    All concentrations use mM except oligo_conc_nM. Mg2+ is the *free*
    magnesium concentration; dNTP chelation is not modelled.
    Errors are returned as {"error": message}, preserving the public API.
    """
    if not isinstance(sequence, str):
        return {"error": "Invalid sequence: must be a non-empty string"}
    seq = "".join(sequence.split()).upper()
    if len(seq) < 2 or any(base not in "ACGT" for base in seq):
        return {"error": "Invalid sequence: at least two ACGT bases required"}
    try:
        ions = (float(na_mM), float(mg_mM), float(oligo_conc_nM))
        if not all(map(math.isfinite, ions)):
            raise ValueError("concentrations must be finite")
        na_mM, mg_mM, oligo_conc_nM = ions
        if na_mM < 0 or mg_mM < 0 or oligo_conc_nM <= 0:
            raise ValueError("Na/Mg must be non-negative and oligo concentration positive")
        if na_mM == 0 and mg_mM == 0:
            raise ValueError("at least one salt concentration must be positive")
        f_gc = (seq.count("G") + seq.count("C")) / len(seq)
        tm_nn = nearest_neighbor_tm(seq, oligo_conc_nM, nn_table)
        tm_na = owczarzy_2004_salt_correction(tm_nn, na_mM / 1000.0, f_gc)
        tm_mg = (owczarzy_2008_mg_correction(
            tm_nn, na_mM / 1000.0, mg_mM / 1000.0, f_gc, len(seq)
        ) if mg_mM else tm_na)
        na_eq = compute_effective_salt_concentration(na_mM, mg_mM, oligo_conc_nM)
    except (ValueError, TypeError, KeyError, OverflowError, ZeroDivisionError) as exc:
        return {"error": str(exc)}
    return {
        "sequence": seq,
        "length": len(seq),
        "gc_fraction": round(f_gc, 4),
        "tm_basic_C": round(2 * (len(seq) - (seq.count("G") + seq.count("C"))) +
                            4 * (seq.count("G") + seq.count("C")), 1),
        "tm_nn_C": round(tm_nn, 2),
        "tm_na_corrected_C": round(tm_na, 2),
        "tm_mg_corrected_C": round(tm_mg, 2),
        "na_mM": na_mM,
        "mg_mM": mg_mM,
        "effective_na_M": round(na_eq, 5),
        "oligo_conc_nM": oligo_conc_nM,
        "method": "SantaLucia nearest-neighbor; Owczarzy 2004/2008 salts",
    }


if __name__ == "__main__":
    import json
    print(json.dumps(calculate_tm_with_mg2("ATCGATCGATCGATCGATCG"), indent=2))
