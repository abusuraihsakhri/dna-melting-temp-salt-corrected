"""Numerical benchmarks independently obtained with Biopython Tm_NN (DNA_NN3).

Reference configuration: dnac1=125, dnac2=125, Na=50 mM,
Mg=0 or 1.5 mM, saltcorr=7, 0 mM dNTPs.
The reference is not a runtime dependency.
"""
import math

import pytest

from mg2_correction import calculate_tm_with_mg2, nearest_neighbor_tm


@pytest.mark.parametrize("seq,mg,expected", [
    ("ATCGATCGATCGATCGATCG", 0.0, 53.09),
    ("ATCGATCGATCGATCGATCG", 1.5, 59.72),
    ("GCGCGCGCGCGCGCGCGCGC", 0.0, 83.51),
    ("GCGCGCGCGCGCGCGCGCGC", 1.5, 86.30),
    ("ATATATATATATATATATAT", 0.0, 22.43),
    ("ATATATATATATATATATAT", 1.5, 31.46),
    ("AGTCAGTCAGTCAGTCAGTC", 1.5, 59.07),
    ("ACTGACTGACTGACTG", 1.5, 52.08),
])
def test_matches_reference_nn(seq, mg, expected):
    result = calculate_tm_with_mg2(seq, na_mM=50.0, mg_mM=mg)
    assert "error" not in result, result
    assert result["tm_mg_corrected_C"] == pytest.approx(expected, abs=0.03)


def test_uncorrected_nn_reference():
    assert nearest_neighbor_tm("ATCGATCGATCGATCGATCG") == pytest.approx(68.51875, abs=0.02)


def test_monovalent_dominant_matches_na_branch():
    result = calculate_tm_with_mg2("AGTCAGTCAGTCAGTCAGTC", na_mM=1000, mg_mM=0.01)
    assert result["tm_mg_corrected_C"] == result["tm_na_corrected_C"]


@pytest.mark.parametrize("kwargs", [
    {"na_mM": float("nan")},
    {"mg_mM": float("inf")},
    {"oligo_conc_nM": -1},
    {"na_mM": 0, "mg_mM": 0},
])
def test_invalid_buffers_are_rejected(kwargs):
    assert "error" in calculate_tm_with_mg2("ATCGATCGATCG", **kwargs)


def test_invalid_oligos_not_silently_accepted():
    assert "error" in calculate_tm_with_mg2("A")
    assert "error" in calculate_tm_with_mg2("ATC N")
    assert calculate_tm_with_mg2("aT cG\nAT CG", mg_mM=0)["length"] == 8
