"""
Tests for mg2_correction module - DNA melting temperature calculations.
"""
import math
import pytest
from mg2_correction import (
    compute_effective_salt_concentration,
    owczarzy_2004_salt_correction,
    owczarzy_2008_mg_correction,
    calculate_tm_with_mg2,
)


class TestComputeEffectiveSaltConcentration:
    def test_no_mg_returns_na_converted(self):
        """Without Mg2+, effective salt should be Na+ converted to M."""
        result = compute_effective_salt_concentration(na_mM=50.0, mg_mM=0.0)
        assert result == pytest.approx(0.05, rel=1e-3)

    def test_with_mg_dominant_regime(self):
        """With Mg2+ present, effective salt should be computed."""
        result = compute_effective_salt_concentration(na_mM=50.0, mg_mM=1.5)
        assert result > 0
        assert isinstance(result, float)

    def test_zero_na_with_mg(self):
        """Zero Na+ with Mg2+ should still produce valid result."""
        result = compute_effective_salt_concentration(na_mM=0.0, mg_mM=1.5)
        assert result >= 0

    def test_negative_mg_treated_as_zero(self):
        """Negative Mg2+ should be treated as zero Mg2+."""
        result = compute_effective_salt_concentration(na_mM=50.0, mg_mM=-1.0)
        assert result == pytest.approx(0.05, rel=1e-3)


class TestOwczarzy2004SaltCorrection:
    def test_zero_na_returns_uncorrected(self):
        """Zero Na+ should return uncorrected Tm."""
        result = owczarzy_2004_salt_correction(tm_nn=60.0, na_M=0.0, f_GC=0.5)
        assert result == 60.0

    def test_positive_na_applies_correction(self):
        """Positive Na+ should apply salt correction."""
        result = owczarzy_2004_salt_correction(tm_nn=60.0, na_M=0.05, f_GC=0.5)
        assert result != 60.0
        assert isinstance(result, float)

    def test_high_gc_increases_correction(self):
        """Higher GC fraction should produce different correction."""
        result_low_gc = owczarzy_2004_salt_correction(tm_nn=60.0, na_M=0.05, f_GC=0.2)
        result_high_gc = owczarzy_2004_salt_correction(tm_nn=60.0, na_M=0.05, f_GC=0.8)
        assert result_low_gc != result_high_gc


class TestOwczarzy2008MgCorrection:
    def test_zero_mg_returns_uncorrected(self):
        """Zero Mg2+ should return uncorrected Tm."""
        result = owczarzy_2008_mg_correction(tm_nn=60.0, na_M=0.05, mg_M=0.0, f_GC=0.5)
        assert result == 60.0

    def test_zero_celsius_tm_is_corrected(self):
        """0 °C is a physical temperature, not a missing-value sentinel."""
        result = owczarzy_2008_mg_correction(tm_nn=0.0, na_M=0.05, mg_M=0.001, f_GC=0.5)
        assert result != 0.0
        assert math.isfinite(result)

    def test_positive_mg_applies_correction(self):
        """Positive Mg2+ should apply correction."""
        result = owczarzy_2008_mg_correction(tm_nn=60.0, na_M=0.05, mg_M=0.001, f_GC=0.5)
        assert result != 60.0
        assert isinstance(result, float)

    def test_nbp_must_be_at_least_2(self):
        """nbp < 2 should raise ValueError."""
        with pytest.raises(ValueError, match="nbp"):
            owczarzy_2008_mg_correction(tm_nn=60.0, na_M=0.05, mg_M=0.001, f_GC=0.5, nbp=1)

    def test_nbp_2_is_valid(self):
        """nbp = 2 should be valid."""
        result = owczarzy_2008_mg_correction(tm_nn=60.0, na_M=0.05, mg_M=0.001, f_GC=0.5, nbp=2)
        assert isinstance(result, float)


class TestCalculateTmWithMg2:
    def test_valid_sequence(self):
        """Valid DNA sequence should return Tm results."""
        result = calculate_tm_with_mg2("ATCGATCGATCGATCGATCG", na_mM=50, mg_mM=1.5)
        assert "error" not in result
        assert "tm_mg_corrected_C" in result
        assert "gc_fraction" in result
        assert result["length"] == 20

    def test_gc_fraction_calculation(self):
        """GC fraction should be calculated correctly."""
        result = calculate_tm_with_mg2("GCGCGCGCGCGCGCGCGCGC", na_mM=50, mg_mM=0)
        assert result["gc_fraction"] == 1.0

    def test_at_only_sequence(self):
        """AT-only sequence should have gc_fraction = 0."""
        result = calculate_tm_with_mg2("ATATATATATATATATATAT", na_mM=50, mg_mM=0)
        assert result["gc_fraction"] == 0.0

    def test_invalid_sequence_returns_error(self):
        """Invalid characters should return error."""
        result = calculate_tm_with_mg2("ATCGXYZ", na_mM=50, mg_mM=1.5)
        assert "error" in result

    def test_empty_sequence_returns_error(self):
        """Empty sequence should return error."""
        result = calculate_tm_with_mg2("", na_mM=50, mg_mM=1.5)
        assert "error" in result

    def test_whitespace_only_returns_error(self):
        """Whitespace-only sequence should return error."""
        result = calculate_tm_with_mg2("   ", na_mM=50, mg_mM=1.5)
        assert "error" in result

    def test_negative_na_returns_error(self):
        """Negative Na+ concentration should return error."""
        result = calculate_tm_with_mg2("ATCGATCG", na_mM=-10, mg_mM=1.5)
        assert "error" in result

    def test_negative_mg_returns_error(self):
        """Negative Mg2+ concentration should return error."""
        result = calculate_tm_with_mg2("ATCGATCG", na_mM=50, mg_mM=-1.0)
        assert "error" in result

    def test_zero_oligo_conc_returns_error(self):
        """Zero oligo concentration should return error."""
        result = calculate_tm_with_mg2("ATCGATCG", na_mM=50, mg_mM=1.5, oligo_conc_nM=0)
        assert "error" in result

    def test_no_mg_skips_mg_correction(self):
        """With mg_mM=0, tm_mg_corrected should equal tm_na_corrected."""
        result = calculate_tm_with_mg2("ATCGATCGATCGATCGATCG", na_mM=50, mg_mM=0)
        assert result["tm_mg_corrected_C"] == result["tm_na_corrected_C"]

    def test_sequence_with_spaces(self):
        """Sequence with spaces should be cleaned."""
        result = calculate_tm_with_mg2("ATCG ATCG ATCG ATCG", na_mM=50, mg_mM=0)
        assert "error" not in result
        assert result["length"] == 16

    def test_sequence_lowercase(self):
        """Lowercase sequence should be handled."""
        result = calculate_tm_with_mg2("atcgatcgatcgatcgatcg", na_mM=50, mg_mM=0)
        assert "error" not in result
        assert result["length"] == 20

    def test_result_contains_expected_fields(self):
        """Result should contain all expected fields."""
        result = calculate_tm_with_mg2("ATCGATCGATCGATCGATCG", na_mM=50, mg_mM=1.5)
        expected_fields = [
            "sequence", "length", "gc_fraction", "tm_basic_C",
            "tm_na_corrected_C", "tm_mg_corrected_C", "na_mM", "mg_mM",
            "effective_na_M", "oligo_conc_nM"
        ]
        for field in expected_fields:
            assert field in result, f"Missing field: {field}"
