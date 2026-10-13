import pytest
from fastapi import HTTPException

from agents.api import TmRequest, api_tm


def test_tm_api_returns_dna_result():
    result = api_tm(TmRequest(sequence="ATCGATCGATCGATCGATCG"))
    assert result["tm_mg_corrected_C"] == 59.72
    assert result["length"] == 20


def test_tm_api_rejects_bad_sequences():
    with pytest.raises(HTTPException) as exc:
        api_tm(TmRequest(sequence="AGCTN"))
    assert exc.value.status_code == 422
