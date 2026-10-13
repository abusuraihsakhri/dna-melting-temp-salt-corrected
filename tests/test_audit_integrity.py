"""Audit integrity checks do not rely on mutable chain links alone."""
from agents.base import AuditTrail


def test_hmac_rechecks_contents_and_chain():
    trail = AuditTrail("test-secret-only-for-unit-tests-0001")
    trail.log("test", "local", "CALC", {"sequence": "ATGC"})
    trail.log("test", "local", "CALC", {"sequence": "CGTA"})
    assert trail.verify_integrity()
    snapshot = trail.get_trail()
    snapshot[0]["actor"] = "modified"
    assert trail.verify_integrity()
    trail.logs[0]["actor"] = "tampered"
    assert not trail.verify_integrity()


def test_modified_signature_fails_even_single_entry():
    trail = AuditTrail("test-secret-only-for-unit-tests-0001")
    trail.log("test", "local", "RUN", {})
    trail.logs[0]["current_hash"] = "0" * 64
    assert not trail.verify_integrity()
