"""
Pytest configuration and shared fixtures for dna-melting-temp-salt-corrected.
"""
import os
import pytest


def pytest_configure(config):
    """Set required environment variables before test collection."""
    # Set a test-only audit secret key if not already configured
    os.environ.setdefault("AUDIT_SECRET_KEY", "test-audit-secret-key-for-pytest-only-32chars")


@pytest.fixture(autouse=True)
def _ensure_audit_key():
    """Ensure AUDIT_SECRET_KEY is available for all tests."""
    os.environ.setdefault("AUDIT_SECRET_KEY", "test-audit-secret-key-for-pytest-only-32chars")
    yield
