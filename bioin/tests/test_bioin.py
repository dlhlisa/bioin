"""
Unit and regression test for the bioin package.
"""

# Import package, test suite, and other packages as needed
import bioin
import pytest
import sys

def test_bioin_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "bioin" in sys.modules
