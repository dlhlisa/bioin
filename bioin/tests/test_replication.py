import bioin
from .. import replication
import pytest
import sys


@pytest.mark.parametrize("pattern, genome, positions", [
    ("ATAT", "GATATATGCATATACTT", [1, 3, 9]),
    ("ACAC", "TTTTACACTTTTTTGTGTAAAAA", [4]),
    ("AAA", "AAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAGCGCACAGACAGATAATAATTACAGAGTACACAACATCCAT", [0, 46, 51, 74]),
    ("TTT", "AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTCACGGATGCCAAGTGCATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTTTATCCAGCGCGTTAACCACGTTCTGTGCCGACTTT", [88, 92, 98, 132]),
    ("ATA", "ATATATA", [0, 2, 4])
    ])
def test_pattern_matching(pattern, genome, positions):
    assert replication.pattern_matching(pattern, genome) == positions


@pytest.mark.parametrize("genome, symbol, array", [
    ("AAAAGGGG", "A", {0: 4, 1: 3, 2: 2, 3: 1, 4: 0, 5: 1, 6: 2, 7: 3}),
    ("AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTCACGGATGCCAAGTGCATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTTTATCCAGCGCGTTAACCACGTTCTGTGCCGACTTT", "CC", {0: 7, 1: 7, 2: 7, 3: 7, 4: 7, 5: 7, 6: 7, 7: 6, 8: 6, 9: 6, 10: 6, 11: 6, 12: 6, 13: 6, 14: 6, 15: 6, 16: 6, 17: 5, 18: 5, 19: 5, 20: 4, 21: 4, 22: 4, 23: 4, 24: 4, 25: 3, 26: 3, 27: 3, 28: 3, 29: 3, 30: 3, 31: 3, 32: 3, 33: 3, 34: 3, 35: 3, 36: 3, 37: 3, 38: 3, 39: 3, 40: 3, 41: 3, 42: 3, 43: 2, 44: 2, 45: 2, 46: 2, 47: 2, 48: 2, 49: 2, 50: 3, 51: 3, 52: 3, 53: 3, 54: 3, 55: 3, 56: 3, 57: 3, 58: 2, 59: 2, 60: 2, 61: 2, 62: 3, 63: 3, 64: 3, 65: 3, 66: 3, 67: 3, 68: 3, 69: 3, 70: 3, 71: 3, 72: 3, 73: 3, 74: 3, 75: 3, 76: 4, 77: 4, 78: 4, 79: 4, 80: 4, 81: 4, 82: 4, 83: 4, 84: 4, 85: 4, 86: 5, 87: 5, 88: 5, 89: 6, 90: 6, 91: 6, 92: 6, 93: 6, 94: 7, 95: 7, 96: 7, 97: 7, 98: 7, 99: 7, 100: 7, 101: 7, 102: 7, 103: 7, 104: 6, 105: 6, 106: 6, 107: 7, 108: 7, 109: 7, 110: 7, 111: 7, 112: 8, 113: 8, 114: 8, 115: 8, 116: 7, 117: 7, 118: 7, 119: 7, 120: 7, 121: 7, 122: 7, 123: 7, 124: 7, 125: 7, 126: 7, 127: 8, 128: 7, 129: 7, 130: 7, 131: 7, 132: 7, 133: 7, 134: 7}),
    ])
def test_symbol_array(genome, symbol, array):
    assert replication.symbol_array(genome, symbol) == array


def test_faster_symbol_array(genome, symbol, array):
    assert replication.symbol_array(genome, symbol) == array

@pytest.mark.parametrize("genome, skew", [
    ("CATGGGCATCGGCCATACGCC", [0, -1, -1, -1, 0, 1, 2, 1, 1, 1, 0, 1, 2, 1, 0, 0, 0, 0, -1, 0, -1, -2]),
    ("AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTCACGGATGCCAAGTGCATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTTTATCCAGCGCGTTAACCACGTTCTGTGCCGACTTT", [0, 0, 1, 0, 1, 1, 2, 1, 0, 1, 1, 1, 1, 1, 1, 1, 2, 1, 0, 1, 0, -1, -1, 0, 0, -1, -2, -2, -1, -2, -2, -1, -2, -1, 0, 0, 1, 2, 1, 0, 0, -1, 0, -1, -2, -1, -1, -2, -2, -2, -3, -3, -4, -3, -2, -2, -2, -1, -2, -3, -3, -3, -2, -2, -1, -2, -2, -2, -2, -1, -1, 0, 1, 1, 1, 2, 1, 2, 2, 3, 2, 2, 2, 2, 3, 4, 4, 5, 6, 6, 6, 6, 5, 5, 5, 5, 4, 5, 4, 4, 4, 4, 4, 4, 3, 2, 2, 3, 2, 3, 2, 3, 3, 3, 3, 3, 2, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 2, 1, 0, 1, 1, 0, 0, 0, 0])
    ])
def test_skew_array(genome, skew):
    assert replication.skew_array(genome) == skew


@pytest.mark.parametrize("p, q, count", [
    ("GGGCCGTTGGT", "GGACCGTTGAC", 3),
    ("AAAA", "TTTT", 4),
    ("ACGTACGT", "TACGTACG", 8),
    ("ACGTACGT", "CCCCCCCC", 6),
    ("ACGTACGT", "TGCATGCA", 8),
    ("GATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACT", "AATAGCAGCTTCTCAACTGGTTACCTCGTATGAGTAAATTAGGTCATTATTGACTCAGGTCACTAACGTCT", 15),
    ("AGAAACAGACCGCTATGTTCAACGATTTGTTTTATCTCGTCACCGGGATATTGCGGCCACTCATCGGTCAGTTGATTACGCAGGGCGTAAATCGCCAGAATCAGGCTG", "AGAAACCCACCGCTAAAAACAACGATTTGCGTAGTCAGGTCACCGGGATATTGCGGCCACTAAGGCCTTGGATGATTACGCAGAACGTATTGACCCAGAATCAGGCTC", 28)
    ])
def test_hamming_distance(p, q, count):
    assert len(p) == len(q)
    assert replication.hamming_distance(p, q) == count


@pytest.mark.parametrize("pattern, text, d, positions", [
    ("ATTCTGGA", "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT", 3, [6, 7, 26, 27]),
    ("AAA", "TTTTTTAAATTTTAAATTTTTT", 2, [4, 5, 6, 7, 8, 11, 12, 13, 14, 15]),
    ("GAGCGCTGG", "GAGCGCTGGGTTAACTCGCTACTTCCCGACGAGCGCTGTGGCGCAAATTGGCGATGAAACTGCAGAGAGAACTGGTCATCCAACTGAATTCTCCCCGCTATCGCATTTTGATGCGCGCCGCGTCGATT", 2, [0, 30, 66]),
    ("AATCCTTTCA", "CCAAATCCCCTCATGGCATGCATTCCCGCAGTATTTAATCCTTTCATTCTGCATATAAGTAGTGAAGGTATAGAAACCCGTTCAAGCCCGCAGCGGTAAAACCGAGAACCATGATGAATGCACGGCGATTGCGCCATAATCCAAACA", 3, [3, 36, 74, 137]),
    ("CCGTCATCC", "CCGTCATCCGTCATCCTCGCCACGTTGGCATGCATTCCGTCATCCCGTCAGGCATACTTCTGCATATAAGTACAAACATCCGTCATGTCAAAGGGAGCCCGCAGCGGTAAAACCGAGAACCATGATGAATGCACGGCGATTGC", 3, [0, 7, 36, 44, 48, 72, 79, 112]),
    ("TTT", "AAAAAA", 3, [0, 1, 2, 3]),
    ("CCA", "CCACCT", 0, [0]),
    ("GTGCCG", "AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTCACGGATGCCAAGTGCATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTTTATCCAGCGCGTTAACCACGTTCTGTGCCGACTTT", 3, [3, 13, 16, 22, 25, 27, 28, 30, 33, 34, 36, 39, 47, 54, 61, 71, 76, 84, 87, 91, 101, 106, 119, 124])
    ])
def test_approximate_pattern_matching(pattern, text, d, positions):
    assert replication.approximate_pattern_matching(pattern, text, d) == positions


@pytest.mark.parametrize("pattern, text, d, count", [
    ("GAGG", "TTTAGAGCCTTCAGAGG", 2, 4),
    ("AA", "AAA", 0, 2),
    ("ATA", "ATA", 1, 1),
    ("GTGCCG", "AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTCACGGATGCCAAGTGCATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTTTATCCAGCGCGTTAACCACGTTCTGTGCCGACTTT", 3, 24)
    ])
def test_approximate_pattern_count(pattern, text, d, count):
    assert replication.approximate_pattern_count(pattern, text, d) == count

