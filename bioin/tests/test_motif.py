from .. import motif
import pytest
import sys


@pytest.mark.parametrize("motifs, count", [
    (["AACGTA", "CCCGTT", "CACCTT", "GGATTA", "TTCCGG"], {'A': [1, 2, 1, 0, 0, 2], 'C': [2, 1, 4, 2, 0, 0], 'G': [1, 1, 0, 2, 1, 1], 'T': [1, 1, 0, 1, 4, 2]}),
    (["GTACAACTGT", "CAACTATGAA", "TCCTACAGGA", "AAGCAAGGGT", "GCGTACGACC", "TCGTCAGCGT", "AACAAGGTCA", "CTCAGGCGTC", "GGATCCAGGT", "GGCAAGTACC"], {'A': [2, 3, 3, 3, 6, 4, 2, 2, 1, 3], 'C': [2, 3, 4, 3, 2, 3, 2, 1, 3, 3], 'G': [4, 2, 3, 0, 1, 3, 4, 5, 5, 0], 'T': [2, 2, 0, 4, 1, 0, 2, 2, 1, 4]})
    ])
def test_count_motif(motifs, count):
    assert motif.count_motif(motifs) == count


@pytest.mark.parametrize("motifs, count", [
    (["AACGTA", "CCCGTT", "CACCTT", "GGATTA", "TTCCGG"], {'A': [0.2, 0.4, 0.2, 0.0, 0.0, 0.4], 'C': [0.4, 0.2, 0.8, 0.4, 0.0, 0.0], 'G': [0.2, 0.2, 0.0, 0.4, 0.2, 0.2], 'T': [0.2, 0.2, 0.0, 0.2, 0.8, 0.4]}),
    (["GTACAACTGT", "CAACTATGAA", "TCCTACAGGA", "AAGCAAGGGT", "GCGTACGACC", "TCGTCAGCGT", "AACAAGGTCA", "CTCAGGCGTC", "GGATCCAGGT", "GGCAAGTACC"], {'A': [0.2, 0.3, 0.3, 0.3, 0.6, 0.4, 0.2, 0.2, 0.1, 0.3], 'C': [0.2, 0.3, 0.4, 0.3, 0.2, 0.3, 0.2, 0.1, 0.3, 0.3], 'G': [0.4, 0.2, 0.3, 0.0, 0.1, 0.3, 0.4, 0.5, 0.5, 0.0], 'T': [0.2, 0.2, 0.0, 0.4, 0.1, 0.0, 0.2, 0.2, 0.1, 0.4]})
    ])
def test_profile_motif(motifs, count):
    assert motif.profile_motif(motifs) == count


@pytest.mark.parametrize("motifs, consensus", [
    (["AACGTA", "CCCGTT", "CACCTT", "GGATTA", "TTCCGG"], "CACCTA"),
    (["GTACAACTGT", "CAACTATGAA", "TCCTACAGGA", "AAGCAAGGGT", "GCGTACGACC", "TCGTCAGCGT", "AACAAGGTCA", "CTCAGGCGTC", "GGATCCAGGT", "GGCAAGTACC"], "GACTAAGGGT")
    ])
def test_consensus_motif(motifs, consensus):
    assert motif.consensus_motif(motifs) == consensus


@pytest.mark.parametrize("motifs, score", [
    (["AACGTA", "CCCGTT", "CACCTT", "GGATTA", "TTCCGG"], 14),
    (["GTACAACTGT", "CAACTATGAA", "TCCTACAGGA", "AAGCAAGGGT", "GCGTACGACC", "TCGTCAGCGT", "AACAAGGTCA", "CTCAGGCGTC", "GGATCCAGGT", "GGCAAGTACC"], 57)
    ])
def test_score_motif(motifs, score):
    assert motif.score_motif(motifs) == score


@pytest.mark.parametrize("text, motif_profile, probability", [
    ("ACGGGGATTACC", {
    'A': [0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0],
    'C': [0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6],
    'G': [0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
    'T': [0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4]
}, 0.000839808),
    ("TCGGGGGCCACC", {
    'A': [0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0],
    'C': [0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6],
    'G': [0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
    'T': [0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4]
}, 3.26592e-05)
])
def test_probability_profile(text, motif_profile, probability):
    assert motif.probability_profile(text, motif_profile) == pytest.approx(probability, abs=1.e-3)


@pytest.mark.parametrize("text, k, motif_profile, kmer", [
    ("ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT", 5, {'A': [0.2, 0.2, 0.3, 0.2, 0.3], 'C': [0.4, 0.3, 0.1, 0.5, 0.1], 'G': [0.3, 0.3, 0.5, 0.2, 0.4], 'T': [0.1, 0.2, 0.1, 0.1, 0.2]}, "CCGAG"),
    ("AGCAGCTTTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATCTGAACTGGTTACCTGCCGTGAGTAAAT", 8, {'A': [0.7, 0.2, 0.1, 0.5, 0.4, 0.3, 0.2, 0.1], 'C': [0.2, 0.2, 0.5, 0.4, 0.2, 0.3, 0.1, 0.6], 'G': [0.1, 0.3, 0.2, 0.1, 0.2, 0.1, 0.4, 0.2], 'T': [0.0, 0.3, 0.2, 0.0, 0.2, 0.3, 0.3, 0.1]}, "AGCAGCTT"),
    ("TTACCATGGGACCGCTGACTGATTTCTGGCGTCAGCGTGATGCTGGTGTGGATGACATTCCGGTGCGCTTTGTAAGCAGAGTTTA", 12, {'A': [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.1, 0.2, 0.3, 0.4, 0.5], 'C': [0.3, 0.2, 0.1, 0.1, 0.2, 0.1, 0.1, 0.4, 0.3, 0.2, 0.2, 0.1], 'G': [0.2, 0.1, 0.4, 0.3, 0.1, 0.1, 0.1, 0.3, 0.1, 0.1, 0.2, 0.1], 'T': [0.3, 0.4, 0.1, 0.1, 0.1, 0.1, 0.0, 0.2, 0.4, 0.4, 0.2, 0.3]}, "AAGCAGAGTTTA"),
    ("AACCGGTT", 3, {'A': [1.0, 1.0, 1.0], 'C': [0.0, 0.0, 0.0], 'G': [0.0, 0.0, 0.0], 'T': [0.0, 0.0, 0.0]}, "AAC"),
    ("TTACCATGGGACCGCTGACTGATTTCTGGCGTCAGCGTGATGCTGGTGTGGATGACATTCCGGTGCGCTTTGTAAGCAGAGTTTA", 5, {'A': [0.2, 0.2, 0.3, 0.2, 0.3], 'C': [0.4, 0.3, 0.1, 0.5, 0.1], 'G': [0.3, 0.3, 0.5, 0.2, 0.4], 'T': [0.1, 0.2, 0.1, 0.1, 0.2]}, "CAGCG")
])
def test_profile_most_probable_kmer(text, k, motif_profile, kmer):
    assert motif.profile_most_probable_kmer(text, k, motif_profile) == kmer


@pytest.mark.parametrize("dna, k, t, best_motifs", [
    (["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"], 3, 5, ["CAG", "CAG", "CAA", "CAA", "CAA"]),
    (["GCCCAA", "GGCCTG", "AACCTA", "TTCCTT"], 3, 4, ["GCC", "GCC", "AAC", "TTC"]),
    (["GAGGCGCACATCATTATCGATAACGATTCGCCGCATTGCC", "TCATCGAATCCGATAACTGACACCTGCTCTGGCACCGCTC", "TCGGCGGTATAGCCAGAAAGCGTAGTGCCAATAATTTCCT", "GAGTCGTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTG", "GACGGCAACTACGGTTACAACGCAGCAACCGAAGAATATT", "TCTGTTGTTGCTAACACCGTTAAAGGCGGCGACGGCAACT", "AAGCGGCCAACGTAGGCGCGGCTTGGCATCTCGGTGTGTG", "AATTGAAAGGCGCATCTTACTCTTTTCGCTTTCAAAAAAA"], 5, 8, ["GAGGC", "TCATC", "TCGGC", "GAGTC", "GCAGC", "GCGGC", "GCGGC", "GCATC"]),
    (["GCAGGTTAATACCGCGGATCAGCTGAGAAACCGGAATGTGCGT", "CCTGCATGCCCGGTTTGAGGAACATCAGCGAAGAACTGTGCGT", "GCGCCAGTAACCCGTGCCAGTCAGGTTAATGGCAGTAACATTT", "AACCCGTGCCAGTCAGGTTAATGGCAGTAACATTTATGCCTTC", "ATGCCTTCCGCGCCAATTGTTCGTATCGTCGCCACTTCGAGTG"], 6, 5, ["GTGCGT", "GTGCGT", "GCGCCA", "GTGCCA", "GCGCCA"]),
    (["GACCTACGGTTACAACGCAGCAACCGAAGAATATTGGCAA", "TCATTATCGATAACGATTCGCCGGAGGCCATTGCCGCACA", "GGAGTCTGGTGAAGTGTGGGTTATGGGGCAGACTGGGAAA", "GAATCCGATAACTGACACCTGCTCTGGCACCGCTCTCATC", "AAGCGCGTAGGCGCGGCTTGGCATCTCGGTGTGTGGCCAA", "AATTGAAAGGCGCATCTTACTCTTTTCGCTTAAAATCAAA", "GGTATAGCCAGAAAGCGTAGTTAATTTCGGCTCCTGCCAA", "TCTGTTGTTGCTAACACCGTTAAAGGCGGCGACGGCAACT"], 5, 8, ["GCAGC", "TCATT", "GGAGT", "TCATC", "GCATC", "GCATC", "GGTAT", "GCAAC"]),
    (["GACCTACGGTTACAACGCAGCAACCGAAGAATATTGGCAA", "TCATTATCGATAACGATTCGCCGGAGGCCATTGCCGCACA",
      "GGAGTCTGGTGAAGTGTGGGTTATGGGGCAGACTGGGAAA", "GAATCCGATAACTGACACCTGCTCTGGCACCGCTCTCATC",
      "AAGCGCGTAGGCGCGGCTTGGCATCTCGGTGTGTGGCCAA", "AATTGAAAGGCGCATCTTACTCTTTTCGCTTAAAATCAAA",
      "GGTATAGCCAGAAAGCGTAGTTAATTTCGGCTCCTGCCAA", "TCTGTTGTTGCTAACACCGTTAAAGGCGGCGACGGCAACT"], 4, 8, ["CGCA", "CGCA", "GGAG", "GGCA", "GGCA", "CGCA", "GGTA", "GGCA"])
])
def test_greedy_motif_search(dna, k, t, best_motifs):
    assert motif.greedy_motif_search(dna, k, t) == best_motifs


@pytest.mark.parametrize("motifs, count", [
    (["AACGTA", "CCCGTT", "CACCTT", "GGATTA", "TTCCGG"], {'A': [2, 3, 2, 1, 1, 3], 'C': [3, 2, 5, 3, 1, 1], 'T': [2, 2, 1, 2, 5, 3], 'G': [2, 2, 1, 3, 2, 2]}),
    (["GTACAACTGT", "CAACTATGAA", "TCCTACAGGA", "AAGCAAGGGT", "GCGTACGACC", "TCGTCAGCGT", "AACAAGGTCA", "CTCAGGCGTC", "GGATCCAGGT", "GGCAAGTACC"], {'A': [3, 4, 4, 4, 7, 5, 3, 3, 2, 4], 'C': [3, 4, 5, 4, 3, 4, 3, 2, 4, 4], 'T': [3, 3, 1, 5, 2, 1, 3, 3, 2, 5], 'G': [5, 3, 4, 1, 2, 4, 5, 6, 6, 1]})
])
def test_count_with_pseudocount(motifs, count):
    assert motif.count_with_pseudocount(motifs) == count


@pytest.mark.parametrize("motifs, profile", [
    (["AACGTA", "CCCGTT", "CACCTT", "GGATTA", "TTCCGG"], {'A': [0.2222222222222222, 0.3333333333333333, 0.2222222222222222, 0.1111111111111111, 0.1111111111111111, 0.3333333333333333], 'C': [0.3333333333333333, 0.2222222222222222, 0.5555555555555556, 0.3333333333333333, 0.1111111111111111, 0.1111111111111111], 'T': [0.2222222222222222, 0.2222222222222222, 0.1111111111111111, 0.2222222222222222, 0.5555555555555556, 0.3333333333333333], 'G': [0.2222222222222222, 0.2222222222222222, 0.1111111111111111, 0.3333333333333333, 0.2222222222222222, 0.2222222222222222]}),
    (["GTACAACTGT", "CAACTATGAA", "TCCTACAGGA", "AAGCAAGGGT", "GCGTACGACC", "TCGTCAGCGT", "AACAAGGTCA", "CTCAGGCGTC", "GGATCCAGGT", "GGCAAGTACC"], {'A': [0.21428571428571427, 0.2857142857142857, 0.2857142857142857, 0.2857142857142857, 0.5, 0.35714285714285715, 0.21428571428571427, 0.21428571428571427, 0.14285714285714285, 0.2857142857142857], 'C': [0.21428571428571427, 0.2857142857142857, 0.35714285714285715, 0.2857142857142857, 0.21428571428571427, 0.2857142857142857, 0.21428571428571427, 0.14285714285714285, 0.2857142857142857, 0.2857142857142857], 'T': [0.21428571428571427, 0.21428571428571427, 0.07142857142857142, 0.35714285714285715, 0.14285714285714285, 0.07142857142857142, 0.21428571428571427, 0.21428571428571427, 0.14285714285714285, 0.35714285714285715], 'G': [0.35714285714285715, 0.21428571428571427, 0.2857142857142857, 0.07142857142857142, 0.14285714285714285, 0.2857142857142857, 0.35714285714285715, 0.42857142857142855, 0.42857142857142855, 0.07142857142857142]})
])
def test_profile_with_pseudocount(motifs, profile):
    assert motif.profile_with_pseudocount(motifs) == profile


@pytest.mark.parametrize("dna, k, t, best_motifs", [
    (["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"], 3, 5, ["TTC", "ATC", "TTC", "ATC", "TTC"]),
    (["AGGCGGCACATCATTATCGATAACGATTCGCCGCATTGCC", "ATCCGTCATCGAATAACTGACACCTGCTCTGGCACCGCTC", "AAGCGTCGGCGGTATAGCCAGATAGTGCCAATAATTTCCT", "AGTCGGTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTG", "AACCGGACGGCAACTACGGTTACAACGCAGCAAGAATATT", "AGGCGTCTGTTGTTGCTAACACCGTTAAGCGACGGCAACT", "AAGCGGCCAACGTAGGCGCGGCTTGGCATCTCGGTGTGTG", "AATTGAAAGGCGCATCTTACTCTTTTCGCTTTCAAAAAAA"], 5, 8, ["AGGCG", "ATCCG", "AAGCG", "AGTCG", "AACCG", "AGGCG", "AGGCG", "AGGCG"]),
    (["GCACATCATTAAACGATTCGCCGCATTGCCTCGATAGGCG", "TCATAACTGACACCTGCTCTGGCACCGCTCATCCGTCGAA", "AAGCGGGTATAGCCAGATAGTGCCAATAATTTCCTTCGGC", "AGTCGGTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTG", "AACCGGACGGCAACTACGGTTACAACGCAGCAAGAATATT", "AGGCGTCTGTTGTTGCTAACACCGTTAAGCGACGGCAACT", "AAGCTTCCAACATCGTCTTGGCATCTCGGTGTGTGAGGCG", "AATTGAACATCTTACTCTTTTCGCTTTCAAAAAAAAGGCG"], 5, 8, ["AGGCG", "TGGCA", "AAGCG", "AGGCA", "CGGCA", "AGGCG", "AGGCG", "AGGCG"]),
    (["GCACATCATTATCGATAACGATTCATTGCCAGGCGGCCGC", "TCATCGAATAACTGACACCTGCTCTGGCTCATCCGACCGC", "TCGGCGGTATAGCCAGATAGTGCCAATAATTTCCTAAGCG", "GTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTGAGTCG", "GACGGCAACTACGGTTACAACGCAGCAAGAATATTAACCG", "TCTGTTGTTGCTAACACCGTTAAGCGACGGCAACTAGGCG", "GCCAACGTAGGCGCGGCTTGGCATCTCGGTGTGTGAAGCG", "AAAGGCGCATCTTACTCTTTTCGCTTTCAAAAAAAAATTG"], 5, 8, ["GGCGG", "GGCTC", "GGCGG", "GGCAG", "GACGG", "GACGG", "GGCGC", "GGCGC"]),
    (["GCACATCATTATCGATAACGATTCATTGCCAGGCGGCCGC", "TCATCGAATAACTGACACCTGCTCTGGCTCATCCGACCGC", "TCGGCGGTATAGCCAGATAGTGCCAATAATTTCCTAAGCG", "GTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTGAGTCG", "GACGGCAACTACGGTTACAACGCAGCAAGAATATTAACCG", "TCTGTTGTTGCTAACACCGTTAAGCGACGGCAACTAGGCG", "GCCAACGTAGGCGCGGCTTGGCATCTCGGTGTGTGAAGCG", "AAAGGCGCATCTTACTCTTTTCGCTTTCAAAAAAAAATTG"], 3, 8, ["GGC", "GGC", "GGC", "GGC", "GGC", "GGC", "GGC", "GGC"])
])
def test_greedy_motif_search_with_pseudocount(dna, k, t, best_motifs):
    assert motif.greedy_motif_search_with_pseudocount(dna, k, t) == best_motifs


@pytest.mark.parametrize("motif_profile, dna, motifs", [
    ({"A" : [0.8, 0.0, 0.0, 0.2],
    "C" : [0.0, 0.6, 0.2, 0.0],
    "G" : [0.2, 0.2, 0.8, 0.0],
    "T" : [0.0, 0.2, 0.0, 0.8]}, ['TTACCTTAAC', 'GATGTCTGTC', 'ACGGCGTTAG', 'CCCTAACGAG', 'CGTCAGAGGT'],
     ['ACCT', 'ATGT', 'GCGT', 'ACGA', 'AGGT']),
    ({ "A" : [0.5, 0.0, 0.2, 0.2],
       "C" : [0.3, 0.6, 0.2, 0.0],
       "G" : [0.2, 0.2, 0.6, 0.0],
       "T" : [0.0, 0.2, 0.0, 0.8]}, ['TTACCTTAAC', 'GATGTCTGTC', 'ACGGCGTTAG', 'CCCTAACGAG', 'CGTCAGAGGT'],
     ['ACCT', 'ATGT', 'GCGT', 'ACGA', 'AGGT'])
])
def test_profile_probable_motifs(motif_profile, dna, motifs):
    assert motif.profile_probable_motifs(motif_profile, dna) == motifs

