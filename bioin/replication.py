from bioin import pattern_count
# import matplotlib.pyplot as plt


def pattern_matching(pattern, genome):
    """
    Find all occurrences of a pattern in a string

    Args:
        pattern: string, pattern in the genome
        genome: string, search space for pattern

    Returns:
        list of int, all starting positions in genome where pattern appears as a substring

    """
    positions = [] # output variable
    for i in range(len(genome)-len(pattern)+1):
        if genome[i:i+len(pattern)] == pattern:
            positions.append(i)
    return positions


# given integers L and t, a k-mer Pattern forms an (L, t)-clump inside a (longer) string Genome
# if there is an interval of Genome of length L in which this k-mer appears at least t times.
# (This definition assumes that the k-mer completely fits within the interval.)
# For example, "TGCA" forms a (25, 3)-clump in the following Genome:
# gatcagcataagggtccCTGCAATGCATGACAAGCCTGCAGTtgttttac


def symbol_array(genome, symbol):
    """
    calculate the symbol count in genome

    Args:
        genome: string, a DNA string as the search space
        symbol: string, the single base to look for

    Returns:
        a dictionary, counts of symbol in genome

    """
    array = {}
    n = len(genome)
    extended_genome = genome + genome[0:n//2]
    for i in range(n):
        array[i] = pattern_count(symbol, extended_genome[i:i+(n//2)])
    return array


def faster_symbol_array(genome, symbol):
    """
    A more faster calculation method for counting a symbol in genome

    Args:
        genome: string, a DNA string as the search space
        symbol: string, the single base to look for

    Returns:
        a dictionary, counts of symbol in genome

    """
    array = {}
    n = len(genome)
    extended_genome = genome + genome[0:n//2]

    # look at the first half of Genome to compute first array value
    array[0] = pattern_count(symbol, genome[0:n//2])

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]

        # the current array value can differ from the previous array value by at most 1
        if extended_genome[i-1] == symbol:
            array[i] = array[i]-1
        if extended_genome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array


# with open('../data/e_coli.txt') as file:
    # e_coli = file.read()

# array = faster_symbol_array(e_coli, "C")

# plt.plot(*zip(*sorted(array.items())))
# plt.show()


def skew_array(genome):
    """
    Compute the skew array of genome as a list

    Args:
        genome: string, a DNA string genome

    Returns:
        list, the i-th element is skew[i], which equals to the number of G minus the number of C,
        in the first i nucleotides of the genome, set skew[0]=0.

    """
    skew = [0]
    # normal way:
    # for i in range(1, len(genome)+1):
    #     g_count = pattern_count('G', genome[0:i])
    #     c_count = pattern_count('C', genome[0:i])
    #     skew.append(g_count-c_count)
    # return skew
    # another more simpler way:
    score = {"A": 0, "T": 0, "C": -1, "G": 1}
    for i in range(1,len(genome)+1):
        skew.append(score[genome[i-1]] + skew[i-1])
    return skew


def minimum_skew(genome):
    """

    Args:
        genome: string, a DNA string genome

    Returns:
        list, the position where the minimum values are in the skew array

    """
    positions = []
    skew = skew_array(genome)
    m = min(skew)
    for i in range(len(skew)):
        if skew[i] == m:
            positions.append(i)
    return positions


def hamming_distance(p, q):
    """
    Calculate the HammingDistance for two given strings

    Args:
        p: string, the first DNA string
        q: string, the second DNA string, p and q should have the same length

    Returns:
        integer, number of different base count between p and q

    """
    count = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            count += 1
    return count


def approximate_pattern_matching(pattern, text, d):
    """

    Args:
        text: string, a DNA string genome
        pattern: string, a substring of DNA genome
        d: integer, at most d mismatches between pattern and a substring in text

    Returns:
        list, all starting positions where pattern appears as a substring of text with at most d mismatches

    """
    positions = []
    for i in range(len(text)-len(pattern)+1):
        if hamming_distance(pattern, text[i:i+len(pattern)]) <= d:
            positions.append(i)
    return positions


def approximate_pattern_count(pattern, text, d):
    """
    Compute the number of occurrences of pattern in text with at most d mismatches

    Args:
        pattern: string, a sub DNA string
        text: string, a DNA string
        d: the number of maximum mismatches

    Returns:
        integer, the number of occurrences of pattern in text with at most d mismatches

    """
    count = 0
    for i in range(len(text) - len(pattern) + 1):
        if hamming_distance(pattern, text[i:i + len(pattern)]) <= d:
            count = count + 1
    return count


print(skew_array("CATTCCAGTACTTCGATGATGGCGTGAAGA"))
print(hamming_distance("TGACCCGTTATGCTCGAGTTCGGTCAGAGCGTCATTGCGAGTAGTCGTTTGCTTTCTCAAACTCC", "GAGCGATTAAGCGTGACAGCCCCAGGGAACCCACAAAACGTGATCGCAGTCCATCCGATCATACA"))

print(symbol_array("AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTCACGGATGCCAAGTGCATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTTTATCCAGCGCGTTAACCACGTTCTGTGCCGACTTT", "C"))