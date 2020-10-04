from bioin import pattern_count
# import matplotlib.pyplot as plt


def pattern_matching(pattern, genome):
    """Find all occurrences of a pattern in a string.

    Args:
        pattern (str): pattern string to search in the genome string.
        genome (str): search space for pattern.

    Returns:
        List, list of int, i.e. all starting positions in genome where pattern appears as a substring.

    Examples:
        Find all the starting positions for a pattern string in the genome string.
        
        1.

        >>> pattern = 'ATAT'
        >>> genome = 'GATATATGCATATACTT'
        >>> positions = pattern_matching(pattern, genome)
        >>> positions
            [1, 3, 9]

        2.

        >>> pattern = 'CTTGATCAT'
        >>> genome = 'CTTGATCATCTTGATCATCTTGATCAT'
        >>> positions = pattern_matching(pattern, genome)
        >>> positions
            [0, 9, 18]
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
    """Calculate the symbol count in genome. Background: Analyzing a genome’s half-strands: Although most bacteria have circular genomes, we have thus far assumed that genomes were linear, a reasonable simplifying assumption because the length of the window is much shorter than the length of the genome. This time, because we are sliding a giant window, we should account for windows that “wrap around” the end of Genome. To do so, we will define a string ExtendedGenome as Genome+Genome[0:n//2]. That is, we copy the first len(Genome)//2 nucleotides of Genome to the end of the string. For example, this genome: 
    
    CTGCTTCGCCCGCCGGACCGGCCTCGTGATGGGGT_CTGCTTCGCCCGCCGGA

    A DNA string Genome (shown before the underscore) containing 35 nucleotides that is extended by its first 17 nucleotides (shown after the underscore) to yield ExtendedGenome (the ExtendeGenome doesn't contain the underscore).

    We will keep track of the total number of occurrences of 'C' that we encounter in each window of ExtendedGenome by using a symbol array. The i-th element of the symbol array is equal to the number of occurrences of the symbol in the window of length len(Genome)//2 starting at position i of ExtendedGenome. For example, array[0] equals the number of A count in the extendedGenome from position index 0 to 0+4=4, i.e. 'AAAAG', there are 4 'A' in it, next from position 1 to 5, i.e. 'AAAGG' there are 3 'A' in it, and so forth. Finally return the key-value pair of all the i: array[i] in a dictionary.

    ExtendedGenome A A A A G G G G A A A A

    xxxxxxxxxxxxxxx i 0 1 2 3 4 5 6 7 

    xxxxxxxxx array[i] 4 3 2 1 0 1 2 3 

    The symbol array for Genome equal to "AAAAGGGG" and symbol equal to "A".

    Args:
        genome (str): a DNA string as the search space.
        symbol (str): the single base to query in the search space.

    Returns:
        Dictionary, a dictionary, position-counts pairs of symbol in each genome sliding window.

    Examples:
        The symbol array for genome equal to "AAAAGGGG" and symbol equal to "A".

        >>> genome = 'AAAAGGGG'
        >>> symbol = 'A'
        >>> position_symbolcount_dict = symbol_array(genome, symbol)
        >>> position_symbolcount_dict
            {0: 4, 1: 3, 2: 2, 3: 1, 4: 0, 5: 1, 6: 2, 7: 3}
    """
    array = {}
    n = len(genome)
    extended_genome = genome + genome[0:n//2]
    for i in range(n):
        array[i] = pattern_count(symbol, extended_genome[i:i+(n//2)])
    return array


def faster_symbol_array(genome, symbol):
    """A faster calculation method for counting a symbol in genome.

    Args:
        genome (str): a DNA string as the search space.
        symbol (str): the single base to query in the search space.

    Returns:
        Dictionary, a dictionary, position-counts pairs of symbol in each genome sliding window.

    Examples:
        The symbol array for genome equal to "AAAAGGGG" and symbol equal to "A".

        >>> genome = 'AAAAGGGG'
        >>> symbol = 'A'
        >>> position_symbolcount_dict = symbol_array(genome, symbol)
        >>> position_symbolcount_dict
            {0: 4, 1: 3, 2: 2, 3: 1, 4: 0, 5: 1, 6: 2, 7: 3}
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
    """Compute the skew array of genome as a list. In the table containing nucleotide counts for T. petrophila (reproduced below), we noted that not just C but also G has peculiar statistics on the forward and reverse half-strands.

    xxxxxxxxxxxxxxxxxxxx   #C          #G          #A          #T

    xxxxx Entire strand  427419          413241          491488          491363

    Reverse half-strand    219518          201634          243963          246641
    
    Forward half-strand    207901          211607          247525          244722
    
    xxxxxxx Difference     +11617          -9973          -3362          +1919
    
    In practice, scientists use a more accurate approach that accounts for both G and C when searching for ori. As the above figure illustrates, the difference between the total amount of guanine and the total amount of cytosine is negative on the reverse half-strand and positive on the forward half-strand.

    We will keep track of the difference between the total number of occurrences of G and the total number of occurrences of C that we have encountered so far in Genome by using a skew array. This array, denoted Skew, is defined by setting Skew[i] equal to the number of occurrences of G minus the number of occurrences of C in the first i nucleotides of Genome (see figure below). We also set Skew[0] equal to zero.

    The array Skew for Genome = "CATGGGCATCGGCCATACGCC".

    xxxxxx i 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21

    Skew[i] 0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2
    
    xxxx Text C A T G G G C A T C G G C C A T A C G C C 

    Given a string Genome, we can form its skew array by setting Skew[0] equal to 0, and then ranging through the genome.  At position i of Genome, if we encounter an A or a T, we set Skew[i+1] equal to Skew[i]; if we encounter a G, we set Skew[i+1] equal to Skew[i]+1; if we encounter a C, we set Skew[i+1] equal to Skew[i]-1.

    Args:
        genome (str): a DNA string genome.

    Returns:
        List, the i-th element is skew[i], which equals to the number of G minus the number of C, in the first i nucleotides of the genome, set skew[0]=0.

    Examples:
        The array Skew for Genome = "CATGGGCATCGGCCATACGCC".
        
        >>> genome = 'CATGGGCATCGGCCATACGCC'
        >>> array_skew = skew_array(genome)
        >>> array_skew
            [0, -1, -1, -1, 0, 1, 2, 1, 1, 1, 0, 1, 2, 1, 0, 0, 0, 0, -1, 0, -1, -2]
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
    """Find a position in a genome where the skew diagram attains a minimum.

    Args:
        genome (str): a DNA string genome.

    Returns:
        List, the position (integers i) where the minimum values (skew[i]) are in the skew array.

    Examples:
        Taking a DNA string Genome as input and returning all integers i minimizing Skew[i] for Genome.

        >>> genome = 'TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT'
        >>> array_skew_min_position = minimum_skew(genome)
        >>> array_skew_min_position = minimum_skew(genome)
        >>> array_skew
            [11, 24]
    """
    positions = []
    skew = skew_array(genome)
    m = min(skew)
    for i in range(len(skew)):
        if skew[i] == m:
            positions.append(i)
    return positions


def hamming_distance(p, q):
    """Calculate the HammingDistance for two strings. We say that position i in k-mers p and q is a mismatch if the symbols at position i of the two strings are not the same. The total number of mismatches between strings p and q is called the Hamming distance between these strings. We will let you implement a function to compute this distance, called HammingDistance(p, q).

    Args:
        p (str): the first DNA string.
        q (str): the second DNA string, p and q of the equal length.

    Returns:
        Integer, number of different base count between p and q, i.e. the Hamming distance between these strings.

    Examples:
        Solving the Hamming distance of two DNA Genomes.

        >>> p = 'GGGCCGTTGGT'
        >>> q = 'GGACCGTTGAC'
        >>> hammingdistance = hamming_distance(p, q)
        >>> hammingdistance
            3
    """
    count = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            count += 1
    return count


def approximate_pattern_matching(pattern, text, d):
    """Find all approximate occurences of a pattern in a string. We say that a k-mer Pattern appears as a substring of Text with at most d mismatches if there is some k-mer substring Pattern' of Text having d or fewer mismatches with Pattern; that is, HammingDistance(Pattern, Pattern') ≤ d. Our observation that a DnaA box may appear with slight variations leads to the following generalization of the Pattern Matching Problem.

    Args:
        text (str): a DNA string genome.
        pattern (str): a substring of DNA genome.
        d (int): at most d mismatches between pattern and a substring in text.

    Returns:
        List, all starting positions where pattern appears as a substring of text with at most d mismatches.

    Examples:
        Positions that a pattern appears as a substring of text with at most d mismatches fulfilled.

        >>> pattern = 'ATTCTGGA'
        >>> text = 'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT'
        >>> d = 3
        >>> positions = approximate_pattern_matching(pattern, text, d)
        >>> positions
            [6, 7, 26, 27]
    """
    positions = []
    for i in range(len(text)-len(pattern)+1):
        if hamming_distance(pattern, text[i:i+len(pattern)]) <= d:
            positions.append(i)
    return positions


def approximate_pattern_count(pattern, text, d):
    """Compute the number of occurrences of pattern in text with at most d mismatches. Given input strings Text and Pattern as well as an integer d, we extend the definition of pattern_count to the function approximate_pattern_count(Pattern, Text, d). This function computes the number of occurrences of Pattern in Text with at most d mismatches. For example, approximate_pattern_count('AAAAA', 'AACAAGCATAAACATTAAAGAG', 1) = 4.

    This is because AAAAA appears four times in this string with at most one mismatch: AACAA, ATAAA, AAACA, and AAAGA. Notice that two of these occurrences overlap.

    Args:
        pattern (str): a sub DNA string.
        text (str): a DNA string.
        d (int): the number of maximum mismatches.

    Returns:
        Integer, the number of occurrences of pattern in text with at most d mismatches.

    Examples:
        The number of times Pattern appears in Text with at most d mismatches.

        >>> pattern = 'GAGG'
        >>> text = 'TTTAGAGCCTTCAGAGG'
        >>> d = 2
        >>> approx_count = approximate_pattern_count(pattern, text, d)
        >>> approx_count
            4
    """
    count = 0
    for i in range(len(text) - len(pattern) + 1):
        if hamming_distance(pattern, text[i:i + len(pattern)]) <= d:
            count = count + 1
    return count


# print(skew_array("CATTCCAGTACTTCGATGATGGCGTGAAGA"))
# print(hamming_distance("TGACCCGTTATGCTCGAGTTCGGTCAGAGCGTCATTGCGAGTAGTCGTTTGCTTTCTCAAACTCC", "GAGCGATTAAGCGTGACAGCCCCAGGGAACCCACAAAACGTGATCGCAGTCCATCCGATCATACA"))

# print(symbol_array("AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTCACGGATGCCAAGTGCATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTTTATCCAGCGCGTTAACCACGTTCTGTGCCGACTTT", "C"))