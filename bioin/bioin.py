"""
bioin.py
Study notes for beginners bioinformatics

Handles the primary functions
"""
# import sys
# lines = sys.stdin.read().splitlines()


def canvas(with_attribution=True):
    """Placeholder function to show example docstring (NumPy format)

    Replace this function and doc string for your own project

    Parameters
    ----------
    with_attribution : bool, Optional, default: True
        Set whether or not to display who the quote is from

    Returns
    -------
    quote : str
        Compiled string including quote and optional attribution
    """

    quote = "The code is but a canvas to our imagination."
    if with_attribution:
        quote += "\n\t- Adapted from Henry David Thoreau"
    return quote


def pattern_count(pattern, text):
    """
    The number of times that a pattern appears as a substring of text.

    Args:
        pattern (str): the substring pattern to find in the given text.
        text (str): the string space for looking.

    Returns:
        String, number of substring pattern that appears in text.

    Examples:
        Count the frequency (overlapping occurrences also counts) of a substring, i.e. pattern in the given string, i.e. text.

        >>> pattern = "GCG" 
        >>> text = "GCGCG"
        >>> count = pattern_count(pattern, text)
        >>> count
            2
    """
    count = 0
    for i in range(len(text) - len(pattern) + 1):
        if text[i:i + len(pattern)] == pattern:
            count = count + 1
    return count


# print(pattern_count(lines[1], lines[0]))


# The Frequent Words Problem
# We say that Pattern is a most frequent k-mer in Text if it maximizes PatternCount(Pattern, Text) among all k-mers.
# You can verify that "ACTAT" is a most frequent 5-mer for Text = "ACAACTATGCATACTATCGGGAACTATCCT",
# and "ATA" is a most frequent 3-mer for Text = "CGATATATCCATAG".


def frequency_map(text, k):
    """
    Find the frequency of all k-mers in a string.

    Args:
        text (str): text.
        k (int): length of the substring (i.e. kmers).

    Returns:
        Dictionary, a dictionary that contains the count of all the k-mers in text.

    Examples:
        Computes the frequency map of a given string (i.e. text) and integer (i.e. k). Return a dictionary of the k-mers and the corresponding frequency for all k-mers that appears in text.
        
        >>> text = "CGATATATCCATAG" 
        >>> k = 3
        >>> kmers_count_map = frequency_map(text, k)
        >>> kmers_count_map
            {'CGA': 1, 'GAT': 1, 'ATA': 3, 'TAT': 2, 'ATC': 1, 'TCC': 1, 'CCA': 1, 'CAT': 1, 'TAG': 1}
    """
    freq = {}
    n = len(text)
    for i in range(n-k+1):
        pattern = text[i:i+k]
        freq[pattern] = 0
        for m in range(n-k+1):
            if text[m:m+k] == pattern:
                freq[pattern] = freq[pattern] + 1
    return freq


def frequent_words(text, k):
    """Find all the most frequent k-mers in text. Depend on function frequency_map.

    Args:
        text (str): text.
        k (int): length of the substring (i.e. kmers).

    Returns:
        List, a list that contains all the most frequent k-mers in text.

    Examples:
        Compare the frequency map of all the k-mers (given string (i.e. text) and integer (i.e. k)), then return a list of the most frequent k-mers.
        
        >>> text = "ACGTTGCATGTCGCATGATGCATGAGAGCT" 
        >>> k = 4
        >>> kmers_list = frequency_map(text, k)
        >>> kmers_list
            ["CATG", "GCAT"]
    """
    words = []
    freq = frequency_map(text, k)
    m = max(freq.values())
    for key in freq:
        # add each key to words whose corresponding frequency value is equal to m
        if freq[key] == m:
            words.append(key)
    return words


def reverse(pattern):
    """Reverse a string sequence. For example, if we reverse the string 'ACGT', we would get 'TGCA'..

    Args:
        pattern (str): a DNA string (i.e. pattern).

    Returns:
        String, a reversed string of the given pattern.

    Examples:
        Reverse a pattern string.
        
        >>> pattern = 'AAAACCCGGT'
        >>> reversed_pattern = reverse(pattern)
        >>> reversed_pattern
            'TGGCCCAAAA'
    """
    reversed_pattern = pattern[::-1]
    return reversed_pattern


def complement(pattern):
    """Compute the complementary string of pattern, with every nucleotide being replaced by its complement.

    Args:
        pattern (str): a DNA string pattern.

    Returns:
        String, a DNA string pattern in complementary with the given pattern.

    Examples:
        Return the complementary string of a pattern string.
        
        >>> pattern = 'AAAACCCGGT'
        >>> complementary_pattern = complement(pattern)
        >>> complementary_pattern
            'TTTTGGGCCA'
    """
    complement_dict = {"A": "T", "C": "G", "G": "C", "T": "A"}
    complement_pattern = "".join(complement_dict[i] for i in pattern)
    return complement_pattern


def reverse_complement(pattern):
    """Find the reverse complement of a DNA string. This is how DNA is replicated.

    Args:
        pattern (str): a DNA string pattern.

    Returns:
        String, the reverse complement string of the given pattern string.

    Examples:
        Return the reversed complementary string of a pattern string.
        
        >>> pattern = 'AAAACCCGGT'
        >>> reversed_complementary_pattern = reverse_complement(pattern)
        >>> reversed_complementary_pattern
            'ACCGGGTTTT'
    """
    pattern = reverse(pattern)  # reverse all letters in a string
    pattern = complement(pattern)  # complement each letter in a string
    #  or simply use return complement(reverse(pattern))
    return pattern


if __name__ == "__main__":
    # Do something if this file is invoked on its own
    print(canvas())
