"""
bioin.py
Study notes for beginners bioinformatics

Handles the primary functions
"""
# import sys
# lines = sys.stdin.read().splitlines()


def canvas(with_attribution=True):
    """
    Placeholder function to show example docstring (NumPy format)

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
    The number of times that a k-mer pattern appears as a substring of text

    Count the frequency (overlapping occurrences also counts) of a substring of pattern in the given text.

    Args:
        pattern: string, the substring pattern to find in the given text.
        text: string, the string space for looking

    Returns:
        Number of substring pattern that appears in text.

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
    Find the most frequent k-mers in a string

    Computes the frequency map of a given string text and integer k.

    Args:
        text: string, text
        k: integer, k

    Returns:
        All most frequent k-mers in text

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
    """
    Find all the most frequent k-mers in text

    A list contains all the most frequent k-mers in text.

    Args:
        text: string, searching space text string
        k: k-mers to find in text.

    Returns:
        A list of all the most frequent k-mers in the given text

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
    """

    Args:
        pattern: string, a DNA string pattern

    Returns:
        string, a reversed string of the given pattern

    """
    reversed_pattern = pattern[::-1]
    return reversed_pattern


def complement(pattern):
    """
    Compute the complementary string of pattern

    With every nucleotide replaced by its complement.

    Args:
        pattern: string, a DNA string pattern

    Returns:
        string, a DNA string pattern in complementary with the given pattern

    """
    complement_dict = {"A": "T", "C": "G", "G": "C", "T": "A"}
    complement_pattern = "".join(complement_dict[i] for i in pattern)
    return complement_pattern


def reverse_complement(pattern):
    """
    Find the reverse complement of a DNA string.


    Args:
        pattern: string, a DNA string pattern

    Returns:
        The reverse complement of pattern

    """
    pattern = reverse(pattern)  # reverse all letters in a string
    pattern = complement(pattern)  # complement each letter in a string
    #  or simply use return complement(reverse(pattern))
    return pattern


if __name__ == "__main__":
    # Do something if this file is invoked on its own
    print(canvas())
