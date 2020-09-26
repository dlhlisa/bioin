import random


def count_motif(motifs):
    """
    Count the number of nucleotides column wise from a motifs matrix

    Args:
        motifs: 2D matrix, matrix of motifs in genome

    Returns:
        dictionary, the count of each nucleotides in each column of the motifs matrix

    """
    count = {}  # initializing the count dictionary
    k = len(motifs[0])  # the length of the first row in motifs is actually the column number of motifs matrix
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
    t = len(motifs)  # t is equal to the row number of motifs matrix
    for i in range(t):
        for j in range(k):
            symbol = motifs[i][j]
            count[symbol][j] += 1
    return count


def profile_motif(motifs):
    """
    The percentage of count number of nucleotides column wise from a motifs matrix

    Args:
        motifs: 2D matrix, matrix of motifs in genome

    Returns:
        dictionary, the percentile count of each nucleotides in each column of the motifs matrix

    """
    t = len(motifs)
    k = len(motifs[0])
    profile = {}
    count = count_motif(motifs)
    for symbol in "ACGT":
        profile[symbol] = []
        for j in range(k):
            profile[symbol].append(count[symbol][j] / t)
    return profile


def consensus_motif(motifs):
    """
    compute the consensus string of a set of k-mers motifs

    Args:
        motifs: list of strings, a set of k-mers motifs

    Returns:
        string, a consensus string of motifs

    """
    consensus = ""
    count = count_motif(motifs)
    k = len(motifs[0])
    for j in range(k):
        m = 0
        frequent_symbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequent_symbol = symbol
        consensus += frequent_symbol
    return consensus


def score_motif(motifs):
    """
    Compute the total different occurrences of nucleotides that consensus have with motifs

    Args:
        motifs: list of strings, a set of k-mers motifs

    Returns:
        integer, number of different nucleotides that consensus have with motifs

    """
    consensus = consensus_motif(motifs)
    count = count_motif(motifs)
    t = len(motifs)
    k = len(motifs[0])
    max_com = 0
    for j in range(k):
        max_com += count[consensus[j]][j]
    score = t * k - max_com
    # another way that using a former function hamming_distance in replication.py
    # score = 0
    #     for i in range(len(motifs)):
    #         score += hamming_distance(motifs[i], consensus_motif(motifs))
    #     return score
    return score


def probability_profile(text, motif_profile):
    """

    Args:
        text:
        motif_profile:

    Returns:
        float, the probability of text calculated from profile

    """
    probability = 1
    for i in range(len(text)):
        probability *= motif_profile[text[i]][i]
    return probability


def profile_most_probable_kmer(text, k, motif_profile):
    """

    Args:
        text:
        k:
        motif_profile:

    Returns:

    """
    probability_kmer = {}
    most_probable_kmer = []
    for i in range(len(text)-k+1):
        probability_kmer[text[i:i+k]] = probability_profile(text[i:i+k], motif_profile)
        # assume that the given motif_profile is the current k-mer, calculate the probability
        # build the dictionary probability_kmer with current kmer string as key, and the probability as value
        # the profile_most_probable_kmer should be the one with the highest probability value in the dictionary
    prob_max = max(probability_kmer.values())
    for key, value in probability_kmer.items():
        if probability_kmer[key] == prob_max:
            most_probable_kmer.append(key)
    # return most_probable_kmer will raise error in the online course
    # (sys.stdout.write(most_prob_kmer + '\n') TypeError: can only concatenate list (not "str") to list)
    # add [0] to only return the first result in the list would solve the problem
    # But in real scenario, I think we should just return most_probable_kmer as a list, not the first element
    return most_probable_kmer[0]


def greedy_motif_search(dna, k, t):
    """
    calculate t k-mers from dna that have the best score (i.e. the most frequently occur in the given dna)

    Args:
        dna: list, a list of dna strings
        k: integer, k-mer
        t: integer, t is the number of k-mers in dna to return

    Returns:
        list, t k-mers

    """
    best_motifs = []
    for i in range(0, t):
        best_motifs.append(dna[i][0:k])
        # best_motifs initialized with all the first k-mer (total number is t) from each dna string

    for m in range(len(dna[0])-k+1):  # loop through the first dna string, motifs[0] = the first k-mer in dna[0]
        motifs = [dna[0][m:m+k]]  # motifs is the m-th k-mer got through looping the first dna string
        for j in range(1, t):
            profile = profile_motif(motifs[0:j])
            # motifs[0:1] is "GGC", the first k-mer, using this to calculate profile will got
            # profile={"A":[0, 0, 0], "C":[0, 0, 1], "G":[1, 1, 0], "T":[0, 0, 0]} in the first iteration
            # when start the second loop, j=2, the motifs already got two elements
            motifs.append(profile_most_probable_kmer(dna[j], k, profile))
            # return the profile_most_probable_kmer from dna[1] and append it to motifs
            # now the list motifs has two element, i.e. the first k-mer from dna[0], and the one returned from dna[1]
            # start the next loop, j = 2, got another k-mer from dna[2], append it to motifs
            # until the inner-loop (j from 1 to t) finished,
            # the motifs list should already have t element (of k-mer) each from different lines of the dna string
        if score_motif(motifs) < score_motif(best_motifs):
            # compare the score between the initialized best_motifs (t first k-mers from each dna lines)
            # and the motifs (t k-mers from the inner-loop)
            # assign the one with the lowest score as best_motifs
            best_motifs = motifs
            # get out of the inner loop, assign the initial motifs as the second k-mer from dna[0], the first dna line
            # repeat the inner loop again and re-assign the best_motifs if true
            # after looping through all the k-mers in dna[0], the first line in dna, return tbe best_motifs
    return best_motifs


"""
# python allows to redefine a function's name just as we can change the value of a variable
def CountWithPseudocounts(Motifs):

def ProfileWithPseudocounts(Motifs):

# Include these lines to have WithPseudocounts versions be used:
Count = CountWithPseudocounts
Profile = ProfileWithPseudocounts

def Consensus(Motifs):

def Score(Motifs):

def Pr(Text, ProfileMatrix):

def ProfileMostProbablePattern(Text, k, Profile):

def GreedyMotifSearch(Dna, k, t):

def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    # Use the WithPseudocounts versions of Count() and Profile()
    # See extra lines added ABOVE
    return GreedyMotifSearch(Dna, k, t)
"""


def count_with_pseudocount(motifs):
    """

    Args:
        motifs:

    Returns:

    """
    t = len(motifs)  # t is equal to the row number of motifs matrix
    k = len(motifs[0])  # the length of the first row in motifs is actually the column number of motifs matrix
    count = {}  # initializing the count dictionary
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(1)
    for i in range(t):
        for j in range(k):
            symbol = motifs[i][j]
            count[symbol][j] += 1
    return count


def profile_with_pseudocount(motifs):
    """

    Args:
        motifs:

    Returns:

    """
    t = len(motifs)
    k = len(motifs[0])
    profile = {}
    count = count_with_pseudocount(motifs)
    for symbol in "ACGT":
        profile[symbol] = []
        for j in range(k):
            profile[symbol].append(count[symbol][j] / (t+4))
    return profile


def consensus_with_pseudocount(motifs):
    """
    compute the consensus string of a set of k-mers motifs

    Args:
        motifs: list of strings, a set of k-mers motifs

    Returns:
        string, a consensus string of motifs

    """
    consensus = ""
    count = count_with_pseudocount(motifs)
    k = len(motifs[0])
    for j in range(k):
        m = 0
        frequent_symbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequent_symbol = symbol
        consensus += frequent_symbol
    return consensus


def score_with_pseudocount(motifs):
    """
    Compute the total different occurrences of nucleotides that consensus have with motifs

    Args:
        motifs: list of strings, a set of k-mers motifs

    Returns:
        integer, number of different nucleotides that consensus have with motifs

    """
    consensus = consensus_with_pseudocount(motifs)
    count = count_with_pseudocount(motifs)
    t = len(motifs)
    k = len(motifs[0])
    max_com = 0
    for j in range(k):
        max_com += count[consensus[j]][j]
    score = t * k - max_com
    # another way that using a former function hamming_distance in replication.py
    # score = 0
    #     for i in range(len(motifs)):
    #         score += hamming_distance(motifs[i], consensus_motif(motifs))
    #     return score
    return score


def probability_with_pseudocount(text, motif_profile):
    """

    Args:
        text:
        motif_profile:

    Returns:
        float, the probability of text calculated from profile

    """
    probability = 1
    for i in range(len(text)):
        probability *= motif_profile[text[i]][i]
    return probability


def profile_most_probable_kmer_with_pseudocount(text, k, motif_profile):
    """

    Args:
        text:
        k:
        motif_profile:

    Returns:

    """
    probability_kmer = {}
    most_probable_kmer = []
    for i in range(len(text)-k+1):
        probability_kmer[text[i:i+k]] = probability_with_pseudocount(text[i:i+k], motif_profile)
        # assume that the given motif_profile is the current k-mer, calculate the probability
        # build the dictionary probability_kmer with current kmer string as key, and the probability as value
        # the profile_most_probable_kmer should be the one with the highest probability value in the dictionary
    prob_max = max(probability_kmer.values())
    for key, value in probability_kmer.items():
        if probability_kmer[key] == prob_max:
            most_probable_kmer.append(key)
    # return most_probable_kmer will raise error in the online course
    # (sys.stdout.write(most_prob_kmer + '\n') TypeError: can only concatenate list (not "str") to list)
    # add [0] to only return the first result in the list would solve the problem
    # But in real scenario, I think we should just return most_probable_kmer as a list, not the first element
    return most_probable_kmer[0]


def greedy_motif_search_with_pseudocount(dna, k, t):
    """

    Args:
        dna:
        k:
        t:

    Returns:

    """
    best_motifs = []
    for i in range(0, t):
        best_motifs.append(dna[i][0:k])
        # best_motifs initialized with all the first k-mer (total number is t) from each dna string

    for m in range(len(dna[0]) - k + 1):  # loop through the first dna string, motifs[0] = the first k-mer in dna[0]
        motifs = [dna[0][m:m + k]]  # motifs is the m-th k-mer got through looping the first dna string
        for j in range(1, t):
            profile = profile_with_pseudocount(motifs[0:j])
            # motifs[0:1] is "GGC", the first k-mer, using this to calculate profile will got
            # profile={"A":[0, 0, 0], "C":[0, 0, 1], "G":[1, 1, 0], "T":[0, 0, 0]} in the first iteration
            # when start the second loop, j=2, the motifs already got two elements
            motifs.append(profile_most_probable_kmer_with_pseudocount(dna[j], k, profile))
            # return the profile_most_probable_kmer from dna[1] and append it to motifs
            # now the list motifs has two element, i.e. the first k-mer from dna[0], and the one returned from dna[1]
            # start the next loop, j = 2, got another k-mer from dna[2], append it to motifs
            # until the inner-loop (j from 1 to t) finished,
            # the motifs list should already have t element (of k-mer) each from different lines of the dna string
        if score_with_pseudocount(motifs) < score_with_pseudocount(best_motifs):
            # compare the score between the initialized best_motifs (t first k-mers from each dna lines)
            # and the motifs (t k-mers from the inner-loop)
            # assign the one with the lowest score as best_motifs
            best_motifs = motifs
            # get out of the inner loop, assign the initial motifs as the second k-mer from dna[0], the first dna line
            # repeat the inner loop again and re-assign the best_motifs if true
            # after looping through all the k-mers in dna[0], the first line in dna, return tbe best_motifs
    return best_motifs


def profile_probable_motifs(motif_profile, dna):
    """

    Args:
        motif_profile:
        dna:

    Returns:

    """
    k = len(motif_profile["A"])
    motifs = []
    for i in range(len(dna)):
        motifs.append(profile_most_probable_kmer(dna[i], k, motif_profile))
    return motifs


def random_motifs(dna, k, t):
    """

    Args:
        dna:
        k:
        t:

    Returns:

    """
    motifs = []
    for i in range(t):
        m = random.randint(0, len(dna[0])-k)
        motifs.append(dna[i][m:m+k])
    return motifs


def randomized_motif_search(dna, k, t):
    m = random_motifs(dna, k, t)
    best_motifs = m
    while True:
        profile = profile_with_pseudocount(m)
        m = profile_probable_motifs(profile, dna)
        if score_with_pseudocount(m) < score_with_pseudocount(best_motifs):
            best_motifs = m
        else:
            return best_motifs


def normalize_probability(probabilities):
    """

    Args:
        probabilities:

    Returns:

    """
    n_probability = {}
    for k, v in probabilities.items():
        n_probability[k] = probabilities[k] / sum(probabilities.values())
    return n_probability


def weighted_die(probabilities):
    """

    Args:
        probabilities:

    Returns:

    """
    n = random.uniform(0, 1)
    for p in probabilities:
        n -= probabilities[p]
        if n <= 0:
            return p


def profile_generated_string(text, motif_profile, k):
    """

    Args:
        text:
        motif_profile:
        k:

    Returns:

    """
    probabilities = {}
    for i in range(0, len(text)-k+1):
        probabilities[text[i:i+k]] = probability_with_pseudocount(text[i:i+k], motif_profile)

    probabilities = normalize_probability(probabilities)
    return weighted_die(probabilities)



a = profile_probable_motifs(profile_with_pseudocount(["TGA", "GTT", "GAA", "TGT"]), ["TGACGTTC", "TAAGAGTT", "GGACGAAA", "CTGTTCGC"])
# print(a)