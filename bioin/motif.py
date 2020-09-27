"""
From motifs to profile matrices and consensus strings:

A computational problem formulation for motif finding would score individual instances of motifs depending on how similar they are to an “ideal” motif (i.e., a transcription factor binding site that binds the best to the transcription factor). However, since the ideal motif is unknown, we attempt to select a k-mer from each string and score these k-mers depending on how similar they are to each other.

To define scoring, consider a list of t DNA strings Dna, where each string has length n, and select a k-mer from each string to form a collection Motifs, which we represent as a t x k motif matrix. In the figure below, which shows the motif matrix for the NF-κB binding sites from the figure below, we indicate the most frequent nucleotide in each column of the motif matrix by upper case letters. If there are multiple most frequent nucleotides in a column, then we arbitrarily select one of them to break the tie. Note that positions 2 and 3 are the most conserved (nucleotide G is completely conserved in these positions), whereas position 10 is the least conserved.

["TCGGGGGTTTTT", "CCGGGTGACTTAC", "ACGGGGATTTTC", "TTGGGGACTTTT", "AAGGGGACTTCC", "TTGGGGACTTCC", "TCGGGGATTCAT", "TCGGGGATTCCT", "TAGGGGAACTAC", "TCGGGTATAACC"]

The NF-κB binding sites form a 10 x 12 motif matrix, with the most frequent nucleotide in each column shown in upper case letters and all other nucleotides shown in lower case letters.

In Python, we will represent a motif matrix as a list of strings Motifs. We can access the i-th string in the motif matrix by calling Motifs[i]; we can access the j-th symbol in this string by calling Motifs[i][j].

By varying the choice of k-mers in each string, we can construct a large number of different motif matrices from a given sample of DNA strings. Our goal is to select k-mers resulting in the most “conserved” motif matrix, meaning the matrix with the most upper case letters (and thus the fewest number of lower case letters). Leaving aside the question of how we select such k-mers, we will first focus on how to score the resulting motif matrices, defining Score(Motifs) as the number of unpopular (lower case) letters in the motif matrix Motifs (see updated figure below). Our goal is to find a collection of k-mers that minimizes this score (for more on motif scoring functions, see DETOUR: Motif Scoring Functions).

["TCGGGGGTTTTT", 

"CCGGGTGACTTAC", 

"ACGGGGATTTTC", 

"TTGGGGACTTTT", 

"AAGGGGACTTCC", 

"TTGGGGACTTCC", 

"TCGGGGATTCAT", 

"TCGGGGATTCCT", 

"TAGGGGAACTAC", 

"TCGGGTATAACC"]

3+4+0+0+1+1+1+5+2+3+6+4 is the Score(Motifs), look at the matrix by column, 3 is the number of nucleotides that is different from the most frequent nucleotides in the first column. 

A: 2  2  0  0  0  0  9  1  1  1  3  0

C: 1  6  0  0  0  0  0  4  1  2  4  6

G: 0  0  10 10 9  9  1  0  0  0  0  0

T: 7  2  0  0  1  1  0  5  8  7  3  4 

The matrix above is the Count(Motifs), we can represent this count matrix in Python as follows:

>>> count = {"A": [2, 2, 0, 0, 0, 0, 9, 1, 1, 1, 3, 0],
         "C": [1, 6, 0, 0, 0, 0, 0, 4, 1, 2, 4, 6],
         "G": [0, 0,10,10, 9, 9, 1, 0, 0, 0, 0, 0],
         "T": [7, 2, 0, 0, 1, 1, 0, 5, 8, 7, 3, 4]
        } 

As shown below, we will further divide all of the elements in the count matrix by t, the number of rows in Motifs. This results in a profile matrix Profile(Motifs) for which element (i,j) is the frequency of the i-th nucleotide in the j-th column of the motif matrix (i.e., the number of occurrences of the i-th nucleotide divided by t, the number of nucleotides in the column). Note that the elements of any column of the profile matrix sum to 1.

A: .2  .2  0  0  0  0  .9  .1  .1  .1  .3  0

C: 1  .6  0  0  0  0  0  .4  .1  .2  .4  .6

G: 0  0  1  1 .9  .9  .1  0  0  0  0  0

T: .7  .2  0  0  .1  .1  0  .5  .8  .7  .3  .4 

Above is the Proifle(motifs) matrix.

Finally, we can form a consensus string, denoted Consensus(Motifs), from the most popular nucleotides in each column of the motif matrix (ties are broken arbitrarily). If we select Motifs correctly from the collection of upstream regions, then Consensus(Motifs) provides a candidate regulatory motif for these regions. For example, as shown below, the consensus string for the NF-κB binding sites is "TCGGGGATTTCC".

The Motif Finding Problem:

Now that we have a good grasp of scoring a collection of k-mers, we are ready to formulate a computational problem for motif finding. 

Motif Finding Problem: Given a collection of strings, find a set of k-mers, one from each string, that minimizes the score of the resulting motif.
 
Input: A collection of strings Dna and an integer k.

Output: A collection Motifs of k-mers, one from each string in Dna, minimizing Score(Motifs) among all possible choices of k-mers.

Brute force search (also known as exhaustive search) is a general problem-solving technique that explores all possible candidate solutions and checks whether each candidate solves the problem. Such algorithms require little effort to design and are guaranteed to produce a correct solution, but they may take an enormous amount of time, and the number of candidates may be too large to check.

A brute force algorithm for the Motif Finding Problem, BruteForceMotifSearch, considers every possible choice of k-mers Motifs from Dna (one k-mer from each string of n nucleotides) and returns the collection Motifs having minimum score.

A k-mer tends to have a higher probability when it is more similar to the consensus string of a profile. For example, for the NF-κB profile matrix (reproduced at bottom) and its consensus string "TCGGGGATTTCC",

Pr("TCGGGGATTTCC", Profile) = .7 · .6 · 1 · 1 · .9 · .9 · .9 · .5 · .8 · .7 · .4 · .6 = 0.0205753 , Pr means probability

which is larger than the value of Pr("ACGGGGATTACC", Profile) = 0.000839808 that we computed on the previous step.

The entropy of the completely conserved third column of the profile matrix in the figure in the first step is 0, which is the minimum possible entropy. On the other hand, a column with equally-likely nucleotides (all probabilities equal to 1/4) has maximum possible entropy 4 * 1/4 * log2 (1/4) = 2. In general, the more conserved the column, the smaller its entropy. Thus, entropy offers an improved method of scoring motif matrices: the entropy of a motif matrix is defined as the sum of the entropies of its columns. In this book, we will continue to use Score(Motifs) for simplicity, but the entropy score is used more often in practice.

What is the probability that the sun will not rise tomorrow?

In 1650, after the Scots proclaimed Charles II as king during the English Civil War, Oliver Cromwell made a famous appeal to the Church of Scotland. Urging them to see the error of their royal alliance, he pleaded,

I beseech you, in the bowels of Christ, think it possible that you may be mistaken.

The Scots rejected the appeal, and Cromwell invaded Scotland in response. His quotation later inspired the statistical maxim called Cromwell’s rule, which states that we should not use probabilities of 0 or 1 unless we are talking about logical statements that can only be true or false. In other words, we should allow a small probability for extremely unlikely events, such as “this book was written by aliens” or “the sun will not rise tomorrow”. We cannot speak to the likelihood of the former event, but in the 18th Century, the French mathematician Pierre-Simon Laplace actually estimated the probability that the sun will not rise tomorrow (1/1826251), given that it has risen every day for the past 5000 years. Although this estimate was ridiculed by his contemporaries, Laplace’s approach to this question now plays an important role in statistics. 

In any observed data set, there is the possibility, especially with low-probability events or small data sets, that an event with nonzero probability does not occur. Its observed frequency is therefore zero; however, setting the empirical probability of the event equal to zero represents an inaccurate oversimplification that may cause problems. By artificially adjusting the probability of rare events, these problems can be mitigated.

Laplace’s Rule of Succession:

Cromwell’s rule is relevant to the calculation of the probability of a string based on a profile matrix. For example, consider the above mentioned Profile

Pr('TCGTGGATTTCC', Profile) = .7 * .6 * 1 * .0 * .9 * .9 * .9 * .5 * .8 * .7 * .4 * .6 = 0, because the column 4 element (differs from the consensus string).

The fourth symbol of "TCGTGGATTTCC" causes Pr("TCGTGGATTTCC", Profile) to be equal to zero. As a result, the entire string is assigned a zero probability, even though "TCGTGGATTTCC" differs from the consensus string at only one position. For that matter, "TCGTGGATTTCC" has the same low probability as "AAATCTTGGAA", which differs from the consensus string at every position.

In order to improve this unfair scoring, bioinformaticians often substitute zeroes with small numbers called pseudocounts. The simplest approach to introducing pseudocounts, called Laplace’s Rule of Succession, is similar to the principle that Laplace used to calculate the probability that the sun will not rise tomorrow. In the case of motifs, pseudocounts often amount to adding 1 (or some other small number) to each element of Count(Motifs). For example, say that we have the following motif, count, and profile matrices:

For example, say that we have the following motif, count, and profile matrices:

Motifs = ['TAAC', 'GTCT', 'ACTA', 'AGGT']

Count(Motifs) = {A:[2, 1, 1, 1], C:[0, 1, 1, 1], G:[1, 1, 1, 0], T:[1, 1, 1, 2]}

Profile(Motifs) = {A:[2/4, 1/4, 1/4, 1/4], C:[0, 1/4, 1/4, 1/4], G:[1/4, 1/4, 1/4, 0], T:[1/4, 1/4, 1/4, 2/4]}

Laplace’s Rule of Succession adds 1 to each element of Count(Motifs), updating the two matrices to the following:

Count(Motifs) = {A:[2+1, 1+1, 1+1, 1+1], C:[0+1, 1+1, 1+1, 1+1], G:[1+1, 1+1, 1+1, 0+1], T:[1+1, 1+1, 1+1, 2+1]}

Profile(Motifs) = {A:[3/8, 2/8, 2/8, 2/8], C:[1/8, 2/8, 2/8, 2/8], G:[2/8, 2/8, 2/8, 1/8], T:[2/8, 2/8, 2/8, 3/8]}

Note that RandomizedMotifSearch may change all t strings in Motifs in a single iteration. This strategy may prove reckless, since some correct motifs (captured in Motifs) may potentially be discarded at the next iteration. GibbsSampler is a more cautious iterative algorithm that discards a single k-mer from the current set of motifs at each iteration and decides to either keep it or replace it with a new one. This algorithm thus moves with more caution in the space of all motifs, as illustrated below.

Like RandomizedMotifSearch, GibbsSampler starts with randomly chosen k-mers in each of t DNA strings, but it makes a random choice at each iteration. It uses a list of t randomly selected k-mers Motifs to come up with another (hopefully better scoring) set of k-mers. In contrast with RandomizedMotifSearch, which defines new motifs as:

Motifs(Profile(Motifs), Dna)

GibbsSampler randomly selects an integer i between 0 and t-1 and then randomly changes a single k-mer Motif[i]. That is, GibbsSampler makes two random choices at each iteration. It uses random.randint(0, t-1) for the first choice (since all t strings in Dna are equally likely), but it does not make sense to use random.randint to choose a k-mer from Motif[i] because some k-mers are more likely than others. Indeed, each k-mer Pattern in Motif[i] may have a different probability Pr(Pattern, Profile) that it was generated by Profile.

We now need to simulate rolling a die so that "ccgG" has probability 4/80, "cgGC" has probability 8/80, and so on. We will do so by generating a random number p between 0 and 1. If p is between 0 and 4/80, then it corresponds to "ccgG". If p is between 4/80 and 4/80 + 8/80, then it corresponds to "cgGC", etc.

How can we generate a random number between 0 and 1? In addition to the randint function, Python’s random module also includes a function called uniform(a, b) that generates a random floating point number (i.e., a decimal) between a and b. We can therefore generate our desired random number p by calling random.uniform(0, 1).


"""
import random


def count_motif(motifs):
    """Count the number of nucleotides (4 types: ACGT) column wise from a motifs matrix.

    Args:
        motifs: 2D matrix, matrix of motifs in genome.

    Returns:
        dictionary, the count of each nucleotides in each column of the motifs matrix.
<<<<<<< HEAD

    Examples:
        Takes a list of strings motifs as input and returns the count matrix of motifs (as a dictionary of lists.)

=======

    Examples:
        Takes a list of strings motifs as input and returns the count matrix of motifs (as a dictionary of lists.)

>>>>>>> bfc59e1c1145bfc70ed604cff789d85f87235a81
        >>> motifs = ['AACGTA', 'CCCGTT', 'CACCTT', 'GGATTA', 'TTCCGG']
        >>> counts_dict = count_motif(motifs)
        >>> counts_dict
            {'A': [1, 2, 1, 0, 0, 2], 'C': [2, 1, 4, 2, 0, 0], 'G': [1, 1, 0, 2, 1, 1], 'T': [1, 1, 0, 1, 4, 2]}
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
    """The percentage of count number of nucleotides column wise from a motifs matrix.

    Args:
        motifs: 2D matrix, matrix of motifs in genome.

    Returns:
        dictionary, the percentile count of each nucleotides in each column of the motifs matrix.

    Examples:
        Takes a list of strings motifs as input and then generate the count_motif(motifs), then divide each element of the count matrix by the number of rows in the count matrix, to obtain the profile_motif matrix (as a dictionary of lists.)

        >>> motifs = ['AACGTA', 'CCCGTT', 'CACCTT', 'GGATTA', 'TTCCGG']
        >>> profile_dict = profile_motif(motifs)
        >>> profile_dict
            {'A': [0.2, 0.4, 0.2, 0.0, 0.0, 0.4], 'C': [0.4, 0.2, 0.8, 0.4, 0.0, 0.0], 'G': [0.2, 0.2, 0.0, 0.4, 0.2, 0.2], 'T': [0.2, 0.2, 0.0, 0.2, 0.8, 0.4]}
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
    """Compute the consensus string of a set of k-mers motifs.

    Args:
        motifs: list of strings, a set of k-mers motifs.

    Returns:
        string, a consensus string of motifs.

    Examples:
        Takes a list of strings motifs as input and then return the string that contain the most frequent element of each column of the matrix.

        >>> motifs = ['AACGTA', 'CCCGTT', 'CACCTT', 'GGATTA', 'TTCCGG']
        >>> consensus = consensus_motif(motifs)
        >>> consensus
            'CACCTA'
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
    """Compute the total different occurrences of nucleotides that consensus have with motifs.

    Args:
        motifs: list of strings, a set of k-mers motifs.

    Returns:
        integer, number of different nucleotides that consensus have with motifs.

    Examples:
        Takes a list of strings motifs as input, compute score_motifs by first ocnstructing consensus_motifs(motifs) and then summing the number of symbols in the j-th column of motifs that do not match the symbol in the position j of the consensus string.

        >>> motifs = ['AACGTA', 'CCCGTT', 'CACCTT', 'GGATTA', 'TTCCGG']
        >>> score = score_motif(motifs)
        >>> score
            14
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
    """Calculate the probability of a string (i.e. text) given a motif profile matrix.

    Args:
        text: string, a text string probability given the profile matrix.
        motif_profile: matrix, represent the probability of occurence of each nucleotide for each column.

    Returns:
        float, the probability of text calculated from profile

    Examples:
        Calculate the probability for the given text string based on the profile matrix.

        >>> text = "ACGGGGATTACC"
        >>> motif_profile = {'A': [0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0], 'C': [0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6], 'G': [0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0], 'T': [0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4]}
        >>> pr = (text, motif_profile) # probability
            0.000839808
    """
    probability = 1
    for i in range(len(text)):
        probability *= motif_profile[text[i]][i]
    return probability


def profile_most_probable_kmer(text, k, motif_profile):
    """Find a profile-most probable k-mer in a string.

    Args:
        text: string, genome string.
        k: int, length of the k-mer.
        motif_profile: a 4 * k matrix profile (i.e. profile matrix).

    Returns:
        a profile-most probable k-mer in text.
<<<<<<< HEAD

    Examples:
        Given a profile matrix (i.e. motif_profile or Profile), we can compute the probability of every k-mer in a string Text and find a profile-most probable k-mer in Text, i.e., a k-mer that was most likely to have been generated by Profile among all k-mers in Text. For the NF-κB profile matrix, "ACGGGGATTACC" is the Profile-most probable 12-mer in "ggtACGGGGATTACCt". Indeed, every other 12-mer in this string has probability 0. In general, if there are multiple Profile-most probable k-mers in Text, then we select the first such k-mer occurring in Text.

=======

    Examples:
        Given a profile matrix (i.e. motif_profile or Profile), we can compute the probability of every k-mer in a string Text and find a profile-most probable k-mer in Text, i.e., a k-mer that was most likely to have been generated by Profile among all k-mers in Text. For the NF-κB profile matrix, "ACGGGGATTACC" is the Profile-most probable 12-mer in "ggtACGGGGATTACCt". Indeed, every other 12-mer in this string has probability 0. In general, if there are multiple Profile-most probable k-mers in Text, then we select the first such k-mer occurring in Text.

>>>>>>> bfc59e1c1145bfc70ed604cff789d85f87235a81
        >>> text = 'ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT'
        >>> k = 5
        >>> motif_profile = {'A': [0.2, 0.2, 0.3, 0.2, 0.3], 'C': [0.4, 0.3, 0.1, 0.5, 0.1], 'G': [0.3, 0.3, 0.5, 0.2, 0.4], 'T': [0.1, 0.2, 0.1, 0.1, 0.2]}
        >>> kmer = profile_most_probable_kmer(text, k, motif_profile)
        >>> kmer
            "CCGAG"
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
    """Calculate t k-mers from dna that have the best score (i.e. the most frequently occur t k-mers in the given dna)

    Args:
        dna: matrix, has t rows. 
        k: integer, k-mer.
        t: integer, t is the number of k-mers in dna to return, also equal to the row number of dna 2D matrix.

    Returns:
        list (or can take it as 2D matrix, a sub matrix of dna), t k-mers.
<<<<<<< HEAD

    Examples:
        GreedyMotifSearch, starts by setting best_motifs equal to the first k-mer from each string in Dna (each row assign a k-mer), then ranges over all possible k-mers in dna[0], the algorithm then builds a profile matrix Profile fro this lone k-mer, and sets Motifs[1] equal to the profile_most_probable k-mer in dna[1]. Then iterates by updating Profile as the profile matrix formed from Motifs[0] and Motifs[1], and sets Motifs[2] equal to the profile_most_probable k-mer in dna[2]. After finding k-mers Motifs in the first i strings of Dna, GreedyMotifSearch constructs Profile(Motifs) and sets Motifs[i] equal to the profile_most_probable k-mer from dna[i] based on this profile matrix.

=======

    Examples:
        GreedyMotifSearch, starts by setting best_motifs equal to the first k-mer from each string in Dna (each row assign a k-mer), then ranges over all possible k-mers in dna[0], the algorithm then builds a profile matrix Profile fro this lone k-mer, and sets Motifs[1] equal to the profile_most_probable k-mer in dna[1]. Then iterates by updating Profile as the profile matrix formed from Motifs[0] and Motifs[1], and sets Motifs[2] equal to the profile_most_probable k-mer in dna[2]. After finding k-mers Motifs in the first i strings of Dna, GreedyMotifSearch constructs Profile(Motifs) and sets Motifs[i] equal to the profile_most_probable k-mer from dna[i] based on this profile matrix.

>>>>>>> bfc59e1c1145bfc70ed604cff789d85f87235a81
        >>> dna = ['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG']
        >>> k = 3
        >>> t = 5
        >>> t_kmers = greedy_motif_search(dna, k, t)
        >>> t_kmers
            ['CAG', 'CAG', 'CAA', 'CAA', 'CAA']
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
    """Count the number of nucleotides (4 types: ACGT) column wise from a motifs matrix, then add 1 to each position, i.e. pseudocount.

    Args:
        motifs: 2D matrix, matrix of motifs in genome.

    Returns:
        dictionary, the count matrix of Motifs with pseudocounts as a dictionary of lists.

    Examples:
        Takes a list of strings motifs as input and returns the pseudocount matrix of motifs (as a dictionary of lists.)

        >>> motifs = ['AACGTA', 'CCCGTT', 'CACCTT', 'GGATTA', 'TTCCGG']
        >>> pseudocounts_dict = count_with_pseudocount(motifs)
        >>> pseudocounts_dict
            {'A': [2, 3, 2, 1, 1, 3], 'C': [3, 2, 5, 3, 1, 1], 'G': [2, 2, 1, 3, 2, 2], 'T': [2, 2, 1, 2, 5, 3]}
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
    """The percentage of pseudocount number of nucleotides column wise from a motifs matrix.

    Args:
        motifs: 2D matrix, matrix of motifs in genome.

    Returns:
        dictionary, the percentile pseudocount of each nucleotides in each column of the motifs matrix.

    Examples:
        Takes a list of strings motifs as input and then generate the count_with_pseudocount(motifs), then divide each element of the pseudocount matrix by the number of rows plus four in the pseudocount matrix, to obtain the profile_with_pseudocount (as a dictionary of lists.)

        >>> motifs = ['AACGTA', 'CCCGTT', 'CACCTT', 'GGATTA', 'TTCCGG']
        >>> profile_pseudo_dict = profile_with_pseudocount(motifs)
        >>> profile_pseudo_dict
            {'A': [0.2222222222222222, 0.3333333333333333, 0.2222222222222222, 0.1111111111111111, 0.1111111111111111, 0.3333333333333333], 'C': [0.3333333333333333, 0.2222222222222222, 0.5555555555555556, 0.3333333333333333, 0.1111111111111111, 0.1111111111111111], 'G': [0.2222222222222222, 0.2222222222222222, 0.1111111111111111, 0.3333333333333333, 0.2222222222222222, 0.2222222222222222], 'T': [0.2222222222222222, 0.2222222222222222, 0.1111111111111111, 0.2222222222222222, 0.5555555555555556, 0.3333333333333333]}
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
    """Compute the consensus string of a set of k-mers motifs. Using pseudocount.

    Args:
        motifs: list of strings, a set of k-mers motifs.

    Returns:
        string, a consensus string of motifs.

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
    """Compute the total different occurrences of nucleotides that consensus_pseudo have with motifs

    Args:
        motifs: list of strings, a set of k-mers motifs

    Returns:
        integer, number of different nucleotides that consensus_pseudo have with motifs

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
    """Calculate the probability of a string (i.e. text) given a motif profile matrix. With pseudocount.

    Args:
        text: string, a text string probability given the profile matrix.
        motif_profile: matrix, represent the probability of occurence of each nucleotide for each column.

    Returns:
        float, the probability of text calculated from profile.

    """
    probability = 1
    for i in range(len(text)):
        probability *= motif_profile[text[i]][i]
    return probability


def profile_most_probable_kmer_with_pseudocount(text, k, motif_profile):
    """Find a profile-most probable k-mer in a string. With pseudocount.

    Args:
        text: string, genome string.
        k: int, length of the k-mer.
        motif_profile: a 4 * k matrix profile (i.e. profile matrix).

    Returns:
        a profile-most probable k-mer in text.

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
    """Calculate t k-mers from dna that have the best score (i.e. the most frequently occur t k-mers in the given dna). With pseudocount.

    Args:
        dna: matrix, has t rows. 
        k: integer, k-mer.
        t: integer, t is the number of k-mers in dna to return, also equal to the row number of dna 2D matrix.

    Returns:
        list (or can take it as 2D matrix, a sub matrix of dna), t k-mers.

    Examples:
        GreedyMotifSearch, starts by setting best_motifs equal to the first k-mer from each string in Dna (each row assign a k-mer), then ranges over all possible k-mers in dna[0], the algorithm then builds a profile matrix Profile fro this lone k-mer, and sets Motifs[1] equal to the profile_most_probable k-mer in dna[1]. Then iterates by updating Profile as the profile matrix formed from Motifs[0] and Motifs[1], and sets Motifs[2] equal to the profile_most_probable k-mer in dna[2]. After finding k-mers Motifs in the first i strings of Dna, GreedyMotifSearch constructs Profile(Motifs) and sets Motifs[i] equal to the profile_most_probable k-mer from dna[i] based on this profile matrix.

        >>> dna = ["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"]
        >>> k = 3 ["TTC", "ATC", "TTC", "ATC", "TTC"]
        >>> t_kmers_pseudo= greedy_motif_search_with_pseudocount(dna, k, t)
        >>> t_kmers_pseudo
            ["TTC", "ATC", "TTC", "ATC", "TTC"]
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
    """Find the profile-most probable k-mers in each string from dna. Without pseudocount.

    Args:
        motif_profile: a profile matrix, Profile.
        dna: a list of strings DNA.

    Returns:
        a list of the Profile-most probable k-mers in each string from DNA.

    Examples:

        >>> motif_profile = {'A': [0.8, 0.0, 0.0, 0.2], 'C': [0.0, 0.6, 0.2, 0.0], 'G': [0.2, 0.2, 0.8, 0.0], 'T': [0.2, 0.2, 0.0, 0.8]}
        >>> dna = ['TTACCTTAAC', 'GATGTCTGTC', 'ACGGCGTTAG', 'CCCTAACGAG', 'CGTCAGAGGT']
        >>> kmers = profile_probable_motifs(motif_profile, dna)
        >>> kmers
            ['ACCT', 'ATGT', 'GCGT', 'ACGA', 'AGGT']
    """
    k = len(motif_profile["A"]) # get the column number from the profile matrix
    motifs = []
    for i in range(len(dna)):
        motifs.append(profile_most_probable_kmer(dna[i], k, motif_profile))
    return motifs


def random_motifs(dna, k, t):
    """Return a list of random k-mers from each of t different strings dna.

    Args:
        dna: matrix, has t rows. 
        k: integer, k-mer.
        t: integer, t is the number of k-mers in dna to return, also equal to the row number of dna 2D matrix.

    Returns:
        a list of t k-mers strings, each k-mer from each of t different strings dna.
<<<<<<< HEAD

    Examples:
        RandomMotifs(Dna, k, t) that uses random.randint to choose a random k-mer from each of t different strings Dna, and returns a list of t strings.

=======

    Examples:
        RandomMotifs(Dna, k, t) that uses random.randint to choose a random k-mer from each of t different strings Dna, and returns a list of t strings.

>>>>>>> bfc59e1c1145bfc70ed604cff789d85f87235a81
        >>> dna = ['TTACCTTAAC', 'GATGTCTGTC', 'ACGGCGTTAG', 'CCCTAACGAG', 'CGTCAGAGGT']
        >>> k = 3
        >>> t = 5
        >>> random_t_kmers = random_motifs(dna, k, t)
        >>> random_t_kmers
            ['TTA', 'GTC', 'ACG', 'ACG', 'GAG']
    """
    motifs = []
    for i in range(t):
        m = random.randint(0, len(dna[0])-k)
        motifs.append(dna[i][m:m+k])
    return motifs


def randomized_motif_search(dna, k, t):
    """Return a list of best k-mers from each of t different strings dna. Compare score_pseudo of the k-mer.

    Args:
        dna: matrix, has t rows. 
        k: integer, k-mer.
        t: integer, t is the number of k-mers in dna to return, also equal to the row number of dna 2D matrix.

    Returns:
        a list of t k-mers strings, each k-mer from each of t different strings dna.

    Examples:
        RandomMotifs(Dna, k, t) that uses random.randint to choose a random k-mer from each of t different strings Dna, and returns a list of t strings.

        >>> dna = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA', 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG', 'TAGTACCGAGACCGAAAGAAGTATACAGGCGT', 'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC', 'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
        >>> k = 8
        >>> t = 5
        >>> random_best_t_kmers = randomized_motif_search(dna, k, t)
        >>> random_best_t_kmers
            ['CGGGGGTG', 'TGTAAGTG', 'TACAGGCG', 'TTCAGGTG', 'TCCACGTG']    
    """
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
    """Rescale a collection of probabilityes so that they sum up to 1.

    Args:
        probabilities: a dictionary of Probabilities, where keys are k-mers and values are the probabilities of these k-mers (which do no necessarily sum up to 1).

    Returns:
        a normalized dictionary where the probability of each k-mer was divided by the sun of all k-mers' probabilities.
<<<<<<< HEAD

    Examples:
        To rescale a collection of probabilities (the sides of the die) so that these probabilities sum to 1, we will write a function called Normalize(Probabilities). This function takes a dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these k-mers (which do not necessarily sum to 1). The function should divide each value in Probabilities by the sum of all values in  Probabilities, then return the resulting dictionary.

=======

    Examples:
        To rescale a collection of probabilities (the sides of the die) so that these probabilities sum to 1, we will write a function called Normalize(Probabilities). This function takes a dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these k-mers (which do not necessarily sum to 1). The function should divide each value in Probabilities by the sum of all values in  Probabilities, then return the resulting dictionary.

>>>>>>> bfc59e1c1145bfc70ed604cff789d85f87235a81
        >>> probabilities = {'A': 0.1, 'C': 0.1, 'G': 0.1, 'T': 0.1}
        >>> normalized_pro = normalize_probability(probabilities)
        >>> normalized_pro
            {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
    """
    n_probability = {}
    for k, v in probabilities.items():
        n_probability[k] = probabilities[k] / sum(probabilities.values())
    return n_probability


def weighted_die(probabilities):
    """Return a randomly chosen k-mer key with respect to the values in Probabilities.

    Args:
        probabilities: dictionary, a dictionary probabilities whose keys are k-mers and whose values are the probabilities of these kmers.

    Returns:
        string, a randomly chosen k-mer with respect to the values in Probabilities.
<<<<<<< HEAD

    Examples:
        Generalize this idea by writing a function WeightedDie(Probabilities). This function takes a dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these k-mers. The function should return a randomly chosen k-mer key with respect to the values in Probabilities.

=======

    Examples:
        Generalize this idea by writing a function WeightedDie(Probabilities). This function takes a dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these k-mers. The function should return a randomly chosen k-mer key with respect to the values in Probabilities.

>>>>>>> bfc59e1c1145bfc70ed604cff789d85f87235a81
        >>> probabilities = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
        >>> weighted_die_kmer = weighted_die(probabilities)
        >>> weighted_die_kmer
            'G'            
    """
    n = random.uniform(0, 1)
    for p in probabilities:
        n -= probabilities[p]
        if n <= 0:
            return p


#tests in which of the intervals defined by list ar the number r lies
def testinterval(ar,r):
    ar.sort()
    if r<= ar[0]:
        return ar[0]

    for i in range(1,len(ar)-1):
        if ar[i-1]<r<=ar[i]:
            return ar[i]

    if ar[len(ar)-2]< r:
        return ar[len(ar)-1]


def WeightedDie(Probabilities):
    sumprob = {}
    s = 0
    for p in Probabilities:
        s += Probabilities[p]
        sumprob[p] = s

    revprob = {}
    for q in sumprob:
        revprob[sumprob[q]] = q
        w = list(sumprob.values())
        r = random.uniform(0,1)
        kmer = revprob[testinterval(w,r)]

    return kmer


def profile_generated_string(text, motif_profile, k):
    """Return a randomly generated k-mer from Text whose probabilities are generated from profile.

    Args:
        text: string, string dna.
        motif_profile: matrix, a profile matrix
        k: int, length of the k-mer

    Returns:
        string, a randomly chosen k-mer key with respect to the values in Probabilities.
<<<<<<< HEAD

    Examples:
        Assemble this code into a function ProfileGeneratedString(Text, profile, k) that takes a string Text, a profile matrix profile , and an integer k as input.  It should then return a randomly generated k-mer from Text whose probabilities are generated from profile, as described above. 

=======

    Examples:
        Assemble this code into a function ProfileGeneratedString(Text, profile, k) that takes a string Text, a profile matrix profile , and an integer k as input.  It should then return a randomly generated k-mer from Text whose probabilities are generated from profile, as described above. 

>>>>>>> bfc59e1c1145bfc70ed604cff789d85f87235a81
        >>> text = 'AAACCCAAACCC'
        >>> motif_profile = {'A': [0.5, 0.1], 'C': [0.3, 0.2], 'G': [0.2, 0.4], 'T': [0.0, 0.3]}
        >>> k = 2
        >>> kmers = profile_generated_string(text, motif_profile, k)
        >>> kmers
            'AC' 
    """
    probabilities = {}
    for i in range(0, len(text)-k+1):
        probabilities[text[i:i+k]] = probability_with_pseudocount(text[i:i+k], motif_profile)

    probabilities = normalize_probability(probabilities)
    return weighted_die(probabilities)


<<<<<<< HEAD
def GibbsSampler(Dna, k, t, N):
    """Using GibbsSampler method to return the best motifs of t k-mers in each of the strings in Dna.
    
    Args:
        Dna: matrix, a collection of strings dna, has t rows. 
        k: int, k-mer.
        t: int, t is the number of k-mers in dna to return, also equal to the row number of dna 2D matrix.
        N:  int, the number of iterations that we plan to run the program.

    Returns：
        matrix, the best motifs, t k-mers of each row of the strings in Dna.

    Examples:
        Although GibbsSampler performs well in many cases, it may converge to a suboptimal solution, particularly for difficult search problems with elusive motifs. A local optimum is a solution that is optimal within a small neighboring set of solutions, which is in contrast to a global optimum, or the optimal solution among all possible solutions. Since GibbsSampler explores just a small subset of solutions, it may “get stuck” in a local optimum. For this reason, similarly to RandomizedMotifSearch, it should be run many times with the hope that one of these runs will produce the best-scoring motifs. Yet convergence to a local optimum is just one of many issues we must consider in motif finding.

        >>> Dna = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA', 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG', 'TAGTACCGAGACCGAAAGAAGTATACAGGCGT', 'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC', 'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
        >>> k = 8
        >>> t = 5
        >>> N = 100
        >>> best_motif_gibs = GibbsSampler(Dna, k, t, N)
        >>> best_motif_gibs
            ['AACGGCCA', 'AAGTGCCA', 'TAGTACCG', 'AAGTTTCA', 'ACGTGCAA']    
    """
    BestMotifs = [] # output variable
    Motifs =  random_motifs(dna, k, t)
    BestMotifs = Motifs

    for j in range(1,N):
        i = random.randint(0,t-1)
        ReducedMotifs = []
        for j in range(0,t):
            if j != i:
                ReducedMotifs.append(Motifs[j])

        Profile = profile_with_pseudocount(ReducedMotifs)
        Motif_i = profile_generated_string(dna[i], Profile, k)
        Motifs[i] = Motif_i
        if score_with_pseudocount(Motifs) < score_with_pseudocount(BestMotifs):
            BestMotifs = Motifs

    return BestMotifs


=======
>>>>>>> bfc59e1c1145bfc70ed604cff789d85f87235a81
# a = profile_probable_motifs(profile_with_pseudocount(["TGA", "GTT", "GAA", "TGT"]), ["TGACGTTC", "TAAGAGTT", "GGACGAAA", "CTGTTCGC"])
# print(a)