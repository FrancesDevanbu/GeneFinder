# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Frances Devanbu

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    if nucleotide == 'A' :
        return 'T'
    elif nucleotide == 'T' :
        return 'A'
    elif nucleotide == 'G' :
        return 'C'
    elif nucleotide == 'C' :
        return 'G'

    # TODO: implemenif nucleotide = 'A' :
        pass



def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    Complement = ''
    x = len(dna)-1 # finds last place in string
    while x >= 0: # goe from last letter to first letter
        Complement = Complement + get_complement(dna[x]) #adds inverse codon in reverse order
        x = x-1 # moves to the next letter (in reverse)
    return Complement

    # TODO: implement this
    pass


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGAGAGAGG") # does it keep the leftovers :)
    'ATGAGAGAGG'
    """

    stopcodons= ['TAG','TGA','TAA'] # the stop codons
    PC = (len(dna)/3)-1 # gives set that reps position of eact codon
    codons = [dna [i:i+3] for i in range(0, len(dna), 3)] # splits into codons
    x = 0
    string =''
    while x <= PC:
        codon = codons[x]
        if codon in stopcodons: #if dna has a stop codon
            return string # return before adding that codon
        string = string + codon # add this codon to the string
        x =x+1
    #leftovers= dna[-(len(dna)%3):]
    #print("leftovers: ", leftovers)
    #print("string + leftovers: ", string + leftovers)
    return dna

    # TODO: implement this
    pass


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """

    codons = [dna [i:i+3] for i in range(0, len(dna), 3)] # splits into codons
    startcodon = 'ATG'
    PC = (len(dna)//3)-1 # gives set that reps position of eact codon
    x= 0
    result = []
    while x <= PC:
        if startcodon == codons[x]:

            result.append(rest_of_ORF(dna[(x)*3:]))
            #print("Appending: ", rest_of_ORF(dna[(x)*3:]))
            x = x + (len(rest_of_ORF(dna[(x)*3:]))//3)
        else:
            x = x+1
    return result


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    dna1 = dna[1:]
    dna2 = dna[2:]
    Frame = find_all_ORFs_oneframe(dna)
    #print("find all orfs Frame: ", Frame)
    Frame1 = find_all_ORFs_oneframe(dna1)
    #print("find all orfs Frame1: ", Frame1)
    Frame2 = find_all_ORFs_oneframe(dna2)
    return Frame + Frame1 + Frame2


    # TODO: implement this
    pass


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    strand1 = dna
    #print("strand1: ", strand1)
    strand2 = get_reverse_complement(dna)
    #print("strand2: ", strand2)
    Frames1 = find_all_ORFs(strand1)
    #print("Frames1: ", Frames1)
    Frames2 = find_all_ORFs(strand2)
    #print("Frames2: ", Frames2)
    return Frames1 + Frames2
    # TODO: implement this
    pass


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    #chooses longest ORF
    return max(find_all_ORFs_both_strands(dna), key=len)

    # TODO: implement this

    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """

    orfshuffled = []
    x = 0
    while  x <= num_trials:
        #making a set of all the longest orfs from each suffle
        orfshuffled =  orfshuffled + [longest_ORF(shuffle_string(dna))]
        x = x + 1
    return len(max(orfshuffled, key=len))


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    codons = [dna [i:i+3] for i in range(0, len(dna), 3)] # splits into codons
    x = 0 # index of condons
    amino_acids = ""
    while x < len(dna)//3:
        amino_acid = aa_table[codons[x]]
        amino_acids = amino_acids + amino_acid
        x = x + 1
    return amino_acids



def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna, 100)
    all_ORFs = find_all_ORFs_both_strands(dna)
    big_aminos = []
    for string in all_ORFs:
        if len(string) > threshold:
            big_aminos = big_aminos + [coding_strand_to_AA(string)]

    return big_aminos

from load import load_seq
dna = load_seq("./data/X73525.fa")



if __name__ == "__main__":
    print (gene_finder(dna) )

    import doctest
    doctest.testmod()
