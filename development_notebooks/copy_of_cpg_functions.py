'''
Functions implemented here for use in my analysis of CpG content in various genomes with a focus on Covid19. 
Meant as an exercise in learning and making cool graphs.
'''

#Bipython imports since we are using Seq objects in these functions
import Bio
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import numpy as np
import matplotlib.pyplot as plt

#useful function from base python
from random import choices, sample
from collections import Counter

def generate_random_sequence(k, weights):
    '''
    Generates a random sequence using random.choices()
    Weights is a list of weightings for each choice. If you assume sum(weights) = 1 
    you can then state that these weights are probabilities for each choice.
    weights list should be in order: ['A','T','C','G']
    Note: this returns a `Seq` object.
    '''
    random_seq_list = choices(['A','T','C','G'], weights=weights, k=k)
    random_sequence = Seq("".join(random_seq_list))
    return random_sequence


def random_shuffled_genome(sequence):
    '''
    This will build a random genome with the exact same base compositions and lengths of the input sequence
    Assumes sequence is a Bio.Seq object.
    Returns as a SeqRecord.
    Based off Biopython cookbook code for this function.
    '''
    raw_shuffled = ''.join(sample(list(sequence), len(sequence)))
    return SeqRecord(Seq(raw_shuffled, sequence.seq.alphabet), id="Shuffled", description="Based on %s" % sequence.id)

def windowed_base_count(sequence, k_window=1, as_ratio=False):
    '''
    This will return the same result as DAMBE for pairs of residues.
    sequence is a Bio.Seq object that you want to know the count of k_sized sub_strings for
    k_window will define how large your sub_strings are.
    For DNA: k_window = 1 gives {A,T,C,G}
    k_window = 2 gives {AT, TT, TA, AA, AG, GG, GT, AC, CC, CT, TC, CA, CG, GA, TG, GC}
    '''
    sequence_sub_strings = []

    for i in range(0,len(sequence)-(k_window-1)):
        sequence_sub_strings.append(str(sequence.seq[i:i + k_window]))

    if as_ratio == False:
        frequencies = Counter(sequence_sub_strings)
    elif as_ratio == True:
        frequencies = Counter(sequence_sub_strings)
        total = sum(frequencies.values())
        for key in frequencies:
            frequencies[key] = frequencies[key]/total

    return frequencies

def calculate_icpg(singles, doubles, diagnostic_print = False):
    '''
    Calculates simple Icpg values given frequency hashes (dict obj) from windowed_base_count 
    where singles is k_window=1 and doubles is k_window=2.
    Returns Icpg as float.
    diagnostic_print == True will return relevant frequency values
    '''
    #calculate frequencies from dict frequency hashes. Currently includes things that are not needed.
    f_CG = doubles['CG']/sum(doubles.values())
    f_GC = doubles['GC']/sum(doubles.values())
    f_G = singles['G']/sum(singles.values())
    f_C = singles['C']/sum(singles.values())
    
    if diagnostic_print == True:
        print('f(CG) = %.3f f(G) = %.3f, f(C) = %.3f'%(f_CG, f_G, f_C))
        print('f(CG) = %.3f f(GC) = %.3f'%(f_CG, f_GC))

    return f_CG/(f_G*f_C)

def symmetrized_Icpg(singles, doubles):
    '''
    Calculates symmetrized Icpg values using only frequencies 
    '''
    f_CG = doubles['CG']/sum(doubles.values())
    f_G = singles['G']/sum(singles.values())
    f_C = singles['C']/sum(singles.values())

    pCG = (2*(f_CG+f_CG))/((f_C+f_G)**2)
    return pCG

def find_string_indicies(search_string, sequence):
    '''
    Returns string indicies of search_string for sequence as a list, index starts at 0, not 1.
    '''
    i = len(sequence)
    results = []
    while True:
        i = sequence.rfind(search_string,0,i)
        if i == -1:
            break
        results.append(i)
    return sorted(results)

def singles_doubles_graph(genome):
    singles = windowed_base_count(genome, k_window=1)
    doubles = windowed_base_count(genome, k_window=2)
    print('Percent GC for ref sequence: %.2f'%((singles['C']+singles['G'])/len(genome)))

    p_CG = calculate_icpg(singles, doubles)
    p_sym_CG = symmetrized_Icpg(singles, doubles)
    print('I_CpG = %.3f\nI*CpG = %.3f'%(p_CG,p_sym_CG))

    #plot this
    fig, ax = plt.subplots(2, figsize=(10,10))

    i = 0
    for results in [singles, doubles]:
        x = list(results.keys())
        heights = list(results.values())
        ax[i].bar(x,heights)
        ax[i].set_xlabel('sequence pairs')
        ax[i].set_ylabel('absolute count')
        i += 1

    #matplotlib does not like recieving hashes so need to replace with list()
    ax[0].plot(list(singles.keys()),[7475 for x in range(0,len(singles.keys()))], 'orange')

    #the E(x) for each doublet of bases is E(single_base)/4 since each doublet event is not independent of all other doublet events.
    ax[1].plot(x,[7475/4 for x in range(0,len(doubles.keys()))], 'orange')

    plt.show()

def find_string_indicies(search_string, sequence):
    '''
    This is a simple search loop to find the indicies for each      search_string in sequence
    '''
    if len(search_string) >= len(sequence):
        raise ValueError("search_string must be shorter length                             than sequence")
    i = len(sequence)
    results = []
    while True:
        i = sequence.rfind(search_string,0,i)
        if i == -1:
            break
        results.append(i)
    return sorted(results)

def barh_substring_graph(sequence, substring):
    '''
    Assumes sequence is a Bio.Seq.seq object
    Makes a broken barh graph using matplotlib of substring locations within sequence. 
    This will NOT show all sustring locations as there is not enough horizontal space for sufficient pixels to show.
    To help with this each width = 10 so that substrings are more visible.
    '''
    substr_indicies = find_string_indicies(substring, sequence)
    barh_xvals = list(
        zip(substr_indicies, [15 for x in range(0, len(substr_indicies))]))
    fig, ax = plt.subplots(figsize=(30,3))
    ax.broken_barh(xranges=barh_xvals, yrange=(1,1))
    plt.show()

def slidingWindowCG(sequence, window=500):
    '''
    Calculates sliding window calculations of frequency of CG% in sequence.
    Returns matplotlib plot.
    window = 500 by default.
    Assumes sequence is Bio.Seq object, don't call .seq attribute!
    '''

    substr_content_windowed = []
    for i in range(0, len(sequence), int(window/10)):
        sub_sequence = sequence[i:i + window]
        base_counts = windowed_base_count(sub_sequence)
        #this ends the function once the window starts sliding off the end of the sequence
        summed_counts = sum(base_counts.values())
        if summed_counts != window:
            break
        C = base_counts['C']/summed_counts
        G = base_counts['G']/summed_counts
        substr_content_windowed.append(C+G)
    fig, ax = plt.subplots(figsize=(30, 7))
    ax.plot(substr_content_windowed)
    plt.show()


#print('pass')
