'''
Functions implemented here for use in my analysis of CpG content in various genomes with a focus on Covid19. 
Meant as an exercise in learning and making cool graphs.
'''

#Bipython imports since we are using Seq objects in these functions
import Bio
from Bio.Seq import Seq

#useful function from base python
from random import choices
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

def windowed_base_count(sequence, k_window=1):
    '''
    This will return the same result as DAMBE for pairs of residues.
    sequence is a Bio.Seq object that you want to know the count of k_sized sub_strings for
    k_window will define how large your sub_strings are.
    For DNA: k_window = 1 gives {A,T,C,G}
    k_window = 2 gives {AT, TT, TA, AA, AG, GG, GT, AC, CC, CT, TC, CA, CG, GA, TG, GC}
    '''
    sequence_sub_strings = []

    for i in range(0,len(sequence)-(k_window-1)):
        sequence_sub_strings.append(str(sequence.seq[i:i+k_window]))

    return Counter(sequence_sub_strings)

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

#print('pass')