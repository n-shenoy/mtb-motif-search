#!/usr/bin/env python

""" Functions for DNA motif search """

__author__ = "Navami Shenoy"

import numpy as np
import random


def readFASTA(filename):
    # parse and read FASTA file
    seqs = []
    with open(filename, 'r') as f:    
        for line in f:
            if line[0] != '>':
                seqs.append(line.rstrip())
    return seqs  



def Hamming_distance(p , t):
    # compute the Hamming distance between two strings 
    #returns the number of mismatches between the text
    #and the pattern at every reading frame
    mismatches = []
    for i in range(0, len(t) - len(p) + 1):
        num_mismatches = 0
        for j in range(0, len(p)):
            if t[i+j] != p[j]:
                num_mismatches += 1

        mismatches.append(num_mismatches)

    return mismatches



def distance_btw_pattern_and_strings(pattern, dna):
    #returns the sum of Hamming distances between pattern
    #and each dna sequence in the list 
    l = len(pattern)
    distance = 0
    for string in dna:
        mismatches = Hamming_distance(pattern, string)
        min_pos = mismatches.index(min(mismatches))
        d = Hamming_distance(pattern, string[min_pos:min_pos+l])
        distance += d[0]
    
    return distance



def all_strings(dna, k):
    #find all possible kmers from each sequence in the list
    kmers = []
    for string in dna:
        for i in range(0, len(string) - k + 1):
            kmers.append(string[i:i+k])

    return set(kmers) 
        

#first algorithm for motif search:
def median_string(dna, k):
    #returns a k-mer in the pattern that minimizes the distance
    #between the pattern and all given sequences
    # among all possible choices of k-mers
    kmers = all_strings(dna, k)
    distances = []
    for kmer in kmers:
        distances.append(distance_btw_pattern_and_strings(kmer, dna))
    
    return list(kmers)[distances.index(min(distances))]




def most_probable_kmer(text, k, matrix):
    #finds the most probable kmer in the text (sequence)
    #based on a probability profile matrix for each
    #each base at each position in the sequence
    bases = {'A' : 0,'C' : 1,'G' : 2,'T' : 3}
    max_prob = 0
    matrix = matrix.transpose()
    most_prob_kmer = ''
    for i in range(0, len(text) - k + 1):
        kmer = text[i:i+k]
        kmer_prob = [] 
        for i in range(0, len(kmer)):
            kmer_prob.append(matrix[i][bases[kmer[i]]])
        prob = np.prod(kmer_prob)
        if prob > max_prob:
            max_prob = prob
            most_prob_kmer = kmer
    return most_prob_kmer




def generate_profile(motifs):
    #given a collection of motifs, return a matrix
    #containing the probability of each base at each
    #position in the motifs 
    k = len(motifs[0]) #length of each motif = kmer length
    t = len(motifs) #number of sequences/motifs
    bases = {'A' : 0,'C' : 1,'G' : 2,'T' : 3}

    matrix_z = np.zeros([4,k])#matrix without Laplace's correction
    matrix = np.ones([4,k]) #matrix with Laplace's correction

    for i in range(0, t):
        motif = motifs[i]
        for j in range(0, k):
            matrix[bases[motif[j]]][j] += 1
            matrix_z[bases[motif[j]]][j] += 1
            
    for i in range(0, len(matrix[0])):
        for j in range(0, len(matrix)):
            matrix[j][i] = matrix[j][i]/(2*t)
            matrix_z[j][i] = matrix_z[j][i]/t
    return matrix, matrix_z




def score(matrix,t):
    #given a profile probability matrix WITHOUT 
    #Laplace's correction, return the score based
    #on the least popular bases at each position
    #higher score = greater mismatch
    scores = []
    for i in range(0, len(matrix[0])):
        column_probs = []
        for j in range(0, len(matrix)):
            column_probs.append(matrix[j][i])
        column_probs.remove(max(column_probs))
        
        for prob in column_probs:
            scores.append(prob*t)
    
    return sum(scores)



#second algorithm for motif search:
def greedy_motif_search(k, t, dna):
    #returns a list of motifs by generating a profile
    #for each kmer from the sequences, finding the most
    #probable kmer for that profile, and returning the kmers
    #with the least scores
    best_motifs = []
    for string in dna:
        best_motifs.append(string[:k])

    for i in range(0,len(dna[0]) - k + 1):
        motif = dna[0][i:i+k]
        motifs = [motif]
        for string in dna[1:]:   
            profile, profile_z = generate_profile(motifs)
            motifs.append(most_probable_kmer(string, k, profile))
        
        __, best_motifs_profile_z = generate_profile(best_motifs)
        __, motifs_profile_z = generate_profile(motifs)
        if score(motifs_profile_z,t) < score(best_motifs_profile_z,len(best_motifs)):
           best_motifs = motifs
    return best_motifs
        
            

    
def consensus(profile_z):
    #given a profile (ideally without Laplace correction,
    #to follow the method in the Bioinformatics course), 
    #returns a consensus string formed from a list of motifs
    bases = {0 : 'A', 1 : 'C', 2 : 'G', 3 : 'T'}
    consensus_string = []
    for i in range(0, len(profile_z[0])):
        probs = []
        for j in range(0, len(profile_z)):
            probs.append(profile_z[j][i])
        
        max_prob = max(probs)
        max_prob_pos = probs.index(max_prob)
        new_probs = probs
        new_probs[max_prob_pos] = 0
        new_max_prob = max(new_probs)
        if new_max_prob == max_prob:
            tie = '(' + bases[max_prob_pos] + '/' + bases[probs.index(new_max_prob)] + ')'
            consensus_string.append(tie)
        else:
            consensus_string.append(bases[max_prob_pos])
    
    return ''.join(consensus_string)



        
def generate_motifs(profile, dna):
    #returns most probable kmers in the given sequences 
    #based on a probability profile matrix
    bases = {0 : 'A', 1 : 'C', 2 : 'G', 3 : 'T'}
    motifs = []
    k = len(profile[0])
    for i in range(0, len(dna)):
        most_prob_kmer = most_probable_kmer(dna[i], k, profile)
        motifs.append(most_prob_kmer)
    
    return motifs


def random_motifs(k,t,dna):
    #generates a list of motifs (kmers) randomly selected
    #from a set of dna sequences
    seq_length = len(dna[0])
    motifs = []
    for i in range(0, t):
        pos = random.randint(0,seq_length-k)
        motifs.append(dna[i][pos:pos+k])
    return motifs


#third motif search algorithm
def randomized_motif_search(k,t,dna):
    #finds a list of motifs that minimizes the score 
    #by randomly generating motifs over and over
    best_score = float('inf')
    best_motifs = []
    i = 0
    while i < 1000:
        m = random_motifs(k,t,dna)
        motifs = m
        while True:
            profile, profile_z = generate_profile(m)
            m = generate_motifs(profile, dna)
            profile2, profile2_z = generate_profile(m)
            if score(profile2_z,t) < score(profile_z,t):
                motifs = m  
            else:
                break
        __, profile_z = generate_profile(motifs)
        current_score = score(profile_z,t)
        if current_score < best_score:
            best_score = current_score
            best_motifs = motifs
        i += 1  
    return best_motifs



#last motif search algorithm
def gibbs_sampler(k,t,N,dna, random_starts):
    #returns a list of best motifs by randomly generating motifs
    #and changing one kmer based on probability 
    bases = {'A' : 0,'C' : 1,'G' : 2,'T' : 3}
    best_motifs = dna
    r = 0
    while r < random_starts:
        m = random_motifs(k,t,dna)
        n = 0 
        while n < N:
            motifs = m[:]
            I = random.randint(0,t-1)
            motifs.remove(motifs[I])
            profile, profile_z = generate_profile(motifs)
            # for each kmer in dna[i], calculate its probability
            # of being generated by the profile:
            removed_seq = dna[I]
            kmers = []
            probs = []
            for i in range(0, len(removed_seq)-k+1):
                kmer = removed_seq[i:i+k]
                prob = 1
                for i in range(0, len(kmer)):
                    prob = prob*profile[bases[kmer[i]]][i]  
                kmers.append(kmer)
                probs.append(prob) 
            
            new_start_pos = random.choices(range(0, len(kmers)), weights=probs)  
            m[I] = kmers[new_start_pos[0]]
            n += 1

        profile_m, __ = generate_profile(m)
        profile_bm, __ = generate_profile(best_motifs)
        if score(profile_m,t) < score(profile_bm,t):
            best_motifs = m[:]
        r+=1
    return best_motifs



 