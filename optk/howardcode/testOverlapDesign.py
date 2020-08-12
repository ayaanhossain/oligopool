## Design oligonucleotide pools containing overlap regions to construct larger DNA fragments

# Test #1
# How fast can bitwise operations calculate the pairwise Hamming distances of N sequence variants? 
# N^2 (XOR + SUM) operations on bitstrings

import random, time
from bitarray import bitarray
from bitarray.util import count_xor,  count_n, strip
from itertools import combinations
from Bio.SeqUtils import MeltingTemp as mt

def toBitArrayList(seqList):
    symbols = {'A' : [False, False],
               'T' : [True, False],
               'C' : [False, True],
               'G' : [True, True]
               }
    
    bitarrayList = []
    for seq in seqList:
        x = []
        for letter in seq:
            x.extend(symbols[letter])
        bitarrayList.append( bitarray(x) )
    return bitarrayList

def calculateHammingArray(a1, a2):
    return a1 ^ a2

def calculateHammingDistance(a1, a2):
    return count_xor(a1, a2)

def calculateTmWindows(seq, WINDOW):
    # WINDOW, nucleotides
    #Accurate, but slow [369 seconds per 10000 sequences]
    #TM_list = [mt.Tm_NN(seq[x:x+WINDOW], Na=50, Tris=10, Mg=1.5, dNTPs=0.6, saltcorr=7) for x in range(len(seq)-WINDOW)]
    
    #Less accurate, but fast! [9.8 seconds per 10000 sequences]
    TM_list = [4.0 * seq[x:x+WINDOW].count('G') + 4.0 * seq[x:x+WINDOW].count('C') + 2.0 * seq[x:x+WINDOW].count('A') + 2.0 * seq[x:x+WINDOW].count('T') for x in range(len(seq)-WINDOW)] 
    return TM_list

def identifyCandidateOverlaps(seq, OVERLAP, MIN_TM = 50, WINDOW = 20):
    # MIN_TM, degC
    # WINDOW, nucleotides
    # L length, nucleotides
    
    #Identify consecutive regions with valid TMs
    candidates = []
    TM_list = calculateTmWindows(seq, WINDOW)
    counter = 0
    currentPos = None
    for (pos, TM) in enumerate(TM_list):
        if TM >= MIN_TM and counter > 0:
            counter += 1
        elif TM >= MIN_TM and counter == 0:
            currentPos = pos
            counter += 1
        elif TM < MIN_TM and counter > 0:
            if counter >= (OVERLAP - WINDOW): candidates.append( (currentPos, counter) )
            counter = 0
        else:
            pass
    if counter >= (OVERLAP - WINDOW): candidates.append( (currentPos, counter) )
    return candidates

def generateTestSequences(N, L):
    seqList = []
    for x in range(N):
        seq = "".join([random.choice(['A','T','C','G']) for n in range(L)])
        seqList.append(seq)
    return seqList

if __name__ == "__main__":

    L = 500
    N = 2000
    MER = 170
    OVERLAP = 40
    HMIN = 10
    
    testSeqList = generateTestSequences(N, L)

    t1 = time.time()
    #O(N) calculation, ~0.05 seconds per sequence
    candidateOverlaps = {}
    for (n,seq) in enumerate(testSeqList):
        candidateOverlaps[n] = identifyCandidateOverlaps(seq, OVERLAP)
    
    oligo_counter = 1
    overlap_range = {}
    for (n,seq) in enumerate(testSeqList):
        overlap_range[n] = []
        for candidate in candidateOverlaps[n]:
            x = oligo_counter*MER
            if  x - candidate[0] > 0: #this candidate overlap fits in the current MER
                overlap_range[n] = (candidate[0], min(x, candidate[0] + candidate[1]))
    
    aList = toBitArrayList(testSeqList)
    combos = list(combinations(range(N), 2))
    
    shifts = {}
    Passed = {}
    
    for n in range(N):
        shifts[n] = 0
    
    for (i1, i2) in combos:
        Passed[(i1,i2)] = False
    
    round = 0    
    while True:
        print "ROUND #%s" % (round + 1)
        for (i1, i2) in combos:
            if not Passed[(i1,i2)]:
                r11 = overlap_range[i1][1]-OVERLAP-shifts[i1]
                r12 = overlap_range[i1][1]-shifts[i1]
                r21 = overlap_range[i2][1]-OVERLAP-shifts[i2]
                r22 = overlap_range[i2][1]-shifts[i2]
        
                y = calculateHammingDistance(aList[i1][r11:r12], aList[i2][r21:r22])
                if y >= HMIN:
                    Passed[(i1,i2)] = True
                else:
                    print "FAILED: ", testSeqList[i1][r11:r12], testSeqList[i2][r21:r22], (i1,i2)
                    shifts[i2] += 5
                    for n in range(N):
                        Passed[(n,i2)] = False
        round += 1
        print "PASSED: ", sum(Passed.values()), " out of ", len(combos)
        if sum(Passed.values()) == len(combos): break
    print shifts
                
    
    #O(N^2) calculation, but super fast per sequence
#    for (i1, i2) in combos:
#        y = calculateHammingArray(x[i1], x[i2])
#        y = calculateHammingDistance(x[i1], x[i2])
    t2 = time.time()
    print "Elapsed time: ", t2 - t1, " seconds."
    
    
        
