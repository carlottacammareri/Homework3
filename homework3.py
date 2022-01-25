#exercise1fib


def Wascally_Wabbits(month, offspringCount):
    if month == 1:
        return 1
    elif month == 2:
        return offspringCount
    oneGen = Wascally_Wabbits(month -1, offspringCount)
    twoGen = Wascally_Wabbits(month -2, offspringCount)
    if month <= 4:
        return(oneGen + twoGen)
    return (oneGen + (twoGen*offspringCount))
print(Wascally_Wabbits(35,2))


#exercise2hamm

s= "GATACTCCCATTAGTAGAGCACGTCGGGGCCCTCCCGGATCCCGCCATCCAAAAGTCGGTACGCCGATAATCTGAATCAGTCCTTCAGTCATTAGTTAGGCGGACTCGGCCCAACGTGATCATTTTATGCTGGTTTACGCACCGCGACTGATCTGTTGTATACTTCCGTTCCCCTCTTCAGAAAGTGTGCGTGGAAATCACGTACTACATGGTACAACGGCATGTTCCAATGTGGGAATTAGTACTTAAAAGTAGTGTGTACAGACCCGGCATTTGGCCACGAACAATCAAATCGCAAGGTCCGTGCTTGGATTCTTCAGTGAGGGTGGCTGCTTTCAACCCTTTCAGCGGCCACGGCAAATATGGTGTAACGATGCTGAGCTAATAGTACTAGCTACTCTGGGGAGTGAATACTTCAGGTGGGGTGCCCTCCAAGTGTTAAGGGCCGTTTGTTGTAGCATGATTTGTTGCTTTAAAAACCACCACACACCTCAGCTTCCACTTTCCCGATATGATGCACATGCGGGGACATCACGTGAGATCCCGTACGGCTCCATAACCAGCTTCTGCCTGCGGCTGCTTTTTTCGCGAGCGGCTCGGGCGTTGGACAGCACCAAGGTGCCTGTTACTAGTCTAATCACTTGTAGGCGACTCGGAGCGCACACTGCATCTAGCCTGCGCTACTACCGGATTACCCCTGTCCAAGTCCTAGCTTTCTTATAGCGATACTGGGGCTGCCTGTGGTTCTCGTAGGGTTAAACTTTCGCTTTCACTGCGATGCCGATGGGCAATATTGCGGGGAGTCCTCGTAATCACCTAAACCACATAGGTGATGTCGCTTACGTCTCAGGGCGTTACTGCAGTTGCCTAATGGCTCGTCTTTCTCAACAGGAATTGCGCCAGTTGTGGAAAATGGTTTAATGCAGAAGGGATAGCACGCAACACGGCGGGTGGACT"
t= "AACGCCCGATTGCCTAGGGCTTATCGGGCCACTTACTGTGCTATCAATCCGTCCTTATGTACGCCGAGGCGCTATAGTGGTGGTCTATTTCGAGGTCTGAAACACTCGACCTGGGGGCATCATACAAGGTCGTTCTGCGTTTGGCGGTGGCTCAGACGTATACTACCCTACTGCTGTGCTGATTGGCAGCGAGGAAGTCAAATCGCTCATTGTCGAAACGCAATACGCCTAGTGTACTACAGTACTGGTTAAAATCGTGCACGTCACTGGCAGAGCATCAAAGATAACCAACTTACTTCGTTGGGGTCCGGTTTGGTAATCCCCGGGCTGTTCATACTAATCGATCCTGTCTCCCTCCGATACTGGTGTATCGCCTTTAAACAAATATGCATAATTGAACCCCGGTGTCTAAGTCACGCGTGTCGGGTGCTCCCAGCTTCATGAGACATCGATAGGTGCTTGACAAGGTCCTGTAAAAGACGCTGTTCTAGTCAGGAGGTACGTAGGGGTTAAGTGCAACAAACGATTAAATCCCGAGTCCAGCCGTATAGGTCCACAGCTTGAGGTTCACTTCGACGTTGTTGACGGCAAACGGAACACGCGCTCCACAGCACGACAGTTCCTCTTACTCGCTTACCGACTGGTAGTGTGAACCTAACCCAAAGGCCATCAAGCCTGTCCCACTACTTCGTTACCCCTGTCCAGCCTCGGGCTTTTTACTCGCGTACCCCCCCCAGCCTGTGGTCCTAGCAAGCATAAAATTACTCTCTCAGTGCGTTTTCAAGTGCTAATATGGACCCAATCCATCCTAAAAAGTCAAACCAGGTTGAAGCGGCCGCCCAAAGGTGAGAGCATTACAGTCCCCACATCGCTGTCTAATCTCCTATGGCTTGATTCATACGATTGTGGATAATCGTAAGCCTCACTAGAGTTTGACGTCATGCCGGACGATAGGAT"

def hamming_distance(s,t):
    if len(s)  != len(t):
        print('no')
    else:
        hamming_dist = 0
        for position in range(len(s)):
            if s[position] != t[position]:
                hamming_dist = hamming_dist +1
        return hamming_dist
print(hamming_distance(s,t))



#exercise3fibd

n = 85                                                                     
m = 17                                                                       
bun = [1, 1]                                                               
months = 2  
count = 0                                                                   
while months < n:                                                              
    if months < m:                                                             
        bun.append(bunnies[-2] + bun[-1])                              
    elif months == m or count == m + 1:                                        
        bun.append(bun[-2] + bun[-1] - 1)                          
    else:                                                                      
        bun.append(bun[-2] + bun[-1] - bun[-(                  
            m + 1)])                                                           
    months += 1                                                                
print(bun[-1])   



#exercise4mrna

codon_map = {
    'UUU': 'F',     'CUU': 'L',     'AUU': 'I',     'GUU': 'V',
    'UUC': 'F',     'CUC': 'L',     'AUC': 'I',     'GUC': 'V',
    'UUA': 'L',     'CUA': 'L',     'AUA': 'I',     'GUA': 'V',
    'UUG': 'L',     'CUG': 'L',     'AUG': 'M',     'GUG': 'V',
    'UCU': 'S',     'CCU': 'P',     'ACU': 'T',     'GCU': 'A',
    'UCC': 'S',     'CCC': 'P',     'ACC': 'T',     'GCC': 'A',
    'UCA': 'S',     'CCA': 'P',     'ACA': 'T',     'GCA': 'A',
    'UCG': 'S',     'CCG': 'P',     'ACG': 'T',     'GCG': 'A',
    'UAU': 'Y',     'CAU': 'H',     'AAU': 'N',     'GAU': 'D',
    'UAC': 'Y',     'CAC': 'H',     'AAC': 'N',     'GAC': 'D',
    'UAA': 'Stop',  'CAA': 'Q',     'AAA': 'K',     'GAA': 'E',
    'UAG': 'Stop',  'CAG': 'Q',     'AAG': 'K',     'GAG': 'E',
    'UGU': 'C',     'CGU': 'R',     'AGU': 'S',     'GGU': 'G',
    'UGC': 'C',     'CGC': 'R',     'AGC': 'S',     'GGC': 'G',
    'UGA': 'Stop',  'CGA': 'R',     'AGA': 'R',     'GGA': 'G',
    'UGG': 'W',     'CGG': 'R',     'AGG': 'R',     'GGG': 'G'
}
 
def codon_freq():
    freq = {}
    for k, v in codon_map.items():
        if v not in freq:
            freq[v] = 0
        freq[v] += 1
    return (freq)
 
def possible_sequences(sequence):
    f = codon_freq()
    n = f['Stop']
    for seq in sequence:
        n *= f[seq]
    return (n % 1000000) 
with open('rosalind_mrna.txt') as file: 
    protein_seq = file.read().strip()
print(possible_sequences(protein_seq))


#exercise5prtm

weights = {'A': 71.03711,             
           'C': 103.00919,            
           'D': 115.02694,            
           'E': 129.04259,            
           'F': 147.06841,            
           'G': 57.02146,             
           'H': 137.05891,            
           'I': 113.08406,            
           'K': 128.09496,            
           'L': 113.08406,            
           'M': 131.04049,            
           'N': 114.04293,            
           'P': 97.05276,             
           'Q': 128.05858,            
           'R': 156.10111,            
           'S': 87.03203,             
           'T': 101.04768,            
           'V': 99.06841,             
           'W': 186.07931,            
           'Y': 163.06333}            

with open('rosalind_prtm.txt', 'r') as f:
    for line in f:                    
        prot_seq = line.strip('\n')   

weight = 0                            
for aa in prot_seq:                   
    weight += weights[aa]             

print('%0.3f' % weight)


#exercise6lcsm

from Bio import SeqIO                      
sequences = []                             
handle = open('rosalind_lcsm.txt', 'r')
for record in SeqIO.parse(handle, 'fasta'):
    sequence = []                          
    seq = ''                               
    for nt in record.seq:                  
        seq += nt                          
    sequences.append(seq)                  
handle.close()
srt_seq = sorted(sequences, key=len)     
short_seq = srt_seq[0]                   
comp_seq = srt_seq[1:]                   
motif = ''                               
for i in range(len(short_seq)):          
    for j in range(i, len(short_seq)):   
        m = short_seq[i:j + 1]           
        found = False                   
        for sequ in comp_seq:            
            if m in sequ:                
                found = True            
            else:                        
                found = False           
                break                   
        if found and len(m) > len(motif):
            motif = m                    
print(motif)


#exercise7perm

import math                                         
import itertools                                    
n = 6                                             
print(math.factorial(n))                            
perm = itertools.permutations(list(range(1, n + 1)))
for i, j in enumerate(list(perm)):                  
    permut = ''                                
    for item in j:                                  
        permut += str(item) + ' '              
    print(permut) 


#exercise8revp

def palindrome(s, k):
	
	pali_pairs = [('A','T'),('T','A'),('G','C'),('C','G')]
	for i in range(len(s)-k+1):
		test = s[i:i+k]
		if all((test[j],test[k-j-1]) in pali_pairs for j in range(k)):
			
			print(i+1, k)

string = 'ATATCGCAAGCGAACCGTGGACTGTGTAGGGATCGCACGGCGCTTGCCTAACCTCAGGGTTACCGCTAATAGGGCAGCAAAATAAATGGCACCTGAGAGTAGGAGGTCGGTTCCACCTTCAGGGGTAAACTGGTACTGTACGGTGCGAGGTCCTCGTACAGGCTACCGAGGCGTAGCTGGTTAGAGCACGGCCGGTCAGTCGGCTCGCCATTCGTGCTGGAAGCTAAAAGGGGGACGTCCCTTCGTTAGGTTAACTAGGCTCGCTTGGAAACGCCTATTGACAAAGAGCTGCCAAACACCCCTTAGATAAAGTGGAACTACGAGTGTCCCAACACCTAACTAAGACCAGCACACACTATCATTTAAGAGATGACATGTGCGAACCAAAATTCTGGCCCGAGCTCGACGCTGTTTGCACTTGAAGCTAAGCGACATTCATTGTCACTGCGATGGCCAAATAGAGTACTAGCAAAACTGTGAACCAGTTACTCAACTAATTAGTTAACCATGAACTCCACATAGTTCTCGTAAACAGCACTGCCGTAATGTTTTGCGGCGAGTGTCCGGCGTGCACTCACCGTTCTTGACCCCCGATACGAACGACAAATATACTCTGTGGCAAAAATGTTTATTCTATTATTAGCAGGCCCGCCGGATCGACAAAAGGTATGGCTACCGATTGACTCCCAAGCGCGGCCAACCTGTTGTGTTTCGCGCACTAATATTTACCTACGAAATCCTCGGACTGTAATGGCTGTGCCTCCGGCATCGTGACTAGTGCCAGGCTCGCGCAGACGTGTTGACGACTCTGGGAGGAGCCCAATAGTCAGAACAGTTCTAATCGGTCTTGACGGGCTCAGTGGATAGTTAAGGCTTCGAACAGCGTGAGAAGGCGTTTTTGATGGGACGAGTTCGACGAGTTATCGGCTCTCTCCGGCTTGTGCGTTAGCATATTGGGTAGATGCGCCTTCTTAACTCTCGAACCTTGCAGCGC'
for i in range(4,14,2):
	palindrome(string, i)



#exercise9lexf

import itertools
 
with open("rosalind_lexf.txt") as f:
    data = f.read().split()
    let = data[:-1]
    n = int(data[-1])
perm = itertools.product(let, repeat = n)
out = []
for i, j in enumerate(list(perm)):
    permutation = ''
    for item in j:
        permutation += str(item)
    out.append(permutation)
out.sort()
for item in out:
    print(item, end="\n")



#exercise10lgis

data = []                                
with open('rosalind_lgis.txt', 'r') as f:   
    for line in f:                       
        for nr in line.split():          
            data.append(int(nr))         
perm = data[1:]


def increasing(seq):
    P = [None] * len(seq)
    M = [None] * len(seq)

    L = 1
    M[0] = 0
    for i in range(1, len(seq)):
        lo = 0
        hi = L
        if seq[M[hi - 1]] < seq[i]:
            j = hi
        else:
            while hi - lo > 1:
                mid = (hi + lo) // 2
                if seq[M[mid - 1]] < seq[i]:
                    lo = mid
                else:
                    hi = mid

            j = lo
        P[i] = M[j - 1]
        if j == L or seq[i] < seq[M[j]]:
            M[j] = i
            L = max(L, j + 1)

    result = []
    pos = M[L - 1]
    for k in range(L):
        result.append(seq[pos])
        pos = P[pos]

    return (result[::-1])


def decreasing(seq):
    P = [None] * len(seq)
    M = [None] * len(seq)

    L = 1
    M[0] = 0
    for i in range(1, len(seq)):
        lo = 0
        hi = L
        if seq[M[hi - 1]] > seq[i]:
            j = hi
        else:
            while hi - lo > 1:
                mid = (hi + lo) // 2
                if seq[M[mid - 1]] > seq[i]:
                    lo = mid
                else:
                    hi = mid

            j = lo
        P[i] = M[j - 1]
        if j == L or seq[i] > seq[M[j]]:
            M[j] = i
            L = max(L, j + 1)

    result = []
    pos = M[L - 1]
    for k in range(L):
        result.append(seq[pos])
        pos = P[pos]

    return (result[::-1])

incr = increasing(perm)
decr = decreasing(perm)

print(*incr)
print(*decr)



