import collections
from collections import OrderedDict
import math
import numpy as np
import pandas as pd
import csv
import os


ent={}
Q=[]
idname=[]
#with open(r"C:\Users\nishanthi\Desktop\ncbi_dataset\DS486008.1.txt", "r") as f:
for filename in os.listdir(r"C:\Users\nishanthi\Desktop\ncbi_dataset\seqout"):
    if filename.endswith('.txt'):
        with open(os.path.join(r"C:\Users\nishanthi\Desktop\ncbi_dataset\seqout" ,filename)) as f:
            for line in f:
                 # move the file pointer to the beginning of the file
                 f.seek(0)

                
                 sequence = f.read()
                 
                 
                 data = "".join(sequence.split("\n")[1:])
                 print(len(data))
                 combinations=["AAA", "AAC",	"AAG",	"AAT",	"ACA",	"ACC",	"ACG",	"ACT",
"AGA",	"AGC",	"AGG",	"AGT",	"ATA",	"ATC",	"ATG",	"ATT",
"CAA",	"CAC",	"CAG",	"CAT",	"CCA",	"CCC",	"CCG",	"CCT",
"CGA",	"CGC",	"CGG",	"CGT",	"CTA",	"CTC",	"CTG",	"CTT",
"GAA",	"GAC",	"GAG",	"GAT",	"GCA",	"GCC",	"GCG",	"GCT",
"GGA",	"GGC",	"GGG",	"GGT",	"GTA",	"GTC",	"GTG",	"GTT",
"TAA",	"TAC",	"TAG",	"TAT",	"TCA",	"TCC",	"TCG",	"TCT",
"TGA",	"TGC",	"TGG",	"TGT",	"TTA",	"TTC",	"TTG",	"TTT",
]

                 
                 def count_kmers(sequence, k):
                     d = collections.defaultdict(int)
                     for i in range(len(sequence)-(k-1)):
                         d[sequence[i:i+k]] +=1

                     for key in list(d):
                          if "N" in key:
                              del d[key]

                     return d
                 def probabilities(kmer_count, k):
                      probabilities = collections.defaultdict(float)
                      N = len(data)
                      for key, value in kmer_count.items():
                          probabilities[key] = float(value) / (N - k + 1)
                      return probabilities

                 def chaos_game_representation(count_kmers, k):
                      array_size = int(math.sqrt(4**k))
                      chaos = []
                      for i in range(array_size):
                           chaos.append([0]*array_size)
                      maxx = array_size
                      maxy = array_size
                      posx = 1
                      posy = 1
                      for key, value in count_kmers.items():
                            for char in key:
                                if char == "T":
                                     posx += maxx // 2
                                elif char == "C":
                                    posy += maxy // 2
                                elif char == "G":
                                    posx += maxx // 2
                                    posy += maxy // 2
                                maxx = maxx // 2
                                maxy //= 2
                            chaos[posy-1][posx-1] = key+"\n"+str(value)
                            maxx = array_size
                            maxy = array_size
                            posx = 1
                            posy = 1

                      return chaos 
                 f3 = count_kmers(data, 1)
                 f4 = count_kmers(data, 2)
                 f5 = count_kmers(data, 3)                 
                 new_d = OrderedDict(sorted(f5.items(), key=lambda a:a[0]))
                 u=dict(new_d)
                 for i in combinations:
                     if i not in u.keys():
                         u.update({i:0})
                 #print(set1)


                 

                 
                 #f6 = count_kmers(data, 4)
                 #f7 = count_kmers(data, 5)
                 ent.update(f5)
                 
                 Z=ent.copy()
                 
                 Q.append(Z)
                 #print(ent)  #dict
                 #print(Q) #list
                 


                
                 f3_prob = probabilities(f3, 1)
                 f4_prob = probabilities(f4, 2)
                 f5_prob = probabilities(f5, 3)
                 #f6_prob = probabilities(f6, 4)
                 #f7_prob = probabilities(f7, 5)
                 
                 '''print(f3_prob)
                 print(f4_prob)
                 print(f5_prob)
                 print(f6_prob)
                 print(f7_prob)
                 
                 print(chaos_game_representation(f3, 1))'''
                 #print(chaos_game_representation(f4, 2))
                 print(chaos_game_representation(f5, 3))
                 #print(chaos_game_representation(f6, 4))
                 #print(chaos_game_representation(f7, 5))
'''
#jaccard similarity               
def jaccard_similarity(set1, set2):
     intersection = set1.keys() & set2.keys()
     #print(intersection)
     union = set1.keys() | set2.keys()
     numerator = sum(min(set1[k], set2[k]) for k in intersection)
     denominator = sum(max(set1.get(k, 0), set2.get(k, 0)) for k in union)
     return numerator / denominator



first_sequence = Q[0]
for i in range(1, len(Q)):
    similarity = jaccard_similarity(first_sequence, Q[i])
    print(f"Jaccard similarity between first sequence and sequence {i+1}: {similarity}")

#jaccard = jaccard_similarity(Q[0], Q[i])
#print(jaccard)'''
for filename in os.listdir(r"C:\Users\nishanthi\Desktop\ncbi_dataset\test"):
    if filename.endswith('.txt'):
        with open(os.path.join(r"C:\Users\nishanthi\Desktop\ncbi_dataset\test" ,filename)) as f:
            for line in f:
                 # move the file pointer to the beginning of the file
                 f.seek(0)

                
                 # read the first line of the file
                 header = f.readline().strip()
                 idname.append(header)
                 sequence = f.read()

with open("kmer.csv", 'a', newline='') as result:
    fieldnames = ['id'] + list(Q[0].keys())
    writer = csv.DictWriter(result, fieldnames=fieldnames)
    writer.writeheader()
    for i, row in enumerate(Q):
        row_values = [idname[i]] + list(row.values())
        writer.writerow(dict(zip(fieldnames, row_values))) 
