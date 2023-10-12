import random
import os
import collections
from collections import OrderedDict
import csv

dic={}
Q=[]
idname=[]

for filename in os.listdir(r"C:\Users\nishanthi\Desktop\ncbi_dataset\seqout"):
    if filename.endswith('.txt'):
        with open(os.path.join(r"C:\Users\nishanthi\Desktop\ncbi_dataset\seqout" ,filename)) as f:
            for line in f:
                 # move the file pointer to the beginning of the file
                 f.seek(0)

                 sequence = f.read()
                 
                 
                 data = "".join(sequence.split("\n")[1:])
                 #print(data)

                 def find_kmers(dna_sequence, k):
                      d = collections.defaultdict(int)
                      for i in range(len(dna_sequence) - k + 1):
                           d[dna_sequence[i:i+k]] +=1
                      for key in list(d):
                          if "N" in key:
                              del d[key]

                      return d


                 def generate_random_weights(k):
    
                    weights = {}
                    for i in k:
                        weights[i] = random.random()
                    return weights

                 def weighted_kmer_count(sequence, weights):
                        
    
                        count = {}
                        for i in range(len(sequence)-k+1):
                            kmer = sequence[i:i+k]
                            if kmer in weights:
                                weight = weights[kmer]
                            if kmer in count:
                                count[kmer] += weight
                            else:
                                count[kmer] = weight
                        return count

                     

                 

                    
                 f5 = find_kmers(data, 3)
                 #print(f5)
                 weights = generate_random_weights(f5)
                 print(weights)
                

                 #sequences = ["ATGCGCTAG", "CCGATAG", "ATAGCGATGCTA"]
                 #weights = [1.0, 2.0, 0.5]
                 k =3 
                 
                 weighted_counts = weighted_kmer_count(data,weights)
                 r = 4
  
                 # loop to iterate for values 
                 wc = dict()
                 for key in weighted_counts:
      
                      # rounding to K using round()
                      wc[key] = round(weighted_counts[key], r)
                  #print(weighted_counts)

                 dic.update(wc)
                 
                 Z=dic.copy()
                 
                 Q.append(Z)


 
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




with open("weightedkmer.csv", 'a', newline='') as result:
    fieldnames = ['id'] + list(Q[0].keys())
    writer = csv.DictWriter(result, fieldnames=fieldnames)
    writer.writeheader()
    for i, row in enumerate(Q):
        row_values = [idname[i]] + list(row.values())
        writer.writerow(dict(zip(fieldnames, row_values)))
