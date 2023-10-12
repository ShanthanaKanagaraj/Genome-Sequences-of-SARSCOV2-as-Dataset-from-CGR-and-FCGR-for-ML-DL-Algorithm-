from collections import defaultdict
from collections import OrderedDict
import os
import math
import csv
import xlwt
from sklearn.preprocessing import MinMaxScaler
import numpy as np

dic={}
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
                 #print(len(data))
                 def normalize_kmer_count(data, k):
                     kmer_counts = defaultdict(int)
                     n = len(data)
                     # Count all k-mers in the sequence
                     for i in range(n - k + 1):
                         kmer = data[i:i+k]
                         kmer_counts[kmer] += 1

                     # Compute the total number of k-mers
                     total_kmers = sum(kmer_counts.values())

                     # Compute the normalized count for each k-mer
                     normalized_counts = {}
                     for kmer, count in kmer_counts.items():
                         normalized_counts[kmer] = count / (total_kmers - k + 1)
                     return normalized_counts
                # Example usage
                 k = 3
                 normalized_counts = normalize_kmer_count(data, k)                 
                 #print(normalized_counts)

                 

                
                 dic.update(normalized_counts)
                 
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

 


with open("Normalizedkmer.csv", 'a', newline='') as result:
    fieldnames = ['id'] + list(Q[0].keys())
    writer = csv.DictWriter(result, fieldnames=fieldnames)
    writer.writeheader()
    for i, row in enumerate(Q):
        row_values = [idname[i]] + list(row.values())
        writer.writerow(dict(zip(fieldnames, row_values)))
     



