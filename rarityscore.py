import pandas as pd
import numpy as np

# Load trimer counts data
df = pd.read_csv('kmer.csv')

# Calculate rarity scores for each trimer
trimer_counts = df.iloc[:, 1:].values
trimer_frequencies = trimer_counts / trimer_counts.sum(axis=1, keepdims=True)
rarity_scores = np.power(trimer_frequencies, 2).sum(axis=1)

# Add rarity score as a new column to the dataframe
df['Rarity Score'] = rarity_scores

# Determine disease status based on rarity score
df['Disease Status'] = np.where(df['Rarity Score'] > 0.01, 1, 0)

# Save updated dataframe as CSV
df.to_csv('trimer_counts_with_target.csv', index=False)
