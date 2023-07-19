
import pandas as pd

#Load datasets
gene_count = pd.read_csv('./datasets/gene_count.csv')
genes = pd.read_csv('./datasets/gene_expression.csv')
proteins = pd.read_csv('./datasets/proteins.csv')
phosphos = pd.read_csv('./datasets/PhosphoSites.csv')

# Preprocess datasets
from preprocess import proteome_preprocess, phosphoproteome_preprocess,transcripteome_preprocess

proteins = proteome_preprocess(proteins)
phosphos = phosphoproteome_preprocess(phosphos)
genes = transcripteome_preprocess(genes, gene_count)

# Preprocess datasets for MOFA
from MOFA_utils import pre_process4MOFA
protein4mofa =  pre_process4MOFA(proteins, 'Majority_protein_IDs', 'Proteome')
print(protein4mofa.head())

phosphos = phosphos.replace('-','_', regex=True)
phosphos4mofa =  pre_process4MOFA(phosphos, 'Phosphosite_ID', 'Phosphoproteome')
phosphos4mofa = phosphos4mofa.replace('_','-', regex=True)
print(phosphos4mofa.head())

genes4mofa =  pre_process4MOFA(genes, 'gene_id', 'Transcriptome')
print(genes4mofa.head())

# Write results to CSV files
write.csv(protein4mofa, file = "protein4mofa.csv")
write.csv(phosphos4mofa, file = "phosphos4mofa.csv")
write.csv(genes4mofa, file = "genes4mofa.csv")
