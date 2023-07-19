# Import nessecery packages
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

# Perform PCA and plot
from PCA_utils import pcaplt

pcaplt(proteins, 'Proteome', con_loc = (0.67,0.56))
pcaplt(phosphos, 'Phosphoproteome',t_loc = 'upper left')
pcaplt(genes, 'Transcriptome', early_time = '6 h',\
            con_loc = 'upper left', t_loc= 'lower center')

