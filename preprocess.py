import numpy as np
import re


def proteome_preprocess(proteins):

    # Replace the space in column names
    proteins.rename(columns= lambda x:x.replace(" ", "_"),inplace=True)
    # Replace the zeros with nan
    proteins.replace(0, np.nan, inplace=True)
    proteins = proteins.drop(columns=['Gene_names'])

    # Define regular express- contains
    multiid = re.compile('.*;.*')
    # Keep the first protein ID if multiple ids exist
    for i in range(proteins.shape[0]): 
        item = proteins.loc[i,'Majority_protein_IDs'] 
        if (re.match(multiid, item) != None): 
            str1 = item.split(';')[0] 
            proteins.loc[i,'Majority_protein_IDs'] = str1 
    return proteins

def phosphoproteome_preprocess(phosphos):

    # Replace the space in column names
    phosphos.rename(columns= lambda x:x.replace(" ", "_"),inplace=True)
    # Replace the zeros with nan
    phosphos.replace(0, np.nan, inplace=True)

    return phosphos

def transcripteome_preprocess(genes, gene_count):

    # Replace the characters in column names
    genes.rename(columns= lambda x:x.replace(" ", "_"),inplace=True)
    genes.rename(columns= lambda x:x.replace(".", ""),inplace=True)
    # Replace the zeros with nan
    genes.replace(0, np.nan, inplace=True)
    genes = genes.drop(columns = ['gene_name'])
    # Keep the genes that with a median greater than 5
    genes = genes[genes['gene_id'].isin(list(gene_count['id']))].reset_index(drop=True)

    return genes
