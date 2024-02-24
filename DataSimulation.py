import numpy as np
import pandas as pd
import random
import json

from typing import List, Tuple, Dict

with open('data.json') as f:
    data = json.load(f)

nts = data['nucleotides']
amino_acids = data['amino_acids']


# Function to simulate data
def random_nt_sequence(length: int) -> str:
    """
    Generate a random nucleotide sequence

    Parameters
    ----------
    length : int, length of the sequence
    """
    prob_nt = {"A": 0.245, "C": 0.245, "G": 0.245, "T": 0.245, "N": 0.02}
    return ''.join(np.random.choice(list(prob_nt.keys()), length, p=list(prob_nt.values())))


def simulate_vj_genes(n_genes: int, len_gene: int, epitope_pool: List) -> Dict:
    """
    Simulate V and J genes

    Parameters
    ----------
    n_genes : int, number of genes to simulate
    len_gene : int, length of the genes
    epitope_pool : List, epitopes
    """
    genes = {}
    for i in range(n_genes):
        gene = random_nt_sequence(len_gene)
        # TODO - A gene can have multiple epitopes
        num_epitopes = random.randint(1, 10)
        gene_epitopes = random.sample(epitope_pool, num_epitopes)
        genes[i] = {'gene': gene, 'epitopes': gene_epitopes}
    return genes


def simulate_epitopes(n_epitopes: int) -> List:
    """
    Simulate epitopes

    Parameters
    ----------
    n_epitopes : int, number of epitopes to simulate
    len_epitope : int, length of the epitopes
    """
    epitopes = []
    amino_acids_list = list(amino_acids.keys())
    len_epitopes = random.randint(8, 11)
    for i in range(n_epitopes):
        epitope = ''.join(random.choice(amino_acids_list) for _ in range(len_epitopes))
        epitopes.append(epitope)
    return epitopes


def simulate_overlap_reads(n_reads: int, len_read: int, v_gene_pool: Dict, j_gene_pool: Dict) -> Dict:
    """
    Simulate reads that overlap with genes

    Parameters
    ----------
    n_reads : int, number of reads to simulate
    len_read : int, length of the reads
    v_gene_pool : Dict, V genes
    j_gene_pool : Dict, J genes
    """
    reads = {}
    count = 0
    while count < n_reads:
        v_gene = random.choice(list(v_gene_pool.keys()))
        j_gene = random.choice(list(j_gene_pool.keys()))
        v_gene_seq = v_gene_pool[v_gene]['gene']
        j_gene_seq = j_gene_pool[j_gene]['gene']
        len_v_overlap = random.randint(0, min(len_read, len(v_gene_seq)))
        len_j_overlap = random.randint(0, min(len_read - len_v_overlap, len(j_gene_seq)))
        len_d = len_read - len_v_overlap - len_j_overlap
        if len_d <= 0 or len_d >= 12: continue
        d_seq = random_nt_sequence(len_d)
        # Take the last v_overlap nucleotides from the V gene, d_seq, and the first j_overlap nucleotides from the J gene
        read_seq = v_gene_seq[-len_v_overlap:] + d_seq + j_gene_seq[:len_j_overlap]
        read_epitopes = list(set(v_gene_pool[v_gene]['epitopes'] + j_gene_pool[j_gene]['epitopes']))
        reads[count] = {'read': read_seq, "v_gene": v_gene_pool[v_gene], "d_gene": d_seq, "j_gene": j_gene_pool[j_gene], "epitopes": read_epitopes}
        count += 1
        assert len(read_seq) == len_read
        assert len(read_seq) == len_v_overlap + len_d + len_j_overlap
        assert read_seq[:len_v_overlap] == v_gene_seq[-len_v_overlap:] if len_v_overlap > 0 else True
        assert read_seq[len_v_overlap:len_v_overlap + len_d] == d_seq if len_d > 0 else True
        assert read_seq[-len_j_overlap:] == j_gene_seq[:len_j_overlap] if len_j_overlap > 0 else True

    return reads


def simulate_random_reads(n_reads: int, len_read: int) -> List:
    """
    Simulate random reads

    Parameters
    ----------
    n_reads : int, number of reads to simulate
    len_read : int, length of the reads
    """
    reads = [random_nt_sequence(len_read) for _ in range(n_reads)]
    return reads



if __name__ == "__main__":
    # Simulate epitopes
    epitopes = simulate_epitopes(10)
    # print(f'Epitopes: {epitopes}')

    # Simulate V and J genes
    v_genes = simulate_vj_genes(10, 75, epitopes)
    j_genes = simulate_vj_genes(10, 75, epitopes)
    # print(f'V genes: {v_genes}')
    # print(f'J genes: {j_genes}')

    # Simulate reads
    overlap_reads = simulate_overlap_reads(10, 75, v_genes, j_genes)
    random_reads = simulate_random_reads(10, 75)
    print(f'Overlap reads: {overlap_reads}')
    # print(f'Random reads: {random_reads}')
