import numpy as np
import pandas as pd
import random
import json
from Read import Read
from Gene import Gene
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


def simulate_vj_genes(n_genes: int, len_gene: int, epitope_pool: List, gene_type: str) -> Dict:
    """
    Simulate V and J genes

    Parameters
    ----------
    n_genes : int, number of genes to simulate
    len_gene : int, length of the genes
    epitope_pool : List, epitopes
    """
    genes = []
    for _ in range(n_genes):
        gene = random_nt_sequence(len_gene)
        num_epitopes = random.randint(1, 10)
        gene_epitopes = random.sample(epitope_pool, num_epitopes)
        gene_obj = Gene(gene, gene_type, gene_epitopes)
        genes.append(gene_obj)
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
    reads = []
    while len(reads) < n_reads:
        v_gene = random.choice(v_gene_pool)
        j_gene = random.choice(j_gene_pool)

        v_gene_seq = v_gene.seq
        j_gene_seq = v_gene.seq

        len_v_overlap = random.randint(0, min(len_read, len(v_gene_seq)))
        len_j_overlap = random.randint(0, min(len_read - len_v_overlap, len(j_gene_seq)))
        len_d = len_read - len_v_overlap - len_j_overlap

        if len_d <= 0 or len_d >= 12: continue

        d_seq = random_nt_sequence(len_d)
        # Take the last v_overlap nucleotides from the V gene, d_seq, and the first j_overlap nucleotides from the J gene
        if len_v_overlap == 0:
            read_seq = d_seq + j_gene_seq[:len_j_overlap]
        elif len_j_overlap == 0:
            read_seq = v_gene_seq[-len_v_overlap:] + d_seq
        else:
            read_seq = v_gene_seq[-len_v_overlap:] + d_seq + j_gene_seq[:len_j_overlap]

        read_epitopes = list(set(v_gene.epitopes + j_gene.epitopes))
        read_obj = Read(read_seq, v_gene_seq, d_seq, j_gene_seq, read_epitopes)
        reads.append(read_obj)

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
    reads_obj = [Read(read, '', '', '', []) for read in reads]
    return reads_obj

def simulate(num_epitopes: int, num_v_genes: int, num_j_genes: int, num_reads: int, len_read: int, json: bool = False):
    """
    Simulate data main function

    Parameters
    ----------
    num_epitopes : int, number of epitopes to simulate
    num_v_genes : int, number of V genes to simulate
    num_j_genes : int, number of J genes to simulate
    num_reads : int, number of reads to simulate
    len_read : int, length of the reads
    """
    epitopes = simulate_epitopes(num_epitopes)
    v_genes = simulate_vj_genes(num_v_genes, len_read, epitopes, 'V')
    j_genes = simulate_vj_genes(num_j_genes, len_read, epitopes, 'J')
    overlap_reads = simulate_overlap_reads(num_reads, len_read, v_genes, j_genes)
    random_reads = simulate_random_reads(num_reads, len_read)
    if json:
        return {
            'epitopes': epitopes,
            'v_genes': [v.__json__() for v in v_genes],
            'j_genes': [j.__json__() for j in j_genes],
            'overlap_reads': [r.__json__() for r in overlap_reads],
            'random_reads': [r.__json__() for r in random_reads]
        }
    else:
        return {
            'epitopes': epitopes,
            'v_genes': v_genes,
            'j_genes': j_genes,
            'overlap_reads': overlap_reads,
            'random_reads': random_reads
        }


if __name__ == "__main__":
    num_epitopes = 10
    num_v_genes = 10
    num_j_genes = 10
    num_reads = 100
    len_read = 75
    data = simulate(num_epitopes, num_v_genes, num_j_genes, num_reads, len_read, True)
    with open('./Simulation/sim_1.json', 'w') as f:
        json.dump(data, f, indent=4)
