import pandas as pd
import numpy as np
import sys
from typing import List, Dict, Iterable, Tuple
sys.setrecursionlimit(100000)


def OverlapVDJAlignment(match_reward: int, mismatch_penalty: int, indel_penalty: int,
                        overlap_match_score: int, overlap_mismatch_score: int,
                        v_gene: Dict, j_gene: Dict, read: Dict, 
                        print_details) -> Tuple[int, str, str]:
    """
    Perform overlap alignment between V, D, and J genes and a read

    Parameters
    ----------
    match_reward : int, reward for matching nucleotides
    mismatch_penalty : int, penalty for mismatching nucleotides
    indel_penalty : int, penalty for indels
    overlap_match_score : int, reward for matching nucleotides in the overlap region
    overlap_mismatch_score : int, penalty for mismatching nucleotides in the overlap region
    v_gene : Dict, V gene
    j_gene : Dict, J gene
    read : Dict, read
    print_details : bool, print details of the alignment

    Returns
    -------
    Tuple[int, str, str, List, str, str, List], score, aligned V tail, aligned read head, aligned read tail, aligned J head
    """
    v_gene, j_gene, read = v_gene['gene'], j_gene['gene'], read['read']
    v_gene_epi, j_gene_epi, read_epi = v_gene['epitopes'], j_gene['epitopes'], read['epitopes']
    score_v_tail, aligned_v_tail, aligned_read_head = OverlapAlignment(match_reward, mismatch_penalty, indel_penalty, v_gene, read, print_details)
    score_j_head, aligned_read_tail, aligned_j_head = OverlapAlignment(match_reward, mismatch_penalty, indel_penalty, read, j_gene, print_details)
    score_overlap_v, score_overlap_j = 0, 0
    if len(aligned_v_tail) > 0 and len(aligned_read_head) > 0:
        score_overlap_v = sum([overlap_match_score if aligned_v_tail[i] == aligned_read_head[i] else overlap_mismatch_score for i in range(len(aligned_v_tail))])
    if len(aligned_read_tail) > 0 and len(aligned_j_head) > 0:
        score_overlap_j = sum([overlap_match_score if aligned_read_tail[i] == aligned_j_head[i] else overlap_mismatch_score for i in range(len(aligned_read_tail))])
    final_score = score_overlap_v + score_overlap_j
    return final_score, aligned_v_tail, aligned_read_head, v_gene_epi, aligned_read_tail, aligned_j_head, j_gene_epi


def OverlapAlignment(match_reward: int, mismatch_penalty: int, indel_penalty: int,
                    s: str, t: str,
                    print_details) -> Tuple[int, str, str]:
    """
    Perform overlap alignment between two sequences
    
    Parameters
    ----------
    match_reward : int, reward for matching nucleotides
    mismatch_penalty : int, penalty for mismatching nucleotides
    indel_penalty : int, penalty for indels
    s : str, first sequence
    t : str, second sequence
    print_details : bool, print details of the alignment
    
    Returns
    -------
    Tuple[int, str, str], score, aligned s, aligned t
    """
    score_matrix, backtrack_matrix = [[0 for j in range(len(t)+1)] for i in range(len(s)+1)], [[0 for j in range(len(t)+1)] for i in range(len(s)+1)]
    max_score, max_i, max_j = float('-inf'), -1, -1
    for i in range(len(s)+1):
        score_matrix[i][0] = 0
    for j in range(len(t)+1):
        score_matrix[0][j] = -j * indel_penalty

    for i in range(1, len(s)+1):
        for j in range(1, len(t)+1):
            score_matrix[i][j] = max(score_matrix[i-1][j] - indel_penalty, 
                                    score_matrix[i][j-1] - indel_penalty, 
                                    score_matrix[i-1][j-1] + (match_reward if s[i-1] == t[j-1] else -mismatch_penalty))
            if score_matrix[i][j] > max_score and i == len(s):
                max_score = score_matrix[i][j]
                max_i, max_j = i, j
            if score_matrix[i][j] == score_matrix[i-1][j-1] + (match_reward if s[i-1] == t[j-1] else -mismatch_penalty):
                backtrack_matrix[i][j] = 3
            elif score_matrix[i][j] == score_matrix[i][j-1] - indel_penalty:
                backtrack_matrix[i][j] = 2
            elif score_matrix[i][j] == score_matrix[i-1][j] - indel_penalty:
                backtrack_matrix[i][j] = 1
    if print_details:
        def print_matrix(matrix):
            for row in matrix:
                print(" ".join(map(str, row)))
        print_matrix(score_matrix)
        print()
        print_matrix(backtrack_matrix)

    s_out, t_out = "", ""
    i, j = max_i, max_j
    while j > 0:
        if backtrack_matrix[i][j] == 1:
            s_out = s[i-1] + s_out
            t_out = "-" + t_out
            i -= 1
        elif backtrack_matrix[i][j] == 2:
            s_out = "-" + s_out
            t_out = t[j-1] + t_out
            j -= 1
        elif backtrack_matrix[i][j] == 3:
            s_out = s[i-1] + s_out
            t_out = t[j-1] + t_out
            i -= 1
            j -= 1
        elif backtrack_matrix[i][j] == 0:
            s_out = "-" + s_out
            t_out = t[j-1] + t_out
            j -= 1
    return max_score, s_out, t_out


if __name__ == "__main__":
    match_reward, mismatch_penalty, indel_penalty = 1, 1, 1
    v_gene = "ZZZZZGCAT"
    read = "ATNNNGC"
    j_gene = "GCATYYYYY"
    overlap_match_score, overlap_mismatch_score = 3, 1
    print_details = False
    print(OverlapVDJAlignment(match_reward, mismatch_penalty, indel_penalty, 
                              overlap_match_score, overlap_mismatch_score, 
                              v_gene, j_gene, read, 
                              print_details))
