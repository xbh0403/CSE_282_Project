import pandas as pd
import numpy as np
import sys
from typing import List, Dict, Iterable, Tuple
sys.setrecursionlimit(100000)


def OverlapVDJAlignment(match_reward: int, mismatch_penalty: int, indel_penalty: int,
                        v_gene: str, j_gene: str, read: str,
                        print_details) -> Tuple[int, str, str]:
    score_v_tail, aligned_v_tail, aligned_read_head = OverlapAlignment(match_reward, mismatch_penalty, indel_penalty, v_gene, read, print_details)
    score_j_head, aligned_read_tail, aligned_j_head = OverlapAlignment(match_reward, mismatch_penalty, indel_penalty, read, j_gene, print_details)


def OverlapAlignment(match_reward: int, mismatch_penalty: int, indel_penalty: int,
                    s: str, t: str,
                    print_details) -> Tuple[int, str, str]:
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