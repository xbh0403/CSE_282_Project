import numpy as np
import sys
from typing import List, Dict, Iterable, Tuple
import json
sys.setrecursionlimit(100000)
import tqdm

from Read import Read
from Gene import Gene
from OverlapAlignment import OverlapVDJAlignment

import concurrent
from concurrent.futures import ProcessPoolExecutor

def JunkReadRecovery(match_reward: int, mismatch_penalty: int, indel_penalty: int,
                     overlap_match_score: int, overlap_mismatch_score: int, threshold: int,
                     data: Dict, save_path: str|None = None, print_progress: bool = False) -> Dict:
    """
    Perform overlap alignment between V, D, and J genes and reads

    Parameters
    ----------
    match_reward : int, reward for matching nucleotides
    mismatch_penalty : int, penalty for mismatching nucleotides
    indel_penalty : int, penalty for indels
    overlap_match_score : int, reward for matching nucleotides in the overlap region
    overlap_mismatch_score : int, penalty for mismatching nucleotides in the overlap region
    threshold : int, threshold for the score
    data : Dict, data
    save_path : str, path to save the results
    print_progress : bool, print progress
    """
    final_json = AlignAllReads(match_reward, mismatch_penalty, indel_penalty, overlap_match_score, overlap_mismatch_score, data, print_progress=print_progress)
    final_json_filtered = KeepHighScoreAlignments(final_json, threshold)
    if save_path is not None:
        final_json_filtered_save = {
            'high_score_overlap': final_json_filtered['high_score_overlap'],
            'high_score_random': final_json_filtered['high_score_random'],
            'all_epitopes': final_json_filtered['all_epitopes'],
            'all_v_genes': [i.__json__() for i in final_json_filtered['all_v_genes']],
            'all_j_genes': [i.__json__() for i in final_json_filtered['all_j_genes']]
        }
        with open(save_path+'_recovered.json', 'w') as f:
            json.dump(final_json_filtered_save, f)

        final_json_save = {
            'overlap': final_json['overlap'],
            'random': final_json['random'],
            'all_epitopes': final_json['all_epitopes'],
            'all_v_genes': [i.__json__() for i in final_json['all_v_genes']],
            'all_j_genes': [i.__json__() for i in final_json['all_j_genes']]
        }
        with open(save_path+'_all.json', 'w') as f:
            json.dump(final_json_save, f)
            
    return final_json_filtered


def KeepHighScoreAlignments(data: Dict, threshold: int) -> Dict:
    """
    Keep reads with high scores

    Parameters
    ----------
    data : Dict, data
    threshold : int, threshold for the score
    """
    overlap_results, random_results = data['overlap'], data['random']
    all_epitopes, all_v_genes, all_j_genes = data['all_epitopes'], data['all_v_genes'], data['all_j_genes']

    high_score_overlap = [result for result in overlap_results if result['final_score'] > threshold]
    high_score_random = [result for result in random_results if result['final_score'] > threshold]

    final_json = {
        'high_score_overlap': high_score_overlap,
        'high_score_random': high_score_random,
        'all_epitopes': all_epitopes,
        'all_v_genes': all_v_genes,
        'all_j_genes': all_j_genes
    }

    return final_json

def AlignAllReads(match_reward: int, mismatch_penalty: int, indel_penalty: int,
                  overlap_match_score: int, overlap_mismatch_score: int,
                  data: dict, print_progress: bool) -> list:
    """
    Perform overlap alignment between V, D, and J genes and reads
    """
    # Unpack data
    all_epitopes, all_v_genes, all_j_genes, overlap_reads, random_reads = \
        data['epitopes'], data['v_genes'], data['j_genes'], data['overlap_reads'], data['random_reads']
    
    # Convert dictionaries to Gene and Read objects if necessary
    all_v_genes = [Gene(**d) if type(d) != Gene else d for d in all_v_genes]
    all_j_genes = [Gene(**d) if type(d) != Gene else d for d in all_j_genes]
    overlap_reads = [Read(**d) if type(d) != Read else d for d in overlap_reads]
    random_reads = [Read(**d) if type(d) != Read else d for d in random_reads]
    
    def align_reads(reads):
        with ProcessPoolExecutor() as executor:
            futures = [executor.submit(AlignOneRead, match_reward, mismatch_penalty, indel_penalty,
                                       overlap_match_score, overlap_mismatch_score, all_v_genes, all_j_genes, read)
                       for read in reads]
            results = []
            for future in tqdm.tqdm(concurrent.futures.as_completed(futures), total=len(futures), disable=not print_progress):
                results.append(future.result())
        return results

    print('Aligning overlap reads')
    results_overlap = align_reads(overlap_reads)

    print('\nAligning random reads')
    results_random = align_reads(random_reads)

    final_json = {
        'overlap': results_overlap,
        'random': results_random,
        'all_epitopes': all_epitopes,
        'all_v_genes': all_v_genes,
        'all_j_genes': all_j_genes
    }

    return final_json


def AlignOneRead(match_reward: int, mismatch_penalty: int, indel_penalty: int,
                 overlap_match_score: int, overlap_mismatch_score: int,
                v_genes: Gene, j_genes: Gene, read: Read) -> Dict:
    """
    Perform overlap alignment between V, D, and J genes and a read

    Parameters
    ----------
    match_reward : int, reward for matching nucleotides
    mismatch_penalty : int, penalty for mismatching nucleotides
    indel_penalty : int, penalty for indels
    overlap_match_score : int, reward for matching nucleotides in the overlap region
    overlap_mismatch_score : int, penalty for mismatching nucleotides in the overlap region
    v_genes : Gene, V genes
    j_genes : Gene, J genes
    read : Read, read
    """
    best_score = -np.inf
    best_result = None
    for v_gene in v_genes:
        for j_gene in j_genes:
            result = OverlapVDJAlignment(match_reward, mismatch_penalty, indel_penalty, overlap_match_score,
                                        overlap_mismatch_score, v_gene, j_gene, read, False)
            if result['final_score'] > best_score:
                best_score = result['final_score']
                best_result = result
    return best_result

if __name__ == "__main__":
    with open('./Simulation/sim_10.json') as f:
        data = json.load(f)
        data = {k: [Gene(**d) for d in v] if k in ['v_genes', 'j_genes'] else [Read(**d) for d in v] if k in ['overlap_reads', 'random_reads'] else v for k, v in data.items()}
    
    match_reward, mismatch_penalty, indel_penalty = 1, 1, 1
    overlap_match_score, overlap_mismatch_score = 1, 1
    threshold = 25
    final_json = JunkReadRecovery(match_reward, mismatch_penalty, indel_penalty, overlap_match_score, overlap_mismatch_score, threshold, data, True, True)

