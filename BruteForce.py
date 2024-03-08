from itertools import combinations
from typing import Dict, Set
import json
from JunkReadRecovery import JunkReadRecovery
from BuildEpiReadDict import build_epitope_reads_dict

def brute_force_max_coverage(epitope_reads_dict: Dict, k: int) -> Set:
    """
    Brute force algorithm to find the set of k epitopes that has the maximum coverage of reads.

    Parameters
    ----------
    epitope_reads_dict : Dict
        A dictionary with epitopes as keys and a list of reads as values (epitopes).
    k : int
        The number of epitopes to select.
    
    Returns
    -------
    Set
        The set of k epitopes that maximizes the coverage.
    """
    max_coverage = 0
    best_combination = set()

    # Generate all possible combinations of k epitopes
    for combination in combinations(epitope_reads_dict.keys(), k):
        # Find the union of reads for the current combination of epitopes
        combined_reads = set()
        for epitope in combination:
            combined_reads.update(epitope_reads_dict[epitope])
        
        # Update max_coverage and best_combination if this combination is better
        if len(combined_reads) > max_coverage:
            max_coverage = len(combined_reads)
            best_combination = set(combination)
    
    return best_combination


if __name__ == "__main__":
    with open('./Simulation/sim_3_7.json') as f:
        data = json.load(f)

    match_reward, mismatch_penalty, indel_penalty = 1, 1, 1
    overlap_match_score, overlap_mismatch_score = 1, 1
    threshold = 25
    final_json = JunkReadRecovery(match_reward, mismatch_penalty, indel_penalty, overlap_match_score, overlap_mismatch_score, threshold, data, save=True, print_progress=True)
    
    epitope_reads_dict = build_epitope_reads_dict(final_json)
    k = 2
    print(brute_force_max_coverage(epitope_reads_dict, k))
