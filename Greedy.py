from typing import Dict, Set
import json
from JunkReadRecovery import JunkReadRecovery
from BuildEpiReadDict import build_epitope_reads_dict


def greedy_max_coverage(epitope_reads_dict: Dict, k: int) -> Set:
    """
    Greedy algorithm to find the set of k epitopes that has the maximum coverage of reads.

    Parameters
    ----------
    epitope_reads_dict : Dict
        A dictionary with epitopes as keys and a list of reads as values.
    k : int
        The number of epitopes to select.

    Returns
    -------
    Set
        The set of selected epitopes that maximizes the coverage of reads.
    """
    selected_epitopes = set()
    uncovered_reads = set(read for reads in epitope_reads_dict.values() for read in reads)
    for _ in range(k):
        if not uncovered_reads:
            break
        best_epitope = None
        best_coverage = 0
        for epitope, reads in epitope_reads_dict.items():
            coverage = len(uncovered_reads.intersection(reads))
            if coverage > best_coverage:
                best_coverage = coverage
                best_epitope = epitope
        if best_epitope:
            selected_epitopes.add(best_epitope)
            uncovered_reads -= set(epitope_reads_dict[best_epitope])
    return selected_epitopes

# Example usage
if __name__ == "__main__":
    with open('./Simulation/sim_3_7.json') as f:
        data = json.load(f)

    match_reward, mismatch_penalty, indel_penalty = 1, 1, 1
    overlap_match_score, overlap_mismatch_score = 1, 1
    threshold = 25
    final_json = JunkReadRecovery(match_reward, mismatch_penalty, indel_penalty, overlap_match_score, overlap_mismatch_score, threshold, data, save=True, print_progress=True)
    
    epitope_reads_dict = build_epitope_reads_dict(final_json)
    k = 2
    print(greedy_max_coverage(epitope_reads_dict, k))