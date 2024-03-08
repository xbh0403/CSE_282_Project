import json
from Gene import Gene
from Read import Read
from JunkReadRecovery import JunkReadRecovery
from collections import defaultdict


def build_epitope_reads_dict(final_json):
    epitope_reads_dict = defaultdict(list)
    read_overlap = final_json['high_score_overlap']
    read_random = final_json['high_score_random']
    for read in read_overlap:
        for epitope in list(set(read['v_gene_epi'] + read['j_gene_epi'])):
            epitope_reads_dict[epitope].append(read['read_id'])
    for read in read_random:
        for epitope in list(set(read['v_gene_epi'] + read['j_gene_epi'])):
            epitope_reads_dict[epitope].append(-read['read_id'])
    return epitope_reads_dict


if __name__ == "__main__":
    with open('./Simulation/sim_3_7.json') as f:
        data = json.load(f)

    match_reward, mismatch_penalty, indel_penalty = 1, 1, 1
    overlap_match_score, overlap_mismatch_score = 1, 1
    threshold = 25
    final_json = JunkReadRecovery(match_reward, mismatch_penalty, indel_penalty, overlap_match_score, overlap_mismatch_score, threshold, data, print_progress=True)
