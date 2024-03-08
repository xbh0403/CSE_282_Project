"""
Grid search on a small simulated dataset to find the best hyperparameters for:
    - match_reward
    - mismatch_penalty
    - indel_penalty
    - overlap_match_score
    - overlap_mismatch_score

The best hyperparameters are then used to align the reads and the results are saved to a JSON file.
"""

import json
import numpy as np
import tqdm
from typing import List, Dict
from Read import Read
from Gene import Gene
from JunkReadRecovery import JunkReadRecovery

from ray import tune


def objective(config):
    match_reward = config["match_reward"]
    mismatch_penalty = config["mismatch_penalty"]
    indel_penalty = config["indel_penalty"]
    overlap_match_score = config["overlap_match_score"]
    overlap_mismatch_score = config["overlap_mismatch_score"]
    threshold = config["threshold"]

    with open("/new-stg/home/banghua/Pavel/Simulation/sim_grid_search.json", "r") as f:
        data = json.load(f)
    
    num_pos = len(data["overlap_reads"])
    num_neg = len(data["random_reads"])

    recovered_reads = JunkReadRecovery(match_reward, mismatch_penalty, indel_penalty, overlap_match_score, overlap_mismatch_score, threshold, data, save=False)
    
    num_pred_pos = len(recovered_reads["high_score_overlap"]) + len(recovered_reads["high_score_random"])
    num_pred_neg = num_pos + num_neg - num_pred_pos

    num_true_pos = len(recovered_reads["high_score_overlap"])
    num_false_pos = len(recovered_reads["high_score_random"])
    num_true_neg = num_neg - num_false_pos
    num_false_neg = num_pos - num_true_pos

    precision = num_true_pos / (num_true_pos + num_false_pos)
    recall = num_true_pos / (num_true_pos + num_false_neg)

    f1 = 2 * (precision * recall) / (precision + recall)
    
    return {"score": f1}


search_space = {
    "match_reward": tune.grid_search(list(range(1, 6))),
    "mismatch_penalty": tune.grid_search(list(range(1, 6))),
    "indel_penalty": tune.grid_search(list(range(1, 6))),
    "overlap_match_score": tune.grid_search(list(range(1, 6))),
    "overlap_mismatch_score": tune.grid_search(list(range(1, 6))),
    "threshold": tune.grid_search(list(range(24, 36)))
}

tuner = tune.Tuner(objective, param_space=search_space)

results = tuner.fit()
print(results.get_best_result(metric="score", mode="max"))

df = results.get_dataframe()
df.to_csv("grid_search_results_big.csv", index=False)
