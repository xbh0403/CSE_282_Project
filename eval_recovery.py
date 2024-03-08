import os
import json
from JunkReadRecovery import JunkReadRecovery
from Gene import Gene
from Read import Read

# {'match_reward': 2, 'mismatch_penalty': 4, 'indel_penalty': 3, 'overlap_match_score': 2, 'overlap_mismatch_score': 3, 'threshold': 24}

all_files = os.listdir("Simulation/algorithm_eval")
all_files = [file for file in all_files if file.endswith(".json")]

for file in all_files:
    with open(f"Simulation/algorithm_eval/{file}", "r") as f:
        data = json.load(f)
    data = {k: [Gene(**d) for d in v] if k in ['v_genes', 'j_genes'] else [Read(**d) for d in v] if k in ['overlap_reads', 'random_reads'] else v for k, v in data.items()}
    file_prefix = file.split(".")[0]
    final_json = JunkReadRecovery(match_reward=2, mismatch_penalty=4, indel_penalty=3,
                                  overlap_match_score=2, overlap_mismatch_score=3, 
                                  threshold=24, data=data, save_path='/new-stg/home/banghua/Pavel/results/'+file_prefix, print_progress=True)
    print(file, "done!")
