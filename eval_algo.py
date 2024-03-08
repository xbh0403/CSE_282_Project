from concurrent.futures import ProcessPoolExecutor
import concurrent
from tqdm import tqdm

from BuildEpiReadDict import build_epitope_reads_dict
from BruteForce import brute_force_max_coverage
from Greedy import greedy_max_coverage

import json
import time
import os
import sys

def eval_algo(final_json_path, k, save_path):
    with open(final_json_path) as f:
        data = json.load(f)

    total_reads = len(data['high_score_overlap']) + len(data['high_score_random'])
    print("Total Reads:", total_reads)
    epitope_reads_dict = build_epitope_reads_dict(data)
    print("Brute Force")
    start = time.time()
    brute_force_result, brute_force_num_covered = brute_force_max_coverage(epitope_reads_dict, k)
    # Time taken in seconds
    print("Time taken:", time.time() - start)
    t_brute_force_result = time.time() - start
    print("Brute Force Num Covered:", brute_force_num_covered)
    print("Brute Force Result:", brute_force_result)
    print("Greedy")
    start = time.time()
    greedy_result, greedy_num_recovered = greedy_max_coverage(epitope_reads_dict, k)
    # Time taken in seconds
    print("Time taken:", time.time() - start)
    t_greedy_result = time.time() - start
    print("Greedy Num Recovered:", greedy_num_recovered)
    print("Greedy Result:", greedy_result)
    result = {
            'brute_force_result': list(brute_force_result),
            'greedy_result': list(greedy_result),
            'brute_force_num_covered': brute_force_num_covered,
            'greedy_num_recovered': greedy_num_recovered,
            'time_brute_force': t_brute_force_result,
            'time_greedy': t_greedy_result,
            'total_reads': total_reads
        }
    with open(save_path, 'w') as f:
        json.dump(result, f)
    
    return result

def process_file(final_json_path, k, results_dir):
    """
    A wrapper function to call eval_algo with the correct save_path based on the file and k.
    """
    save_path = os.path.join(results_dir, f"{os.path.splitext(final_json_path)[0]}_k{k}_result.json")
    eval_algo(os.path.join('./results/Alignment', final_json_path), k, save_path)

if __name__ == "__main__":
    k_range = range(1, 11)
    all_final_jsons = os.listdir('./results/Alignment')
    all_final_jsons = [f for f in all_final_jsons if "recovered" in f]

    results_dir = './results/Algorithm'
    os.makedirs(results_dir, exist_ok=True)

    # Prepare a list of tasks with all combinations of k and JSON files
    tasks = [(final_json, k) for final_json in all_final_jsons for k in k_range]

    # Use ProcessPoolExecutor to execute the tasks in parallel
    with ProcessPoolExecutor() as executor:
        # Using list comprehension to create and immediately start all tasks
        # tqdm is used to display progress
        futures = [executor.submit(process_file, final_json, k, results_dir) for final_json, k in tasks]
        # Ensure all tasks are complete before moving on
        results = []
        for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures)):
            results.append(future.result())
