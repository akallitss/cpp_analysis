#!/usr/bin/env python3
# -- coding: utf-8 --
"""
Created on March 24 13:53 2025
Created in PyCharm
Created as cpp_analysis/check_condor_run.py

@author: Alexandra Kallitsopoulou, akallits
"""

import os
import re

def main():


    #add txt file that contains the list of files submitted to condor
    submitted_files = find_all_files('.', 'processed', '.txt')


    if len(submitted_files) == 0:
        print('Error: No submitted files found. Exiting...')
        return
    elif len(submitted_files) > 1:
        print('Warning: More than one submitted files found. Using the first one...')
    run_list_file = submitted_files[0]

    #split the submitted files into pairs of Run and Pool numbers
    run_pool_numbers = []
    with open(run_list_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            parts = line.strip().split(',')
            if len(parts) >= 2:
                run, pool = parts[:2]
                run_pool_numbers.append([int(run), int(pool)])
    expected_jobs = len(run_pool_numbers)
    print(f"Extracted {expected_jobs} Run-Pool numbers:")

    # Find all .out files in 'out/' directory
    out_files_condor = find_all_files('out/', 'job', '.out')


    #split the out files into lines and print the lines that contain the following strings
    finished_files = []
    finished_run_pools, unfinished_run_pools = [], []
    for out_file in out_files_condor:
        with open(out_file, 'r') as f:
            parts = out_file.split('_')
            if len(parts) == 3:
                run = int(parts[1].replace('run', ''))
                pool = int(parts[2].replace('pool', '').replace('.out', ''))
                run_pool_list = [run, pool]

            finished = False
            lines = f.readlines()
            for line in lines:
                if 'End of script!' in line:

                    finished_files.append(out_file)
                    finished = True
            if finished:
                finished_run_pools.append(run_pool_list)
            else:
                unfinished_run_pools.append(run_pool_list)

    print(f'Finished Run-Pool Pairs ({len(finished_run_pools)}): {sorted(finished_run_pools)}')
    print(f'Unfinished Run-Pool Pairs ({len(unfinished_run_pools)}): {sorted(unfinished_run_pools)}')

    # Get the run-pool pairs that are not in either finished or unfinished lists
    missing_run_pools = [run_pool for run_pool in run_pool_numbers if run_pool not in finished_run_pools and run_pool not in unfinished_run_pools]
    print(f'Missing Run-Pool Pairs ({len(missing_run_pools)}): {sorted(missing_run_pools)}')

    print()
    print(f'Number of expected jobs: {expected_jobs}')
    print(f'Number of total files in out: {len(out_files_condor)}')
    print(f'Number of finished files in out: {len(finished_files)}')
    print(f'Fraction of jobs missing: {len(missing_run_pools)/expected_jobs * 100:.2f}%')
    print(f'Fraction of finished jobs: {len(finished_run_pools)/expected_jobs * 100:.2f}%')
    print(f'Fraction of finished files in out: {len(finished_files)/len(out_files_condor) * 100:.2f}%')

    print('bonzo')

def find_all_files(directory, name_contains, extension):
    matching_files = [
        os.path.join(directory, file)  # Use os.path.join for proper path handling
        for file in os.listdir(directory)
        if file.endswith(extension) and name_contains in file
    ]
    return matching_files


if __name__ == '__main__':
    main()
