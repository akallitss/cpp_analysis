#!/usr/bin/env python3
# -- coding: utf-8 --
"""
Created on March 23 16:01 2025
Created in PyCharm
Created as cpp_analysis/make_condor_input_list.py

@author: Alexandra Kallitsopoulou, akallits
"""

import os
import re


def main():

    # Generate run-pool list file
    output_file_name = 'processedTrees_run_list.txt'
    myfunctions_h_path = find_file('../', 'MyFunctions', '.h')

    outdirname_path = get_directory_path_from_myfunctions_h(myfunctions_h_path, 'OUTDIRNAME')
    print(outdirname_path)
    run_pool_numbers = extract_run_pool_numbers_from_dir(outdirname_path)
    with open(output_file_name, 'w') as f:
        for run_pool in run_pool_numbers:
            f.write(f'{run_pool[0]},{run_pool[1]}\n')


    # Generate condor bash script
    myfunctions_c_path = find_file('../', 'MyFunctions', '.C')
    analysetreepicosec_c_path = find_file('../', 'AnalyseTreePicosec', '.C')
    codedirname_path = get_directory_path_from_myfunctions_h(myfunctions_h_path, 'CODEDIR')

    myfunctions_c_name = os.path.basename(myfunctions_c_path)
    analysetreepicosec_c_name = os.path.basename(analysetreepicosec_c_path)

    make_bash_script(codedirname_path, myfunctions_c_name, analysetreepicosec_c_name)


    print('bonzo')

def find_file(directory, name_contains, extension):
    for file in os.listdir(directory):
        if file.endswith(extension):
            if name_contains in file:
                return f'{directory}{file}'
    return None

def get_directory_path_from_myfunctions_h(myfunctions_h_path, path_var_name):
    with open(myfunctions_h_path, 'r') as f:
        lines = f.readlines()

    for line in lines:
        if path_var_name in line:
            if '//' not in line:
                path = line.split('"')[1]
                return path
    return None

def extract_run_pool_numbers_from_dir(dir_path):
    run_pool_numbers = []
    for file in os.listdir(dir_path):
        if file.endswith('.root'):
            if 'DEBUG' not in file:
                match = re.search(r'Run(\d+)-Pool(\d+)', file)
                if match:
                    run_number = match.group(1)
                    pool_number = match.group(2)
                    run_pool_numbers.append([run_number, pool_number])
    return run_pool_numbers

def make_bash_script(code_dir_path, myfunctions_c_name, analysetreepicosec_c_name):
    bash_script = f"""#!/bin/bash

cd {code_dir_path} || exit 1  # Ensure script exits if cd fails
rm -r *pcm *.so *d
echo $1
echo $2
root -l -b -q "{myfunctions_c_name}++"
root -l -b -q "{analysetreepicosec_c_name}++($1,$2)"
"""
    with open('run_code.sh', 'w') as f:
        f.write(bash_script)
    os.chmod('run_code.sh', 0o755)

if __name__ == '__main__':
    main()
