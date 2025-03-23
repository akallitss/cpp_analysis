#!/bin/bash

python3 make_condor_input_list.py
chmod 755 processedTrees_run_list.txt
condor_submit job.submit
