#!/bin/bash

python3 make_condor_input_list.py
condor_submit job.submit
