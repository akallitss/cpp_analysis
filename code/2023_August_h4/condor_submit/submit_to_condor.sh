#!/bin/bash

# Can clean and compile ../ directory and recompile the code, if needed
rm -r *pcm *.so *d *_ACLiC_*
root -l -b -q "../MyFunctions_2023_August.C++"
root -l -b -q "../AnalyseTreePicosec_2023_August.C++"
python3 make_condor_input_list.py
chmod 755 processedTrees_run_list.txt
condor_submit job.submit
