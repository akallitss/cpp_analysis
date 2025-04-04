#!/bin/bash

# Can clean and compile ../ directory and recompile the code, if needed
rm -r *pcm *.so *d *_ACLiC_*
rm err/*
rm out/*
rm log/*
#root -l -b -q "../MyFunctions_2023_April.C++" #not needed when transfering files to condor node
#root -l -b -q "../AnalyseTreePicosec_2023_April.C++" #not needed when transfering files to condor node
#python3 make_condor_input_list.py
#chmod 755 processedTrees_run_list.txt
condor_submit job.submit
