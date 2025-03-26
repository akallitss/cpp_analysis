#!/bin/bash

# Can clean and compile ../ directory and recompile the code, if needed
rm -r *pcm *.so *d *_ACLiC_*
rm err_make_tree/*
rm out_make_tree/*
rm log_make_tree/*
root -l -b -q "../MyFunctions_2023_April.C++"
root -l -b -q "../MakeTreefromRawTreePicosecApril23.C++"
condor_submit job_make_tree.submit
