#!/bin/bash

# load an environment (see https://lcginfo.cern.ch), ideally the same as in .bash_profile
#source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_105  x86_64-el9-gcc13-opt 
export LD_LIBRARY_PATH=/usr/lib64/root/libCore.so.6.30:$LD_LIBRARY_PATH

# change to working directory, can be EOS or AFS
cd /eos/project-p/picosec/analysis/Saclay/code/Bin2Tree/
make 
## execute commands. Change with your commands
./bin2tree 302 2