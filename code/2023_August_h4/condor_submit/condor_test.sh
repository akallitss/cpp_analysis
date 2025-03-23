#!/bin/bash


cd /afs/cern.ch/user/a/akallits/PicoAnalysis/cpp_analysis/code || exit 1  # Ensure script exits if cd fails
echo "Current directory after cd: $(pwd)"
echo "ROOT is located at: $(which root)"
root --version
root -q -e '.! ls'
root -q -e '.! pwd'
root -l -b -q 'AnalyseTreePicosec.C++(224,2)'
