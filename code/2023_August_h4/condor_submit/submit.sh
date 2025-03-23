#!/bin/bash


cd /afs/cern.ch/user/a/akallits/PicoAnalysis/cpp_analysis/code/2023_April_h4 || exit 1  # Ensure script exits if cd fails
rm -r *pcm *.so *d
root -l -b -q 'MyFunctions_2023_April.C++'
root -l -b -q 'AnalyseTreePicosec_2023_April.C++(224,2)'
