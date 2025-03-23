#!/bin/bash

cd /afs/cern.ch/user/a/akallits/PicoAnalysis/cpp_analysis/code/2024_June_h4 || exit 1  # Ensure script exits if cd fails
rm -r *pcm *.so *d
root -l -b -q 'MyFunctions_2024_June.C++'
root -l -b -q 'AnalyseTreePicosec_2024_June.C++($1,$2)'
