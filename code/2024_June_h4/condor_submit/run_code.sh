#!/bin/bash

cd /afs/cern.ch/user/a/akallits/PicoAnalysis/cpp_analysis/code/2023_June_h4 || exit 1  # Ensure script exits if cd fails
root -l -b -q "AnalyseTreePicosec_2023_June.C+($1,$2)"
