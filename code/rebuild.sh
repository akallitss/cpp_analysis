rm -r *pcm *.so *d
root -l -b -q MyFunctions.C+
root -l 'AnalyseTreePicosec.C++(224, 2, 1)'