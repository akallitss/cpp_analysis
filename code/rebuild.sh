rm -r *pcm *.so *d *dict
root -l -b -q MyFunctions.C+
root -l 'AnalyseTreePicosec.C++(224, 2, 1)'