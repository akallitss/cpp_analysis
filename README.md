# cpp_analysis
Picosec Analysis for extracting analysis code parameters CPP

==============================================================
README FILE
This analysis code was developped from CEA/IRFU//DEDIP Group
for the TestBeam data analysis of PICOSEC Micromegas Detector

The main code is contained in the code/ folder,
while on the data/ you will find, organised by TestBeam period,
the data used and produced in the analysis procedure.


==============================================================
HOW TO RUN
----In the /paht_to/code/directory you will find the header MyFunctions.h.
------MyFunctions.h-----------------
First you need to define the directories of the raw-binary data, tracking data
and all the other directories you may need for your analysis.
In this version we use the following structure:


const char *CODEDIR="/sw/akallits/PicoAnalysis/Saclay_Analysis/code";
const char *BASEDIRNAME="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2022_October_h4";
const char *WORKDIR="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2022_October_h4/wdir";
const char *PLOTDIR="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2022_October_h4/plots";
const char *DATADIRNAME="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2022_October_h4/dataTrees";
const char *TRACKDIRNAME="/sw/akallits/PicoAnalysis/TestBeams/2022_October_h4/tracking/";
const char *OUTDIRNAME="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2022_October_h4/processedTrees";
const char *PARAMDIRNAME="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2022_October_h4/processedTrees/ParameterTrees";
const char *DATA_PATH_NAME="/sw/akallits/PicoAnalysis/TestBeams/2022_October_h4";


All those directories will be created by simply running:

    % root -l
    % .x ConstructDirTree.C
----Having your environment set up, you need to move to the Bin2Tree/ folder


bin2tree.cxx is tha main code which together with the makefile
will read the LEcroy software and binary files

    % make
    % ./bin2tree runNo poolNo


!!!!!!!!!!!!!!!!!!!!!!!!ATTENTION for GDD scope you use poolNo == 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
This will result to a root file containing all trc waveforms in a tree format.
You can find this in the directory
/path_to/data/2022_October_h4/dataTrees/
with the name RUNNO-POOLNO_TESTBEAMraw_tree.root
this contatins all the waveform points for each channel and all the oscilloscope information during data taking
i.e, gain, offset, npoints, dt, etc.

----For the next step you will start to process the tree with the RAW data.
First you have to create a logbook txt file (see the example) with RUN parameter information,
named OsciloscopeSetup.txt in the data/TESTBEAMPERIODDATA/ directory
In the code/ directory, you will find the codes: makeTree.cxx and MakeTreefromRawTreePicosec.C
In the MakeTreefromRawTreePicosec.C you read both the RUNNO-POOLNO_TESTBEAMraw_tree.root raw dataTree file
and the OsciloscopeSetup information of the individual run.
Those txt run informations will then be stored to the output root file of the data.
At this point the tracking information is added and the SRS decoding is processed.

    % root -l  
    % .x MakeTreefromRawTreePicosec.C++(RUNNO,POOLNO)
This will result to the outputfile /path_to/data/2022_October_h4/processedTrees/Run224-Pool2_TESTBEAM_tree.root
Which contains to trees :
OsciloscopeSetup with all the information read from the txt file for the run parameters
RawDataTree with all waveforms per channel, the tracking information and the SRS number from decoding

----Last step is to follow the offline analysis procedure based on the CFD method
This will use the AnalyseTreePicosec.C and.....to be continued
