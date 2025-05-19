
==============================================================
		README FILE
==============================================================

This analysis code was developped from CEA/IRFU//DEDIP Group 
for the TestBeam data analysis of PICOSEC Micromegas Detector

-------------------------------------------------------------

The main code is contained in the code/ folder, 
while on the data/ you will find, organised by TestBeam period, 
the data used and produced in the analysis procedure. 

-------------------------------------------------------------

==============================================================
		HOW TO RUN 
==============================================================


----In the /paht_to/code/directory you will find the header MyFunctions.h. 

------MyFunctions.h-----------------

First you need to define the directories of the raw-binary data, tracking data
and all the other directories you may need for your analysis. 
In this version we use the following structure:  

const char *CODEDIR="/path/to/code";
const char *BASEDIRNAME="/path/to/data/TestBeamPeriod";
const char *WORKDIR="/path/to/data/TestBeamPeriod/wdir";
const char *PLOTDIR="/path/to/data/TestBeamPeriod/plots";
const char *DATADIRNAME="/path/to/data/TestBeamPeriod/dataTrees";
const char *TRACKDIRNAME="/path/to/TestBeamPeriod/tracking/";
const char *OUTDIRNAME="/path/to/data/TestBeamPeriod/processedTrees";
const char *PARAMDIRNAME="/path/topath/to/data/TestBeamPeriod/processedTrees/ParameterTrees";
const char *DATA_PATH_NAME="/path/to/TestBeamPeriod";


All those directories will be created by simply running: 

        % root -l
	% .x ConstructDirTree.C

----Having your environment set up, you need to move to the Bin2Tree/ folder 

- bin2tree.cxx is tha main code which together with the makefile
will read the LEcroy software and binary files 

	 % make 
	 % ./bin2tree runNo poolNo 

!!!!!!!!!!!!!!!!!!!!!!!!ATTENTION for GDD scope you use poolNo == 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

This will result to a root file containing all trc waveforms in a tree format. 
You can find this in the directory 

/path_to/data/TestBeamPeriod/dataTrees/

with the name RUNNO-POOLNO_TESTBEAMraw_tree.root 
this contatins all the waveform points for each channel and all the oscilloscope information during data taking
i.e, gain, offset, npoints, dt, etc. 

--------------------------------------------------------------


----For the next step you will start to process the tree with the RAW data.

First you have to create a logbook txt file (see the example) with RUN parameter information, 
named OsciloscopeSetup.txt in the data/TestBeamPeriodData/ directory 
 
In the code/ directory, you will find the codes: makeTree.cxx and MakeTreefromRawTreePicosec.C

In the MakeTreefromRawTreePicosec.C you read both the RUNNO-POOLNO_TESTBEAMraw_tree.root raw dataTree file
and the OsciloscopeSetup information of the individual run.

Those txt run informations will then be stored to the output root file of the data.

At this point the tracking information is added and the SRS decoding is processed. 

	% root -l  
	% .x MakeTreefromRawTreePicosec.C++(RUNNO,POOLNO)

This will result to the outputfile /path_to/data/TestBeamPeriod/processedTrees/RunNo-PoolNo_TESTBEAM_tree.root
Which contains to trees : 
	OsciloscopeSetup with all the information read from the txt file for the run parameters
	RawDataTree with all waveforms per channel, the tracking information and the SRS number from decoding

------------------------------------------------------------

----Last step is to follow the offline analysis procedure based on the CFD method 

This will use the AnalyseTreePicosec.C by running : 


	% root -l  
	% .x AnalyseTreePicosec.C++(RUNNO,POOLNO) % optionally you can add an extra argument 1 to have event by event display monitoring

 In this part of the code we are using the TestBeam_tree as input and we follow the Constant Fraction Discrimination Method for the timing analysis of the PICOSEC and reference signals. 

 Tmining analysis, amplitude, charge and all the important parameters of the signal processing are saved in a Class called PEAKPARAM 

 All the information is saved per channel (maximum 4) as a separate brance in the outputfile /path_to/data/TestBeamPeriod/processedTrees/RunNo-PoolNo_treeParam.root
