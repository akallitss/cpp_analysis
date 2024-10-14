# cpp_analysis

Picosec Analysis for extracting analysis code parameters in C++.

## Table of Contents
1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Examples](#examples)
5. [Contributing](#contributing)
6. [License](#license)
7. [Contact](#contact)

## Introduction
This analysis code was developed by the CEA/IRFU/DEDIP Group for the TestBeam data analysis of the PICOSEC Micromegas Detector.

## Installation
1. Clone the repository:
    ```bash
    git clone https://github.com/akallitss/cpp_analysis.git
    cd cpp_analysis
    ```
2. Ensure you have the required dependencies (e.g., ROOT).

## Usage
1. Define the directories in `code/MyFunctions.h`:
    ```cpp
    const char *CODEDIR = "/path/to/code";
    const char *BASEDIRNAME = "/path/to/data";

    const char *CODEDIR="/sw/akallits/PicoAnalysis/Saclay_Analysis/code";
    const char *BASEDIRNAME="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2022_October_h4";
    const char *WORKDIR="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2022_October_h4/wdir";
    const char *PLOTDIR="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2022_October_h4/plots";
    const char *DATADIRNAME="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2022_October_h4/dataTrees";
    const char *TRACKDIRNAME="/sw/akallits/PicoAnalysis/TestBeams/2022_October_h4/tracking/";
    const char *OUTDIRNAME="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2022_October_h4/processedTrees";
    const char *PARAMDIRNAME="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2022_October_h4/processedTrees/ParameterTrees";
    const char *DATA_PATH_NAME="/sw/akallits/PicoAnalysis/TestBeams/2022_October_h4";

    ...
    ```
2. Create the directory structure:
    ```bash
    root -l
    .x ConstructDirTree.C
    ```
3. Compile and run the first analysis part code:
    ```bash
    cd Bin2Tree
    make
    ./bin2tree runNo poolNo
    ```
This will convert from the binary input data an output Tree with the data or each channel of the oscilloscope

## Examples
- Example command to generate root file:
    ```bash
    ./bin2tree 302 1
    ```
Now having the data in ROOT format, you can start the analysis by matching the tracking data information ( already in a reconstructed data file) 
Make a file named OscilloscopeSetUp.txt, where it contains the information about the run, which detectors were connected to which channel, 
and their description (photocathode, anode/cathode voltages, position on the tracking, amplifier connected) 

- Example command to process raw data:
    ```bash
    root -l
    .x MakeTreefromRawTreePicosec.C++(302, 1)
    ```
At this point the output is again a ROOT file that matches the desired tracking data on the specific detector, as well as the event number from the tracking
and after decoding the 16bit waveform from the reference oscilloscope channel
This will make the output file /path_to/data/processedTrees/ 
- Example command to process the data with the right tracking info
    ```bash
    root -l
    .x AnalyseTreePicosec.C++(302, 1)
Here we follow the offline analysis procedure based on the CFD method
and we generate the final important parameters for the analysis, 

The results are saved in the /path_to/data/processedTrees/ParameterTrees

*TO BE CONTINUED* 

## Contact
For questions or support, please contact [alexandra.kallitsopoulou@cea.fr](alexandra.kallitsopoulou@cea.fr)
