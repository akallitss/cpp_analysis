//
// Created by ak271430 on 27/11/24.
//

#include <TROOT.h>
#include <TSystem.h>
#include <TMath.h>
#include <TH1.h>
#include <TTree.h>
#include <TFile.h>
 #include <TH1F.h>
#include <TF1.h>
#include <TGaxis.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TLatex.h>
#include <TString.h>
#include <TTimeStamp.h>
#include <TExec.h>
#include <TLine.h>
#include <TObjString.h>
#include<TMultiGraph.h>
#include<TClonesArray.h>
#include "MyFunctions.C"

#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>

void testAnalyseLongPulseCiv() {

//    gROOT->LoadMacro("MyFunctions.h"); // Load your macro
    gROOT->LoadMacro("MyFunctions.C"); // Load your macro
  // Parameters
    const int points = 10002; // Number of points in the waveform
//    const double dt = 0.1;   // Time step in microseconds
    const double threshold = -0.05; // Threshold for the trigger
    const double sig_shift = 2.0;   // Signal shift for analysis
    const int tshift = 0;           // Time shift
    double dt = 0.0;
    double rms = 0.0;
    double bsl = 0.0;
    double data[points];
    double drv[points];

    gStyle->SetLabelSize(0.045,"X");
    gStyle->SetLabelSize(0.045,"Y");
    gStyle->SetLabelFont(132,"X");
    gStyle->SetLabelFont(132,"Y");
    gStyle->SetTitleSize(0.045,"X");
    gStyle->SetTitleSize(0.045,"Y");
    gStyle->SetTitleFont(22,"X");
    gStyle->SetTitleFont(22,"Y");
    gStyle->SetOptStat(1001111);
    gStyle->SetNdivisions(507);
    gStyle->SetMarkerSize(100);


    // Open the file for reading
//    std::ifstream inFile("/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2022_October_h4/plots/Run224/Pool2/moreplots/waveform_data.txt");

    std::ifstream inFile("/home/akallits/Documents/PicoAnalysis/Saclay_Analysis/data/2022_October_h4/plots/Run224/Pool2/moreplots/waveform_data.txt");

    if (!inFile) {
        std::cerr << "Error: Could not open file for reading!" << std::endl;
        return;
    }

    // Parse the file line by line
    std::string line;
    while (std::getline(inFile, line)) {
        if (line.find("dt:") != std::string::npos) {
            // Parse the dt value
            std::istringstream iss(line.substr(line.find(":") + 1));
            iss >> dt;
        } else if (line.find("RMS: ") != std::string::npos) {
            // Parse the RMS value
            std::istringstream iss(line.substr(line.find(":") + 1));
            iss >> rms;
            std::cout << "RMS: " << rms << std::endl;
        } else if (line.find("BSL: ") != std::string::npos) {
            // Parse the BSL value
            std::istringstream iss(line.substr(line.find(":") + 1));
            iss >> bsl;
            std::cout << "BSL: " << bsl << std::endl;
        } else if (line.find("double data") != std::string::npos) {
            // Parse the data array
            line = line.substr(line.find("{") + 1);
            line.pop_back(); // Remove the closing brace
            std::istringstream iss(line);
            for (int i = 0; i < points && iss; ++i) {
                std::string value;
                if (std::getline(iss, value, ',')) {
                    data[i] = std::stod(value);
                }
            }
        } else if (line.find("double drv") != std::string::npos) {
            // Parse the drv array
            line = line.substr(line.find("{") + 1);
            line.pop_back(); // Remove the closing brace
            std::istringstream iss(line);
            for (int i = 0; i < points && iss; ++i) {
                std::string value;
                if (std::getline(iss, value, ',')) {
                    drv[i] = std::stod(value);
                }
            }
        }
    }

    // Close the file
    inFile.close();

    // Verify the data (optional)
    std::cout << "Read dt = " << dt << std::endl;
    std::cout << "First few data points: ";
    for (int i = 0; i < 5; ++i) {
        std::cout << data[i] << " ";
    }
    std::cout << "\nFirst few drv points: ";
    for (int i = 0; i < 5; ++i) {
        std::cout << drv[i] << " ";
    }
    std::cout << std::endl;


    // Prepare PEAKPARAM structure
    PEAKPARAM par;
    par.Reset();
    par.rms = rms;
    par.bsl = bsl;

    // Call the AnalyseLongPulseCiv function
    int evNo = 27; // Example event number
    int result = AnalyseLongPulseCiv(points, evNo, data, dt, drv, &par, threshold, sig_shift, tshift);

    // Check result
    if (result < 0) {
        std::cout << "Analysis failed with result: " << result << std::endl;
        return;
    }
    std::cout << "Analysis succeeded, result: " << result << std::endl;

}
