//
// Created by dylan on 10/20/24.
//

#include <TFile.h>
#include <TTree.h>
#include <TRandom.h>
#include "WaveformFit.h"
#include <vector>
#include <TInterpreter.h>
#include <TROOT.h>

ClassImp(WaveformFit)

void write_waveform_fits() {
    gROOT->ProcessLine(".L WaveformFit.h+");  // Load the header and compile with ROOT

    // Create the ROOT file and TTree
    TFile *file = new TFile("waveform_fits.root", "RECREATE");
    TTree *tree = new TTree("eventTree", "Tree of events with waveform fits");

    const int maxFits = 10;  // Max number of fits we expect per event
    WaveformFit fits[maxFits];  // Fixed-size array for storing WaveformFit objects
    int nFits;  // To store the number of fits in the current event

    // Set up branches
    tree->Branch("nFits", &nFits, "nFits/I");
    tree->Branch("fits", fits, "amplitude[nFits]/F:time[nFits]/F:chi2[nFits]/F");

    TRandom randGen;
    for (int i = 0; i < 100; ++i) {  // Generate 100 events
        nFits = randGen.Poisson(2);  // Poisson-distributed number of fits (1-3 fits per event)

        for (int j = 0; j < nFits; ++j) {
            fits[j].amplitude = randGen.Gaus(10, 2);  // Mock fit values
            fits[j].time = randGen.Gaus(5, 1);
            fits[j].chi2 = randGen.Uniform(0, 1);
        }

        tree->Fill();
    }

    tree->Write();
    file->Close();
}
