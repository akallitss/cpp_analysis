//
// Created by dylan on 10/20/24.
//


#include <TFile.h>
#include <TTree.h>
#include <vector>
#include "WaveformFit.h"

void read_waveform_fits() {
    TFile *file = TFile::Open("waveform_fits.root");
    TTree *tree = (TTree*)file->Get("eventTree");

    std::vector<WaveformFit> *fits = nullptr;
    tree->SetBranchAddress("fits", &fits);

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        std::cout << "Event " << i << " has " << fits->size() << " fits." << std::endl;

        for (const auto &fit : *fits) {
            std::cout << "  Fit: amplitude = " << fit.amplitude
                      << ", time = " << fit.time
                      << ", chi2 = " << fit.chi2 << std::endl;
        }
    }

    file->Close();
}
