//
// Created by Dylan on 3/15/2025.
//

#include "rms_baseline_calculator.h"

#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <vector>


int rms_baseline_calculator() {
    std::string base_path = "C:/Users/Dylan/Desktop/picosec/";
    std::string file_name = "Run224-Pool2_TESTBEAM_tree.root";
    std::string file_path = base_path + file_name;
    std::string tree_name = "RawDataTree";

    std::cout << "Opening file: " << file_path << std::endl;

    // Open the ROOT file
    TFile *file = TFile::Open(file_path.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file " << file_path << std::endl;
        return 1;
    }

    // Get the TTree
    TTree *tree = (TTree*) file->Get(tree_name.c_str());
    if (!tree) {
        std::cerr << "Error: Could not find tree " << tree_name << " in file " << file_path << std::endl;
        file->Close();
        return 1;
    }

    // Variables to store branch data
    ULong64_t epoch;
    Double_t ampl;

    // Set branch addresses
    tree->SetBranchAddress("epoch", &epoch);
    tree->SetBranchAddress("amplC4", &ampl);

    // Read entries
    int event_start = 0;
    int event_end = std::min(1000, (int)tree->GetEntries());

    std::vector<float> epoch_data, ampl_data;

    for (int i = event_start; i < event_end; i++) {
        tree->GetEntry(i);
        epoch_data.push_back(epoch);
        ampl_data.push_back(ampl);
    }

    // Print some values for verification
    std::cout << "Read " << epoch_data.size() << " entries.\n";
    for (size_t i = 0; i < std::min(size_t(10), epoch_data.size()); i++) {
        std::cout << "Entry " << i << ": epoch = " << epoch_data[i] << ", ampl = " << ampl_data[i] << std::endl;
    }

    std::cout << "donzo" << std::endl;

    file->Close();
    return 0;
}
