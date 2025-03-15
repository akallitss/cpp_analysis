//
// Created by dylan on 3/15/25.
//

#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TMath.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <map>

class RMSBaselineCalculator {
public:
    RMSBaselineCalculator(const std::string& filePath, const std::string& treeName, int channel, int nPoints = 1000, double sigmaCut = 3.0, int minEventsPerEpoch = 5)
        : filePath(filePath), treeName(treeName), channel(channel), nPoints(nPoints), sigmaCut(sigmaCut), minEventsPerEpoch(minEventsPerEpoch) {}

    void Process();

private:
    std::string filePath;
    std::string treeName;
    int channel;
    int nPoints;
    double sigmaCut;
    int minEventsPerEpoch;
};

void RMSBaselineCalculator::Process() {
    TFile file(filePath.c_str(), "READ");
    if (!file.IsOpen()) {
        std::cerr << "Error opening file: " << filePath << std::endl;
        return;
    }

    TTree* tree = (TTree*)file.Get(treeName.c_str());
    if (!tree) {
        std::cerr << "Error: Tree " << treeName << " not found in file." << std::endl;
        return;
    }

    Double_t waveform[20000];
    ULong64_t epoch;
    int eventNo;
    tree->SetBranchAddress(Form("amplC%d", channel), waveform);
    tree->SetBranchAddress("epoch", &epoch);
    tree->SetBranchAddress("eventNo", &eventNo);

    std::vector<double> rmses;
    std::vector<int> eventNumbers;
    std::vector<ULong64_t> epochs;
    std::vector<double> waveforms;

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        double rms = TMath::Abs(waveform[0]);
        rmses.push_back(rms);
        eventNumbers.push_back(eventNo);
        epochs.push_back(epoch);
        waveforms.push_back(waveform[0]);
    }

    double medianRMS = TMath::Median(rmses.size(), rmses.data());
    double rmsRMS = TMath::RMS(rmses.size(), rmses.data());

    std::map<ULong64_t, std::pair<double, int>> epochSums;
    std::map<ULong64_t, int> epochCounts;

    for (size_t i = 0; i < rmses.size(); ++i) {
        if (rmses[i] <= medianRMS + sigmaCut * rmsRMS) {
            epochSums[epochs[i]].first += waveforms[i];
            epochSums[epochs[i]].second++;
            epochCounts[epochs[i]]++;
        }
    }

    std::map<ULong64_t, std::pair<double, int>> mergedEpochs;
    ULong64_t lastEpoch = 0;
    for (const auto& pair : epochSums) {
        ULong64_t epoch = pair.first;
        double sumBaseline = pair.second.first;
        int count = pair.second.second;

        if (lastEpoch != 0 && mergedEpochs[lastEpoch].second < minEventsPerEpoch) {
            mergedEpochs[lastEpoch].first += sumBaseline;
            mergedEpochs[lastEpoch].second += count;
        } else {
            mergedEpochs[epoch] = {sumBaseline, count};
            lastEpoch = epoch;
        }
    }

    for (const auto& pair : mergedEpochs) {
        double avgBaseline = pair.second.first / pair.second.second;
        std::cout << "Epoch: " << pair.first << " Avg Baseline: " << avgBaseline << " Events: " << pair.second.second << std::endl;
    }

    file.Close();
}

int calculate_rms_baseline() {
    std::string file_path = "/home/dylan/Desktop/picosec/data/Run224-Pool2_TESTBEAM_tree.root";
    RMSBaselineCalculator calculator(file_path, "RawDataTree", 2);
    calculator.Process();
    return 0;
}