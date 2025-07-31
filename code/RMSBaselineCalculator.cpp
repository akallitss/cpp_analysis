#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TMath.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <map>
#include <chrono>
#include "MyFunctions.C"

class RMSBaselineCalculator {
public:
    RMSBaselineCalculator() = default; // Default constructor

    RMSBaselineCalculator(const std::string& filePath, const std::string& treeName, int channel, double tRMSCut = 100, double sigmaCut = 3.0, int minEventsPerEpoch = 5, double integrationTime = 5.0)
        : filePath(filePath), treeName(treeName), channel(channel + 1), tRMSCut(tRMSCut), sigmaCut(sigmaCut), minEventsPerEpoch(minEventsPerEpoch), integrationTime(integrationTime){}

    void Process();
    void Plot();

    double get_epoch_baseline(ULong64_t epoch) {
        if (epoch_baselines.find(epoch) == epoch_baselines.end()) {
            std::cerr << "Error: Epoch " << epoch << " not found in get_epoch_baseline, for channel: " << channel <<std::endl;
            return -1;
        }
        return epoch_baselines[epoch];
    }
    double get_epoch_integral_rms(ULong64_t epoch) {
        if (epoch_integral_rmses.find(epoch) == epoch_integral_rmses.end()) {
            std::cerr << "Error: Epoch " << epoch << " not found in get_epoch_integral_rms, for channel: " << channel <<std::endl;
            return -1;
        }
        return epoch_integral_rmses[epoch];
    }
    double get_epoch_rms(ULong64_t epoch) {
        if (epoch_rmses.find(epoch) == epoch_rmses.end()) {
            std::cerr << "Error: Epoch " << epoch << " not found in get_epoch_rms, for channel: " << channel << std::endl;
            return -1;
        }
        return epoch_rmses[epoch];
    }

private:
    std::string filePath;
    std::string treeName;
    int channel;
    double tRMSCut;
    double sigmaCut;
    int minEventsPerEpoch;
    double integrationTime;

    map<ULong64_t, double> epoch_baselines;
    map<ULong64_t, double> epoch_rmses;
    map<ULong64_t, double> epoch_integral_rmses;

    TGraph* CreateGraph(const std::map<ULong64_t, double>& data, const char* xTitle, const char* yTitle) {
        int n = data.size();
        double* x = new double[n];
        double* y = new double[n];

        int i = 0;
        for (const auto& [key, value] : data) {
            x[i] = key;
            y[i] = value;
            i++;
        }

        TGraph* graph = new TGraph(n, x, y);
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(1.2);
        graph->SetTitle("");
        graph->GetXaxis()->SetTitle(xTitle);
        graph->GetYaxis()->SetTitle(yTitle);

        delete[] x;
        delete[] y;

        return graph;
    }
};

void RMSBaselineCalculator::Process() {
  	auto start_time = std::chrono::high_resolution_clock::now();
//    cout << "Processing file: " << filePath << endl;
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

    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus(Form("amplC%d", channel), 1);
    tree->SetBranchStatus("epoch", 1);
    tree->SetBranchStatus("eventNo", 1);
    tree->SetBranchStatus("sumpoints", 1);
    tree->SetBranchStatus("dt", 1);

//    cout << "Getting array size" << endl;
    const int ARRAYSIZE = tree->GetMaximum("sumpoints");
    double dt = tree->GetMaximum("dt");
    if (dt != tree->GetMinimum("dt")) {
        std::cerr << "Error: dt values are not consistent." << std::endl;
        return;
    }

//    cout << "Creating time values" << endl;
    double timeValues[ARRAYSIZE];
    for (int i = 0; i < ARRAYSIZE; ++i) {
        timeValues[i] = i * dt;
    }

    int integrationPoints = TMath::Floor(integrationTime / dt);
    int iRMSCut = convert_x_to_index(timeValues, ARRAYSIZE, tRMSCut);

    Double_t waveform[ARRAYSIZE];
    ULong64_t epoch;
    int eventNo;
    tree->SetBranchAddress(Form("amplC%d", channel), waveform);
    tree->SetBranchAddress("epoch", &epoch);
    tree->SetBranchAddress("eventNo", &eventNo);

    std::vector<double> rmses;
    std::vector<double> integral_rmses;
    std::vector<double> baselines;
    std::vector<int> eventNumbers;
    std::vector<ULong64_t> epochs;
    std::vector<std::vector<double>> waveforms;

//    cout << "Getting event by event stats" << endl;
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        std::vector<double> t_values(timeValues, timeValues + iRMSCut);
        std::vector<double> y_values(waveform, waveform + iRMSCut);

        auto [x_int, y_int] = IntegratePulse_std(t_values, y_values, integrationPoints);

        // Get RMS of y_values from 0 to index iRMSIntegrationCut
        double integral_rms = TMath::RMS(y_int.begin(), y_int.end());
        double rms = TMath::RMS(y_values.begin(), y_values.end());
        double baseline = TMath::Mean(y_values.begin(), y_values.end());

        baselines.push_back(baseline);
        integral_rmses.push_back(integral_rms);
        rmses.push_back(rms);
        waveforms.push_back(y_values);
        eventNumbers.push_back(eventNo);
        epochs.push_back(epoch);
    }

//    cout << "Calculating medians" << endl;
    double medianIntegralRMS = TMath::Median(integral_rmses.size(), integral_rmses.data());
    double rmsIntegralRMS = TMath::RMS(integral_rmses.begin(), integral_rmses.end());

    std::map<ULong64_t, vector<double>> epochBaselinesAll;
    std::map<ULong64_t, vector<double>> epochIntegralRMSesAll;
    std::map<ULong64_t, vector<double>> epochRMSesAll;
    std::map<ULong64_t, int> epochCounts;

//    cout << "Filtering events" << endl;
    for (size_t i = 0; i < rmses.size(); ++i) {
//        if (rmses[i] <= medianRMS + sigmaCut * rmsRMS) {
        if (integral_rmses[i] <= medianIntegralRMS + sigmaCut * rmsIntegralRMS) {
			epochIntegralRMSesAll[epochs[i]].push_back(integral_rmses[i]);
            epochRMSesAll[epochs[i]].push_back(rmses[i]);
            epochBaselinesAll[epochs[i]].push_back(baselines[i]);
            epochCounts[epochs[i]]++;
        }
    }

//    cout << "Calculating medians for epochs" << endl;
    ULong64_t lastEpoch = epochCounts.begin()->first;
    for (auto epoch_it = epochCounts.begin(); epoch_it != epochCounts.end(); ++epoch_it) {
        ULong64_t epoch = epoch_it->first;
        if (epoch_it->second >= minEventsPerEpoch) {
          	epoch_baselines[epoch] = TMath::Median(epochBaselinesAll[epoch].size(), epochBaselinesAll[epoch].data());
          	epoch_integral_rmses[epoch] = TMath::Median(epochIntegralRMSesAll[epoch].size(), epochIntegralRMSesAll[epoch].data());
          	epoch_rmses[epoch] = TMath::Median(epochRMSesAll[epoch].size(), epochRMSesAll[epoch].data());
        } else {
            epoch_baselines[epoch] = epoch_baselines[lastEpoch];
            epoch_integral_rmses[epoch] = epoch_integral_rmses[lastEpoch];
            epoch_rmses[epoch] = epoch_rmses[lastEpoch];
        }
    }

    // Calculate for the first epoch and use the second epoch if not enough events
    ULong64_t firstEpoch = epochCounts.begin()->first;
    if (epochCounts.begin()->second < minEventsPerEpoch) {
        epoch_baselines[firstEpoch] = epoch_baselines[(++epochCounts.begin())->first];
        epoch_integral_rmses[firstEpoch] = epoch_integral_rmses[(++epochCounts.begin())->first];
        epoch_rmses[firstEpoch] = epoch_rmses[(++epochCounts.begin())->first];
    } else {
        epoch_baselines[firstEpoch] = TMath::Median(epochBaselinesAll[firstEpoch].size(), epochBaselinesAll[firstEpoch].data());
        epoch_integral_rmses[firstEpoch] = TMath::Median(epochIntegralRMSesAll[firstEpoch].size(), epochIntegralRMSesAll[firstEpoch].data());
        epoch_rmses[firstEpoch] = TMath::Median(epochRMSesAll[firstEpoch].size(), epochRMSesAll[firstEpoch].data());
    }

    // Print the results
    // for (auto& [epoch, baseline] : epoch_baselines) {
        // std::cout << "Epoch: " << epoch << " Baseline: " << baseline << " Integral RMS: " << epoch_integral_rmses[epoch] << " RMS: " << epoch_rmses[epoch] << std::endl;
    // }

    auto stop_time = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop_time - start_time);

    std::cout << "Channel " << channel << " took " << duration.count() << " seconds" << std::endl;

    file.Close();
}

void RMSBaselineCalculator::Plot() {
        // Creating three canvases
        TCanvas* c1 = new TCanvas(("c1_channel_" + std::to_string(channel)).c_str(), "Baseline Plot", 800, 600);
        TCanvas* c2 = new TCanvas(("c2_channel_" + std::to_string(channel)).c_str(), "Integral RMS Plot", 800, 600);
        TCanvas* c3 = new TCanvas(("c3_channel_" + std::to_string(channel)).c_str(), "RMS Plot", 800, 600);

        // Creating graphs
        TGraph* graph_baseline = CreateGraph(epoch_baselines, "Epoch", "Baseline Value");
        TGraph* graph_integral_rms = CreateGraph(epoch_integral_rmses, "Epoch", "Integral RMS Value");
        TGraph* graph_rms = CreateGraph(epoch_rmses, "Epoch", "RMS Value");

        // Drawing graphs on respective canvases
        c1->cd();
        graph_baseline->Draw("APL");

        c2->cd();
        graph_integral_rms->Draw("APL");

        c3->cd();
        graph_rms->Draw("APL");
}

int calculate_rms_baseline() {
    std::string file_path = "/home/akallits/Documents/PicoAnalysis/Saclay_Analysis/data/2022_October_h4/processedTrees/Run224-Pool2_TESTBEAM_tree.root";
    RMSBaselineCalculator calculator_c4(file_path, "RawDataTree", 3);
    calculator_c4.Process();
    calculator_c4.Plot();
    RMSBaselineCalculator calculator_c2(file_path, "RawDataTree", 1);
    calculator_c2.Process();
    calculator_c2.Plot();
    return 0;
}