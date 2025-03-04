//
// Created by akallits on 12/12/24.
//
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMath.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <numeric>
#include <algorithm>
#include <cmath>
#include "MyFunctions.C"

using namespace std;

// Function to read data from the text file
bool ReadWaveformData(const string& filename, double& dt, vector<double>& dataValues, double& threshold, double& bsl, double& rms, double& total_bkg_rejection_probability, double& ion_tail_end_point_threshold_fraction, double& CIVIDEC_PULSE_DURATION) {
    ifstream infile(filename);
    if (!infile.is_open()) {
        cerr << "Error: Unable to open file " << filename << endl;
        return false;
    }

    string line;
    while (getline(infile, line)) {
        // Parse the lines based on their content
        if (line.find("dt:") != string::npos) {
            // Extract dt value
            istringstream iss(line.substr(line.find(":") + 1));
            iss >> dt;
        } else if (line.find("double data") != string::npos) {
            // Extract the data array
            size_t start = line.find("{");
            size_t end = line.find("}");
            if (start != string::npos && end != string::npos) {
                string dataStr = line.substr(start + 1, end - start - 1);
                istringstream dataStream(dataStr);
                double value;
                while (dataStream >> value) {
                    dataValues.push_back(value);
                    if (dataStream.peek() == ',') dataStream.ignore();
                }
            }
        } else if (line.find("Threshold: ") != std::string::npos) {
            // Parse the threshold value
            std::istringstream iss(line.substr(line.find(":") + 1));
            iss >> threshold;
            std::cout << "Threshold: " << threshold << std::endl;
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
        }
        else if (line.find("INTEGRATION_TIME_TRIGGER: ") != std::string::npos) {
            // Parse the INTEGRATION_TIME_TRIG value
            std::istringstream iss(line.substr(line.find(":") + 1));
            iss >> INTEGRATION_TIME_TRIG;
            std::cout << "INTEGRATION_TIME_TRIGGER: " << bsl << std::endl;
        }
        else if (line.find("total_bkg_rejection_probability: ") != std::string::npos) {
            // Parse the total_bkg_rejection_probability value
            std::istringstream iss(line.substr(line.find(":") + 1));
            iss >> total_bkg_rejection_probability;
            std::cout << "total_bkg_rejection_probability: " << bsl << std::endl;
        }
        else if (line.find("ion_tail_end_point_threshold_fraction: ") != std::string::npos) {
            // Parse the ion_tail_end_point_threshold_fraction value
            std::istringstream iss(line.substr(line.find(":") + 1));
            iss >> ion_tail_end_point_threshold_fraction;
            std::cout << "ion_tail_end_point_threshold_fraction: " << bsl << std::endl;
        }
        else if (line.find("CIVIDEC_PULSE_DURATION: ") != std::string::npos) {
            // Parse the CIVIDEC_PULSE_DURATION value
            std::istringstream iss(line.substr(line.find(":") + 1));
            iss >> CIVIDEC_PULSE_DURATION;
            std::cout << "CIVIDEC_PULSE_DURATION: " << bsl << std::endl;
        }
    }



    infile.close();
    return true;
}

int IntegratePulse() {
    gROOT->LoadMacro("MyFunctions.C");
//    // Define parameters
//    int npoints = 100;
//    double dt = 1.0; // time step
//    double data[npoints];
//
//    // Generate mock data (a sine wave with noise)
//    for (int i = 0; i < npoints; i++) {
//        double time = i * dt;
//        data[i] = 10 * TMath::Sin(0.1 * time) + 2 * gRandom->Gaus();
//    }
    // File containing waveform data
     string filename ="/home/akallits/Documents/PicoAnalysis/Saclay_Analysis/data/2022_October_h4/plots/Run224/Pool2/moreplots/waveform_data.txt" ;

    //string filename ="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2022_October_h4/plots/Run224/Pool2/moreplots/waveform_data.txt" ;
    // Variables to hold dt and data values
    double dt = 0.0;
    double threshold = 0.0;
    double bsl = 0.0;
    double rms = 0.0;
    double epeak_width = 5.0; // Width of the peak in ns
    double ion_tail_width = 100.0; // Width of the ion tail in ns
    // double ion_tail_end_point_threshold = 0.01; // Fraction of the ion tail width to integrate to
    double total_bkg_rejection_probability = 0.0; // for the total bkg waveform points to be rejected
    double ion_tail_end_point_threshold_fraction = 0.0; // Set ion tail end point fraction
    int maxpoints = 10002;


    vector<double> dataValues;

    // Read the waveform data
    if (!ReadWaveformData(filename, dt, dataValues, threshold, bsl, rms,
        total_bkg_rejection_probability, ion_tail_end_point_threshold_fraction, CIVIDEC_PEAK_DURATION)) {
        return -1; // Exit if file reading fails
    }

    // Convert dataValues vector to array for compatibility with ROOT
    int npoints = dataValues.size();
	//double data[npoints];     // Populate with data values
	//double integral[npoints]; // To store integral values
    double* data = &dataValues[0];
    threshold = -6 * rms; // Convert threshold to units of RMS
    double ion_tail_end_point_threshold = 0.2 * threshold; // Set ion tail end point fraction to 10% of threshold


    // double tintValues[] = {1.0, 2.0, 5.0, 10.0};
    // double tintValues[] = {1.0, 2.0, 5.0, 10.0, 100, 150, 200};
    double tintValues[] = {epeak_width}; // Width of the peak in ns (5.0)
    int nTint = sizeof(tintValues) / sizeof(tintValues[0]);

    // Colors for the different tint plots
    int colors[] = {kRed, kBlue, kGreen, kMagenta, kOrange};

    // Calculate data min and max for dynamic axis ranges
    double minData = *min_element(data, data + npoints);
    double maxData = *max_element(data, data + npoints);
	double xValues[npoints];
	for (int i = 0; i < npoints; i++) {
    		xValues[i] = i * dt; // Replace with your time data if reading from a file
	}

    // Create a canvas with two pads
    TCanvas *c1 = new TCanvas("c1", "Integration and Original Data", 800, 800);
    c1->Divide(1, 2);


    // Pad 1: Original Data
    c1->cd(1);
    TGraph *graphData = new TGraph(npoints, xValues, data);
    for (int i = 0; i < npoints; i++) {
        graphData->SetPoint(i, i * dt, data[i]);
    }
    graphData->SetTitle("Original Data;Time [ns];Amplitude [V]");
    graphData->SetLineColor(kBlack);
    graphData->SetLineWidth(2);
    graphData->Draw("AL");
    graphData->GetXaxis()->SetRangeUser(0, npoints * dt);
    graphData->GetYaxis()->SetRangeUser(minData - 1, maxData + 1);


    // Pad 2: Integral Curves
    c1->cd(2);

    // Legend for integral plots
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->SetHeader("Tint Values", "C");

    // Track min and max for integral data to adjust axes
    double minIntegral = 1e9, maxIntegral = -1e9;

    for (int t = 0; t < nTint; t++) {
    //for (int t = 0; t < 1; t++) {
        double tint = tintValues[t];
        // double integral[npoints];
		int npt = (int)TMath::FloorNint(tint / dt);
        double threshold_npt = threshold * sqrt(npt);
        double ion_tail_end_point_threshold_npt = ion_tail_end_point_threshold * sqrt(npt);
        // Perform integration
                cout<<"npoints = "<< npt<<endl;
        // IntegratePulse(npoints, data, integral, dt, tint);
        // vector<double> t_values(npoints);
        // for (int i = 0; i < npoints; i++) {
        //     t_values[i] = i * dt;
        // }

        //convert time vector to array
        double t_values[npoints];
        for (int i = 0; i < npoints; i++) {
            t_values[i] = i * dt;
        }

        vector<pair<double, double>> trigger_windows = GetTriggerWindows(t_values, npoints, data, dt,
            threshold);

        // vector<double> data_vec = {data, data + npoints};
        // auto [x_int, y_int] = IntegratePulse_std(t_values, data_vec, npt);
        //
        // // Find pulse bounds
        // cout << "threshold_npt = " << threshold_npt << endl;
        // cout << "ion_tail_end_point_fraction_npt = " << ion_tail_end_point_threshold_npt << endl;
        // vector<pair<double, double>> pulse_bounds = find_pulse_bounds(x_int, y_int, threshold_npt, ion_tail_width,
        // ion_tail_end_point_threshold_npt);
        //
        // //For each bound, add npoints/2*dt to the left and subtract npoints/2*dt from the right to get the full pulse
        // adjust_pulse_bounds(pulse_bounds, npt, dt);

        if (tint == 5.0) {
            c1->cd(1);
            // Draw vertical lines to indicate the pulse bounds
            for (size_t i = 0; i < pulse_bounds.size(); i++) {
                double x_left = pulse_bounds[i].first;
                double x_right = pulse_bounds[i].second;
                TLine *line_left = new TLine(x_left, minData, x_left, maxData);
                TLine *line_right = new TLine(x_right, minData, x_right, maxData);
                line_left->SetLineStyle(2);
                line_right->SetLineStyle(2);
                line_left->SetLineColor(colors[t % 5]);
                line_right->SetLineColor(colors[t % 5]);
                line_left->Draw();
                line_right->Draw();
            }
            c1->cd(2);
        }


        double localMin = *min_element(y_int.begin(), y_int.end());
        double localMax = *max_element(y_int.begin(), y_int.end());
        minIntegral = min(minIntegral, localMin);
        maxIntegral = max(maxIntegral, localMax);




        TGraph *graphIntegral = new TGraph(x_int.size(), &x_int[0], &y_int[0]);
        for (int i = 0; i < x_int.size(); i++) {
            //graphIntegral->SetPoint(i, x_int[i], y_int[i]);
            graphIntegral->SetPoint(i, x_int[i], y_int[i]);
        }

        // Set graph style
        graphIntegral->SetLineColor(colors[t % 5]);
        graphIntegral->SetLineWidth(2);
        graphIntegral->SetTitle(";Time [ns];Integral [V*ns]");

        // Add graph to canvas
        if (t == 0) {
            graphIntegral->Draw("AL"); // Draw the first graph with axes
        } else {
            graphIntegral->Draw("L"); // Overlay subsequent graphs
        }

        // Add entry to legend
        legend->AddEntry(graphIntegral, Form("Tint = %.1f ns", tint), "l");

        //Add vertical lines to indicate the pulse bounds
        for (size_t i = 0; i < pulse_bounds.size(); i++) {
            double x_left = pulse_bounds[i].first;
            double x_right = pulse_bounds[i].second;
            TLine *line_left = new TLine(x_left, minIntegral, x_left, maxIntegral);
            TLine *line_right = new TLine(x_right, minIntegral, x_right, maxIntegral);
            line_left->SetLineStyle(2);
            line_right->SetLineStyle(2);
            line_left->SetLineColor(colors[t % 5]);
            line_right->SetLineColor(colors[t % 5]);
            line_left->Draw();
            line_right->Draw();
        }



    }

    // Adjust axis ranges for integral plot
    c1->cd(2);
    gPad->Modified();
    gPad->Update();
    TGraph *g = (TGraph*) gPad->FindObject("Graph"); // Grab the first graph
    if (g) {
        g->GetXaxis()->SetRangeUser(0, npoints * dt);
        g->GetYaxis()->SetRangeUser(minIntegral - 1, maxIntegral + 1);
    }


    // Draw legend
    legend->Draw();

    // Save canvas
    c1->SaveAs("IntegrationAndOriginalData.pdf");

    return 0;
}
