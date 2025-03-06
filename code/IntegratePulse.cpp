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
bool ReadWaveformData(const string& filename, double& dt, vector<double>& dataValues,vector<double>& timeValues, vector<double>& drvValues,
                      double& threshold, double& bsl, double& rms, double INTEGRATION_TIME_TRIG,
                      double& total_bkg_rejection_probability, double& ion_tail_end_point_threshold_fraction,
                      double CIVIDEC_PULSE_DURATION) {
    ifstream infile(filename);
    if (!infile.is_open()) {
        cerr << "Error: Unable to open file " << filename << endl;
        return false;
    }

    string line;
    while (getline(infile, line)) {
        if (line.find("dt:") != string::npos) {
            istringstream iss(line.substr(line.find(":") + 1));
            iss >> dt;
        } else if (line.find("double data") != string::npos) {
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
        }else if (line.find("double drv") != string::npos) {
            size_t start = line.find("{");
            size_t end = line.find("}");
            if (start != string::npos && end != string::npos) {
                string drvStr = line.substr(start + 1, end - start - 1);
                istringstream drvStream(drvStr);
                double value;
                while (drvStream >> value) {
                    drvValues.push_back(value);
                    if (drvStream.peek() == ',') drvStream.ignore();
                }
            }
        }else if (line.find("double time") != string::npos) {
            size_t start = line.find("{");
            size_t end = line.find("}");
            if (start != string::npos && end != string::npos) {
                string timeStr = line.substr(start + 1, end - start - 1);
                istringstream timeStream(timeStr);
                double value;
                while (timeStream >> value) {
                    timeValues.push_back(value);
                    if (timeStream.peek() == ',') timeStream.ignore();
                }
            }
        }else if (line.find("Threshold:") != string::npos) {
            istringstream iss(line.substr(line.find(":") + 1));
            iss >> threshold;
            cout << "Threshold: " << threshold << endl;
        } else if (line.find("RMS:") != string::npos) {
            istringstream iss(line.substr(line.find(":") + 1));
            iss >> rms;
            cout << "RMS: " << rms << endl;
        } else if (line.find("BSL:") != string::npos) {
            istringstream iss(line.substr(line.find(":") + 1));
            iss >> bsl;
            cout << "BSL: " << bsl << endl;
        } else if (line.find("INTEGRATION_TIME_TRIGGER:") != string::npos) {
            istringstream iss(line.substr(line.find(":") + 1));
            iss >> INTEGRATION_TIME_TRIG;
            cout << "INTEGRATION_TIME_TRIGGER: " << INTEGRATION_TIME_TRIG << endl;
        } else if (line.find("total_bkg_rejection_probability:") != string::npos) {
            istringstream iss(line.substr(line.find(":") + 1));
            iss >> total_bkg_rejection_probability;
            cout << "total_bkg_rejection_probability: " << total_bkg_rejection_probability << endl;
        } else if (line.find("ion_tail_end_point_threshold_fraction:") != string::npos) {
            istringstream iss(line.substr(line.find(":") + 1));
            iss >> ion_tail_end_point_threshold_fraction;
            cout << "ion_tail_end_point_threshold_fraction: " << ion_tail_end_point_threshold_fraction << endl;
        } else if (line.find("CIVIDEC_PULSE_DURATION:") != string::npos) {
            istringstream iss(line.substr(line.find(":") + 1));
            iss >> CIVIDEC_PULSE_DURATION;
            cout << "CIVIDEC_PULSE_DURATION: " << CIVIDEC_PULSE_DURATION << endl;
        }
    }

    infile.close();
    return true;
}

int IntegratePulse() {
    gROOT->LoadMacro("MyFunctions.C");

    // File containing waveform data
     string filename ="/home/akallits/Documents/PicoAnalysis/Saclay_Analysis/data/2022_October_h4/plots/Run224/Pool2/moreplots/waveform_data.txt" ;
    //string filename ="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2022_October_h4/plots/Run224/Pool2/moreplots/waveform_data.txt" ;


    // Variables to hold dt and data values
    double dt = 0.0;
    double threshold = 0.0;
    double bsl = 0.0;
    double rms = 0.0;
    double epeak_width = 5.0; // Width of the peak in ns
    //double ion_tail_width = 100.0; // Width of the ion tail in ns
    // double ion_tail_end_point_threshold = 0.01; // Fraction of the ion tail width to integrate to
    double total_bkg_rejection_probability = 0.0; // for the total bkg waveform points to be rejected
    double ion_tail_end_point_threshold_fraction = 0.0; // Set ion tail end point fraction
    int maxpoints = 10002;


    vector<double> dataValues;
    vector<double> drvValues;
    vector<double> timeValues;

    // Read the waveform data
    if (!ReadWaveformData(filename, dt, dataValues, timeValues, drvValues, threshold, bsl, rms, INTEGRATION_TIME_TRIG, total_bkg_rejection_probability, ion_tail_end_point_threshold_fraction, CIVIDEC_PULSE_DURATION)) {
        return -1; // Exit if file reading fails
    }

    // Convert dataValues vector to array for compatibility with ROOT
    int npoints = dataValues.size();
    double* data = &dataValues[0];
    double* drv = &drvValues[0];

    double tintValues[] = {epeak_width}; // Width of the peak in ns (5.0)
    int nTint = sizeof(tintValues) / sizeof(tintValues[0]);

    // Colors for the different tint plots
    int colors[] = {kRed, kBlue, kGreen, kMagenta, kOrange};

    // Calculate data min and max for dynamic axis ranges
    double minData = *min_element(data, data + npoints);
    double maxData = *max_element(data, data + npoints);

    // Determine min and max for y-axis range
    double minDrv = *min_element(drv, drv + npoints);
    double maxDrv = *max_element(drv, drv + npoints);

    double globalMin = min(minData, minDrv) - 1;
    double globalMax = max(maxData, maxDrv) + 1;



	//Replace xCoordinates with timeValues
    double* xValues = &timeValues[0];

    // Create a canvas with two pads
    TCanvas *c1 = new TCanvas("c1", "Integration and Original Data", 800, 800);
    c1->Divide(1, 2);


    // Pad 1: Original Data
    c1->cd(1);
    TGraph *graphData = new TGraph(npoints, xValues, data);
    for (int i = 0; i < npoints; i++) {
        graphData->SetPoint(i, i * dt, data[i]);
    }
    //graphData->SetTitle("Original Data;Time [ns];Time [ns]");
    graphData->SetLineColor(kBlack);
    graphData->SetLineWidth(3);
    //graphData->Draw("AL");
    // graphData->GetXaxis()->SetRangeUser(0, npoints * dt);
    // graphData->GetYaxis()->SetRangeUser(minData - 1, maxData + 1);

    //add the derivative of the data in the same graph as the original data
    TGraph *graphDrv = new TGraph(npoints, xValues, drv);
    for (int i = 0; i < npoints; i++) {
        graphDrv->SetPoint(i, i * dt, drv[i]);
    }
    graphDrv->SetTitle("Derivative Data;Time [ns];Amplitude [V/ns]");
    graphDrv->SetLineColor(kRed);
    graphDrv->SetLineWidth(2);
    //graphDrv->Draw("SAME");
    // graphDrv->GetXaxis()->SetRangeUser(0, npoints * dt);


    // graphData->GetYaxis()->SetRangeUser(globalMin, globalMax); // Adjust Y-axis
    // graphData->Draw("AL");  // Original data
    // graphDrv->Draw("L SAME");  // Derivative data

    //calculate the cumulative distribution of the waveform data

    vector<double> cdfValues(npoints, 0.0);
    cdfValues[0] = dataValues[0];


    for (int i = 1; i < npoints; i++) {
        cdfValues[i] = cdfValues[i - 1] + dataValues[i];
    }

    // Normalize CDF
    for (int i = 0; i < npoints; i++) {
        cdfValues[i] /= -fabs(cdfValues[npoints - 1]);  // Normalize to range [0,1] -- this makes the pulse positive
    }
    double* cdf = &cdfValues[0];

    // Determine min and max for y-axis range
    double minCdf = *min_element(cdf, cdf + npoints);
    double maxCdf = *max_element(cdf, cdf + npoints);

    double globalMinCDF = min(minData, minCdf) - 1;
    double globalMaxCDF = max(maxData, maxCdf) + 1;

    TGraph *graphCDF = new TGraph(npoints, xValues, &cdfValues[0]);
    graphCDF->SetTitle("Cumulative Distribution Function;Time [ns];CDF");
    graphCDF->SetLineColor(kBlue);
    graphCDF->SetLineWidth(2);

    //Derivative of the CDF
    vector<double> drvCDFValues(npoints, 0.0);
    //convert drvCDFValues vector to array


    for (int i = 1; i < npoints; i++) {
        drvCDFValues[i] = (cdfValues[i] - cdfValues[i - 1]) / dt;
    }
    double* drvCDF = &drvCDFValues[0];
    double drvCFD2[npoints];
    double smoothCFD[npoints];
    double smoothdata[npoints];
    double smoothdata2[npoints];

    SmoothArray(cdf, smoothCFD, npoints,11, 1);
    // SmoothArray(data, smoothdata, npoints,50, 1);
    // SmoothArray(smoothdata, smoothdata2, npoints,173, 1);
    // DerivateArray(cdf, drvCFD2, npoints, dt,3);
    DerivateArray(smoothCFD, drvCFD2, npoints, dt,5);
    for (int i=0; i<npoints; i++) {
        drvCFD2[i]*=-10;
    }
    // Determine min and max for y-axis range
    double mindrvCDF = *min_element(drvCDF, drvCDF + npoints);
    double maxdrvCDF = *max_element(drvCDF, drvCDF + npoints);

    double globalMindrvCDF  = min(minData, mindrvCDF) - 1;
    double globalMaxdrvCDF  = max(maxData, maxdrvCDF) + 1;


    //plot the derivative of the CDF
    TGraph *graphDrvCDF = new TGraph(npoints, xValues, drvCDF);
    graphDrvCDF->SetTitle("Derivative of the Cumulative Distribution Function;Time [ns];dCDF/dt");
    graphDrvCDF->SetLineColor(kGreen);
    graphDrvCDF->SetLineWidth(2);

    TGraph *graphDrvCDF2 = new TGraph(npoints, xValues, drvCFD2);
    graphDrvCDF2->SetTitle("Derivative of the Cumulative Distribution Function;Time [ns];dCDF/dt");
    graphDrvCDF2->SetLineColor(kMagenta);
    graphDrvCDF2->SetLineWidth(2);

    c1->cd(1);
    graphData->GetYaxis()->SetRangeUser(globalMin, globalMaxCDF); // Adjust Y-axis
    graphData->Draw("AL");  // Original data
    // graphCDF->Draw("AL SAME");  // CDF data
    // graphDrvCDF->Draw("AL SAME"); // Derivative of the CDF data
    // graphDrv->Draw("AL"); // Derivative data
    // graphCDF->Draw("AL");  // CDF data
    graphDrvCDF2->Draw("L SAME"); // Derivative of the CDF data
    // graphDrv->Draw("L SAME"); // Derivative data



    // Add a legend
    TLegend *legend1 = new TLegend(0.7, 0.7, 0.9, 0.85);
    legend1->AddEntry(graphData, "Original Data", "l");
    legend1->AddEntry(graphDrv, "Derivative Data", "l");
    legend1->AddEntry(graphCDF, "CDF Data", "l");
    legend1->AddEntry(graphDrvCDF, "Derivative of CDF Data", "l");
    legend1->AddEntry(graphDrvCDF2, "Derivative of Smoothed CDF Data", "l");
    legend1->Draw();

    // Pad 2: Integral Curves
    c1->cd(2);

    // Legend for integral plots
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->SetHeader("Tint Values", "C");

    // Track min and max for integral data to adjust axes
    double minIntegral = 1e9, maxIntegral = -1e9;

    for (int t = 0; t < nTint; t++) {
        double tint = tintValues[t];
		int npt = (int)TMath::FloorNint(tint / dt);

        //double threshold_npt = threshold * sqrt(npt);
        //double ion_tail_end_point_threshold_npt = ion_tail_end_point_threshold * sqrt(npt);
        // Perform integration
        cout<<"Integration npoints = "<< npt<<endl;


        vector<pair<double, double>> trigger_windows = GetTriggerWindows(xValues, npoints, data, dt,
            threshold);

        // Extract integral data from GetTriggerWindows
        vector<double> integrated_data_vec(data, data + npoints);
        //print the integrated data vector values
        for (int i; i<integrated_data_vec.size(); i++){
            cout << integrated_data_vec[i] << endl;
        }

        // cout << "Size of integrated_data_vec: " << integrated_data_vec.size() << endl;
        // cin.get();
        auto [x_int, y_int] = IntegratePulse_std(timeValues, integrated_data_vec, npt);

        //Call function to plot integral with bounds
        PlotIntegralWithBounds(x_int, y_int, trigger_windows, minIntegral, maxIntegral, t, c1, legend, std::vector<int>(colors, colors + sizeof(colors) / sizeof(colors[0])), tint);

        //plot trigger windows on the original data
        c1->cd(1);
        for (size_t i = 0; i < trigger_windows.size(); i++) {
            double x_left = trigger_windows[i].first;
            double x_right = trigger_windows[i].second;
            TLine *line_left = new TLine(x_left, minData, x_left, maxData);
            TLine *line_right = new TLine(x_right, minData, x_right, maxData);
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

    vector<double> integrated_data_vec(data, data + npoints);
    //print the integrated data vector values
    for (int i; i<integrated_data_vec.size(); i++){
        cout << integrated_data_vec[i] << endl;
    }
    // Calculate CDF of integrated data

    vector<double> cdfValues_int(npoints, 0.0);
    cdfValues_int[0] = integrated_data_vec[0];
    for (int i = 1; i < npoints; i++) {
        cdfValues_int[i] = cdfValues_int[i - 1] + integrated_data_vec[i];
    }

    // Normalize CDF of integrated data
    for (int i = 0; i < npoints; i++) {
        cdfValues_int[i] /= -fabs(cdfValues_int[npoints - 1]);  // Normalize to range [0,1] -- this makes the pulse positive
    }
    double* cdf_int = &cdfValues_int[0];
    double* cdf_int_smoothed = &cdfValues_int[0];

    SmoothArray(cdf_int,cdf_int_smoothed, npoints,173, 1);

    TGraph *graphCDF_int = new TGraph(npoints, xValues, &cdfValues_int[0]);
    graphCDF_int->SetTitle("Cumulative Distribution Function;Time [ns];CDF");
    graphCDF_int->SetLineColor(kRed);
    graphCDF_int->SetLineWidth(1);

    TGraph *graphCDF_int_smoothed = new TGraph(npoints, xValues, &cdf_int_smoothed[0]);
    graphCDF_int_smoothed->SetTitle("Cumulative Distribution Function;Time [ns];CDF");
    graphCDF_int_smoothed->SetLineColor(kOrange);
    graphCDF_int_smoothed->SetLineWidth(2);


    //graphCDF_int->Draw("AL");
    // Draw legend
    legend->Draw();

    // Save canvas
    c1->SaveAs("IntegrationAndOriginalData.pdf");

    TCanvas* c2 = new TCanvas("c2", "mpla mpla", 800, 800);
    graphCDF->Draw("AL");  // CDF data
    graphCDF_int->Draw("L SAME");  // CDF data
    graphCDF_int_smoothed->Draw("L SAME");  // CDF int data smoothed
    c2->Update();

    return 0;
}

