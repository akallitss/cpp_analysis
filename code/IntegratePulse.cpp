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
                      double CIVIDEC_PULSE_DURATION, double CIVIDEC_PEAK_DURATION) {
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
        }else if (line.find("CIVIDEC_PEAK_DURATION:") != string::npos) {
            istringstream iss(line.substr(line.find(":") + 1));
            iss >> CIVIDEC_PEAK_DURATION;
            cout << "CIVIDEC_PEAK_DURATION: " << CIVIDEC_PEAK_DURATION << endl;
        }
    }

    infile.close();
    return true;
}

int IntegratePulse() {
    gROOT->LoadMacro("MyFunctions.C");

    // File containing waveform data
     string filename ="/home/akallits/Documents/PicoAnalysis/Saclay_Analysis/data/2022_October_h4/plots/Run224/Pool2/moreplots/waveform_data.txt" ;
     // string filename ="/home/akallits/Documents/PicoAnalysis/Saclay_Analysis/data/2023_July_h4/plots/Run057/Pool4/moreplots/waveform_data.txt" ;
    //string filename ="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2022_October_h4/plots/Run224/Pool2/moreplots/waveform_data.txt" ;


    // Variables to hold dt and data values
    double dt = 0.0, threshold = 0.0, bsl = 0.0, rms = 0.0;
    double epeak_width = 5.0; // Width of the peak in ns
    double ion_tail_width = 150.0; // Width of the ion tail in ns
    double total_bkg_rejection_probability = 0.0; // for the total bkg waveform points to be rejected
    double ion_tail_end_point_threshold_fraction = 0.0; // Set ion tail end point fraction


    vector<double> dataValues, drvValues, timeValues, cdfValues, cdfValues_smooth150, cdfValues_int;


    // Read the waveform data
    if (!ReadWaveformData(filename, dt, dataValues, timeValues, drvValues, threshold, bsl, rms, INTEGRATION_TIME_TRIG, total_bkg_rejection_probability, ion_tail_end_point_threshold_fraction, CIVIDEC_PULSE_DURATION, CIVIDEC_PEAK_DURATION)) {
        return -1; // Exit if file reading fails
    }

    // Convert dataValues vector to array for compatibility with ROOT
    int npoints = dataValues.size();

    //convert from vector to array
    double* data = &dataValues[0];
    double* drv = &drvValues[0];
    double* xValues = &timeValues[0];     //Replace xCoordinates with timeValues
    double* cdf;
    double* cdf_smooth150;
    double* cdf_int;

    double drvCFD[npoints];
    double drvsmoothCFD_11_5[npoints];
    double smoothCFD_11[npoints];
    double smoothdata[npoints];
    double smoothdata_150[npoints];
    double cdf_int_smoothed_173[npoints];


    //convert from array to vector

    vector<double> smoothdataVec = vector<double>(smoothdata, smoothdata + npoints);
    vector<double> smoothdata_150Vec = vector<double>(smoothdata_150, smoothdata_150 + npoints);
    vector<double> drvCFDVec = vector<double>(drvCFD, drvCFD + npoints);
    vector<double> drvCFD2Vec = vector<double>(drvsmoothCFD_11_5, drvsmoothCFD_11_5 + npoints);
    vector<double> smoothCFDVec = vector<double>(smoothCFD_11, smoothCFD_11 + npoints);


    // double tintValues[] = {epeak_width, ion_tail_width}; // Width of the peak in ns (5.0)
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


    // Create a canvas with two pads
    TCanvas *c1 = new TCanvas("c1", "Integration and Original Data", 800, 800);
    c1->Divide(1, 2);

    // Pad 1: Original Data
    c1->cd(1);
    TGraph *graphData = new TGraph(npoints, xValues, data);
    for (int i = 0; i < npoints; i++) {
        graphData->SetPoint(i, i * dt, data[i]);
    }
    graphData->SetLineColor(kBlack);
    graphData->SetLineWidth(2);

    //recalculate the baseline here using the first 80ns of the data
    auto baseline_region_end = std::find_if(timeValues.begin(), timeValues.end(), [](double t) { return t > 80.0; });
    size_t baseline_region_end_point = std::distance(timeValues.begin(), baseline_region_end);
    auto [baseline_rms_test, baseline_level_test] = FindBaselineLevel(data, baseline_region_end_point);

    cout<<"Baseline RMS = "<<baseline_rms_test<<endl;
    cout<<"Baseline Level = "<<baseline_level_test<<endl;
    // return 222222;

    // subtract the baseline from the data
    for (int i = 0; i < npoints; i++) {
        data[i] -= baseline_level_test;
    }

    //add the derivative of the data in the same graph as the original data
    TGraph *graphDrv = new TGraph(npoints, xValues, drv);
    for (int i = 0; i < npoints; i++) {
        graphDrv->SetPoint(i, i * dt, drv[i]);
    }
    graphDrv->SetTitle("Derivative Data;Time [ns];Amplitude [V/ns]");
    graphDrv->SetLineColor(kGreen);
    graphDrv->SetLineWidth(2);

    //calculate the cumulative distribution of the waveform data
    cdfValues = CDF(dataValues);
    SmoothArray(data, smoothdata, npoints,1500, 1);
    cdfValues_smooth150 = CDF(smoothdataVec);
    cdf = &cdfValues[0];
    cdf_smooth150 = &cdfValues_smooth150[0];

    SmoothArray(cdf, smoothCFD_11, npoints,11, 1);
    DerivateArray(smoothCFD_11, drvsmoothCFD_11_5, npoints, dt,5);
    for (int i=0; i<npoints; i++) {
        drvsmoothCFD_11_5[i]*=-10;
    }

    // SmoothArray(data, smoothdata_50, npoints,50, 1);
    // SmoothArray(smoothdata_50, smoothdata_150, npoints,150, 1);

    DerivateArray(cdf, drvCFD, npoints, dt,3); //Derivative of the CDF
    vector<double> drvCDFValues(npoints, 0.0);

    //convert drvCDFValues vector to array
    double* cdfDrv = &drvCDFValues[0];

    // Determine min and max for y-axis range
    double maxCdf = *max_element(cdf, cdf + npoints);

    double globalMaxCDF = max(maxData, maxCdf) + 1;

    TGraph *graphCDF = new TGraph(npoints, xValues, &cdfValues[0]);
    graphCDF->SetTitle("Data CDF;Time [ns];CDF");
    graphCDF->SetLineColor(kBlue);
    graphCDF->SetLineWidth(5);

    TGraph *graphCDF_smooth150 = new TGraph(npoints, xValues, &cdfValues_smooth150[0]);
    graphCDF_smooth150->SetTitle(" Smooth 150ns Data CDF ;Time [ns];CDF");
    graphCDF_smooth150->SetLineColor(kMagenta+3);
    graphCDF_smooth150->SetLineWidth(2);

    //plot the derivative of the CDF
    TGraph *graphDrvCDF = new TGraph(npoints, xValues, cdfDrv);
    graphDrvCDF->SetTitle("Derivative of the Cumulative Distribution Function;Time [ns];dCDF/dt");
    graphDrvCDF->SetLineColor(kGreen);
    graphDrvCDF->SetLineWidth(2);

    TGraph *graphSmoothData = new TGraph(npoints, xValues, smoothdata);
    graphSmoothData->SetTitle("Smooth Data;Time [ns];dCDF/dt");
    graphSmoothData->SetLineColor(kViolet+1);
    graphSmoothData->SetLineWidth(2);

    TGraph *graphDrvCDF2 = new TGraph(npoints, xValues, drvsmoothCFD_11_5);
    graphDrvCDF2->SetTitle("Derivative of the Data CDF Smoothed 11 pts;Time [ns];dCDF/dt");
    graphDrvCDF2->SetLineColor(kMagenta);
    graphDrvCDF2->SetLineWidth(2);

    c1->cd(1);
    graphData->GetYaxis()->SetRangeUser(globalMin, globalMaxCDF); // Adjust Y-axis
    graphData->Draw("AL");  // Original data
    // graphCDF->Draw("AL SAME");  // CDF data
    // graphDrvCDF->Draw("AL SAME"); // Derivative of the CDF data
    // graphDrv->Draw("AL"); // Derivative data
    graphCDF->Draw("L SAME");  // CDF data
    graphDrvCDF2->Draw("L SAME"); // Derivative of the CDF data
    graphDrvCDF->Draw("L SAME"); // Derivative of the CDF data
    graphCDF_smooth150->Draw("L SAME");  // Smooth of the CDF data

    // graphSmoothData->Draw("L SAME"); // Smooth of the data
    //graphDrv->Draw("L SAME"); // Derivative data


    // Pad 2: Integral Curves
    c1->cd(2);


    // Track min and max for integral data to adjust axes
    double minIntegral = 1e9, maxIntegral = -1e9;

    for (int t = 0; t < nTint; t++) {
        double tint = tintValues[t];
		int npt = (int)TMath::FloorNint(tint / dt);
        int npt_ion_tail = (int)TMath::FloorNint(ion_tail_width / dt);

        //double threshold_npt = threshold * sqrt(npt);
        //double ion_tail_end_point_threshold_npt = ion_tail_end_point_threshold * sqrt(npt);
        // Perform integration
        cout<<"Integration npoints = "<< npt<<endl;


        vector<pair<double, double>> trigger_windows = GetTriggerWindows(xValues, npoints, data, dt,
            threshold);

        // Extract integral data from GetTriggerWindows
        vector<double> integrated_data_vec(data, data + npoints);
        vector<double> integrated_data_ion_tail_vec(data, data + npoints);

        auto [x_int, y_int] = IntegratePulse_std(timeValues, integrated_data_vec, npt);
        auto [x_int_ion_tail, y_int_ion_tail] = IntegratePulse_std(timeValues, integrated_data_ion_tail_vec, npt_ion_tail);

        //extract integration threshold
        //Call function to plot integral with bounds
        PlotIntegralWithBounds(x_int, y_int, trigger_windows, minIntegral, maxIntegral, t, c1, std::vector<int>(colors, colors + sizeof(colors) / sizeof(colors[0])), tint);

        //plot trigger windows on the original data
        c1->cd(1);
        for(size_t i = 0; i < trigger_windows.size(); i++) {
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
    int npt = (int)TMath::FloorNint(epeak_width/ dt);

    auto [x_int, y_int] = IntegratePulse_std(timeValues, integrated_data_vec, npt);
    //access the x,y_int vector

    vector<double> integrated_data_ion_tail_vec(data, data + npoints);
    int npt_ion_tail = (int)TMath::FloorNint(ion_tail_width / dt);

    auto [x_int_ion_tail, y_int_ion_tail] = IntegratePulse_std(timeValues, integrated_data_ion_tail_vec, npt_ion_tail);


    int int_secondary_points = static_cast<int>(CIVIDEC_PEAK_DURATION/ dt);
    //integrate to peak width strictly to get secondary pulses
    auto [x_int_sec, y_int_sec] = IntegratePulse_std(timeValues, dataValues, int_secondary_points);
    //derivate the integral to get the derivative
    auto[x_der_sec, y_der_sec] = DerivatePulse_std(x_int_sec, y_int_sec);
    //smooth the derivative to remove noise
    // auto [x_der_sec_smooth, y_der_sec_smooth] = SmoothPulse_std(x_der_sec, y_der_sec, int_secondary_points);
    auto [x_der_sec_int, y_der_sec_int] = IntegratePulse_std(x_der_sec, y_der_sec, int_secondary_points);

    // Plot original data
    TGraph *graph_data = new TGraph(npoints, xValues, data);
    graph_data->SetTitle("Waveform Data;Time [ns];Amplitude [V]");
    graph_data->SetLineColor(kBlack);
    graph_data->SetLineWidth(2);
    graph_data->Draw("L SAME");

    //plot the integral of the data with ion tail integration
    TGraph *graph_int_ion_tail = new TGraph(x_int_ion_tail.size(), x_int_ion_tail.data(), y_int_ion_tail.data());
    graph_int_ion_tail->SetTitle("Integral of the Data with Ion Tail;Time [ns];Integral [V]");
    graph_int_ion_tail->SetLineColor(kViolet);
    graph_int_ion_tail->SetLineWidth(2);
    graph_int_ion_tail->Draw("L SAME");

    // Make x_der_sec_int secondary waveform graph and plot it
    TGraph *graph_int_sec = new TGraph(x_der_sec_int.size(), x_der_sec_int.data(), y_der_sec_int.data());
    graph_int_sec->SetTitle("Secondary Integral of the Data;Time [ns];Integral [V]");
    graph_int_sec->SetLineColor(kOrange);
    graph_int_sec->SetLineWidth(2);
    // graph_int_sec->Draw("L SAME");

    double secondary_pulse_threshold_fraction = 0.5; // Threshold for secondary pulses
    // Find the minimum of y_der_sec_int and plot a line at that value times the threshold fraction
    double secondary_pulse_threshold = *min_element(y_der_sec_int.begin(), y_der_sec_int.end()) * secondary_pulse_threshold_fraction;
    cout << "Secondary Pulse Threshold: " << secondary_pulse_threshold << endl;

    // cdfValues_int = CDF(integrated_data_vec);
    cdfValues_int = CDF(y_int);
    cdf_int = &cdfValues_int[0];

    SmoothArray(cdf_int,cdf_int_smoothed_173, npoints,173, 1);
    double integration_threshold = threshold * sqrt(epeak_width / dt);

    TGraph *graphCDF_int = new TGraph(npoints, xValues, &cdfValues_int[0]);
    graphCDF_int->SetTitle("CDF of int 5ns ;Time [ns];CDF");
    graphCDF_int->SetLineColor(kRed);
    graphCDF_int->SetLineWidth(2);

    TGraph *graphCDF_int_smoothed = new TGraph(npoints, xValues, &cdf_int_smoothed_173[0]);
    graphCDF_int_smoothed->SetTitle(" Smoothed 173 pts CDF of integral 5ns;Time [ns];CDF");
    graphCDF_int_smoothed->SetLineColor(kOrange);
    graphCDF_int_smoothed->SetLineWidth(2);

    double end_thresh_ion_tail = threshold * sqrt(ion_tail_width/dt)* ion_tail_end_point_threshold_fraction;

    gPad->BuildLegend(0.7, 0.7, 0.9, 0.9);
    graphDrv->Draw("L SAME"); // Derivative data
    //draw a line for the threshold
    TLine *y_line_thres = new TLine(0, threshold, npoints * dt, threshold);
    TLine *y_integration_threshold = new TLine(0, integration_threshold, npoints * dt, integration_threshold);
    TLine *y_int_ion_tail_threshold = new TLine(0, end_thresh_ion_tail, npoints * dt, end_thresh_ion_tail);
    TLine *y_line_sec_thres = new TLine(0, secondary_pulse_threshold, npoints * dt, secondary_pulse_threshold);
    y_line_thres->SetLineColor(kBlue);
    y_integration_threshold->SetLineColor(kViolet);
    y_int_ion_tail_threshold->SetLineColor(kGreen);
    y_line_sec_thres->SetLineColor(kOrange);
    y_line_thres->SetLineStyle(2);
    y_integration_threshold->SetLineStyle(2);
    y_line_sec_thres->SetLineStyle(2);
    y_line_thres->Draw("SAME");
    y_integration_threshold->Draw("SAME");
    y_int_ion_tail_threshold->Draw("SAME");
    y_line_sec_thres->Draw("SAME");
// Save canvas
    c1->SaveAs("IntegrationAndOriginalData.pdf");

    TCanvas* c2 = new TCanvas("c2", "CDFs", 800, 800);
    c2->Divide(2, 1);

    c2->cd(1);
    graphCDF->Draw("AL");  // CDF data
    graphCDF_int->Draw("L SAME");  // CDF integrated data
    graphCDF_int_smoothed->Draw("L SAME");  // CDF int data smoothed
    graphCDF_smooth150->Draw("L SAME");  // Smooth of the CDF data
    gPad->BuildLegend(0.7, 0.7, 0.9, 0.9);


    c2->cd(2);
    graphCDF->Draw("AL");  // CDF data
    graphCDF_int->Draw("L SAME");  // CDF integrated data
    graphCDF_int_smoothed->Draw("L SAME");  // CDF int data smoothed 170pts

    gPad->BuildLegend(0.7, 0.7, 0.9, 0.9);

    c2->Update();

    c1->cd(1);
    graphCDF_int_smoothed->Draw("L SAME");  // CDF int data smoothed 170pts
    gPad->BuildLegend(0.7, 0.7, 0.9, 0.9);
    c1->Update();
    return 0;
}

