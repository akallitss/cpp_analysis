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
	//double data[npoints];     // Populate with data values
	//double integral[npoints]; // To store integral values
    double* data = &dataValues[0];
    double* drv = &drvValues[0];

    //print drv values
    // for (int i=0; i<drvValues.size(); i++){
    //     cout << "drv[" << i << "]: " << drvValues[i] << endl;
    // }
    //
    // cin.get();


    //threshold = -6 * rms; // Convert threshold to units of RMS
    //double ion_tail_end_point_threshold = 0.2 * threshold; // Set ion tail end point fraction to 10% of threshold

    //print the threshold
    // cout<<"Threshold = "<<threshold<<endl;
    // cin.get();

    //// double tintValues[] = {1.0, 2.0, 5.0, 10.0};
    // double tintValues[] = {1.0, 2.0, 5.0, 10.0, 100, 150, 200};
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

	// double xValues[npoints];
	// for (int i = 0; i < npoints; i++) {
    // 		xValues[i] = i * dt; // Replace with your time data if reading from a file
	// }

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
    graphData->SetTitle("Original Data;Time [ns];Amplitude [V]");
    graphData->SetLineColor(kBlack);
    graphData->SetLineWidth(2);
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


    graphData->GetYaxis()->SetRangeUser(globalMin, globalMax); // Adjust Y-axis
    graphData->Draw("AL");  // Original data
    graphDrv->Draw("L SAME");  // Derivative data


    // Add a legend
    TLegend *legend1 = new TLegend(0.7, 0.7, 0.9, 0.85);
    legend1->AddEntry(graphData, "Original Data", "l");
    legend1->AddEntry(graphDrv, "Derivative Data", "l");
    legend1->Draw();

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

        //double threshold_npt = threshold * sqrt(npt);
        //double ion_tail_end_point_threshold_npt = ion_tail_end_point_threshold * sqrt(npt);
        // Perform integration
        cout<<"npoints = "<< npt<<endl;
        // IntegratePulse(npoints, data, integral, dt, tint);
        // vector<double> t_values(npoints);
        // for (int i = 0; i < npoints; i++) {
        //     t_values[i] = i * dt;
        // }

        //convert time vector to array
        // double t_values[npoints];
        // for (int i = 0; i < npoints; i++) {
        //     t_values[i] = i * dt;
        // }

        //convert xValues to timeValues
        //double* t_values = &timeValues[0];

        // cout << "Size of dataValues: " << dataValues.size() << endl;
        // cout << "Size of timeValues: " << timeValues.size() << endl;

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



        //Print the x_int and y_int vectors
        // cout << "Size of x_int: " << x_int.size() << endl;
        // cout << "Size of y_int: " << y_int.size() << endl;

        // for(int i =0; i < 10; i++){
            // cout << "x_int[" << i << "]: " << x_int[i] << ", y_int[" << i << "]: " << y_int[i] << endl;
        // }
        // cin.get();
        // Plot Integral with Bounds

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

    // cout << "minIntegral: " << minIntegral << endl;
    // cout << "maxIntegral: " << maxIntegral << endl;
    // cin.get(); //Pause to see the values

    // Draw legend
    legend->Draw();

    // Save canvas
    c1->SaveAs("IntegrationAndOriginalData.pdf");




    return 0;
}
/////backup from the for loop ////////
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

        // // Set grap

        // if (tint == 5.0) {
        //     c1->cd(1);
        //     // Draw vertical lines to indicate the pulse bounds
        //     for (size_t i = 0; i < pulse_bounds.size(); i++) {
        //         double x_left = pulse_bounds[i].first;
        //         double x_right = pulse_bounds[i].second;
        //         TLine *line_left = new TLine(x_left, minData, x_left, maxData);
        //         TLine *line_right = new TLine(x_right, minData, x_right, maxData);
        //         line_left->SetLineStyle(2);
        //         line_right->SetLineStyle(2);
        //         line_left->SetLineColor(colors[t % 5]);
        //         line_right->SetLineColor(colors[t % 5]);
        //         line_left->Draw();
        //         line_right->Draw();
        //     }
        //     c1->cd(2);
        // }


    //     double localMin = *min_element(y_int.begin(), y_int.end());
    //     double localMax = *max_element(y_int.begin(), y_int.end());
    //     minIntegral = min(minIntegral, localMin);
    //     maxIntegral = max(maxIntegral, localMax);
    //
    //
    //
    //
    //     TGraph *graphIntegral = new TGraph(x_int.size(), &x_int[0], &y_int[0]);
    //     for (int i = 0; i < x_int.size(); i++) {
    //         //graphIntegral->SetPoint(i, x_int[i], y_int[i]);
    //         graphIntegral->SetPoint(i, x_int[i], y_int[i]);
    //     }
    //
    //     // Set graph style
    //     graphIntegral->SetLineColor(colors[t % 5]);
    //     graphIntegral->SetLineWidth(2);
    //     graphIntegral->SetTitle(";Time [ns];Integral [V*ns]");
    //
    //     // Add graph to canvas
    //     if (t == 0) {
    //         graphIntegral->Draw("AL"); // Draw the first graph with axes
    //     } else {
    //         graphIntegral->Draw("L"); // Overlay subsequent graphs
    //     }
    //
    //     // Add entry to legend
    //     legend->AddEntry(graphIntegral, Form("Tint = %.1f ns", tint), "l");
    //
    //     //Add vertical lines to indicate the pulse bounds
    //     for (size_t i = 0; i < pulse_bounds.size(); i++) {
    //         double x_left = pulse_bounds[i].first;
    //         double x_right = pulse_bounds[i].second;
    //         TLine *line_left = new TLine(x_left, minIntegral, x_left, maxIntegral);
    //         TLine *line_right = new TLine(x_right, minIntegral, x_right, maxIntegral);
    //         line_left->SetLineStyle(2);
    //         line_right->SetLineStyle(2);
    //         line_left->SetLineColor(colors[t % 5]);
    //         line_right->SetLineColor(colors[t % 5]);
    //         line_left->Draw();
    //         line_right->Draw();
    //     }
    //
    //
////
// // Adjust axis ranges for integral plot
// c1->cd(2);
// gPad->Modified();
// gPad->Update();
// TGraph *g = (TGraph*) gPad->FindObject("Graph"); // Grab the first graph
// if (g) {
//     g->GetXaxis()->SetRangeUser(0, npoints * dt);
//     g->GetYaxis()->SetRangeUser(minIntegral - 1, maxIntegral + 1);
// }
//
//
// // Draw legend
// legend->Draw();
//
// // Save canvas
// c1->SaveAs("IntegrationAndOriginalData.pdf");
