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

using namespace std;

// IntegratePulse function
double IntegratePulse(int npoints, const double* data, double* integral, double dt, double tint) {
    double sum = data[0] * dt;

    if (npoints < 2) return sum;

    int nint = (int) TMath::FloorNint(tint / dt); // number of points in the time window
    integral[0] = data[0] * dt;
    if (nint < 2) nint = 2;

    for (int i = 1; i < nint; i++) {
        integral[i] = data[i] * dt + integral[i - 1];
        sum += data[i] * dt;
    }


   for (int i = nint; i < npoints; i++) {
        integral[i] = data[i] * dt + integral[i - 1] - data[i - nint] * dt;
        sum += data[i] * dt;
    }

    return sum;

}

// Function to read data from the text file
bool ReadWaveformData(const string& filename, double& dt, vector<double>& dataValues) {
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
        }
    }



    infile.close();
    return true;
}

int IntegratePulse() {
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
    // string filename ="/home/akallits/Documents/PicoAnalysis/Saclay_Analysis/data/2022_October_h4/plots/Run224/Pool2/moreplots/waveform_data.txt" ;

    string filename ="/sw/akallits/PicoAnalysis/Saclay_Analysis/data/2022_October_h4/plots/Run224/Pool2/moreplots/waveform_data.txt" ;
    // Variables to hold dt and data values
    double dt = 0.0;

    vector<double> dataValues;

    // Read the waveform data
    if (!ReadWaveformData(filename, dt, dataValues)) {
        return -1; // Exit if file reading fails
    }

    // Convert dataValues vector to array for compatibility with ROOT
    int npoints = dataValues.size();
	//double data[npoints];     // Populate with data values
	//double integral[npoints]; // To store integral values
    double* data = &dataValues[0];
	double x_values[npoints];
    double x_Avg[npoints];
    // Array of tint values
    // double tintValues[] = {1.0,20.0, 50.0, 100.0, 150.0};
    double tintValues[] = {1.0,20.0, 50.0, 100.0, 150.0};
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
        double integral[npoints];
		int npt = (int)TMath::FloorNint(tint / dt);
        // Perform integration
                cout<<"npoints = "<< npt<<endl;
        IntegratePulse(npoints, data, integral, dt, tint);

        // Update min and max for integral data
        double localMin = *min_element(integral, integral + npoints);
        double localMax = *max_element(integral, integral + npoints);
        minIntegral = min(minIntegral, localMin);
        maxIntegral = max(maxIntegral, localMax);

        // Create a graph for this tint
		//double xAvg[npoints];

        //ComputeAveragedX(npoints, xValues, xAvg, tint);
		//TGraph *graphIntegral = new TGraph(npoints, xAvg, integral);

        TGraph *graphIntegral = new TGraph(npoints-npt, &x_values[npt], integral);
        for (int i = 0; i < npoints-npt; i++) {
            graphIntegral->SetPoint(i, i * dt, integral[i+npt]);
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
