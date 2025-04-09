//
// Created by akallits on 06/04/25.
//

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "RMSBaselineCalculator_plotting.cpp"


int rms_vs_integration_points() {
    std::string file_path = "/home/akallits/Documents/PicoAnalysis/Saclay_Analysis/data/2022_October_h4/processedTrees/Run224-Pool2_TESTBEAM_tree.root";

    vector<double> integration_times = {0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0};
    double dt = 0.1; //ns
     vector<double>  integration_points;
    for (double& integration_time : integration_times) {
        integration_points.push_back(integration_time / dt); // Convert to ns
    }
    vector<double> integration_rmses_c4;
    vector<double> integration_rmses_c2;
    vector<double> rmses_c4;
    vector<double> rmses_c2;

    ULong64_t secondEpoch_c4;
    ULong64_t secondEpoch_c2;
    for (double integration_time : integration_times) {
         cout<<"Integration time: " << integration_time << " ns" << endl;
        RMSBaselineCalculator calculator_c4(file_path, "RawDataTree", 3, 100, 3.0, 5, integration_time);
        calculator_c4.Process();
        RMSBaselineCalculator calculator_c2(file_path, "RawDataTree", 1, 100, 3.0, 5, integration_time);
        calculator_c2.Process();
        vector<ULong64_t> epochs_c4 = calculator_c4.get_epochs();
        vector<ULong64_t> epochs_c2 = calculator_c2.get_epochs();
        if (epochs_c4.size() > 1) {
            secondEpoch_c4 = epochs_c4[1];
            secondEpoch_c2 = epochs_c2[1];
            std::cout << "Second epoch C4: " << secondEpoch_c4 << std::endl;
            std::cout << "Second epoch C2: " << secondEpoch_c2 << std::endl;
        }

        double rms_integral_c4 = calculator_c4.get_epoch_integral_rms(secondEpoch_c4);
        double rms_integral_c2 = calculator_c2.get_epoch_integral_rms(secondEpoch_c2);
        double rms_c4 = calculator_c4.get_epoch_rms(secondEpoch_c4);
        double rms_c2 = calculator_c2.get_epoch_rms(secondEpoch_c2);
        integration_rmses_c4.push_back(rms_integral_c4);
        integration_rmses_c2.push_back(rms_integral_c2);
        rmses_c4.push_back(rms_c4);
        rmses_c2.push_back(rms_c2);
    }

    // Plot integration rmses vs integration times.
    TGraph* graph = new TGraph(integration_points.size(), &integration_points[0], &integration_rmses_c4[0]);
    graph->SetTitle("Integration RMS vs Integration Time");
    graph->GetXaxis()->SetTitle("Integration Points");
    graph->GetYaxis()->SetTitle("RMS");
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.2);
    graph->Draw("AP");

    TGraph* graph_c2 = new TGraph(integration_points.size(), &integration_points[0], &integration_rmses_c2[0]);
    graph_c2->SetMarkerStyle(20);
    graph_c2->SetMarkerSize(1.2);
    graph_c2->SetMarkerColor(kBlue);
    graph_c2->Draw("P SAME");


    std::vector<double> sqrt_points_c4;
    std::vector<double> sqrt_points_c2;
    for (double pt : integration_points) {
        sqrt_points_c4.push_back(std::sqrt(pt)*rmses_c4[0]);
        sqrt_points_c2.push_back(std::sqrt(pt)*rmses_c2[0]);
    }

    TGraph* sqrtGraph = new TGraph(integration_points.size(), &integration_points[0], &sqrt_points_c4[0]);
    sqrtGraph->SetLineColor(kRed);
    sqrtGraph->SetLineWidth(2);
    sqrtGraph->SetLineStyle(2); // dashed line
    sqrtGraph->Draw("L SAME"); // Draw line on same canvas

    TGraph* sqrtGraph_c2 = new TGraph(integration_points.size(), &integration_points[0], &sqrt_points_c2[0]);
    sqrtGraph_c2->SetLineColor(kBlue);
    sqrtGraph_c2->SetLineWidth(2);
    sqrtGraph_c2->SetLineStyle(2); // dashed line
    sqrtGraph_c2->Draw("L SAME"); // Draw line on same canvas





    return 0;
}