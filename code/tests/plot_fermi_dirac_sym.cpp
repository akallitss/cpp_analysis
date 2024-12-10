//
// Created by akallits on 27/11/24.
//
#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"

double fermi_dirac_sym_double(double *x, double *par) {
    // Function definition
    //    double fdreturn = (par[0] / (1. + TMath::Exp(-par[2] * (x[0] - par[1]))) + par[3]) *
    //                      (1. / (1. + TMath::Exp(-par[5] * (x[0] - par[4]))));
    double fdreturn = (par[0] / (1. + TMath::Exp(-par[2] * (x[0] - par[1])))) *
                      (1. / (1. + TMath::Exp(-par[5] * (x[0] - par[4])))) + par[3];
    return fdreturn;
}

void plot_fermi_dirac_sym() {
    // Create a canvas
    TCanvas *c1 = new TCanvas("c1", "Fermi-Dirac Symmetric Function", 800, 600);

    // Define the function and its parameter ranges
    TF1 *fd_func = new TF1("fd_func", fermi_dirac_sym_double, 220, 245, 6);

    // Set parameter names and initial values
    fd_func->SetParName(0, "A");
    fd_func->SetParName(1, "x0_1");
    fd_func->SetParName(2, "beta_1");
    fd_func->SetParName(3, "C");
    fd_func->SetParName(4, "x0_2");
    fd_func->SetParName(5, "beta_2");

    // Hardcoded parameter values
    fd_func->SetParameter(0, -4174.0);  // A
    fd_func->SetParameter(1, 226.4);  // x0_1
    fd_func->SetParameter(2, 4.99);  // beta_1
    fd_func->SetParameter(3, 48.02);  // C
    fd_func->SetParameter(4, 209.1);  // x0_2
    fd_func->SetParameter(5, -0.5946);  // beta_2

    // Draw the function
    fd_func->SetTitle("Fermi-Dirac Symmetric Double Function; x; f(x)");
    fd_func->SetLineColor(kBlue);
    fd_func->SetLineWidth(2);
    fd_func->Draw();

    // Update canvas
    c1->Update();
}
