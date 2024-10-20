//
// Created by Dylan on 10/20/2024.
//


#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include "Event.h"

#define NUM_EVENTS_PER_ENTRY 10  // Number of events per entry

void read_class_arrays_from_tree() {
    // Open the ROOT file
    TFile *file = TFile::Open("events_class_arrays.root");

    // Retrieve the TTree
    TTree *tree = (TTree*)file->Get("eventTree");

    // Define an array of Event objects to hold the data
    Event events[NUM_EVENTS_PER_ENTRY];

    // Set the branch address
    tree->SetBranchAddress("events", events);

    // Create histograms for the different variables
    TH1F *hAmplitude = new TH1F("hAmplitude", "Amplitude", 100, 0, 10);
    TH1F *hStartTime = new TH1F("hStartTime", "Start Time", 100, 0, 100);
    TH1F *hCharge = new TH1F("hCharge", "Charge", 50, 0, 10);

    // Loop over the entries in the tree
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);  // Get the entry data

        // Loop over each event in the array
        for (int j = 0; j < NUM_EVENTS_PER_ENTRY; ++j) {
            hAmplitude->Fill(events[j].amplitude);
            hStartTime->Fill(events[j].start_time);
            hCharge->Fill(events[j].charge);
        }
    }

    // Create a canvas to display the histograms
    TCanvas *canvas = new TCanvas("canvas", "Event Histograms", 1200, 400);
    canvas->Divide(3, 1);  // Divide into 3 pads

    // Draw the amplitude histogram
    canvas->cd(1);
    hAmplitude->Draw();

    // Draw the start time histogram
    canvas->cd(2);
    hStartTime->Draw();

    // Draw the charge histogram
    canvas->cd(3);
    hCharge->Draw();

    // Save the canvas as an image
    canvas->SaveAs("event_histograms.png");

    // Clean up
    delete hAmplitude;
    delete hStartTime;
    delete hCharge;
    delete canvas;
    file->Close();
    delete file;
}