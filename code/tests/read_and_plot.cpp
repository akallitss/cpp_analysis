#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>

// Define the same C struct
struct Event {
    float amplitude;
    float start_time;
    float charge;
};

void read_and_plot() {
    // Open the ROOT file containing the tree
    TFile *file = new TFile("events.root", "READ");

    // Get the tree from the file
    TTree *tree = (TTree*)file->Get("eventTree");

    // Create an instance of the struct to hold the data
    Event event;

    // Set the branch address to point to the struct
    tree->SetBranchAddress("event", &event);

    // Create histograms for each variable
    TH1F *h_amplitude = new TH1F("h_amplitude", "Amplitude;Amplitude;Counts", 50, 0, 10);
    TH1F *h_start_time = new TH1F("h_start_time", "Start Time;Start Time;Counts", 50, 0, 100);
    TH1F *h_charge = new TH1F("h_charge", "Charge;Charge;Counts", 50, 0, 10);

    // Loop over all entries in the tree
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

	cout << event.amplitude << endl;

        // Fill histograms with data from the struct
        h_amplitude->Fill(event.amplitude);
        h_start_time->Fill(event.start_time);
        h_charge->Fill(event.charge);
    }

    // Create a canvas to draw the histograms
    TCanvas *c1 = new TCanvas("c1", "Histograms", 900, 600);
    c1->Divide(2, 2);  // Create a 2x2 grid for histograms

    // Draw each histogram in its own pad
    c1->cd(1);
    h_amplitude->Draw();
    
    c1->cd(2);
    h_start_time->Draw();

    c1->cd(3);
    h_charge->Draw();

    // Save the canvas as a PNG
    c1->SaveAs("histograms.png");

    // Clean up
    //delete h_amplitude;
    //delete h_start_time;
    //delete h_charge;
    //file->Close();
    //delete file;
}
