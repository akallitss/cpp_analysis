#include <TFile.h>
#include <TTree.h>
#include <TObject.h>
#include "Event.h"

#define EventDef(name,id)

void save_class_to_tree() {
//    gROOT->Load("Event.h");
    // Create a ROOT file to store the tree
    TFile *file = new TFile("events_class.root", "RECREATE");

    // Create a TTree
    TTree *tree = new TTree("eventTree", "Tree storing Event struct");

    // Create an instance of the class
    Event event;

    // Create a branch for the whole struct
    tree->Branch("event", &event);

    // Fill the tree with mock events
    for (int i = 0; i < 100; ++i) {
        event.amplitude = static_cast<float>(rand() % 100) / 10.0;  // Random amplitude
        event.start_time = static_cast<float>(rand() % 1000) / 10.0; // Random start time
        event.charge = static_cast<float>(rand() % 50) / 5.0;        // Random charge

        // Fill the tree with the event
        tree->Fill();
    }

    // Write the tree to the file
    tree->Write();

    // Close the file
    file->Close();

    // Clean up
    delete file;
}
