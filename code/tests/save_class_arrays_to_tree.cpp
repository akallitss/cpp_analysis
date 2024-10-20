#include <TFile.h>
#include <TTree.h>
#include <TObject.h>
#include "Event.h"

#define EventDef(name,id)

#define NUM_EVENTS_PER_ENTRY 10  // Number of events per entry

void save_class_arrays_to_tree() {
//    gROOT->Load("Event.h");
    // Create a ROOT file to store the tree
    TFile *file = new TFile("events_class_arrays.root", "RECREATE");

    // Create a TTree
    TTree *tree = new TTree("eventTree", "Tree storing arrays of Event objects");

    // Define an array of Event objects
    Event events[NUM_EVENTS_PER_ENTRY];

    // Create a branch for the array of Event objects
    tree->Branch("events", &events, Form("events[%d]/F", NUM_EVENTS_PER_ENTRY));

    // Fill the tree with mock events
    for (int i = 0; i < 100; ++i) {
        for (int j = 0; j < NUM_EVENTS_PER_ENTRY; ++j) {
            events[j].amplitude = static_cast<float>(rand() % 100) / 10.0;  // Random amplitude
            events[j].start_time = static_cast<float>(rand() % 1000) / 10.0; // Random start time
            events[j].charge = static_cast<float>(rand() % 50) / 5.0;        // Random charge
        }

        // Fill the tree with the array of events
        tree->Fill();
    }

    // Write the tree to the file
    tree->Write();

    // Close the file
    file->Close();

    // Clean up
    delete file;
}