#include <TFile.h>
#include <TTree.h>
#include "Event.h"

void write_class_to_tree() {
    // Create a ROOT file to store the tree
    TFile *file = new TFile("events_class.root", "RECREATE");

    // Create a TTree
    TTree *tree = new TTree("eventTree", "Tree storing Event class objects");

    // Create an instance of the class
    Event *event = new Event();

    // Create a branch for the class object
    tree->Branch("event", &event);

    // Fill the tree with mock events
    for (int i = 0; i < 100; ++i) {
        event->SetValues(static_cast<float>(rand() % 100) / 10.0,    // Random amplitude
                         static_cast<float>(rand() % 1000) / 10.0,   // Random start time
                         static_cast<float>(rand() % 50) / 5.0);     // Random charge

        // Fill the tree with the event object
        tree->Fill();
    }

    // Write the tree to the file
    tree->Write();

    // Close the file
    file->Close();

    // Clean up
    delete file;
    delete event;
}

