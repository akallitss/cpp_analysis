#include <TFile.h>
#include <TTree.h>

// Define the C struct
struct Event {
    float amplitude;
    float start_time;
    float charge;
};

void save_struct_to_tree() {
    // Create a ROOT file to store the tree
    TFile *file = new TFile("events.root", "RECREATE");

    // Create a TTree
    TTree *tree = new TTree("eventTree", "Tree storing Event struct");

    // Create an instance of the struct
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
