#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <iostream>
#include <vector>

void make_LL_dir(const char* input_dir) {
    // Get a list of files in the directory
    void* dirp = gSystem->OpenDirectory(input_dir);
    const char* file_name;
    
    // Loop over files in the directory
    while ((file_name = gSystem->GetDirEntry(dirp))) {
        // Only process files matching the pattern "PAIReD_*11.root"
        if (TString(file_name).BeginsWith("PAIReD") && TString(file_name).EndsWith(".root")) {
            TString input_file_path = TString::Format("%s/%s", input_dir, file_name);
            TString output_file_path = TString::Format("%s/LL_%s", input_dir, file_name);
            
            // Open the original ROOT file
            TFile* input_file = TFile::Open(input_file_path, "READ");
            if (!input_file || input_file->IsZombie()) {
                std::cerr << "Error: could not open " << input_file_path << std::endl;
                continue;
            }
            
            // Get the TTree named "tree" from the file
            TTree* tree = (TTree*)input_file->Get("tree");
            if (!tree) {
                std::cerr << "Error: could not find tree in " << input_file_path << std::endl;
                input_file->Close();
                continue;
            }

            // Create a new ROOT file to store the modified tree
            TFile* output_file = new TFile(output_file_path, "RECREATE");
            
            // Clone the original tree structure
            TTree* new_tree = tree->CloneTree(0);
            
            Char_t label_BB, label_CC;
            tree->SetBranchAddress("label_BB", &label_BB);
            tree->SetBranchAddress("label_CC", &label_CC);
            
            // Create a new branch for label_LL
            Char_t label_LL;
            TBranch* b_label_LL = new_tree->Branch("label_LL", &label_LL, "label_LL/B");
            
            // Loop over all events and fill the new tree
            Long64_t n_entries = tree->GetEntries();
            for (Long64_t i = 0; i < n_entries; ++i) {
                tree->GetEntry(i);
                label_LL = !(label_BB || label_CC);
                new_tree->Fill();  // Automatically fills label_LL with value 1
            }

            // Write the new tree to the output file
            new_tree->Write();
            output_file->Close();
            input_file->Close();
            
            std::cout << "Processed and saved: " << output_file_path << std::endl;
        }
    }

    gSystem->FreeDirectory(dirp);  // Clean up the directory pointer
}

void make_LL() {
    make_LL_dir("XtoHH/MX-500-1000");
    make_LL_dir("test/XtoHH/MX-500-1000");
    make_LL_dir("XtoHH/MX-1000-4000");
    make_LL_dir("test/XtoHH/MX-1000-4000");
}
