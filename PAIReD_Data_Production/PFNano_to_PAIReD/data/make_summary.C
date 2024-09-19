void make_summary() {
    // Create output file for histograms
    TFile *outputFile = new TFile("summary.root", "RECREATE");

    // Directories and prefixes
    std::vector<std::string> directories = {"DY", "TT"};
    std::string combinedDirPrefix = "XtoHH";  // Prefix for combined XtoHH directories

    // Variables to create histograms for
    std::vector<std::string> variables = {
        "MC_higgs_pt", "MC_higgs_eta", "MC_higgs_phi", "MC_higgs_mass", "MC_higgs_flav", 
        "MC_vector_flav", "MC_lepton_channel", "MC_gendijet_pt", "MC_gendijet_eta", 
        "MC_gendijet_phi", "MC_gendijet_mass", "MC_genjet1_flav", "MC_genjet2_flav", 
        "MC_genjet1_matched", "MC_genjet2_matched", "MC_drqq", "MC_n_c"
    };

    // Labels for additional cuts
    std::vector<std::string> labels = {
        "label_BB", "label_CC", "label_bb", "label_bx", "label_cx", "label_ll"
    };

    // Loop through DY and TT directories
    for (const std::string &dir : directories) {
        TChain chain("tree");

        // Add all files using wildcard
        chain.Add(Form("%s/PAIReD_*.root", dir.c_str()));

        std::string prefix = dir;

        // Loop over variables and create histograms
        for (const std::string &var : variables) {
            std::cout<<dir<<"  "<<var<<std::endl;
            // With zeros
            chain.Draw((var + ">>" + prefix + "_" + var + "_with_zeros").c_str());

            // Without zeros, only MC_higgs_valid == 1
            chain.Draw((var + ">>" + prefix + "_" + var).c_str(), ("MC_higgs_valid == 1 && " + var + " > 0").c_str());

            // Create histograms for each label
            for (const std::string &label : labels) {
                // With zeros for each label
                chain.Draw((var + ">>" + prefix + "_" + var + "_with_zeros_" + label).c_str(), (label + " == 1").c_str());

                // Without zeros, MC_higgs_valid == 1 and label == 1
                chain.Draw((var + ">>" + prefix + "_" + var + "_" + label).c_str(), 
                            ("MC_higgs_valid == 1 && " + var + " > 0 && " + label + " == 1").c_str());
            }
        }
    }

    // Handle combined XtoHH_500-1000 and XtoHH_1000-4000 directories
    TChain chainXtoHH("tree");
    chainXtoHH.Add("XtoHH_500-1000/PAIReD_*.root");
    chainXtoHH.Add("XtoHH_1000-4000/PAIReD_*.root");

    // Loop over variables and create histograms for combined XtoHH
    for (const std::string &var : variables) {
        // With zeros
        chainXtoHH.Draw((var + ">>" + combinedDirPrefix + "_" + var + "_with_zeros").c_str());

        // Without zeros, only MC_higgs_valid == 1
        chainXtoHH.Draw((var + ">>" + combinedDirPrefix + "_" + var).c_str(), ("MC_higgs_valid == 1 && " + var + " > 0").c_str());

        // Create histograms for each label
        for (const std::string &label : labels) {
            // With zeros for each label
            chainXtoHH.Draw((var + ">>" + combinedDirPrefix + "_" + var + "_with_zeros_" + label).c_str(), (label + " == 1").c_str());

            // Without zeros, MC_higgs_valid == 1 and label == 1
            chainXtoHH.Draw((var + ">>" + combinedDirPrefix + "_" + var + "_" + label).c_str(), 
                            ("MC_higgs_valid == 1 && " + var + " > 0 && " + label + " == 1").c_str());
        }
    }

    // Write all histograms to the output file
    outputFile->Write();
    outputFile->Close();
}

