import os

# Define the directory containing the files
directories = ["/isilon/hadoop/store/group/common/bruxhcc/ZH_Hto2C_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v2_BTV_Run3_2023_Comm_MINIAODv4/250506_010433/0000/",
               "/isilon/hadoop/store/group/common/bruxhcc/ZH_Hto2C_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v2_BTV_Run3_2023_Comm_MINIAODv4/250507_220629/0000/",
               "/isilon/hadoop/store/group/common/bruxhcc/ZH_Hto2C_Zto2Nu_M-125_TuneCP5_13p6TeV_powheg-minlo-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v15-v2_BTV_Run3_2023_Comm_MINIAODv4/250506_010443/0000/",
               "/isilon/hadoop/store/group/common/bruxhcc/ZH_Hto2C_Zto2Nu_M-125_TuneCP5_13p6TeV_powheg-minlo-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v15-v2_BTV_Run3_2023_Comm_MINIAODv4/250507_220639/0000/"]
prefixes = ["L_A_", "L_B_", "Nu_A_", "Nu_B_"]

directories = ["/isilon/hadoop/store/group/common/bruxhcc/ZH_Hto2B_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v2_BTV_Run3_2023_Comm_MINIAODv4/250506_010504/0000/",
               "/isilon/hadoop/store/group/common/bruxhcc/ZH_Hto2B_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v2_BTV_Run3_2023_Comm_MINIAODv4/250507_220659/0000/",
               "/isilon/hadoop/store/group/common/bruxhcc/ZH_Hto2B_Zto2Nu_M-125_TuneCP5_13p6TeV_powheg-minlo-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v15-v2_BTV_Run3_2023_Comm_MINIAODv4/250506_010515/0000/",
               "/isilon/hadoop/store/group/common/bruxhcc/ZH_Hto2B_Zto2Nu_M-125_TuneCP5_13p6TeV_powheg-minlo-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v15-v2_BTV_Run3_2023_Comm_MINIAODv4/250507_220709/0000/"]

directories = ["/isilon/data/common/smondal5/VH_gen/ZHbb_PFNano/",
               "/isilon/data/common/smondal5/VH_gen/ZHcc_PFNano/"
              ]
# Open the output file in append mode
with open("inputs_ZHvar.txt", "w") as output_file:
    # Loop through each file in the directory
    for i, directory in enumerate(directories):
        print(directory)
        for filename in os.listdir(directory):
            # Check if the file matches the pattern
            if filename.startswith("ZH") and filename.endswith(".root"):
                # Extract the part of the filename that corresponds to X
                #X = prefixes[i]+filename.replace("MC_defaultAK4_2023_", "").replace(".root", "")
                X = filename.replace(".root", "")
                # Create the line in the format you described using .format()
                line = "{} PAIReD_{}.root PAIReD_{}".format(directory+filename, X, X)
                # Write the line to the output file
                output_file.write(line + "\n")
