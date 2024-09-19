import os

# Define the directory containing the files
directory = "/isilon/hadoop/store/group/common/bruxhcc/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v2_BTV_Run3_2023_Comm_MINIAODv4/240916_150242/0000/"

# Open the output file in append mode
with open("inputs_TT.txt", "a") as output_file:
    # Loop through each file in the directory
    for filename in os.listdir(directory):
        # Check if the file matches the pattern
        if filename.startswith("MC_defaultAK4_2023_") and filename.endswith(".root"):
            # Extract the part of the filename that corresponds to X
            X = filename.replace("MC_defaultAK4_2023_", "").replace(".root", "")
            # Create the line in the format you described using .format()
            line = "{} PAIReD_{}.root PAIReD_{}".format(filename, X, X)
            # Write the line to the output file
            output_file.write(line + "\n")
