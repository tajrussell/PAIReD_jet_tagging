import uproot
import awkward as ak
import numpy as np
import os

# Directory containing the ROOT files
directory = "DY"

# Loop over all ROOT files in the directory
for filename in os.listdir(directory):
    if filename.endswith(".root"):
        file_path = os.path.join(directory, filename)
        
        with uproot.open(file_path) as file:
            tree = file["tree"]  # Assuming the TTree is named 'tree'
            
            # Load the branches as awkward arrays
            part_px = tree["part_px"].array(library="ak")
            part_py = tree["part_py"].array(library="ak")
            part_pz = tree["part_pz"].array(library="ak")
            part_energy = tree["part_energy"].array(library="ak")
            part_mass = tree["part_mass"].array(library="ak")
            
            # Calculate recm^2 and part_mass^2 using awkward arrays
            recm2 = part_energy**2 - part_px**2 - part_py**2 - part_pz**2
            part_mass2 = part_mass**2
            
            # Find entries where recm^2 is negative
            negative_mask = recm2 < 0
            negative_recm2 = recm2[negative_mask]
            negative_mass2 = part_mass2[negative_mask]
            
            # Calculate and print the portion of negative recm^2 entries
            portion_negative = ak.count(negative_recm2) / ak.count(recm2) if ak.count(recm2) > 0 else 0
            print("File: %s" % filename)
            print("Portion of entries with negative recm^2: %.4f" % portion_negative)
            
            # Print the arrays for the negative recm^2 entries
            print("Negative recm^2 values:")
            print(negative_recm2)
            print("Corresponding part_mass^2 values:")
            print(negative_mass2)

