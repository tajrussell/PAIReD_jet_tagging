import uproot
import numpy as np
import glob

# Define directories
directories = ["DY", "XtoHH/MX-500-1000", "XtoHH/MX-1000-4000"]

# Initialize an empty list to collect dijet_pt values
dijet_pt_values = []

# Loop through each directory and read matching files
for directory in directories:
    # Construct the file pattern
    file_pattern = f"{directory}/PAIReD_*.root"
    files = glob.glob(file_pattern)

    for file in files:
        with uproot.open(file) as root_file:
            # Access the tree, assuming the tree name is known, e.g., "tree"
            tree = root_file["tree"]

            # Read the "dijet_pt" branch
            dijet_pt = tree["dijet_pt"].array(library="np")

            # Collect the values
            dijet_pt_values.extend(dijet_pt)

# Convert to a numpy array
dijet_pt_values = np.array(dijet_pt_values)

# Calculate the 10 quantiles
quantiles = np.percentile(dijet_pt_values, np.linspace(0, 100, 11))

print("10 quantiles of dijet_pt:")
print(quantiles)
