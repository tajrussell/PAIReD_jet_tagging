import uproot
import numpy as np
import sys
import glob
import os

def check_nans(file_numbers=None):
    # Define the directories to search for ROOT files
    directories = [
        "DY/PAIReD_%s.root",
        "TT/PAIReD_%s.root",
        "XtoHH/MX-500-1000/PAIReD_%s.root",
        "XtoHH/MX-1000-4000/PAIReD_%s.root"
    ]

    # Iterate over each directory pattern
    for dir_pattern in directories:
        print("Checking files in directory: %s" % os.path.dirname(dir_pattern))
        
        # Create a list of file paths based on provided numbers or gather all files
        if file_numbers:
            file_paths = set(dir_pattern % number for number in file_numbers)
            file_paths = [fp for fp in file_paths if os.path.exists(fp)]
        else:
            file_paths = glob.glob(dir_pattern % "*")

        # Initialize totals for each directory
        total_nan_count = 0
        total_negative_count = 0
        events_with_conditions = 0
        
        # Dictionary to keep track of NaN and negative counts by branch name
        branch_stats = {}

        # Iterate over all specified or found files
        for file_path in file_paths:
            print("  Checking file: %s" % file_path)

            # Open the ROOT file and get the tree
            try:
                file = uproot.open(file_path)
                tree = file["tree"]
            except Exception as e:
                print("    Could not open file %s: %s" % (file_path, str(e)))
                continue

            # Iterate over all branches in the tree
            for branch_name, branch in tree.items():
                try:
                    # Read the branch data as a NumPy array
                    data = branch.array(library="np")

                    # Check if the branch contains vectors (arrays of arrays)
                    if isinstance(data[0], np.ndarray):
                        # Flatten the vector data
                        flattened_data = np.concatenate(data)
                    else:
                        # Use data as-is for scalar branches
                        flattened_data = data

                    # Count NaNs and negatives
                    nan_count = np.sum(np.isnan(flattened_data))
                    negative_count = np.sum(flattened_data < 0)

                    # Update totals
                    total_nan_count += nan_count
                    total_negative_count += negative_count

                    # Store the counts in the dictionary
                    branch_stats[branch_name] = branch_stats.get(branch_name, {'nans': 0, 'negatives': 0})
                    branch_stats[branch_name]['nans'] += nan_count
                    branch_stats[branch_name]['negatives'] += negative_count

                except Exception as e:
                    print("    Could not check branch %s: %s" % (branch_name, str(e)))

            # Check for the specific condition for MC_higgs_mass, label_BB, and label_CC
            if "MC_higgs_mass" in tree and "label_BB" in tree and "label_CC" in tree:
                mc_higgs_mass = tree["MC_higgs_mass"].array(library="np")
                label_bb = tree["label_BB"].array(library="np")
                label_cc = tree["label_CC"].array(library="np")

                # Count events where MC_higgs_mass is 0 and either label_BB or label_CC is true
                condition_mask = (mc_higgs_mass == 0) & (label_bb | label_cc)
                events_with_conditions += np.sum(condition_mask)

        # Print branch statistics for the directory
        for branch_name, stats in branch_stats.items():
            print("  Branch: %s - Total NaN values: %d, Total Negative values: %d" %
                  (branch_name, stats['nans'], stats['negatives']))

        # Print totals for the directory
        print("Directory: %s - Total NaN values: %d, Total Negative values: %d" %
              (os.path.dirname(dir_pattern), total_nan_count, total_negative_count))

        # Print total counts for the specific condition for the directory
        print("Total events with MC_higgs_mass == 0 and (label_BB or label_CC) is true: %d" %
              events_with_conditions)

# Entry point to run the function
if __name__ == "__main__":
    # Check for command line arguments for file numbers
    if len(sys.argv) > 1:
        file_numbers = sys.argv[1:]
    else:
        file_numbers = None

    # Convert input numbers to string format
    if file_numbers:
        file_numbers = list(set(str(num) for num in file_numbers))

    check_nans(file_numbers)

