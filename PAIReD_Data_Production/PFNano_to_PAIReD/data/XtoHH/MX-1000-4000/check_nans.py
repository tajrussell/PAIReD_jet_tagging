import uproot
import numpy as np

def check_nans():
    # Open the ROOT file and get the tree
    file = uproot.open("PAIReD_72.root")
    tree = file["tree"]

    # Iterate over all branches in the tree
    for branch_name, branch in tree.items():
        try:
            # Read the branch data as a NumPy array
            data = branch.array(library="np")

            # Check if the branch contains vectors (arrays of arrays)
            if isinstance(data[0], np.ndarray):
                # If the branch contains vectors, flatten them to handle NaNs
                flattened_data = np.concatenate(data)
                nan_count = np.sum(np.isnan(flattened_data))
            else:
                # For scalar branches
                nan_count = np.sum(np.isnan(data))

            print("Branch: %s has %d NaN values." % (branch_name, nan_count))

        except Exception as e:
            print("Could not check branch %s: %s" % (branch_name, str(e)))

# Run the function
check_nans()

