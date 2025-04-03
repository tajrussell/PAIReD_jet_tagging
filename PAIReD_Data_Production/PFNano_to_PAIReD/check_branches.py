import glob
import uproot

def check_branch_exists(file_path, branch_name):
    """Check if a given branch exists in a ROOT file."""
    try:
        with uproot.open(file_path) as file:
            tree = file[file.keys()[0]]  # Assuming the first tree is the relevant one
            return branch_name in tree.keys()
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return False

def main():
    branch_name = "dijet_ak4_mass"
    missing_files = []
    
    # Define paths
    paths = [
        "data/test/XtoHH/*/PAIReD_*.root",
        "data/DY/PAIReD_*.root"
    ]
    
    for path in paths:
        files = glob.glob(path)
        print(f"Checking {len(files)} files in {path}...")
        
        for file_path in files:
            if not check_branch_exists(file_path, branch_name):
                print(f"Missing branch in: {file_path}")
                missing_files.append(file_path)
    
    if missing_files:
        print("\nFiles missing the branch:")
        for f in missing_files:
            print(f)
    else:
        print("\nAll files contain the branch!")

if __name__ == "__main__":
    main()
