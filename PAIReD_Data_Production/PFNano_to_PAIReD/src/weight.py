import uproot
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt
from collections import defaultdict
import glob

#dijet_pt_bins = np.array([0.0, 195.8, 381.5, 99999.9])
#dijet_pt_bins = np.array([0.0, 77.0, 137.9, 195.8, 251.4, 310.4, 381.5, 496.1, 715.7, 1080.7, 99999.9])
#dijet_pt_bins = np.array([0.0, 137.9, 251.4, 381.5, 715.7, 99999.9])

mc_mass_bins = np.array([25, 35, 45, 55, 65, 75, 85, 95, 105, 115, 125, 135, 145, 155, 165, 175, 185, 195, 205, 215, 225, 235, 245])
reweight_classes = ['label_BB', 'label_CC', 'label_bb', 'label_bx', 'label_cx', 'label_ll']
class_weights = {'label_BB': 1.0, 'label_CC': 1.0, 'label_bb': 0.08, 'label_bx': 0.18, 'label_cx': 0.28, 'label_ll': 0.7}
class_counts = {'label_BB': 0, 'label_CC': 0, 'label_bb': 0, 'label_bx': 0, 'label_cx': 0, 'label_ll': 0}
reweight_threshold = 0.1

directories = ["../data/clustered_1GeV/XtoHH/MX-500-1000/", "../data/clustered_1GeV/XtoHH/MX-1000-4000/", "../data/clustered_1GeV/DY", "../data/clustered_1GeV/TT", "../data/clustered_1GeV/ZHvar"]
event_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))

def find_quantiles(directories):
    cc_dijet_pt = []
    cc_gendijet_mass = []
    num_jets = dict()
    for directory in directories:
        files = glob.glob(f"{directory}/*.root")
        num_jets[directory] = 0
        for file_name in files:
            print(f"Processing file: {file_name}")
            with uproot.open(file_name) as f:
                tree = f["tree"]
                arrays = tree.arrays(["dijet_pt", "MC_gendijet_mass", "genweight", "label_CC"], library="np") 
                dijet_pt = arrays["dijet_pt"]
                num_jets[directory] += len(dijet_pt)/0.9
                gendijet_mass = arrays["MC_gendijet_mass"]
                cc_mask = arrays["label_CC"] == 1
                cc_dijet_pt.extend(dijet_pt[cc_mask])
                cc_gendijet_mass.extend(gendijet_mass[cc_mask])
        print(directory, num_jets[directory], "jets total")
    dijet_pt_quartiles = np.percentile(np.array(cc_dijet_pt), [25, 50, 75])
    gendijet_mass_quartiles = np.percentile(np.array(cc_gendijet_mass), [10,20,30,40,50,60,70,80,90,100])
    print("HERE ARE THE QUANTILES")
    print(gendijet_mass_quartiles)
    print("HERE ARE THE JET COUNTS")
    for key, val in num_jets.items():
        print(key, val)
    return np.concatenate(([0], dijet_pt_quartiles, [9999999]))#, np.concatenate(([0], gendijet_mass_quartiles, [9999999]))

dijet_pt_bins = find_quantiles(directories)
print(dijet_pt_bins)
#print(mc_mass_bins)


def process_files(directory):
    files = glob.glob(f"{directory}/*.root")
    for file_name in files:
        print(f"Processing file: {file_name}")
        with uproot.open(file_name) as f:
            tree = f["tree"]
            arrays = tree.arrays(["dijet_pt", "MC_gendijet_mass", "genweight"] + reweight_classes, library="np")
            dijet_pt = arrays["dijet_pt"]
            mc_mass = arrays["MC_gendijet_mass"]

            pt_indices = np.digitize(dijet_pt, dijet_pt_bins) - 1
            mass_indices = np.digitize(mc_mass, mc_mass_bins) - 1

            for label in reweight_classes:
                mask = arrays[label] == 1  # Select events for this class
                for pt_idx, mass_idx in zip(pt_indices[mask], mass_indices[mask]):
                    if pt_idx >= 0 and mass_idx >= 0 and pt_idx < len(dijet_pt_bins)-1 and mass_idx < len(mc_mass_bins)-1: event_counts[pt_idx][mass_idx][label] += 1
                    class_counts[label] += 1
    for label, count in class_counts.items():
        print(label, count)
    

def compute_weights():
    weights = {}
    for pt_idx in event_counts:
        for mass_idx in event_counts[pt_idx]:
            bin_weights = {}
            for label in reweight_classes:
                bin_weights[label] = class_weights[label] * class_counts[label] / event_counts[pt_idx][mass_idx][label]
            weights[(pt_idx, mass_idx)] = bin_weights
    return weights

def plot_weights(weights, final=False):
    pt_vals = sorted(set(pt_idx for pt_idx, mass_idx in weights.keys()))
    mass_vals = sorted(set(mass_idx for pt_idx, mass_idx in weights.keys()))
    for label in reweight_classes:
        label_weights = np.zeros((len(pt_vals), len(mass_vals)))
        for (pt_idx, mass_idx), bin_weights in weights.items():
            pt_idx_pos = pt_vals.index(pt_idx)
            mass_idx_pos = mass_vals.index(mass_idx)
            label_weights[pt_idx_pos, mass_idx_pos] = bin_weights.get(label, 0)
        fig, ax = plt.subplots(figsize=(6, 6))
        cax = ax.imshow(label_weights, cmap='viridis', aspect='auto', origin='lower', interpolation='nearest')
        ax.set_title(f'Weights for {label}')
        ax.set_xlabel('Mass Index')
        ax.set_ylabel('Pt Index')
        fig.colorbar(cax, ax=ax, label='Weight')
        fig.tight_layout()
        if final: fig.savefig(f'{label}_weights_final.png')
        else: fig.savefig(f'{label}_weights.png')
        plt.close(fig)

def adjust_weights(weights):
    max_weight = 0
    for pt_idx in event_counts:
        for mass_idx in event_counts[pt_idx]:
            for label in reweight_classes:
                if weights[(pt_idx, mass_idx)][label] > max_weight: max_weight = weights[(pt_idx, mass_idx)][label]
    for pt_idx in event_counts:
        for mass_idx in event_counts[pt_idx]:
            for label in reweight_classes:
                weights[(pt_idx, mass_idx)][label] = max(weights[(pt_idx, mass_idx)][label]/max_weight, reweight_threshold)
    return weights

for directory in directories:
    process_files(directory)

weights = compute_weights()
plot_weights(weights)
weights = adjust_weights(weights)
plot_weights(weights, final=True)
sorted_weights = sorted(weights.items())
print(dijet_pt_bins)
print(sorted_weights)
output_array = ak.Array([{
    'pt_idx': pt_idx,
    'mass_idx': mass_idx,
    **bin_weights
} for (pt_idx, mass_idx), bin_weights in sorted_weights])

with uproot.recreate("reweight.root") as file:
    file["tree"] = {"weights": output_array}

print("Weights saved to reweight.root")


