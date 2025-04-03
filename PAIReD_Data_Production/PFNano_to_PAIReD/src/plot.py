import uproot
import matplotlib.pyplot as plt
import awkward as ak
import numpy as np
# Define file names and labels
files = {
    "hex_ellipse.root": "Ellipse",
    "hex_cluster.root": "Clustered",
    "hex_500MeV.root": "Clustered + >500MeV"
}

# Initialize dictionaries for storing data
npart_data = {}
part_pt_data = {}
part_eta_data = {}
# Load data from ROOT files
for file, label in files.items():
    with uproot.open(file) as f:
        tree = f["tree"]  # Adjust the tree name if needed
        part_pt = tree["part_pt"].array()
        part_eta = tree["part_eta"].array()
        part_pt_data[label] = ak.flatten(part_pt)
        part_eta_data[label] = ak.flatten(part_eta)
        npart_data[label] = ak.num(part_pt)

        if label == "Ellipse":
            part_not_clustered = np.logical_not(np.logical_or(tree["part_in_jet1"].array(), tree["part_in_jet2"].array()))
            part_pt_data["Not Clustered"] = ak.flatten(part_pt[part_not_clustered])
            part_eta_data["Not Clustered"] = ak.flatten(part_eta[part_not_clustered])
            part_pt_data["No Jet"] = ak.flatten(tree["part_jetindex"].array())
            print(ak.sum(part_pt_data["No Jet"] == -1))
            part_eta_data["No Jet"] = ak.flatten(part_eta[(tree["part_jetindex"].array() == -1)])
# Plot npart distribution
plt.figure(figsize=(10, 5))
bins=np.linspace(0,100,50)
for label, data in npart_data.items():
    print(data)
    plt.hist(ak.to_numpy(data), bins=bins, alpha=0.6, label=label, histtype='step', linewidth=2)
plt.xlabel("npart")
plt.ylabel("Counts")
plt.title("npart Distribution")
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("npart.png")

# Plot part_pt distribution
plt.figure(figsize=(10, 5))
bins=np.linspace(0,1500,50)
for label, data in part_pt_data.items():
    if label != "No Jet": continue
    else:
        plt.hist(ak.to_numpy(data), bins=np.linspace(-1.5, 16.5, 19), alpha=0.6, label=label, histtype='step', linewidth=2)
        continue
    plt.hist(ak.to_numpy(data), bins=bins, alpha=0.6, label=label, histtype='step', linewidth=2, density=True)
plt.xlabel("part_pt (MeV)")
plt.ylabel("Fraction of PFCands with Greater PT")
plt.yscale('log')
plt.title("part_pt Distribution")
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("part_pt.png")
plt.clf()

'''
# Plot part_p distribution
plt.figure(figsize=(10, 5))
bins=np.linspace(0,1500,50)
for label, data in part_pt_data.items():
    p = ak.to_numpy(data) * np.cosh(ak.to_numpy(part_eta_data[label]))
    if label == "Ellipse": continue
    else:
        ellipse_p = ak.to_numpy(part_pt_data["Ellipse"]) * np.cosh(ak.to_numpy(part_eta_data["Ellipse"]))
        ellipse_hist, _ = np.histogram(ellipse_p, bins=bins)
        p_hist, _ = np.histogram(p, bins=bins)
        ratio = p_hist / ellipse_hist
    plt.plot((bins[:-1] + bins[1:])/2, ratio, label=label)
    #plt.hist(ratio, bins=bins, alpha=0.6, label=label, histtype='step', linewidth=2, density=True)
plt.xlabel("part_p (MeV)")
plt.ylabel("Counts")
plt.title("part_pt Distribution")
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("part_p.png")
plt.clf()

plt.figure(figsize=(5, 5))
with uproot.open("hex_ellipse.root") as f:
    tree = f["tree"] 
    jet1_matched = tree["MC_genjet1_matched"].array()
    jet2_matched = tree["MC_genjet2_matched"].array()
    jets_matched = jet1_matched + jet2_matched
    extra_jets = tree["extra_jets"].array()

plt.hist(jets_matched, bins=[-0.5,0.5,1.5,2.5])
plt.xlabel("Number of jets matched to gen jets")
plt.ylabel("Number of dijets")
plt.yscale('log')
plt.tight_layout()
plt.savefig("matched.png")
plt.clf()

plt.figure(figsize=(5, 5))
plt.hist(extra_jets, bins=[-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5])
plt.xlabel("Number of extra jets inside ellipse")
plt.ylabel("Number of dijets")
plt.yscale('log')
plt.tight_layout()
plt.savefig("extras.png")
plt.clf()'''