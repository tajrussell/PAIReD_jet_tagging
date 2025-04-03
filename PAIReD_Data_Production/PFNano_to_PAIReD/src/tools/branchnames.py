import awkward as ak

# branch names for the PAIReD file format
BranchNames = {
    "part":  ["eta", "phi", "pt", "charge", "pid", "d0val", 
              "d0err", "dzval", "dzerr", "mass", "puppiweight",
              "px", "py", "pz", "energy",
              "deta1", "deta2", "dphi1", "dphi2",
              "jetindex", "in_jet1", "in_jet2","jetindex_cut"],
    "sv":    ["charge", "chi2", "dlen", "dlenSig", 
              "dxy", "dxySig", "eta", "mass", "ndof",
              "ntracks", "pAngle", "phi", "pt",
              "x", "y", "z",
              "deta1", "deta2", "dphi1", "dphi2"],
    "jet":   ["eta", "phi", "pt", "nparticles", "energy", "rawfactor",
              "index"],
    "dijet": ["eta", "phi", "pt", "mass", "nparticles",
              "index", "ak4_mass", "pnet_mass"],
    "singleValues": ['event', 'genweight','run', 'pv_n', 'pv_ngood', 'fgrfccpu', 'fgrfcc', 'Pileup_nPU',
                     'MC_higgs_pt','MC_higgs_eta','MC_higgs_phi','MC_higgs_mass',
                     'MC_gendijet_pt','MC_gendijet_eta','MC_gendijet_phi','MC_gendijet_mass',
                     'MC_genjet1_flav','MC_genjet2_flav','MC_genjet1_matched','MC_genjet2_matched',
                     'MC_n_c','MC_vector_flav', 'MC_lepton_channel',
                     'MC_drqq','MC_higgs_flav', 'MC_higgs_valid',
                     'MC_physics_process', 'extra_jets',
                     'label_bb','label_cx','label_ll','label_bx', 'label_CC', 'label_BB', 'label_LL']
}


# small helpers function that converts the branch names into the ak-type of the tree
def branchName2akType(branch):
    """
    Input: branch (str) - can be "part", "sv", "jet" or "dijet"
    """
    if branch == "jet" or branch == "dijet":
        s = "{"
    else:
        s = "var * {"

    for name in BranchNames[branch]:
        # add the subbranches
        s += name + ": float32, "
    
    # remove ", " at the end of the s string
    s = s[:-2]

    s += "}"

    return ak.types.from_datashape(s, highlevel=False)



outputTreeType =  {
    'event': "uint64",
    'genweight': "float32",
    'run': "uint8",
    'Pileup_nPU': "uint8",
    'pv_n': "uint8",
    'pv_ngood': "uint8",
    'fgrfccpu': "float32", 
    'fgrfcc': "float32",
    'jet1': branchName2akType("jet"),
    'jet2': branchName2akType("jet"),
    'dijet': branchName2akType("dijet"),
    'part': branchName2akType("part"),
    'sv': branchName2akType("sv"),
    'MC_higgs_pt': "float32",
    'MC_higgs_eta': "float32",
    'MC_higgs_phi': "float32",
    'MC_higgs_mass': "float32",
    'MC_gendijet_pt': "float32",
    'MC_gendijet_eta': "float32",
    'MC_gendijet_phi': "float32",
    'MC_gendijet_mass': "float32",
    'MC_genjet1_flav': "int8",
    'MC_genjet2_flav': "int8",
    'MC_genjet1_matched': "int8",
    'MC_genjet2_matched': "int8",
    'MC_drqq': "float32",
    'MC_n_c': "uint8",
    'MC_higgs_flav': "uint8",
    'MC_higgs_valid': "uint8",
    'MC_vector_flav': "uint8",
    'MC_lepton_channel': "int8",
    'MC_physics_process': "uint8",
    'label_BB': "int8",
    'label_CC': "int8",
    'label_bb': "int8",
    'label_bx': "int8",
    'label_cx': "int8",
    'label_ll': "int8",
    'label_LL': "int8",
    'extra_jets': "int8",
    }
