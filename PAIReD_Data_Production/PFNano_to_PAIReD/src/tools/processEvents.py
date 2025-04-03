"""
 * Function: processEvents
 * -----------------------
Process the events to produce PAIReD jets with all relevant quantities.

Parameters
----------
Events : event tree object
    Object containing all information of multiple simulated events.
ttbar : bool
    True if the data file processed is a ttbar Event.


Returns
-------
DataPAIReD : dict
    A dictionary with the content summarized in /notes/PAIReD_ROOT_file_content.md
"""

# Include dependencies
import awkward as ak
import numpy as np
import vector
import uproot
from tools.branchnames import BranchNames
from tools.helpers import isoLeptonCut, deltaPhi, getJetClusterIndex, _is_rootcompat, isClustered, isHighPt, getJetClusterIndexCut
from tools.getMCInfo import getMCInfo
from tools.PAIReD_geometries.ellipse import isInPAIReD
import fastjet
import sys
#sys.path.append('/home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/src/')
#import fjcontrib



def processEvents(Events, physics_process=0, PAIReD_geometry="Ellipse"):

    # import the PAIReD jet geometry
    if PAIReD_geometry=="Ellipse":
        from tools.PAIReD_geometries.ellipse import isInPAIReD
    elif PAIReD_geometry=="AK4":
        from tools.PAIReD_geometries.ak4 import isInPAIReD
    elif PAIReD_geometry=="Dumbbell":
        from tools.PAIReD_geometries.dumbbell import isInPAIReD
    elif PAIReD_geometry=="Clustered":
        from tools.PAIReD_geometries.ellipse import isInPAIReD
    elif PAIReD_geometry=="Clustered_HighPT":
        from tools.PAIReD_geometries.ellipse import isInPAIReD
    elif PAIReD_geometry=="AKX":
        from tools.PAIReD_geometries.ellipse import isInPAIReD

    # check if events have two or more jets that survive the selection cut
    Jetcut = (((Events.Jet_pt * (1-Events.Jet_rawFactor)) > 20) & 
                (abs(Events.Jet_eta) < 2.5) & 
                (isoLeptonCut(Events)) &
                (Events.Jet_jetId > 4)
            )
    Eventcut = (ak.sum(Jetcut, axis=1) > 1)

    # cut all events that do not meet this criterion
    Events = Events[Eventcut]

    # check if there are events with two or more jets
    if sum(Eventcut)==0:
        print("   !!! No event in the batch had more than one jet surviving the selection cut !!!")
        return False

    # ATTENTION!!! defined arrays are copies of the data!!!

    # sort the particles according to their pt
    s = ak.argsort(Events.PFCands_pt, ascending=False, axis=1)
    # and only include particles not considered pileup
    s = s[Events.PFCands_puppiWeight[s] != 0]
    Part = {"phi": Events.PFCands_phi[s], "eta": Events.PFCands_eta[s],
            "pt": Events.PFCands_pt[s], "charge": Events.PFCands_charge[s],
            "pid": Events.PFCands_pdgId[s], "d0val": Events.PFCands_d0[s],
            "d0err": Events.PFCands_d0Err[s], "dzval": Events.PFCands_dz[s],
            "dzerr": Events.PFCands_dzErr[s], "mass": Events.PFCands_mass[s],
            "puppiweight": Events.PFCands_puppiWeight[s],
            "jetindex": getJetClusterIndex(Events)[s],
            "jetindex_cut": getJetClusterIndexCut(Events, Jetcut)[s]}
    
    # sort the secondary vertices (SV) according to their dlenSig
    s = ak.argsort(Events.SV_dlenSig, ascending=False, axis=1)
    SV = {k: Events["SV_"+k][s] for k in BranchNames["sv"][:16]}
    
    # define 4-vector for jets in order to calculate the energy
    Jetcut = (((Events.Jet_pt * (1-Events.Jet_rawFactor)) > 20) & 
            (abs(Events.Jet_eta) < 2.5) & 
            (isoLeptonCut(Events)) &
            (Events.Jet_jetId > 4)
            )
    Jet4 = vector.zip({"eta": Events.Jet_eta[Jetcut], "phi": Events.Jet_phi[Jetcut],
            "pt": Events.Jet_pt[Jetcut], "mass": Events.Jet_mass[Jetcut]})
    # create jet object
    Jet = ak.zip({"phi": Events.Jet_phi[Jetcut], "eta": Events.Jet_eta[Jetcut],
            "pt": Events.Jet_pt[Jetcut], "energy": Jet4.E,
            "nparticles": Events.Jet_nConstituents[Jetcut],
            "rawfactor": Events.Jet_rawFactor[Jetcut],
            "index": ak.local_index(Events.Jet_pt, axis=1)[Jetcut],
            "PNetRegPtRawCorr": Events.Jet_PNetRegPtRawCorr[Jetcut]})

    # get all jet pair combinations
    Jet = ak.combinations(Jet, 2, axis=1, fields=["j1", "j2"],
            replacement=False)
    # get also the combinations for the genJetIdx
    Jet_genJetIdx = ak.combinations(Events.Jet_genJetIdx[Jetcut], 2, 
            axis=1, fields=["j1", "j2"], replacement=False)

    isInJets = isInPAIReD(Jet.j1.eta, Jet.j2.eta, Jet.j1.phi, Jet.j2.phi, Events.Jet_eta[Jetcut], Events.Jet_phi[Jetcut])
    #print(isInJets.tolist())
    ExtraJets = ak.sum(isInJets, axis=2) - 2

    # get particles inside the PAIReD jet
    if "Cluster" in PAIReD_geometry:
        isInPart, isJ1, isJ2 = isClustered(Jet.j1.index, Jet.j2.index, Part["jetindex"])
        if "HighPT" in PAIReD_geometry:
            isInEllipse = isInPAIReD(Jet.j1.eta, Jet.j2.eta, Jet.j1.phi, Jet.j2.phi, Part["eta"], Part["phi"])           
            isHighPtInEllipse = isHighPt(isInEllipse, Part["pt"])

            #total_parts = ak.sum((isInPart) > -1)
            total_clustered = ak.sum(isInPart)
            #total_ellipse = ak.sum(isInEllipse)
            isClusteredInEllipse = (isInEllipse) & (isInPart)
            #isJ1InEllipse = (isJ1) & (isInEllipse)
            #isJ2InEllipse = (isJ2) & (isInEllipse)
            #total_j1_ellipse = ak.sum(isJ1InEllipse)
            #total_j2_ellipse = ak.sum(isJ2InEllipse)
            total_clustered_ellipse = ak.sum(isClusteredInEllipse)
            
            #print("Ratio of particles that are clustered:", total_clustered/total_parts)
            #print("Ratio of particles that are in ellipses:", total_ellipse/total_parts)
            #print("Ratio of particles in ellipse that are clustered:", total_clustered_ellipse/total_ellipse)
            #print("Ratio of clustered particles that are in the ellipse:", total_clustered_ellipse/total_clustered)
            #print("Ratio of particles in ellipse clustered to jet 1:", total_j1_ellipse/total_ellipse)
            #print("Ratio of particles in ellipse clustered to jet 2:", total_j2_ellipse/total_ellipse)
            
            isInPart = (isInPart) | (isHighPtInEllipse)
            #print(ak.sum(isInEllipse, axis=-1))
            #print(ak.sum(isInPart, axis=-1))
    elif PAIReD_geometry == "AKX":
        isInPart, major_axes = isInPAIReD(Jet.j1.eta, Jet.j2.eta, Jet.j1.phi, Jet.j2.phi, Part["eta"], Part["phi"], return_semimajor=True)
    else:
        isInPart = isInPAIReD(Jet.j1.eta, Jet.j2.eta, Jet.j1.phi, Jet.j2.phi, Part["eta"], Part["phi"])
        #print(ak.sum(isInPart, axis=-1))
    # get the generated particles inside the PAIReD jet
    isInGenPart = isInPAIReD(Jet.j1.eta, Jet.j2.eta, Jet.j1.phi, Jet.j2.phi,
                Events.GenPart_eta, Events.GenPart_phi)

    # get the secondary vertices inside the PAIReD jet
    isInSV = isInPAIReD(Jet.j1.eta, Jet.j2.eta, Jet.j1.phi, Jet.j2.phi,
                SV["eta"], SV["phi"])

    # make part arrays broadcastable with isInPart
    # and drop all particles that are not in the PAIReD jet
    for name in Part.keys():
        Part[name] = ak.unflatten(Part[name], 1, axis=0)
        Part[name], isInPart = ak.broadcast_arrays(Part[name], isInPart)
        #print(ak.sum(isInPart, axis=-1))
        Part[name] = Part[name][isInPart]
        #print(name)
        #print(Part[name])
    # Run AK clustering for AKX geometry and redo Part
    if PAIReD_geometry == "AKX":
        vector.register_awkward()
        major_axes = ak.flatten(major_axes)
        part_4vecs = vector.zip({"eta": Part["eta"], "phi": Part["phi"], "pt": Part["pt"], "mass": Part["mass"]})
        jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 3.14)
        cluster = fastjet.ClusterSequence(part_4vecs, jetdef)
        soft_drop_groomer = fastjet.contrib.SoftDrop(0.1, 0.0)
        jets = cluster.inclusive_jets()
        softdrop_jets = [soft_drop_groomer.result(jet) for jet in jets]
        print(softdrop_jets)
        #print(jets[0])
        #print(cluster.constituent_index()[0])
        #softdrop_remain = [jet[0].soft_drop(0.1, 0).constituent_index() for jet in jets]
        #print(softdrop_remain)
                #print(part_dijet_idx.to_list())
                #new_isInPart = (part_dijet_idx == 0)
                #for name in Part.keys():
                #    Part[name][i][j] = Part[name][i][j][new_isInPart]


    # make SV arrays broadcastable with isInSV
    # and drop all SVs that are not in the PAIReD jet
    for name in SV.keys():
        SV[name] = ak.unflatten(SV[name], 1, axis=0)
        SV[name], isInSV = ak.broadcast_arrays(SV[name], isInSV)
        SV[name] = SV[name][isInSV]

    # flags indicating if particles are clustered in one of seed jets
    Part["in_jet1"] = Part["jetindex"] == Jet.j1.index
    Part["in_jet2"] = Part["jetindex"] == Jet.j2.index

    # calculate remaining variables (px,py,pz,E) from 4-vector
    Part4 = vector.zip({"eta": Part["eta"], "phi": Part["phi"],
                "pt": Part["pt"], "mass": Part["mass"]})
    Part["px"] = Part4.px
    Part["py"] = Part4.py
    Part["pz"] = Part4.pz
    Part["energy"] = Part4.E
    #print(Part["energy"])
    # calculate deta and dphi with respect to the jets
    # (part_eta - jet1_eta) * (jet1_eta > 0 ? 1 : -1)
    Part["deta1"] = (Part["eta"] - Jet.j1.eta) * (-1)**(Jet.j1.eta < 0)
    # (part_eta - jet2_eta) * (jet2_eta > 0 ? 1 : -1)
    Part["deta2"] = (Part["eta"] - Jet.j2.eta) * (-1)**(Jet.j2.eta < 0)
    # dphi
    Part["dphi1"] = deltaPhi(Part["phi"], Jet.j1.phi)
    Part["dphi2"] = deltaPhi(Part["phi"], Jet.j2.phi)

    # calculate the deta and dphi for the secondary vertices
    SV["deta1"] = (SV["eta"] - Jet.j1.eta) * (-1)**(Jet.j1.eta < 0)
    SV["deta2"] = (SV["eta"] - Jet.j2.eta) * (-1)**(Jet.j2.eta < 0)
    SV["dphi1"] = deltaPhi(SV["phi"], Jet.j1.phi)
    SV["dphi2"] = deltaPhi(SV["phi"], Jet.j2.phi)

    # calculate dijets from dijet 4-vector
    Dijet4 = vector.zip({"px": ak.sum(Part4.px, axis=2),
                    "py": ak.sum(Part4.py, axis=2),
                    "pz": ak.sum(Part4.pz, axis=2),
                    "E":  ak.sum(Part4.E, axis=2)})
    #print(ak.num(Part["pt"], axis=-1))       
    if "Cluster" not in PAIReD_geometry:
        # get indices of the particles in one of the AK4 jets
        s = ak.local_index(Part["pt"], axis=2)[(Part["in_jet1"]) | (Part["in_jet2"])]
        # get indices of the particles in no jet
        s = ak.concatenate([s, ak.local_index(Part["pt"], axis=2)[(~Part["in_jet1"]) & (~Part["in_jet2"])]], axis=2)
        # bring particles clustered to the seed jets to the front of the list
        for name in Part.keys():
            Part[name] = Part[name][s]
    #print(ak.num(Part["pt"], axis=-1))
    # prepare jet and particle objects for the tree
    part = ak.zip(Part)
    #print("part", part)
    sv = ak.zip(SV)
    jet1 = Jet.j1
    jet2 = Jet.j2
    #print(jet1)
    jet1_4vec = vector.zip({"eta": jet1.eta, "phi": jet1.phi, "pt": jet1.pt, "E": jet1.energy})
    jet2_4vec = vector.zip({"eta": jet2.eta, "phi": jet2.phi, "pt": jet2.pt, "E": jet2.energy})
    print(1-jet1.rawfactor)
    jet1_4vec_pnet = vector.zip({"eta": jet1.eta, "phi": jet1.phi, "pt": jet1.pt * jet1.PNetRegPtRawCorr / (1-jet1.rawfactor), "E": jet1.energy})
    jet2_4vec_pnet = vector.zip({"eta": jet2.eta, "phi": jet2.phi, "pt": jet2.pt * jet2.PNetRegPtRawCorr / (1-jet2.rawfactor), "E": jet2.energy})
    dijet = ak.zip({"eta": Dijet4.eta, "phi": Dijet4.phi,
                    "pt": Dijet4.pt, "mass": Dijet4.mass,
                    "nparticles": ak.num(part, axis=2),
                    "index": ak.local_index(Jet.j1, axis=1),
                    "ak4_mass": (jet1_4vec+jet2_4vec).mass,
                    "pnet_mass": (jet1_4vec_pnet+jet2_4vec_pnet).mass})

    # prepare other single value branches
    ones = ak.ones_like(Jet.j1.phi)
    event = Events.event * ones
    genweight = Events.genWeight * ones
    run = Events.run * ones
    Pileup_nPU = Events.Pileup_nPU * ones
    pv_n = Events.PV_npvs * ones
    pv_ngood = Events.PV_npvsGood * ones
    fgrfcc = Events.Rho_fixedGridRhoFastjetCentralCalo * ones
    fgrfccpu = Events.Rho_fixedGridRhoFastjetCentralChargedPileUp * ones
    MC_physics_process = physics_process * ones
    #print("part.pt")
    #print(ak.num(part.pt, axis=-1))
    # get summarized information from the MC simulation
    MCInfo = getMCInfo(Events, isInGenPart, Jet_genJetIdx, Jetcut, physics_process)

    if MCInfo == False:
        return False

    DataPAIReD = {
        "event" : event,
        "genweight": genweight,
        "run" : run,
        "Pileup_nPU" : Pileup_nPU,
        "pv_n" : pv_n,       
        "pv_ngood" : pv_ngood,
        "fgrfcc" : fgrfcc,
        "fgrfccpu" : fgrfccpu,
        "jet1" : jet1,
        "jet2" : jet2,
        "dijet" : dijet,
        "part" : part,
        "sv" : sv,
        "MC_physics_process" : MC_physics_process,
        "extra_jets": ExtraJets
    }

    for name in MCInfo.keys():
        DataPAIReD[name] = MCInfo[name]
    higgs_cut = ak.sum(DataPAIReD["MC_higgs_valid"], axis=1) > 0
    nonempty_cut = ak.num(DataPAIReD["part"]["pt"], axis=1) > 0
    total_cut = higgs_cut & nonempty_cut
    flattened_Data = dict()
    # remove the dimension of N_events
    for name in DataPAIReD.keys():
        flat_precut = ak.copy(ak.flatten(DataPAIReD[name], axis=1))
        if DataPAIReD[name].fields:
            postcut_dict = dict()
            for i,n in enumerate(DataPAIReD[name].fields):
                precut = ak.unzip(DataPAIReD[name], highlevel=True)[i]
                postcut_dict[n] = precut[total_cut]
                #if name == "part" and n == "pt":
                #    print(precut)
                #    print(ak.num(precut, axis=-1))
                #    print(postcut_dict[n])
                #    print(ak.num(postcut_dict[n], axis=-1))
            postcut = ak.zip(postcut_dict)
        else:
            postcut = ak.copy(DataPAIReD[name][total_cut])
        flattened = ak.copy(ak.flatten(postcut, axis=1))
        flattened_Data[name] = flattened
    '''# add reweighting
    reweights = dict()
    with uproot.open("tools/reweight.root") as file:
        reweights["label_BB"] = file["tree"]["weights_label_BB"].array()
        reweights["label_CC"] = file["tree"]["weights_label_CC"].array()
        reweights["label_bx"] = file["tree"]["weights_label_bx"].array()
        reweights["label_cx"] = file["tree"]["weights_label_cx"].array()
        reweights["label_ll"] = file["tree"]["weights_label_ll"].array()
    #dijet_pt_bins = np.array([0.0, 195.8, 381.5, 99999.9])
    dijet_pt_bins = np.array([0.00000000e+00, 3.03419350e+02, 4.66926422e+02, 9.08118271e+02, 9.99999900e+06])
    mc_mass_bins = np.array([25, 35, 45, 55, 65, 75, 85, 95, 105, 115, 125, 135, 145, 155, 165, 175, 185, 195, 205, 215, 225, 235, 245])
    pt_idx = np.digitize(flattened_Data["dijet"]["pt"], dijet_pt_bins) - 1
    mass_idx = np.digitize(flattened_Data["MC_gendijet_mass"], mc_mass_bins) - 1
    weight_mask = (pt_idx >= 0) & (pt_idx < len(dijet_pt_bins)-1) & (mass_idx >= 0) & (mass_idx < len(mc_mass_bins)-1)
    weight_idx = pt_idx*(len(mc_mass_bins)-1) + mass_idx
    weight_idx = np.where(weight_mask, weight_idx, 0)
    training_weights = ak.zeros_like(weight_idx, dtype=float)
    for label in ["label_BB", "label_CC", "label_bx", "label_cx", "label_ll"]:
        label_mask = ak.Array(weight_mask) & ak.fill_none(flattened_Data[label], False)
        label_reweights = ak.Array(reweights[label][weight_idx])
        training_weights = training_weights + ak.where(label_mask, label_reweights, 0)
    flattened_Data["training_weight"] = training_weights'''

    num_jets = len(flattened_Data["event"])
    np.random.seed(0)
    rand_idx = np.random.permutation(num_jets)
    split_idx = int(0.9 * num_jets)
    split_data, split_data_test = dict(), dict()
    #print(ak.num(flattened_Data["part"].pt))
    for name in flattened_Data.keys():
        if flattened_Data[name].fields:
            split_data_dict = dict()
            split_data_test_dict = dict()
            for i,n in enumerate(flattened_Data[name].fields):
                presplit = ak.unzip(flattened_Data[name], highlevel=True)[i]
                split_data_test_dict[n] = presplit[rand_idx[split_idx:]]
                split_data_dict[n] = presplit[rand_idx[:split_idx]]
            split_data_test[name] = ak.zip(split_data_test_dict)
            split_data[name] = ak.zip(split_data_dict)
        else:
            split_data_test[name] = flattened_Data[name][rand_idx[split_idx:]]
            split_data[name] = flattened_Data[name][rand_idx[:split_idx]]
    #print(ak.num(split_data["part"].pt))
    return split_data, split_data_test
