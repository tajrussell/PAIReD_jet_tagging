def descript():
    """
    * Function: getMCInfo
    * ----------------------
    Get the relevant true MC information of the event.

    Parameters
    ----------
    Events : event tree object
        Object containing all information of multiple simulated events.
    isInGenPart : bool (array)
        Boolean holding the information whether the generated particle is inside the
        PAIReD jet or not.
    Jet_genJetIdx : array / ak-zip type object
        Indices of the matched generated jets to the reconstructed AK4 jets
    Jetcut : bool (array)
        Boolean holding the information whether the AK4 jet is accepted or not.


    Returns
    -------
    MCInfo : dict
        A dictionary with the following content:
            "MC_higgs_pt" : float (array)
                Higgs transverse momentum (length of the array = number of jet pairs).
            "MC_higgs_eta" : float (array)
                Higgs pseudorapidity (length of the array = number of jet pairs).
            "MC_higgs_phi" : float (array)
                Higgs phi angle (length of the array = number of jet pairs).
            "MC_higgs_mass" : float (array)
                Higgs mass (length of the array = number of jet pairs).
            "MC_higgs_flav" : int (array)
                Generated flavor number in Higgs decay:
                - digit represents the flavor of the daughter particle 
                (e.g., 1 for d, 3 for s, or 5 for b)
                - if there is no vector boson, = 0
            "MC_vector_flav" : int (array)
                Generated flavor number in vector boson decay:
                - digit represents the flavor of the daughter particles 
                (11 for e-like,  13 for mu-like, 15 for tau-like)
                - if there is no vector boson, = 0
            "MC_lepton channel" : int (array)
                number of charged leptons the vector boson decays into (0L, 1L, 2L).
            "MC_gendijet_pt" : float (array)
                Combined GenJet transverse momentum (length of the array = number of jet pairs).
            "MC_gendijet_eta" : float (array)
                Combined GenJet pseudorapidity (length of the array = number of jet pairs).
            "MC_gendijet_phi" : float (array)
                Combined GenJet phi angle (length of the array = number of jet pairs).
            "MC_gendijet_mass" : float (array)
                Combined GenJet mass (length of the array = number of jet pairs).
            "MC_genjet1_flav" : float (array)
                Parton flavour of the GenJet matched to the first AK4 jet. 
                In case of no match = 0 (length of the array = number of jet pairs).
            "MC_genjet2_flav" : float (array)
                Parton flavour of the GenJet matched to the second AK4 jet. 
                In case of no match = 0 (length of the array = number of jet pairs).
            "MC_genjet1_matched" : bool (array)
                Boolean whether first AK4 jet has a matched GenJet 
                (length of the array = number of jet pairs).
            "MC_genjet2_matched" : bool (array)
                Boolean whether second AK4 jet has a matched GenJet 
                (length of the array = number of jet pairs).
            "MC_drqq" : float (array)
                Delta R between decay products of the Higgs (length of array = number of jet
                pairs).
            "MC_n_c" : int (array)
                Number of c jets in the event (length of the array = number of jet pairs).
            "label_BB" : bool (array)
                True label for BB PAIReD jets.
            "label_CC" : bool (array)
                True label for CC PAIReD jets.
            "label_LL" : bool (array)
                True label for LL PAIReD jets.
            "label_ll" : bool (array)
                True label for ll PAIReD jets.
            "label_bb" : bool (array)
                True label for bb PAIReD jets.
            "label_bl" : bool (array)
                True label for bl PAIReD jets.
            "label_cc" : bool (array)
                True label for cc PAIReD jets.
            "label_cl" : bool (array)
                True label for cl PAIReD jets.    
    """
    return

# Include dependencies
import awkward as ak
import numpy as np
import vector
from tools.helpers import deltaR

def processHiggs_HH(Events, isInGenPart, higgs_idx, MCInfo, Jetcut):
    # process fake higgs
    ones = ak.ones_like(isInGenPart[:,:,0])
    singleInfo = dict()
    for i in range(2):
        singleInfo["MC_higgs_pt_"+str(i)] = ones * Events.GenPart_pt[higgs_idx][:,i]
        singleInfo["MC_higgs_eta_"+str(i)] = ones * Events.GenPart_eta[higgs_idx][:,i]
        singleInfo["MC_higgs_phi_"+str(i)] = ones * Events.GenPart_phi[higgs_idx][:,i]
        singleInfo["MC_higgs_mass_"+str(i)] = ones * Events.GenPart_mass[higgs_idx][:,i]
        d1 = ak.argmax(1*(Events.GenPart_genPartIdxMother == higgs_idx[:,i]), axis=1, keepdims=True)
        d2 = ak.argmax(ak.local_index(Events.GenPart_pdgId, axis=1)*(Events.GenPart_genPartIdxMother == higgs_idx[:,i]), axis=1, keepdims=True)
        # kinematics of the daughter particles
        d1_eta = Events.GenPart_eta[d1]
        d2_eta = Events.GenPart_eta[d2]
        d1_phi = Events.GenPart_phi[d1]
        d2_phi = Events.GenPart_phi[d2]

        jet_shape = ak.ones_like(Events.Jet_eta[Jetcut])
        d1_jet_dR = deltaR(ak.unflatten(d1_eta, 1)*jet_shape, Events.Jet_eta[Jetcut], ak.unflatten(d1_phi, 1)*jet_shape, Events.Jet_phi[Jetcut])
        d2_jet_dR = deltaR(ak.unflatten(d2_eta, 1)*jet_shape, Events.Jet_eta[Jetcut], ak.unflatten(d2_phi, 1)*jet_shape, Events.Jet_phi[Jetcut])
        d1_jet_idx = ak.argmin(d1_jet_dR, axis=1)
        d2_jet_idx = ak.argmin(d2_jet_dR, axis=1)
        d1_jet = vector.zip({"eta": Events.Jet_eta[Jetcut][d1_jet_idx], "phi": Events.Jet_phi[Jetcut][d1_jet_idx], "pt": Events.Jet_pt[Jetcut][d1_jet_idx], "mass": Events.Jet_mass[Jetcut][d1_jet_idx]})
        d2_jet = vector.zip({"eta": Events.Jet_eta[Jetcut][d2_jet_idx], "phi": Events.Jet_phi[Jetcut][d2_jet_idx], "pt": Events.Jet_pt[Jetcut][d2_jet_idx], "mass": Events.Jet_mass[Jetcut][d2_jet_idx]})
        d1_jet_pnet = vector.zip({"eta": Events.Jet_eta[Jetcut][d1_jet_idx], "phi": Events.Jet_phi[Jetcut][d1_jet_idx], "pt": Events.Jet_pt[Jetcut][d1_jet_idx]*Events.Jet_PNetRegPtRawCorr[Jetcut][d1_jet_idx]*(1-Events.Jet_rawFactor[Jetcut][d1_jet_idx]), "mass": Events.Jet_mass[Jetcut][d1_jet_idx]})
        d2_jet_pnet = vector.zip({"eta": Events.Jet_eta[Jetcut][d2_jet_idx], "phi": Events.Jet_phi[Jetcut][d2_jet_idx], "pt": Events.Jet_pt[Jetcut][d2_jet_idx]*Events.Jet_PNetRegPtRawCorr[Jetcut][d2_jet_idx]*(1-Events.Jet_rawFactor[Jetcut][d2_jet_idx]), "mass": Events.Jet_mass[Jetcut][d2_jet_idx]})
        singleInfo["dijet_hmass_"+str(i)] = ak.flatten(ak.where(d1_jet_idx != d2_jet_idx, (d1_jet+d2_jet).mass, d1_jet.mass))
        singleInfo["dijet_hmass_pnet_"+str(i)] = ak.flatten(ak.where(d1_jet_idx != d2_jet_idx, (d1_jet_pnet+d2_jet_pnet).mass, d1_jet_pnet.mass))
        singleInfo["d1_jet_idx_"+str(i)] = ak.flatten(d1_jet_idx)
        singleInfo["d2_jet_idx_"+str(i)] = ak.flatten(d2_jet_idx)
        
        d1_pid = abs(Events.GenPart_pdgId[d1])
        d2_pid = abs(Events.GenPart_pdgId[d2])
        # ΔR_qq between the two daughter particles
        singleInfo["MC_drqq_"+str(i)] = ones * deltaR(d1_eta, d2_eta, d1_phi, d2_phi)[:,0]
        # check if daughter particles are b or c
        d1_is_b = ak.unflatten(d1_pid==5, 1, axis=0)
        d2_is_b = ak.unflatten(d2_pid==5, 1, axis=0)
        d1_is_c = ak.unflatten(d1_pid==4, 1, axis=0)
        d2_is_c = ak.unflatten(d2_pid==4, 1, axis=0)
        d1_is_s = ak.unflatten(d1_pid==3, 1, axis=0)
        d2_is_s = ak.unflatten(d2_pid==3, 1, axis=0)
        singleInfo["h_hadronic_"+str(i)] = ((d1_is_b & d2_is_b) | (d1_is_c & d2_is_c) | (d1_is_s & d2_is_s)) * ones
        # check if daughter particles are in PAIReD jet
        ones_ = ak.unflatten(ones, 1, axis=-1)
        d1_isIn = isInGenPart[ones_*d1]
        d2_isIn = isInGenPart[ones_*d2]
        # count how many lie inside (0, 1 or 2?)
        singleInfo["h_daughters_in_"+str(i)] = (1*d1_isIn + d2_isIn)
        # label the jet pairs:
        # - BB if both daughters are b quarks and inside the PAIReD jet
        # - CC if both daughters are c quarks and inside the PAIReD jet
        # - LL else
        singleInfo["label_BB_"+str(i)] = ( d1_is_b * d2_is_b )[:,:,0]
        singleInfo["label_CC_"+str(i)] = ( d1_is_c * d2_is_c )[:,:,0]
        singleInfo["MC_higgs_flav_"+str(i)] = ones * d1_pid[:,0]

    use_h_0 = ((singleInfo["h_daughters_in_0"] == 2) * (singleInfo["h_daughters_in_1"] == 0) * (singleInfo["h_hadronic_0"]))[:,:,0]
    use_h_1 = ((singleInfo["h_daughters_in_1"] == 2) * (singleInfo["h_daughters_in_0"] == 0) * (singleInfo["h_hadronic_1"]))[:,:,0]
    singleInfo["label_BB_0"] = singleInfo["label_BB_0"] * use_h_0
    singleInfo["label_BB_1"] = singleInfo["label_BB_1"] * use_h_1
    singleInfo["label_CC_0"] = singleInfo["label_CC_0"] * use_h_0
    singleInfo["label_CC_1"] = singleInfo["label_CC_1"] * use_h_1

    MCInfo["dijet_hmass"] = (singleInfo["dijet_hmass_0"] * use_h_0) + (singleInfo["dijet_hmass_1"] * use_h_1)
    MCInfo["dijet_hmass_pnet"] = (singleInfo["dijet_hmass_pnet_0"] * use_h_0) + (singleInfo["dijet_hmass_pnet_1"] * use_h_1)
    MCInfo["d1_jet_idx"] = (singleInfo["d1_jet_idx_0"] * use_h_0) + (singleInfo["d1_jet_idx_1"] * use_h_1)
    MCInfo["d2_jet_idx"] = (singleInfo["d2_jet_idx_0"] * use_h_0) + (singleInfo["d2_jet_idx_1"] * use_h_1)

    MCInfo["label_BB"] = singleInfo["label_BB_0"] | singleInfo["label_BB_1"]
    MCInfo["label_CC"] = (~MCInfo["label_BB"]) & (singleInfo["label_CC_0"] | singleInfo["label_CC_1"])

    MCInfo["MC_higgs_valid"] = ((singleInfo["h_daughters_in_0"] + singleInfo["h_daughters_in_1"]) < 3)[:,:,0]
    MCInfo["MC_higgs_valid"] = MCInfo["MC_higgs_valid"] * ((singleInfo["h_daughters_in_0"] == 0) | (singleInfo["h_hadronic_0"] > 0))[:,:,0]
    MCInfo["MC_higgs_valid"] = MCInfo["MC_higgs_valid"] * ((singleInfo["h_daughters_in_1"] == 0) | (singleInfo["h_hadronic_1"] > 0))[:,:,0]

    MCInfo["MC_higgs_mass"] = ((use_h_0) * singleInfo["MC_higgs_mass_0"]) + ((use_h_1) * singleInfo["MC_higgs_mass_1"])
    MCInfo["MC_higgs_pt"] = ((use_h_0) * singleInfo["MC_higgs_pt_0"]) + ((use_h_1) * singleInfo["MC_higgs_pt_1"])
    MCInfo["MC_higgs_eta"] = ((use_h_0) * singleInfo["MC_higgs_eta_0"]) + ((use_h_1) * singleInfo["MC_higgs_eta_1"])
    MCInfo["MC_higgs_phi"] = ((use_h_0) * singleInfo["MC_higgs_phi_0"]) + ((use_h_1) * singleInfo["MC_higgs_phi_1"])
    MCInfo["MC_drqq"] = ((use_h_0) * singleInfo["MC_drqq_0"]) + ((use_h_1) * singleInfo["MC_drqq_1"])
    MCInfo["MC_higgs_flav"] = ((use_h_0) * singleInfo["MC_higgs_flav_0"]) + ((use_h_1) * singleInfo["MC_higgs_flav_1"])

def processHiggs_VH(Events, isInGenPart, higgs_idx, MCInfo, Jetcut):
    # process fake higgs
    ones = ak.ones_like(isInGenPart[:,:,0])

    MCInfo["MC_higgs_pt"] = ones * Events.GenPart_pt[higgs_idx][:,0]
    MCInfo["MC_higgs_eta"] = ones * Events.GenPart_eta[higgs_idx][:,0]
    MCInfo["MC_higgs_phi"] = ones * Events.GenPart_phi[higgs_idx][:,0]
    MCInfo["MC_higgs_mass"] = ones * Events.GenPart_mass[higgs_idx][:,0]
    d1 = ak.argmax(1*(Events.GenPart_genPartIdxMother == higgs_idx[:,0]), axis=1, keepdims=True)
    d2 = ak.argmax(ak.local_index(Events.GenPart_pdgId, axis=1)*(Events.GenPart_genPartIdxMother == higgs_idx[:,0]), axis=1, keepdims=True)
    # kinematics of the daughter particles
    d1_eta = Events.GenPart_eta[d1]
    d2_eta = Events.GenPart_eta[d2]
    d1_phi = Events.GenPart_phi[d1]
    d2_phi = Events.GenPart_phi[d2]

    jet_shape = ak.ones_like(Events.Jet_eta[Jetcut])
    d1_jet_dR = deltaR(ak.unflatten(d1_eta, 1)*jet_shape, Events.Jet_eta[Jetcut], ak.unflatten(d1_phi, 1)*jet_shape, Events.Jet_phi[Jetcut])
    d2_jet_dR = deltaR(ak.unflatten(d2_eta, 1)*jet_shape, Events.Jet_eta[Jetcut], ak.unflatten(d2_phi, 1)*jet_shape, Events.Jet_phi[Jetcut])
    MCInfo["d1_jet_idx"] = ak.argmin(d1_jet_dR, axis=1)
    MCInfo["d2_jet_idx"] = ak.argmin(d2_jet_dR, axis=1)
    d1_jet = vector.zip({"eta": Events.Jet_eta[Jetcut][MCInfo["d1_jet_idx"]], 
                         "phi": Events.Jet_phi[Jetcut][MCInfo["d1_jet_idx"]], 
                         "pt": Events.Jet_pt[Jetcut][MCInfo["d1_jet_idx"]], 
                         "mass": Events.Jet_mass[Jetcut][MCInfo["d1_jet_idx"]]
                        })
    d2_jet = vector.zip({"eta": Events.Jet_eta[Jetcut][MCInfo["d2_jet_idx"]], 
                         "phi": Events.Jet_phi[Jetcut][MCInfo["d2_jet_idx"]], 
                         "pt": Events.Jet_pt[Jetcut][MCInfo["d2_jet_idx"]], 
                         "mass": Events.Jet_mass[Jetcut][MCInfo["d2_jet_idx"]]
                        })
    d1_jet_pnet = vector.zip({"eta": Events.Jet_eta[Jetcut][MCInfo["d1_jet_idx"]], 
                              "phi": Events.Jet_phi[Jetcut][MCInfo["d1_jet_idx"]], 
                              "pt": Events.Jet_pt[Jetcut][MCInfo["d1_jet_idx"]]*Events.Jet_PNetRegPtRawCorr[Jetcut][MCInfo["d1_jet_idx"]]*(1-Events.Jet_rawFactor[Jetcut][MCInfo["d1_jet_idx"]]), 
                              "mass": Events.Jet_mass[Jetcut][MCInfo["d1_jet_idx"]]
                            })
    d2_jet_pnet = vector.zip({"eta": Events.Jet_eta[Jetcut][MCInfo["d2_jet_idx"]], 
                              "phi": Events.Jet_phi[Jetcut][MCInfo["d2_jet_idx"]], 
                              "pt": Events.Jet_pt[Jetcut][MCInfo["d2_jet_idx"]]*Events.Jet_PNetRegPtRawCorr[Jetcut][MCInfo["d2_jet_idx"]]*(1-Events.Jet_rawFactor[Jetcut][MCInfo["d2_jet_idx"]]), 
                              "mass": Events.Jet_mass[Jetcut][MCInfo["d2_jet_idx"]]
                            })
    MCInfo["dijet_hmass"] = ak.where(MCInfo["d1_jet_idx"] != MCInfo["d2_jet_idx"], (d1_jet+d2_jet).mass, d1_jet.mass)
    MCInfo["dijet_hmass_pnet"] = ak.where(MCInfo["d1_jet_idx"] != MCInfo["d2_jet_idx"], (d1_jet_pnet+d2_jet_pnet).mass, d1_jet_pnet.mass)
    
    d1_pid = abs(Events.GenPart_pdgId[d1])
    d2_pid = abs(Events.GenPart_pdgId[d2])
    # ΔR_qq between the two daughter particles
    MCInfo["MC_drqq"] = ones * deltaR(d1_eta, d2_eta, d1_phi, d2_phi)[:,0]
    # check if daughter particles are b or c
    d1_is_b = ak.unflatten(d1_pid==5, 1, axis=0)
    d2_is_b = ak.unflatten(d2_pid==5, 1, axis=0)
    d1_is_c = ak.unflatten(d1_pid==4, 1, axis=0)
    d2_is_c = ak.unflatten(d2_pid==4, 1, axis=0)
    d1_is_s = ak.unflatten(d1_pid==3, 1, axis=0)
    d2_is_s = ak.unflatten(d2_pid==3, 1, axis=0)
    h_hadronic = ((d1_is_b & d2_is_b) | (d1_is_c & d2_is_c) | (d1_is_s & d2_is_s)) * ones
    # check if daughter particles are in PAIReD jet
    ones_ = ak.unflatten(ones, 1, axis=-1)
    d1_isIn = isInGenPart[ones_*d1]
    d2_isIn = isInGenPart[ones_*d2]
    # count how many lie inside (0, 1 or 2?)
    h_daughters_in = (1*d1_isIn + d2_isIn)
    use_h = ((h_daughters_in == 2) & (h_hadronic))[:,:,0]
    # label the jet pairs:
    MCInfo["label_BB"] = ( d1_is_b * d2_is_b )[:,:,0] * use_h
    MCInfo["label_CC"] = ( d1_is_c * d2_is_c )[:,:,0] * use_h
    MCInfo["MC_higgs_flav"] = ones * d1_pid[:,0]
    MCInfo["MC_higgs_valid"] = ((h_daughters_in == 0) | (h_hadronic > 0))[:,:,0]
    MCInfo["dijet_hmass"] = ak.flatten(MCInfo["dijet_hmass"]) * use_h
    MCInfo["dijet_hmass_pnet"] = ak.flatten(MCInfo["dijet_hmass_pnet"]) * use_h
    MCInfo["d1_jet_idx"] = ak.flatten(MCInfo["d1_jet_idx"]) * use_h
    MCInfo["d2_jet_idx"] = ak.flatten(MCInfo["d2_jet_idx"]) * use_h

    #print(MCInfo["label_BB"])
    #print(MCInfo["MC_higgs_valid"])
    lepton_indices = ak.local_index(Events.GenPart_pdgId, axis=1)[((abs(Events.GenPart_pdgId) == 11) | (abs(Events.GenPart_pdgId) == 13) | (abs(Events.GenPart_pdgId) == 15)) * (Events.GenPart_statusFlags%2 == 1) * ((Events.GenPart_statusFlags >> 12)%2 == 1)]
    if ak.all(ak.num(lepton_indices) == 2):
        #print("ITS WORKING DDONT WORRY")
        for i in range(2):
            lep_isIn = isInGenPart[ones_*lepton_indices[:,i]]
            # count how many lie inside (0, 1 or 2?)
            MCInfo["MC_higgs_valid"] = MCInfo["MC_higgs_valid"] * (~lep_isIn[:,:,0])

def processHiggs_DY(Events, isInGenPart, MCInfo):
    ones = ak.ones_like(isInGenPart[:,:,0])
    ones_ = ak.unflatten(ones, 1, axis=-1)
    MCInfo["MC_higgs_pt"] = ak.zeros_like(ones)
    MCInfo["MC_higgs_eta"] = ak.zeros_like(ones)
    MCInfo["MC_higgs_phi"] = ak.zeros_like(ones)
    MCInfo["MC_higgs_mass"] = ak.zeros_like(ones)
    MCInfo["MC_higgs_flav"] = ak.zeros_like(ones)
    MCInfo["MC_higgs_valid"] = ak.ones_like(ones)
    MCInfo["MC_drqq"] = ak.zeros_like(ones)
    MCInfo["label_BB"] = ak.zeros_like(ones)
    MCInfo["label_CC"] = ak.zeros_like(ones)
    MCInfo["dijet_hmass"] = ak.zeros_like(ones)
    MCInfo["dijet_hmass_pnet"] = ak.zeros_like(ones)
    MCInfo["d1_jet_idx"] = ak.zeros_like(ones)
    MCInfo["d2_jet_idx"] = ak.zeros_like(ones)
    # get indices of vector bosons
    lepton_indices = ak.local_index(Events.GenPart_pdgId, axis=1)[((abs(Events.GenPart_pdgId) == 11) | (abs(Events.GenPart_pdgId) == 13) | (abs(Events.GenPart_pdgId) == 15)) * (Events.GenPart_statusFlags%2 == 1) * ((Events.GenPart_statusFlags >> 12)%2 == 1)]
    #print('leptons', lepton_indices)
    #print(lepton_indices[:,0])
    #print(ones_*lepton_indices[:,0])
    for i in range(2):
        lep_isIn = isInGenPart[ones_*lepton_indices[:,i]]
        # count how many lie inside (0, 1 or 2?)
        MCInfo["MC_higgs_valid"] = MCInfo["MC_higgs_valid"] * (~lep_isIn[:,:,0])

def processHiggs_TT(Events, isInGenPart, MCInfo):
    ones = ak.ones_like(isInGenPart[:,:,0])
    MCInfo["MC_higgs_pt"] = ak.zeros_like(ones)
    MCInfo["MC_higgs_eta"] = ak.zeros_like(ones)
    MCInfo["MC_higgs_phi"] = ak.zeros_like(ones)
    MCInfo["MC_higgs_mass"] = ak.zeros_like(ones)
    MCInfo["MC_higgs_flav"] = ak.zeros_like(ones)
    MCInfo["MC_higgs_valid"] = ak.ones_like(ones)
    MCInfo["MC_drqq"] = ak.zeros_like(ones)
    MCInfo["label_BB"] = ak.zeros_like(ones)
    MCInfo["label_CC"] = ak.zeros_like(ones)
    MCInfo["dijet_hmass"] = ak.zeros_like(ones)
    MCInfo["dijet_hmass_pnet"] = ak.zeros_like(ones)
    MCInfo["d1_jet_idx"] = ak.zeros_like(ones)
    MCInfo["d2_jet_idx"] = ak.zeros_like(ones)

def getMCInfo(Events, isInGenPart, Jet_genJetIdx, Jetcut, physics_process):
    # check if there is a Higgs in the event
    isHiggs = ak.any(Events.GenPart_pdgId == 25, axis=1)
    if ak.all(isHiggs):
        local_idx = ak.local_index(Events.GenPart_pdgId, axis=1)
        hmother_idx = Events.GenPart_genPartIdxMother * (Events.GenPart_pdgId == 25)
        hdaughter_idx = ak.any(ak.broadcast_arrays(local_idx[:, None], hmother_idx)[0] == hmother_idx, axis=1)
        higgs_indices = ak.local_index(Events.GenPart_pdgId, axis=1)[(Events.GenPart_pdgId == 25)*(~hdaughter_idx)]
    elif ak.any(isHiggs):
        print("   !!! Skip batch as there are events with Higgs and events without Higgs !!!")
        return False
    # define information dictionary and set defaults
    ones = ak.ones_like(isInGenPart[:,:,0])
    MCInfo = {
        "MC_higgs_pt" : None,
        "MC_higgs_eta" : None,
        "MC_higgs_phi" : None,
        "MC_higgs_mass" : None,
        "MC_higgs_flav" : None,
        "MC_higgs_valid" : None,
        "MC_vector_flav" : ak.ones_like(ones),
        "MC_lepton_channel" : ak.ones_like(ones),
        "MC_gendijet_pt" : None,
        "MC_gendijet_eta" : None,
        "MC_gendijet_phi" : None,
        "MC_gendijet_mass" : None,
        "MC_genjet1_flav" : None,
        "MC_genjet2_flav" : None,
        "MC_genjet1_matched" : None,
        "MC_genjet2_matched" : None,
        "MC_drqq" : ak.zeros_like(ones),
        "MC_n_c" : ak.ones_like(ones),
        "label_BB" : None,
        "label_CC" : None,
        "label_bb" : ak.full_like(ones, False),
        "label_bx" : None,
        "label_cx" : None,
        "label_ll" : None,
        "label_LL": None,
        "dijet_hmass": None,
        "dijet_hmass_pnet": None,
        "d1_jet_idx": None,
        "d2_jet_idx": None
    }
    # get the higgs information
    if physics_process == 25: processHiggs_HH(Events, isInGenPart, higgs_indices, MCInfo, Jetcut)
    elif physics_process == 26: processHiggs_VH(Events, isInGenPart, higgs_indices, MCInfo, Jetcut)
    elif physics_process == 23: processHiggs_DY(Events, isInGenPart, MCInfo)
    elif physics_process == 66: processHiggs_TT(Events, isInGenPart, MCInfo)
    #print(ak.num(MCInfo["dijet_hmass"], axis=1))
    # get c/b quarks in GenParts
    # ensure that statusFlags bit 13 is True (isLastCopy)
    # reject quarks without mother particle (mostly pile-up)
    bit13 = (1 << 13)
    charms = ((abs(Events.GenPart_pdgId) == 4) & 
            ((Events.GenPart_statusFlags & bit13) == bit13) &
            (Events.GenPart_genPartIdxMother != -1))
    bottoms = ((abs(Events.GenPart_pdgId) == 5) & 
            ((Events.GenPart_statusFlags & bit13) == bit13) &
            (Events.GenPart_genPartIdxMother != -1))

    # make charms and bottoms broadcastable and count number of b/c quarks in PAIReD jet
    charms = ak.unflatten(charms, 1, axis=0)
    charms, isInGenPart = ak.broadcast_arrays(charms, isInGenPart)
    Ncharms = ak.sum(isInGenPart[charms], axis=2)

    bottoms = ak.unflatten(bottoms, 1, axis=0)
    bottoms, isInGenPart = ak.broadcast_arrays(bottoms, isInGenPart)
    Nbottoms = ak.sum(isInGenPart[bottoms], axis=2)

    # label the PAIReD jets according to the number of b/c quarks inside them:
    #  - label bb: if two or more b quarks inside PAIReD jet
    #  - label bl: if one b quark inside PAIReD jet
    #  - label cc: if no b quarks but two or more c quarks
    #  - label cl: if no b quarks but one or more c quarks
    #  - label ll: if no b and no c quarks
    
    MCInfo["label_bb"] = (~MCInfo["label_CC"]) & (~MCInfo["label_BB"]) & (Nbottoms > 1)
    MCInfo["label_LL"] = (~MCInfo["label_CC"]) & (~MCInfo["label_BB"])
    MCInfo["label_bx"] = (~MCInfo["label_CC"]) & (~MCInfo["label_BB"]) & (~MCInfo["label_bb"]) & (Nbottoms > 0)
    MCInfo["label_cx"] = (~MCInfo["label_CC"]) & (~MCInfo["label_BB"]) & (Nbottoms == 0) & (Ncharms > 0)
    MCInfo["label_ll"] = (~MCInfo["label_CC"]) & (~MCInfo["label_BB"]) & (Nbottoms == 0) & (Ncharms == 0)

    # check if there is a vector boson in the event
    hasVector = ak.any((abs(Events.GenPart_pdgId) == 23) | (abs(Events.GenPart_pdgId) == 24), axis=1)

    # index of the vector boson (W/Z)
    vector_idx = ak.argmax(ak.local_index(Events.GenPart_pdgId, axis=1)*((abs(Events.GenPart_pdgId) == 23) | (abs(Events.GenPart_pdgId) == 24)), axis=1, keepdims=True)

    # indices of the daughter particles of the vector boson
    vector_d1 = ak.argmax(1*(Events.GenPart_genPartIdxMother == vector_idx), axis=1, keepdims=True)
    vector_d2 = ak.argmax(ak.local_index(Events.GenPart_pdgId, axis=1)*(Events.GenPart_genPartIdxMother == vector_idx), axis=1, keepdims=True)

    # PID of the daughters of the vector boson
    vector_d1_pid = abs(Events.GenPart_pdgId[vector_d1])[:,0]
    vector_d2_pid = abs(Events.GenPart_pdgId[vector_d2])[:,0]

    # generated vector flavor:
    # - digit represents the flavor of the daughter particles
    #   (11 for e-like,  13 for mu-like, 15 for tau-like)
    # - if there is no vector boson, = 0
    vector_flavs = (11 * ((vector_d1_pid==11) | (vector_d1_pid==12)) +
                   13 * ((vector_d1_pid==13) | (vector_d1_pid==14)) +
                   15 * ((vector_d1_pid==15) | (vector_d1_pid==16))) * hasVector
    MCInfo["MC_vector_flav"] = MCInfo["MC_vector_flav"] * vector_flavs

    # lepton channel:
    # - number of charged leptons the vector boson decays into (0L, 1L, 2L)
    # - if there is no vector boson, = -1
    lepton_channels = (1* ((vector_d1_pid==11) | (vector_d1_pid==13) | (vector_d1_pid==15)) + 
                      1* ((vector_d2_pid==11) | (vector_d2_pid==13) | (vector_d2_pid==15))) * hasVector - 1*(~hasVector)
    MCInfo["MC_lepton_channel"] = MCInfo["MC_lepton_channel"] * lepton_channels


    # add the four-vectors of the two generated jets of the ellipse 
    # (matched to the two underlying AK4 jets)
    # and save the combined invariant mass and kinematics
    # 4-vector of the two generated jets
    GenJet1 = vector.zip({
        "eta":  Events.GenJet_eta[Jet_genJetIdx.j1], 
        "phi":  Events.GenJet_phi[Jet_genJetIdx.j1],
        "pt":   Events.GenJet_pt[Jet_genJetIdx.j1], 
        "mass": Events.GenJet_mass[Jet_genJetIdx.j1]
    })
    GenJet2 = vector.zip({
        "eta":  Events.GenJet_eta[Jet_genJetIdx.j2], 
        "phi":  Events.GenJet_phi[Jet_genJetIdx.j2],
        "pt":   Events.GenJet_pt[Jet_genJetIdx.j2], 
        "mass": Events.GenJet_mass[Jet_genJetIdx.j2]
    })

    # combined 4-vector
    GenDijet4 = GenJet1 + GenJet2

    # check if jets have matched GenJet
    MCInfo["MC_genjet1_matched"] = (Jet_genJetIdx.j1 != -1)
    MCInfo["MC_genjet2_matched"] = (Jet_genJetIdx.j2 != -1)

    # save kinematics and set them to zero for PAIReD jets
    # in which at least one Jet is not matched to a GenJet
    bothJetHaveGenJet = MCInfo["MC_genjet1_matched"] * MCInfo["MC_genjet2_matched"]
    MCInfo["MC_gendijet_eta"]  = GenDijet4.eta  * bothJetHaveGenJet
    MCInfo["MC_gendijet_phi"]  = GenDijet4.phi  * bothJetHaveGenJet
    MCInfo["MC_gendijet_pt"]   = GenDijet4.pt   * bothJetHaveGenJet
    MCInfo["MC_gendijet_mass"] = GenDijet4.mass * bothJetHaveGenJet

    # get parton flavour of matched GenJet
    MCInfo["MC_genjet1_flav"] = Events.GenJet_partonFlavour[Jet_genJetIdx.j1] * MCInfo["MC_genjet1_matched"]
    MCInfo["MC_genjet2_flav"] = Events.GenJet_partonFlavour[Jet_genJetIdx.j2] * MCInfo["MC_genjet2_matched"]

    # number of c jets in event
    num_cjets = ak.sum(abs(Events.Jet_partonFlavour[Jetcut])==4, axis=1)  # take selection cut into account
    MCInfo["MC_n_c"] = MCInfo["MC_n_c"] * num_cjets

    # vectorize the scalars for the number jet pairs
    # and define dictionary with all the outputs
    
    # return the summarized MC information
    return MCInfo


