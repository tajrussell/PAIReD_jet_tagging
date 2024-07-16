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

# Include dependencies
import awkward as ak
import numpy as np
import vector
from tools.helpers import deltaR

def getMCInfo(Events, isInGenPart, Jet_genJetIdx, Jetcut):

    # check if there is a Higgs in the event
    isHiggs = ak.any(Events.GenPart_pdgId == 25, axis=1)

    # define array of ones that matches the wanted ouput shape
    ones = ak.ones_like(isInGenPart[:,:,0])

    # if Higgs event
    if ak.all(isHiggs):
        # index of the Higgs boson
        higgs_idx = ak.argmax(ak.local_index(Events.GenPart_pdgId, axis=1)*(Events.GenPart_pdgId == 25), axis=1, keepdims=True)

        # Higgs properties
        higgs_pt   = Events.GenPart_pt[higgs_idx][:,0]
        higgs_eta  = Events.GenPart_eta[higgs_idx][:,0]
        higgs_phi  = Events.GenPart_phi[higgs_idx][:,0]
        higgs_mass = Events.GenPart_mass[higgs_idx][:,0]

        # indices of the daughter particles of the Higgs
        d1 = ak.argmax(1*(Events.GenPart_genPartIdxMother == higgs_idx), axis=1, keepdims=True)
        d2 = ak.argmax(ak.local_index(Events.GenPart_pdgId, axis=1)*(Events.GenPart_genPartIdxMother == higgs_idx), axis=1, keepdims=True)

        # kinematics of the daughter particles
        d1_eta = Events.GenPart_eta[d1]
        d2_eta = Events.GenPart_eta[d2]
        d1_phi = Events.GenPart_phi[d1]
        d2_phi = Events.GenPart_phi[d2]
        d1_pid = abs(Events.GenPart_pdgId[d1])
        d2_pid = abs(Events.GenPart_pdgId[d2])

        # ΔR_qq between the two daughter particles
        drqq = deltaR(d1_eta, d2_eta, d1_phi, d2_phi)[:,0]

        # check if daughter particles are b or c
        d1_is_b = ak.unflatten(d1_pid==5, 1, axis=0)
        d2_is_b = ak.unflatten(d2_pid==5, 1, axis=0)
        d1_is_c = ak.unflatten(d1_pid==4, 1, axis=0)
        d2_is_c = ak.unflatten(d2_pid==4, 1, axis=0)

        # check if daughter particles are in PAIReD jet
        ones_ = ak.unflatten(ones, 1, axis=-1)
        d1_isIn = isInGenPart[ones_*d1]
        d2_isIn = isInGenPart[ones_*d2]

        # count how many lie inside (0, 1 or 2?)
        Ndaughters_isIn = (1*d1_isIn + d2_isIn)

        # label the jet pairs:
        # - BB if both daughters are b quarks and inside the PAIReD jet
        # - CC if both daughters are c quarks and inside the PAIReD jet
        # - LL else
        label_BB = ( d1_is_b * d2_is_b * (Ndaughters_isIn == 2) )[:,:,0]
        label_CC = ( d1_is_c * d2_is_c * (Ndaughters_isIn == 2) )[:,:,0]
        label_LL = (~label_CC) & (~label_BB)

        # generated higgs flavor:
        # - digit represents the flavor of the daughter particle 
        #   (e.g., 1 for d, 3 for s, or 5 for b)
        higgs_flav = d1_pid[:,0]

    # if Higgs but not all
    elif ak.any(isHiggs):
        print("   !!! Skip batch as there are events with Higgs and events without Higgs !!!")
        return False

    # if not a Higgs event
    else:
        # get original partons and gluons from the collision (no mother particle
        # and PDG Id of 21 or less)
        partons = (Events.LHEPart_status==1) & ((abs(Events.LHEPart_pdgId)<=5) | (abs(Events.LHEPart_pdgId)==21))
        
        # check which events have less than two initial partons
        less_than_two_partons = ak.sum(partons, axis=1) < 2

        # get two partons of highest pT
        s = ak.argsort(Events.LHEPart_pt[partons], axis=1, ascending=False)[:, :2]
        partons4 = vector.zip({"phi" : Events.LHEPart_phi[partons][s],
                                "eta" : Events.LHEPart_eta[partons][s],
                                "pt"  : Events.LHEPart_pt[partons][s],
                                "mass": Events.LHEPart_mass[partons][s]})

        # combine two highest pT partons to a "fake Higgs" / proxy Higgs
        fakeHiggs4 = vector.zip({"px": ak.sum(partons4.px, axis=1),
                                "py": ak.sum(partons4.py, axis=1),
                                "pz": ak.sum(partons4.pz, axis=1),
                                "E":  ak.sum(partons4.E, axis=1)})

        # set higgs_mass to -1 for events with less than two initial partons
        higgs_pt = ak.nan_to_num(fakeHiggs4.pt)*(~less_than_two_partons) - 1*less_than_two_partons
        higgs_eta = ak.nan_to_num(fakeHiggs4.eta)*(~less_than_two_partons)
        higgs_phi = ak.nan_to_num(fakeHiggs4.phi)*(~less_than_two_partons)
        higgs_mass = ak.nan_to_num(fakeHiggs4.mass)*(~less_than_two_partons) - 1*less_than_two_partons

        # set drqq and higgs_flav to default=0
        drqq = ones*0
        higgs_flav = ones*0

        # set signal labels to False and background label to True
        label_LL = ones == 1
        label_CC = ~label_LL
        label_BB = ~label_LL


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
    label_bb = (label_LL) & (Nbottoms >= 2)
    label_bl = (label_LL) & (Nbottoms == 1)
    label_cc = (label_LL) & (Nbottoms == 0) & (Ncharms >= 2)
    label_cl = (label_LL) & (Nbottoms == 0) & (Ncharms == 1)
    label_ll = (label_LL) & (Nbottoms == 0) & (Ncharms == 0)

    
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
    vector_flav = (11 * ((vector_d1_pid==11) | (vector_d1_pid==12)) +
                   13 * ((vector_d1_pid==13) | (vector_d1_pid==14)) +
                   15 * ((vector_d1_pid==15) | (vector_d1_pid==16))) * hasVector

    # lepton channel:
    # - number of charged leptons the vector boson decays into (0L, 1L, 2L)
    # - if there is no vector boson, = -1
    lepton_channel = (1* ((vector_d1_pid==11) | (vector_d1_pid==13) | (vector_d1_pid==15)) + 
                      1* ((vector_d2_pid==11) | (vector_d2_pid==13) | (vector_d2_pid==15))) * hasVector - 1*(~hasVector)


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
    genjet1_matched = (Jet_genJetIdx.j1 != -1)
    genjet2_matched = (Jet_genJetIdx.j2 != -1)

    # save kinematics and set them to zero for PAIReD jets
    # in which at least one Jet is not matched to a GenJet
    bothJetHaveGenJet = genjet1_matched * genjet2_matched
    gendijet_eta  = GenDijet4.eta  * bothJetHaveGenJet
    gendijet_phi  = GenDijet4.phi  * bothJetHaveGenJet
    gendijet_pt   = GenDijet4.pt   * bothJetHaveGenJet
    gendijet_mass = GenDijet4.mass * bothJetHaveGenJet

    # get parton flavour of matched GenJet
    genjet1_flav = Events.GenJet_partonFlavour[Jet_genJetIdx.j1] * genjet1_matched
    genjet2_flav = Events.GenJet_partonFlavour[Jet_genJetIdx.j2] * genjet2_matched

    # number of c jets in event
    n_c = ak.sum(abs(Events.Jet_partonFlavour[Jetcut])==4, axis=1)  # take selection cut into account

    # vectorize the scalars for the number jet pairs
    # and define dictionary with all the outputs
    MCInfo = {
        "MC_higgs_pt" : higgs_pt * ones,
        "MC_higgs_eta" : higgs_eta * ones,
        "MC_higgs_phi" : higgs_phi * ones,
        "MC_higgs_mass" : higgs_mass * ones,
        "MC_higgs_flav" : higgs_flav * ones,
        "MC_vector_flav" : vector_flav * ones,
        "MC_lepton_channel" : lepton_channel * ones,
        "MC_gendijet_pt" : gendijet_pt,
        "MC_gendijet_eta" : gendijet_eta,
        "MC_gendijet_phi" : gendijet_phi,
        "MC_gendijet_mass" : gendijet_mass,
        "MC_genjet1_flav" : genjet1_flav,
        "MC_genjet2_flav" : genjet2_flav,
        "MC_genjet1_matched" : genjet1_matched,
        "MC_genjet2_matched" : genjet2_matched,
        "MC_drqq" : drqq * ones,
        "MC_n_c" : n_c * ones,
        "label_BB" : label_BB,
        "label_CC" : label_CC,
        "label_LL" : label_LL,
        "label_bb" : label_bb,
        "label_bl" : label_bl,
        "label_cc" : label_cc,
        "label_cl" : label_cl,
        "label_ll" : label_ll
    }

    # return the summarized MC information
    return MCInfo