"""helpers.py

This python file contains the following helpers functions for PFNano_to_PAIReD:

    * ak_arange
    * getJetCombinationIndices
    * getAllJetCombinationIndices
    * Phi_mpi_pi
    * deltaPhi
    * deltaR
    * shiftPhi

It requires the following modules to be installed in your python environment:

    * numpy
    * vector
    * awkward
    * itertools

"""

# Dependencies
import numpy as np
import vector
import awkward as ak
import itertools
from tools.branchnames import BranchNames

# ******************************************************************************

"""
 * Function: ak_arange
 * -------------------
Returns an awkward array with indices from 0 to a specified number for each row.
E.g.:
>>> ak_arange([1,3,0,2])
    [[0],
     [0, 1, 2],
     [],
     [0, 1]]

Parameters
----------
lens : list of int
    Number of indices per row

Returns
-------
IndicesList : array
    List of indices
"""

def ak_arange(lens):
    # rectangular array with indices from 0 to max(lens)-1
    arr = np.ones((len(lens),max(lens))) * np.arange(max(lens))
    # mask for selected indices in each row
    mask = arr < ak.unflatten(lens, 1)
    # return only selected ones in form of an awkward array
    return ak.drop_none(ak.mask(arr, mask))

# ******************************************************************************


"""
 * Function: getJetCombinationIndices
 * ----------------------------------
Returns a list of index pairs containing all possible jet combinations. E.g. for
N=3 jets, it will return [[0, 1], [0, 2], [1, 2]]. For N=0 or 1, it returns an
empty list [].

Parameters
----------
N_Jet : int
    Number of jets in the event

Returns
-------
IndicesList : list
    List of possible index pairs
"""

def getJetCombinationIndices(N_Jet):
    # create list of indice tuples:
    # e.g. for N_Jet=3: [(0, 1), (0, 2), (1, 2)]
    IndicesList = list(itertools.combinations(np.arange(N_Jet), 2))
    # return the same list as a list of lists:
    # e.g. [[0, 1], [0, 2], [1, 2]]
    return [list(ele) for ele in IndicesList]

# ******************************************************************************


"""
 * Function: getAllJetCombinationIndices
 * -------------------------------------
Returns a list of lists of index pairs containing all possible jet combinations
for all possible numbers of jets up to Nmax_Jet.
E.g. for Nmax=3 jets, it will return
    [ [],                          # 0 jets
      [],                          # 1 jets
      [[0, 1]],                    # 2 jets
      [[0, 1], [0, 2], [1, 2]]     # 3 jets = Nmax jets
    ]
The idea is to create a look-up table at the beginning, instead of running the
getJetCombinationIndices function for each event individually.

Parameters
----------
Nmax_Jet : int
    Maximum number of jets in one event accross all events in the tree

Returns
-------
IndicesList : list
    List of lists of possible index pairs
"""

def getAllJetCombinationIndices(Nmax_Jet):
    # create empty list
    IndicesList = []
    # loop over all possible numbers of jets
    for n in range(Nmax_Jet + 1):
        IndicesList.append(getJetCombinationIndices(n))
    return IndicesList

# ******************************************************************************


"""
 * Function: Phi_mpi_pi
 * --------------------
Outputs the angle in the interval [-pi, pi).

Parameters
----------
angle : float or array of floats
    Angle [in radians]

Returns
-------
angle : float or array of floats
    Angle in [-pi, pi)
"""

def Phi_mpi_pi(angle):
    # angle in between 0 and 2 pi:
    angle_0_2pi = angle % (2*np.pi)
    # return angle in between -pi and pi
    # therefore, subtract 2 pi if angle is in [pi, 2pi)
    return angle_0_2pi - 2*np.pi * (angle_0_2pi >= np.pi)

# ******************************************************************************


"""
 * Function: deltaPhi
 * ------------------
Computes the angular difference between two angles in radians.
The result lies in the interval [-pi, pi).

Parameters
----------
phi1 : float or array of floats
    First angle [in radians]
phi2 : float or array of floats
    Second angle [in radians]

Returns
-------
deltaPhi : float or array of floats
    Angular difference of the two angles
"""

def deltaPhi(phi1, phi2):
    return Phi_mpi_pi(phi1 - phi2)

# ******************************************************************************


"""
 * Function: deltaR
 * ----------------
Computes the distance between two jets in the eta-phi space.

Parameters
----------
eta1 : float or array of floats
    Pseudorapidity of the first jet
eta2 : float or array of floats
    Pseudorapidity of the second jet
phi1 : float or array of floats
    Phi angle of the first jet
phi2 : float or array of floats
    Phi angle of the second jet

Returns
-------
deltaR : float or array of floats
    Distance between two jets in the eta-phi space:
        deltaR = sqrt((eta1 - eta2)^2 + (phi1 - phi2)^2)
"""

def deltaR(eta1, eta2, phi1, phi2):
    deta = eta1 - eta2
    dphi = deltaPhi(phi1, phi2)
    return np.hypot(deta, dphi)

# ******************************************************************************


"""
 * Function: isoLeptonCut
 * ----------------------
Defines a selection cut for all jets in events according to the following 
condition: It finds all isolated leptons (muons and electrons) in each event, 
and checks which jets coincide with those isoLeptons. Since those jets likely
originated from these isoLeptons, we discard them later on.

Parameters
----------
Events : event tree object
    Object containing all information of multiple simulated events.

Returns
-------
cut : array of bools
    Boolean for each jet saying whether it came from an isoLepton
"""

def isoLeptonCut(Events):
    # define collection of isolated leptons
    isoMuons = (Events.Muon_pt > 12) & (abs(Events.Muon_eta) < 2.4) & (Events.Muon_tightId == 1) & (Events.Muon_pfIsoId >= 4)
    isoElectrons = (Events.Electron_pt > 20) & (abs(Events.Electron_eta) < 2.5) & (Events.Electron_mvaIso_WP90 == 1)

    # check which jets originate from these isoLeptons and discard them
    # cut jets from isoMuons
    cut = ak.all(deltaR(ak.unflatten(Events.Jet_eta, 1, axis=1),
                        ak.unflatten(Events.Muon_eta[isoMuons], 1, axis=0),
                        ak.unflatten(Events.Jet_phi, 1, axis=1),
                        ak.unflatten(Events.Muon_phi[isoMuons], 1, axis=0)) > 0.4, axis=2)
    # cut jets from isoElectrons
    cut = cut & ak.all(deltaR(ak.unflatten(Events.Jet_eta, 1, axis=1),
                            ak.unflatten(Events.Electron_eta[isoElectrons], 1, axis=0),
                            ak.unflatten(Events.Jet_phi, 1, axis=1),
                            ak.unflatten(Events.Electron_phi[isoElectrons], 1, axis=0)) > 0.4, axis=2)
    
    # return True for isNotIsoLep and False for isIsoLep
    return cut

# ******************************************************************************


"""
 * Function: shiftPhi
 * ------------------
Shifts the input angle based on some conditions:
    if shift == 0:
        return phi
    elif shift > 0:
        if phi < lim:
            return phi+shift
    else:
        if phi > lim:
            return phi+shift
    return phi

Parameters
----------
phi : float or array of floats
    Angle in radians
shift : float or array of floats
    Potential shift dPhi of the angle in radians
lim : float or array of floats
    Limit for the condition of the shift

Returns
-------
phi : float or array of floats
    Shifted angle
"""

def shiftPhi(phi, shift, lim):
    return phi + shift * (shift>0) * (phi<lim) + shift * (shift<0) * (phi>lim)

# ******************************************************************************


