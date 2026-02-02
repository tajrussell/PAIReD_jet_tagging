"""
 * Function: isInPAIReD
 * ---------------------
Determines which of the input particles are inside the PAIReD jet spanned by the
given jet pair.

Parameters
----------
eta_j1 : float or array of floats
    Pseudorapidity of the first jet, or jets if array is given
eta_j2 : float or array of floats
    Pseudorapidity of the second jet, or jets if array is given
phi_j1 : float or array of floats
    Phi angle of the first jet, or jets if array is given
phi_j2 : float or array of floats
    Phi angle of the second jet, or jets if array is given
eta_p : float or array of floats
    Pseudorapidity of the particle, or particles if array is given
phi_p : float or array of floats
    Phi angle of the particle, or particles if array is given

Returns
-------
isInPAIReD : bool (array)
    Boolean holding the information whether the particle is inside the PAIReD jet
    or not. Shape: (N_events, N_jetpairs, N_particles)
"""

# Include dependencies
import awkward as ak
import numpy as np
from tools.helpers import deltaR, deltaPhi, Phi_mpi_pi, shiftPhi

def isInPAIReD(eta_j1, eta_j2, phi_j1, phi_j2, eta_p, phi_p):

    # in order to be able to broadcast the arrays with the particle arrays to
    # the format
    # (N_events, N_jetpairs, N_particles),
    # we have to add one dimension to the jet related data. So, we have to add a
    # one:
    # (N_events, N_jetpairs) -> (N_events, N_jetpairs, 1)
    eta_j1 = ak.unflatten(eta_j1, 1, axis=-1)
    eta_j2 = ak.unflatten(eta_j2, 1, axis=-1)
    phi_j1 = ak.unflatten(phi_j1, 1, axis=-1)
    phi_j2 = ak.unflatten(phi_j2, 1, axis=-1)
    # and add one to the particle related data:
    # (N_events, N_particles) -> (N_events, 1, N_particles)
    eta_p = ak.unflatten(eta_p, 1, axis=0)
    phi_p = ak.unflatten(phi_p, 1, axis=0)

    # center point
    phi_c = Phi_mpi_pi(phi_j1 + Phi_mpi_pi(phi_j2-phi_j1) / 2)

    # in order to treat the jump in the phi space from phi=pi to phi=-pi
    # correctly, the phi values of the particles might have to be shifted:
    #if phi_c >= 0:
    #    shift = 2*np.pi
    #    lim = phi_c - np.pi
    #else:
    #    shift = -2*np.pi
    #    lim = phi_c + np.pi

    # define shift and lim according to the upper if statement
    shift = 2*np.pi * (phi_c>=0) - 2*np.pi * (phi_c<0)
    lim = phi_c - np.pi * (phi_c>=0) + np.pi * (phi_c<0)

    # eventually shift phi_p
    phi_p = shiftPhi(phi_p, shift, lim)

    # distance in between the particle and jet 1
    d_1p = deltaR(eta_j1, eta_p, phi_j1, phi_p)

    # distance in between the particle and jet 2
    d_2p = deltaR(eta_j2, eta_p, phi_j2, phi_p)


    radius = 0.4

    # condition 1 for circle around jet 1
    cond_1 = d_1p <= radius

    # condition 2 for circle around jet 2
    cond_2 = d_2p <= radius

    # is in AKJets if either condition 1 or condition 2 is True
    selected = cond_1 | cond_2

    return selected