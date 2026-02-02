"""
 * Function: isInPAIReD
 * ---------------------
Determines which of the input particles are inside the ellipse spanned by the
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
    Boolean holding the information whether the particle is inside the ellipse
    or not. Shape: (N_events, N_jetpairs, N_particles)
"""

# Include dependencies
import awkward as ak
import numpy as np
from tools.helpers import deltaR, deltaPhi, Phi_mpi_pi, shiftPhi

def isInPAIReD(eta_j1, eta_j2, phi_j1, phi_j2, eta_p, phi_p, return_semimajor=False):
    
    semimajoradd = 1.

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

    # distance in between the jets 1 and 2
    d_12 = deltaR(eta_j1, eta_j2, phi_j1, phi_j2)

    # semiminor and semimajor axis of the enveloping ellipse
    semiminor = 1.5
    semimajor = np.maximum(semiminor, d_12/2 + semimajoradd)

    # distance of focus points to center
    focus = np.sqrt(semimajor**2 - semiminor**2)

    # center point
    eta_c = (eta_j1 + eta_j2) / 2
    phi_c = Phi_mpi_pi(phi_j1 + Phi_mpi_pi(phi_j2-phi_j1) / 2)

    # focus 1
    eta_f1 = eta_c + focus / (d_12/2) * (eta_j1-eta_c)
    phi_f1 = Phi_mpi_pi(phi_c + focus/(d_12/2) * deltaPhi(phi_j1, phi_c))

    # focus 2
    eta_f2 = eta_c + focus / (d_12/2) * (eta_j2-eta_c)
    phi_f2 = Phi_mpi_pi(phi_c + focus/(d_12/2) * deltaPhi(phi_j2, phi_c))

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

    # Distances to the focus points
    d_f1 = np.hypot(eta_p-eta_f1, phi_p-shiftPhi(phi_f1, shift, lim))
    d_f2 = np.hypot(eta_p-eta_f2, phi_p-shiftPhi(phi_f2, shift, lim))
    d_sum = d_f1 + d_f2

    # if particle is in the ellipse,
    # the sum of the distances will be less than 2*semimajor
    if return_semimajor: return d_sum < 2*semimajor, semimajor*2#+d_12/2
    else: return d_sum < 2*semimajor