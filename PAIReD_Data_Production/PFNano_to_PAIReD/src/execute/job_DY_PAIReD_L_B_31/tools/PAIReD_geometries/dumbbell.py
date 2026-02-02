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
    eta_c = (eta_j1 + eta_j2) / 2
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
    phi_p_shifted = shiftPhi(phi_p, shift, lim)
    phi_j1_shifted = shiftPhi(phi_j1, shift, lim)
    phi_j2_shifted = shiftPhi(phi_j2, shift, lim)


    # translate data, so that midpoint between jets lies in the origin
    eta_trans = (eta_j1 + eta_j2)/2
    phi_trans = (phi_j1_shifted + phi_j2_shifted)/2

    eta_p_trans = eta_p - eta_trans
    phi_p_trans = phi_p_shifted - phi_trans

    eta_j1_trans = eta_j1 - eta_trans
    phi_j1_trans = phi_j1_shifted - phi_trans

    eta_j2_trans = eta_j2 - eta_trans
    phi_j2_trans = phi_j2_shifted - phi_trans

    # there are 4 cases for rotations, so that jet 1 is left from jet 2
    # it depends on which value of which jet is greater 

    # rotation 1: eta_j1 > eta_j2, phi_j1 > phi_j2
    cond_rot_1_1 = eta_j1 >= eta_j2
    cond_rot_1_2 = phi_j1_shifted >= phi_j2_shifted
    cond_rot_1 = cond_rot_1_1 & cond_rot_1_2

    # rotation 2: eta_j1 > eta_j2, phi_j1 < phi_j2
    cond_rot_2_1 = eta_j1 >= eta_j2
    cond_rot_2_2 = phi_j1_shifted <= phi_j2_shifted
    cond_rot_2 = cond_rot_2_1 & cond_rot_2_2

    # rotation 3: eta_j1 < eta_j2, phi_j1 > phi_j2
    cond_rot_3_1 = eta_j1 <= eta_j2
    cond_rot_3_2 = phi_j1_shifted >= phi_j2_shifted
    cond_rot_3 = cond_rot_3_1 & cond_rot_3_2

    # rotation 4: eta_j1 < eta_j2, phi_j1 < phi_j2
    cond_rot_4_1 = eta_j1 <= eta_j2
    cond_rot_4_2 = phi_j1_shifted <= phi_j2_shifted
    cond_rot_4 = cond_rot_4_1 & cond_rot_4_2

    # overall angle
    alpha = np.where(eta_j1 != eta_j2, np.arctan((phi_j1_shifted - phi_j2_shifted)/(eta_j1 - eta_j2)), 1)
    
    # angle for each rotation
    angle_rad = cond_rot_1 * (np.pi - alpha) + cond_rot_2 * (-1) * (np.pi + alpha) + cond_rot_3 * (-1) * alpha + cond_rot_4 * (2*np.pi - alpha)
    
    # particle related data
    eta_p_rot = eta_p_trans * np.cos(angle_rad) - phi_p_trans * np.sin(angle_rad)
    phi_p_rot = eta_p_trans * np.sin(angle_rad) + phi_p_trans * np.cos(angle_rad)

    # both jets
    eta_j1_rot = eta_j1_trans * np.cos(angle_rad) - phi_j1_trans * np.sin(angle_rad)
    phi_j1_rot = eta_j1_trans * np.sin(angle_rad) + phi_j1_trans * np.cos(angle_rad)

    eta_j2_rot = eta_j2_trans * np.cos(angle_rad) - phi_j2_trans * np.sin(angle_rad)
    phi_j2_rot = eta_j2_trans * np.sin(angle_rad) + phi_j2_trans * np.cos(angle_rad)


    # shift eta in order to be achieve -pi <= eta <= pi
    # if eta <= -pi
    eta_p_rot = eta_p_rot + 2 * np.pi * (eta_p_rot <= -np.pi)
    eta_j1_rot = eta_j1_rot + 2 * np.pi * (eta_j1_rot <= -np.pi)
    eta_j2_rot = eta_j2_rot + 2 * np.pi * (eta_j2_rot <= -np.pi)
    # if eta >= pi
    eta_p_rot = eta_p_rot - 2 * np.pi * (eta_p_rot >= np.pi)
    eta_j1_rot = eta_j1_rot - 2 * np.pi * (eta_j1_rot >= np.pi)
    eta_j2_rot = eta_j2_rot - 2 * np.pi * (eta_j2_rot >= np.pi)


    # distance in between the jets 1 and 2
    d_12 = deltaR(eta_j1_rot, eta_j2_rot, phi_j1_rot, phi_j2_rot)

    # distance in between the particle and jet 1
    d_1p = deltaR(eta_j1_rot, eta_p_rot, phi_j1_rot, phi_p_rot)

    # distance in between the particle and jet 1
    d_2p = deltaR(eta_j2_rot, eta_p_rot, phi_j2_rot, phi_p_rot)


    # parameter for circles and rectangle
    # x is the shortest distance from jet to rectangle
    radius = 0.4 # d_12 / 4
    height = 0.8 # d_12 / 5
    theta = np.arcsin(height/(2*radius))
    x = radius * np.cos(theta)
    width = d_12 - 2 * x


    # condition 1 (circle around jet 1)
    cond_1_1 = eta_p_rot <= eta_j1_rot + x
    cond_1_2 = d_1p <= radius
    
    cond_1 = cond_1_1 & cond_1_2
    # print(cond_1)

    # condition 2 (rectangle between jets)
    cond_2_1 = eta_p_rot >= eta_j1_rot + x
    cond_2_2 = eta_p_rot <= eta_j2_rot - x
    cond_2_3 = phi_p_rot <= phi_j1_rot + height/2
    cond_2_4 = phi_p_rot >= phi_j1_rot - height/2

    cond_2 = cond_2_1 & cond_2_2 & cond_2_3 & cond_2_4
    # print(cond_2)
    
    # condition 3 (circle around jet 2)
    cond_3_1 = eta_p_rot >= eta_j2_rot - x
    cond_3_2 = d_2p <= radius
    
    cond_3 = cond_3_1 & cond_3_2
    # print(cond_3)

    selected = cond_1 | cond_2 | cond_3

    return selected