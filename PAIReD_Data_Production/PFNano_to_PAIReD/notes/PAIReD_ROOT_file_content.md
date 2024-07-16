# Overview of the file structure and content of the PAIReD output files
This note gives you a short overview of the structure (ROOT trees) and content (ROOT branches) of the produced output files of the [`makeNtuplesPAIReDjointMC()`](../src/processFileToPAIReD.py#L76) function. Those output files are meant to be used for training of the PAIReD tagger and therefore consist of both observables and true MC information.

There is only one tree in the ROOT file:
* `tree`: containing all the PAIReD jets with all the relevant information for training.

## Some notes on PAIReD jets

### Used seed jets
Note that the jets which are used as seeds for the PAIReD jet construction are AK4 jets meet the following criteria:
* raw tranverse momentum $p_T^{raw} = (1-rawFactor)\cdot p_T > 20~\text{GeV}$
* $|\eta| < 2.5$
* $JetId > 4$
* jet does not come from an isolated lepton (electron/muon)

### Labeling scheme
The PAIReD jets are assigned a truth label by the following scheme according to their content.
|              | Original  PAIReD 3 labels | Adapted PAIReD 6  labels | The PAIReDEllipse jet covers...                                   | Used in training |
|--------------|:-:|:-:|-------------------------------------------------------------------|:-:|
| Higgs decay | CC                        | CC                       | ...both c quarks from the decay of a Higgs $H\to c\bar{c}$       | x                |
| Higgs decay             | BB                        | BB                       | ...both b quarks from the decay of a Higgs $H\to b\bar{b}$       | x                |
| Background   | LL                        | bb                       | ...no full Higgs decay but two or more $b$ quarks                  | (x)              |
| Background             | LL                          | bl                       | ...no full Higgs decay but one $b$ quark                           | x                |
| Background             | LL                          | cc                       | ...no full Higgs decay and no $b$ quarks but two or more $c$ quarks |                  |
| Background             | LL                          | cl                       | ...no full Higgs decay and no $b$ quarks but one $c$ quark           | x                |
| Background             | LL                          | ll                       | ...neither $b$ nor $c$ quarks                                         | x                |

### Physics process indices meaning
The PAIReD jets are assigned a truth MC flag called `MC_physics_process`, an integer which indicates the simulated physics process of the event. The dictionary for the index to event type conversion is the following.
| MC_physics_process (integer) | Physics process | Short tag |
|--:|:--|:--|
| 1 | $W^+ H(H\to c\bar{c})\quad$    with  $(W^+\to l^+\nu)$ | `WpHcc` |
| 2 | $W^- H(H\to c\bar{c})\quad$    with  $(W^-\to l^-\bar{\nu})$ | `WmHcc` |
| 3 | $ZH(H\to c\bar{c})\quad$    with  $(Z\to l^+l^-)$ | `ZHcc` |
| 11 | $W^+ H(H\to b\bar{b})\quad$    with  $(W^+\to l^+\nu)$ | `WpHbb` |
| 12 | $W^- H(H\to b\bar{b})\quad$    with  $(W^-\to l^-\bar{\nu})$ | `WmHbb` |
| 13 | $ZH(H\to b\bar{b})\quad$    with  $(Z\to l^+l^-)$ | `ZHbb` |
| 21 | QCD ($HT=200-400$)  | `QCD` |
| 22 | $W$ + Jets  | `W` |
| 23 | DY + Jets  | `DY` |
| 24 | $t\bar{t}\to l\nu qq$  | `TT` |


## The `tree` tree
Note that the particles `part` are sorted accoding to their transverse momentum $p_T$ while the secondary vertices `sv` are sorted according to their decay length significance.
| Object | Type | Description |
|:--|:--|:--|
| `event` | uint64_t | Event index of the MC event from the original simulation file |
| `genweight` | float | Generator weight of the event |
| `run` | uint8_t | Run/i from the original simulation file |
| `Pileup_nPU` | uint8_t | the number of pileup interactions that have been added to the event in the current bunch crossing |
| `pv_n` | uint8_t | total number of reconstructed primary vertices |
| `pv_ngood`| uint8_t | number of good reconstructed primary vertices. selection: `!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2` |
| `fgrfcc` | float | `Rho_fixedGridRhoFastjetCentralCalo`: rho from calo towers with `abs(eta) < 2.5`, used e.g. egamma PFCluster isolation |
| `fgrfccpu` | float | `Rho_fixedGridRhoFastjetCentralChargedPileUp`: rho from charged PF Candidates for central region, used e.g. for JECs |
| `jet1_eta` | float | Pseudo-rapidity of the first jet |
| `jet1_phi` | float | Phi angle of the first jet |
| `jet1_pt` | float | Transverse momentum of the first jet |
| `jet1_nparticles` | float | Number of particles in the first jet |
| `jet1_energy` | float | Energy of the first jet |
| `jet1_index` | float | Jet ID of the first jet |
| `jet2_eta` | float | Pseudorapidity of the second jet |
| `jet2_phi` | float | Phi angle of the second jet |
| `jet2_pt` | float | Transverse momentum of the second jet |
| `jet2_nparticles` | float | Number of particles in the second jet |
| `jet2_energy` | float | Energy of the second jet |
| `jet2_index` | float | Jet ID of the second jet |
| `dijet_eta` | float | Pseudo-rapidity of the dijet (sum of all particles in the constructed PAIReD jet ellipse) |
| `dijet_phi` | float | Phi angle of the dijet (sum of all particles in the constructed PAIReD jet ellipse) |
| `dijet_pt` | float | Transverse momentum of the dijet (sum of all particles in the constructed PAIReD jet ellipse) |
| `dijet_mass` | float | Invariant mass of the dijet (sum of all particles in the constructed PAIReD jet ellipse) |
| `dijet_nparticles` | float | Number of particles in the dijet (inside the constructed PAIReD jet ellipse) |
| `dijet_index` | float | Index of the dijet (with respect to the `event`, therefore each `event` has a PAIReD jet `1`) |
| `npart` | int32_t | Number of particles in the dijet (essentially, the same as `dijet_nparticles`). Maybe just remove the other one! |
| `part_eta` | float[] | Array of the pseudo-rapidities of all particles in the PAIReD jet |
| `part_phi` | float[] | Array of the phi angles of all particles in the PAIReD jet |
| `part_pt` | float[] | Array of the transverse momenta of all particles in the PAIReD jet |
| `part_charge` | float[] | Array of the charges of all particles in the PAIReD jet |
| `part_pid` | float[] | Array of the PDG particle IDs of all particles in the PAIReD jet |
| `part_d0val` | float[] | Array of the transverse impact parameter values of all particles in the PAIReD jet |
| `part_d0err` | float[] | Array of the d0 uncertainties of all particles in the PAIReD jet |
| `part_dzval` | float[] | Array of the longitudinal impact parameter values (dz with sign wrt first PV in cm) of all particles in the PAIReD jet |
| `part_dzerr` | float[] | Array of the dz uncertainties in cm of all particles in the PAIReD jet |
| `part_mass` | float[] | Array of the masses of all particles in the PAIReD jet |
| `part_px` | float[] | Array of the x components of the momenta of all particles in the PAIReD jet |
| `part_py` | float[] | Array of the y components of the momenta of all particles in the PAIReD jet |
| `part_pz` | float[] | Array of the z components of the momenta of all particles in the PAIReD jet |
| `part_energy` | float[] | Array of the energies of all particles in the PAIReD jet |
| `part_deta1` | float[] | Array of the pseudo-rapidity differences of all particles in the PAIReD jet to the first jet (multiplied by (-1) if eta of jet 1 is negative) |
| `part_deta2` | float[] | Array of the pseudo-rapidity differences of all particles in the PAIReD jet to the second jet (multiplied by (-1) if eta of jet 2 is negative) |
| `part_dphi1` | float[] | Array of the phi angle differences of all particles in the PAIReD jet to the first jet |
| `part_dphi2` | float[] | Array of the phi angle differences of all particles in the PAIReD jet to the second jet |
| `part_puppiweight` | float[] | Array of the pile-up weights of all particles in the PAIReD jet |
| `nsv` | int32_t | Number of secondary vertices |
| `sv_charge` | float[] | Array of the sums of the charge of the SV tracks of all SV in the PAIReD jet |
| `sv_chi2` | float[] | Array of the reduced chi2, i.e. chi/ndof, of all SV in the PAIReD jet |
| `sv_dlen` | float[] | Array of the decay lengths in cm of all SV in the PAIReD jet |
| `sv_dlenSig` | float[] | Array of the decay length significances of all SV in the PAIReD jet |
| `sv_dxy` | float[] | Array of the 2D decay lengths in cm of all SV in the PAIReD jet |
| `sv_dxySig` | float[] | Array of the 2D decay length significances of all SV in the PAIReD jet |
| `sv_eta` | float[] | Array of the etas of all SV in the PAIReD jet |
| `sv_mass` | float[] | Array of the masses of all SV in the PAIReD jet |
| `sv_ndof` | float[] | Array of the numbers of degrees of freedom of all SV in the PAIReD jet |
| `sv_ntracks` | float[] | Array of numbers of tracks of all SV in the PAIReD jet |
| `sv_pAngle` | float[] | Array of the pointing angles, i.e. acos(p_SV * (SV - PV)), of all SV in the PAIReD jet |
| `sv_phi` | float[] | Array of the phis of all SV in the PAIReD jet |
| `sv_pt` | float[] | Array of the pts of all SV in the PAIReD jet |
| `sv_x` | float[] | Array of the secondary vertex X positions in cm of all SV in the PAIReD jet |
| `sv_y` | float[] | Array of the secondary vertex Y positions in cm of all SV in the PAIReD jet |
| `sv_z` | float[] | Array of the secondary vertex Z positions in cm of all SV in the PAIReD jet |
| `sv_deta1` | float[] | Array of the unsigned pseudo-rapidity differences of all SV in the PAIReD jet to the first jet |
| `sv_deta2` | float[] | Array of the unsigned pseudo-rapidity differences of all SV in the PAIReD jet to the second jet |
| `sv_dphi1` | float[] | Array of the phi angle differences of all SV in the PAIReD jet to the first jet |
| `sv_dphi2` | float[] | Array of the phi angle differences of all SV in the PAIReD jet to the second jet |
| `MC_physics_process` | uint8_t | An index/integer that indicates the physics process simulated in the event; See dictionary above |
| `MC_higgs_pt` | float | For Higgs events: transverse momentum of Higgs; For non-Higgs events: tranverse momentum of constructed diparton (two highest-$p_T$ partons combined) |
| `MC_higgs_eta` | float | For Higgs events: pseudo-rapidity of Higgs; For non-Higgs events: pseudo-rapidity of constructed diparton (two highest-$p_T$ partons combined) |
| `MC_higgs_phi` | float | For Higgs events: phi angle of Higgs; For non-Higgs events: phi angle of constructed diparton (two highest-$p_T$ partons combined) |
| `MC_higgs_mass` | float | For Higgs events: mass of Higgs; For non-Higgs events: mass of constructed diparton (two highest-$p_T$ partons combined) |
| `MC_higgs_flav` | uint8_t | For Higgs events: represents the flavor of the daughter particles of the Higgs boson (e.g., 1 for d, 4 for c, or 5 for b); For non-Higgs events: 0 |
| `MC_vector_flav` | uint8_t | Represents the flavor of the daughter particles of the vector boson if one is present (e.g., 11 for e or nu_e, 13 for mu or nu_mu, or 15 for tau or nu_tau); If no vector boson found in GenPart: 0 |
| `MC_lepton_channel` | int8_t | Number of charged leptons the vector boson decays into (0, 1 or 2); If no vector boson found in GenPart: -1 (Note: due to unwanted type uint, some earlier produced files may contain garbage values instead of -1) |
| `MC_gendijet_pt` | float | Transverse momentum of the di-jet (sum of the two GenJets that are matched to the seed RecoJets used for the PAIReD jet); If no matched GenJets: 0 |
| `MC_gendijet_eta` | float | Pseudo-rapidity of the di-jet (sum of the two GenJets that are matched to the seed RecoJets used for the PAIReD jet); If no matched GenJets: 0 |
| `MC_gendijet_phi` | float | Phi angle of the di-jet (sum of the two GenJets that are matched to the seed RecoJets used for the PAIReD jet); If no matched GenJets: 0 |
| `MC_gendijet_mass` | float | Mass of the di-jet (sum of the two GenJets that are matched to the seed RecoJets used for the PAIReD jet); If no matched GenJets: 0 |
| `MC_genjet1_flav` | int8_t | Flavor of the GenJet that was matched to the first seed RecoJet used for the PAIReD jet; If no matched GenJet: 0 |
| `MC_genjet2_flav` | int8_t | Flavor of the GenJet that was matched to the second seed RecoJet used for the PAIReD jet; If no matched GenJet: 0 |
| `MC_genjet1_matched` | int8_t | Boolean saying whether the first RecoJet used has a matched GenJet |
| `MC_genjet2_matched` | int8_t | Boolean saying whether the second RecoJet used has a matched GenJet |
| `MC_drqq` | float | For Higgs events: Delta R between the two decay products of the Higgs; For non-Higgs events: 0 |
| `MC_n_c` | uint8_t | Number of c jets in the event |
| `label_BB` | int8_t | True label for the PAIReD jet regarding BB |
| `label_CC` | int8_t | True label for the PAIReD jet regarding CC |
| `label_LL` | int8_t | True label for the PAIReD jet regarding LL |
| `label_ll` | int8_t | True label for the PAIReD jet regarding ll |
| `label_cl` | int8_t | True label for the PAIReD jet regarding cl |
| `label_cc` | int8_t | True label for the PAIReD jet regarding cc |
| `label_bl` | int8_t | True label for the PAIReD jet regarding bl |
| `label_bb` | int8_t | True label for the PAIReD jet regarding bb |
