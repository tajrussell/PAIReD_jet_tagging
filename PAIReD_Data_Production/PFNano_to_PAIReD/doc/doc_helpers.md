# Documentation for the `helpers` file

Below you can find the documentation for the file `src/tools/helpers.py`.

This module contains the following functions which are used for the data conversion done in `PFNano_to_PAIReD`:
* [`isoLeptonCut()`](#isoLeptonCut)
* [`Phi_mpi_pi()`](#Phi_mpi_pi)
* [`deltaPhi()`](#deltaPhi)
* [`deltaR()`](#deltaR)
* [`shiftPhi()`](#shiftPhi)
* [`getJetCombinationIndices()`](#getJetCombinationIndices): not used in the current version
* [`getAllJetCombinationIndices()`](#getAllJetCombinationIndices): not used in the current version

The module requires the following python packages:
* `numpy`
* `vector`
* `awkward`
* `itertools`


## Function documentation
<hr class="solid">

### `processEvents` <a name="processEvents">  </a>

**`processEvents(Events)`**

Process the events to produce PAIReD jets with all relevant quantities.

<dl>
  <dt>Parameters:</dt>
  <dd><dl>
    <dt>Events : event tree object</dt>
    <dd>
      <p>Object containing all information of multiple simulated events.</p>
    </dd>
    </dl>
  </dd>
  <dt>Returns:</dt>
  <dd><dl>
    <dt>DataPAIReD : dict</dt>
    <dd>
      <p>A dictionary with the content summarized in <a href="/notes/PAIReD_ROOT_file_content.md">PAIReD_ROOT_file_content.md</a></p>
    </dd>
    </dl>
  </dd>
</dl>

<hr class="solid">

### `isInEllipse` <a name="isInEllipse">  </a>

**`isInEllipse(eta_j1, eta_j2, phi_j1, phi_j2, eta_p, phi_p, semimajoradd = 1.)`**

Determines which of the input particles are inside the ellipse spanned by the given jet pair. The output array has the shape $(N_e, N_{jp}, N_p)$, where $N_e$ is the number of events, $N_{jp}$ the number of jet pairs (length of the jet input arrays) and $N_p$ the number of particles (length of the particle input arrays). Note that $N_{jp}$ and $N_p$ are event dependend which is why the lengths of these axes vary for the different entries along the event axis.

<dl>
  <dt>Parameters:</dt>
  <dd><dl>
    <dt>eta_j1 : float or array of floats</dt>
    <dd>
      <p>Pseudorapidity of the first jet, or jets if array is given.</p>
    </dd>
    <dt>eta_j2 : float or array of floats</dt>
    <dd>
      <p>Pseudorapidity of the second jet, or jets if array is given.</p>
    </dd>
    <dt>phi_j1 : float or array of floats</dt>
    <dd>
      <p>Phi angle of the first jet, or jets if array is given.</p>
    </dd>
    <dt>phi_j2 : float or array of floats</dt>
    <dd>
      <p>Phi angle of the second jet, or jets if array is given.</p>
    </dd>
    <dt>eta_p : float or array of floats</dt>
    <dd>
      <p>Pseudorapidity of the particle, or particles if array is given.</p>
    </dd>
    <dt>phi_p : float or array of floats</dt>
    <dd>
      <p>Phi angle of the particle, or particles if array is given.</p>
    </dd>
    <dt>semimajoradd : float, optional (default is 1.0)</dt>
    <dd>
      <p>Distance between the individual jet centers and the tip of the ellipse. The semimajor axis is therefore:
        $$semimajor = \frac{jetdistance}{2} + semimajoradd$$</p>
    </dd>
    </dl>
  </dd>
  <dt>Returns:</dt>
  <dd><dl>
    <dt>isInEllipse : bool or array of bools</dt>
    <dd>
      <p>Boolean holding the information whether the particle is inside the ellipse or not. The with index pair $(i,j)$ holds the information whether particle $j$ is in jet pair $i$.</p>
    </dd>
    </dl>
  </dd>
</dl>

#### Examples

```python
eta_j1 = ak.Array([[0, 1, 2],
                   [],
                   [-2]])
phi_j1 = ak.Array([[1, 2,-3],
                   [],
                   [-1]])

eta_j2 = ak.Array([[1, 2,-6],
                   [],
                   [1]])
phi_j2 = ak.Array([[2,-2, 0],
                   [],
                   [0]])

eta_p = ak.Array([[2, 0,-5, 8],
                  [2,1,3],
                  [7,2,-2]])
phi_p = ak.Array([[2, 2,-1, 3],
                  [-2,-1,3],
                  [1,1,-3]])

isInEllipse(eta_j1, eta_j2, phi_j1, phi_j2, eta_p, phi_p)
```
```
Out[]:   [[[True, True, False, False], [True, ...], [False, False, True, False]],
         [],
         [[False, False, False]]]
        ------------------------------------------------------------------------
        type: 3 * var * var * bool
```

<hr class="solid">

### `getMCInfo` <a name="getMCInfo">  </a>

**`getMCInfo(Events, isInGenPart, Jet_genJetIdx, Jetcut)`**

Get the relevant true MC information of the event.

<dl>
  <dt>Parameters:</dt>
  <dd><dl>
    <dt>Event : event tree object</dt>
    <dd>
      <p>Object containing all information of a simulated event.</p>
    </dd>
    <dt>isInGenPart : bool (array)</dt>
    <dd>
      <p>Boolean holding the information whether the particle is inside the PAIReD jet
        or not. (Output of <code>isInEllipse()</code>)</p>
    </dd>
    <dt>Jet_genJetIdx : array / ak-zip type object</dt>
    <dd>
      <p>Indices of the matched generated jets to the reconstructed AK4 jets.</p>
    </dd>
    <dt>Jetcut : bool (array)</dt>
    <dd>
      <p>Boolean holding the information whether the AK4 jet is accepted or not.</p>
    </dd>
    </dl>
  </dd>
  <dt>Returns:</dt>
  <dd><dl>
    <dt>higgs_pt : float (array)</dt>
    <dd>
      <p>Higgs transverse momentum (length of the array = number of jet pairs).</p>
    </dd>
    <dt>higgs_eta : float (array)</dt>
    <dd>
      <p>Higgs pseudorapidity (length of the array = number of jet pairs).</p>
    </dd>
    <dt>higgs_phi : float (array)</dt>
    <dd>
      <p>Higgs phi angle (length of the array = number of jet pairs).</p>
    </dd>
    <dt>higgs_mass : float (array)</dt>
    <dd>
      <p>Higgs mass (length of the array = number of jet pairs).</p>
    </dd>
    <dt>drqq : float (array)</dt>
    <dd>
      <p>Delta R between decay products of the Higgs (length of array = number of jet
        pairs).</p>
    </dd>
    <dt>gen_flav : float (array)</dt>
    <dd>
      <p>Generated flavor number:
      <ul>
        <li>first two digits represent the mother particle (e.g., 25 for Higgs,
          or 23 for Z)</li>
        <li>last digit represents the flavor of the daughter particle (e.g., 1 for d,
          3 for s, or 5 for b)</li>
      </ul></p>
    </dd>
    <dt>n_c : int (array)</dt>
    <dd>
      <p>Number of c jets in the event (length of the array = number of jet pairs).</p>
    </dd>
    <dt>label_bb : bool (array)</dt>
    <dd>
      <p>True label for BB PAIReD jets. See labeling scheme <a href="../notes/PAIReD_ROOT_file_content.md">here</a>.</p>
    </dd>
    <dt>label_cc : bool (array)</dt>
    <dd>
      <p>True label for CC PAIReD jets. See labeling scheme <a href="../notes/PAIReD_ROOT_file_content.md">here</a>.</p>
    </dd>
    <dt>label_ll : bool (array)</dt>
    <dd>
      <p>True label for LL PAIReD jets. See labeling scheme <a href="../notes/PAIReD_ROOT_file_content.md">here</a>.</p>
    </dd>
    <dt>label_bl : bool (array)</dt>
    <dd>
      <p>True label for BL PAIReD jets. See labeling scheme <a href="../notes/PAIReD_ROOT_file_content.md">here</a>.</p>
    </dd>
    <dt>label_cl : bool (array)</dt>
    <dd>
      <p>True label for CL PAIReD jets. See labeling scheme <a href="../notes/PAIReD_ROOT_file_content.md">here</a>.</p>
    </dd>
    </dl>
  </dd>
</dl>

<hr class="solid">

### `isoLeptonCut` <a name="isoLeptonCut"> </a>

**`isoLeptonCut(Events)`**

Defines a selection cut for all jets in events according to the following condition: It finds all isolated leptons (muons and electrons) in each event, and checks which jets coincide with those isoLeptons. Since those jets likely originated from these isoLeptons, we discard them later on.

<dl>
  <dt>Parameters:</dt>
  <dd><dl>
    <dt>Events : event tree object</dt>
    <dd>
      <p>Object containing all information of multiple simulated events.</p>
    </dd>
    </dl>
  </dd>
  <dt>Returns:</dt>
  <dd><dl>
    <dt>cut : array of bools</dt>
    <dd>
      <p>Boolean for each jet saying whether it came from an isoLepton.</p>
    </dd>
    </dl>
  </dd>
</dl>

<hr class="solid">

### `Phi_mpi_pi` <a name="Phi_mpi_pi"> </a>

**`Phi_mpi_pi(angle)`**

Outputs the angle in the interval $[-\pi, \pi)$.

<dl>
  <dt>Parameters:</dt>
  <dd><dl>
    <dt>angle : float or array of floats</dt>
    <dd>
      <p>Angle [in radians].</p>
    </dd>
    </dl>
  </dd>
  <dt>Returns:</dt>
  <dd><dl>
    <dt>angle : float or array of floats</dt>
    <dd>
      <p>Angle in $[-\pi, \pi)$. Output has the same shape as the input.</p>
    </dd>
    </dl>
  </dd>
</dl>

<hr class="solid">

### `deltaPhi` <a name="deltaPhi"> </a>

**`deltaPhi(phi1, phi2)`**

Computes the angular difference between two angles in radians. The result lies in the interval $[-\pi, \pi)$.

<dl>
  <dt>Parameters:</dt>
  <dd><dl>
    <dt>phi1 : float or array of floats</dt>
    <dd>
      <p>First angle [in radians].</p>
    </dd>
    <dt>phi2 : float or array of floats</dt>
    <dd>
      <p>Second angle [in radians].</p>
    </dd>
    </dl>
  </dd>
  <dt>Returns:</dt>
  <dd><dl>
    <dt>deltaPhi : float or array of floats</dt>
    <dd>
      <p>Angular difference of the two angles in $[-\pi, \pi)$.</p>
    </dd>
    </dl>
  </dd>
</dl>

<hr class="solid">

### `deltaR` <a name="deltaR"> </a>

**`deltaR(eta1, eta2, phi1, phi2)`**

Computes the distance between two jets in the $(\eta, \phi)$ space.

<dl>
  <dt>Parameters:</dt>
  <dd><dl>
    <dt>eta1 : float or array of floats</dt>
    <dd>
      <p>Pseudorapidity of the first jet.</p>
    </dd>
    <dt>eta2 : float or array of floats</dt>
    <dd>
      <p>Pseudorapidity of the second jet.</p>
    </dd>
    <dt>phi1 : float or array of floats</dt>
    <dd>
      <p>Phi angle of the first jet.</p>
    </dd>
    <dt>phi2 : float or array of floats</dt>
    <dd>
      <p>Phi angle of the second jet.</p>
    </dd>
    </dl>
  </dd>
  <dt>Returns:</dt>
  <dd><dl>
    <dt>deltaR : float or array of floats</dt>
    <dd>
      <p>Distance between two jets in the $(\eta,\phi)$ space:
        $$\Delta R = \sqrt{(\eta_1 - \eta_2)^2 + (\phi_1 - \phi_2)^2}$$</p>
    </dd>
    </dl>
  </dd>
</dl>

<hr class="solid">

### `shiftPhi` <a name="shiftPhi"> </a>

**`shiftPhi(phi, shift, lim)`**

Shifts the input angle based on some conditions:
```python
if shift == 0:
    return phi
elif shift > 0:
    if phi < lim:
        return phi+shift
else:
    if phi > lim:
        return phi+shift
return phi
```
This function is necessary to treat the discontinuity in the $\phi$ dimension at $\pi \widehat{=} -\pi$.

**Note:** if the inputs are arrays, they have to be broadcastable.

<dl>
  <dt>Parameters:</dt>
  <dd><dl>
    <dt>phi : float or array of floats</dt>
    <dd>
      <p>Angle in radians.</p>
    </dd>
    <dt>shift : float or array of floats</dt>
    <dd>
      <p>Potential shift dPhi of the angle in radians.</p>
    </dd>
    </dl>
    <dt>lim : float or array of floats</dt>
    <dd>
      <p>Limit for the condition of the shift.</p>
    </dd>
    </dl>
  </dd>
  <dt>Returns:</dt>
  <dd><dl>
    <dt>phi : float or array of floats</dt>
    <dd>
      <p>Shifted angle.</p>
    </dd>
    </dl>
  </dd>
</dl>

<hr class="solid">

### `getJetCombinationIndices` <a name="getJetCombinationIndices"> </a>

**`getJetCombinationIndices(N_Jet)`**

Returns a list of index pairs containing all possible jet combinations. E.g. for
$N=3$ jets, it will return `[[0, 1], [0, 2], [1, 2]]`. For $N=0$ or $1$, it returns an
empty list `[]`.

<dl>
  <dt>Parameters:</dt>
  <dd><dl>
    <dt>N_Jet : int</dt>
    <dd>
      <p>Number of jets in the event.</p>
    </dd>
    </dl>
  </dd>
  <dt>Returns:</dt>
  <dd><dl>
    <dt>IndicesList : list</dt>
    <dd>
      <p>List of possible index pairs.</p>
    </dd>
    </dl>
  </dd>
</dl>

<hr class="solid">

### `getAllJetCombinationIndices` <a name="getAllJetCombinationIndices"> </a>

**`getAllJetCombinationIndices(Nmax_Jet)`**

Returns a list of lists of index pairs containing all possible jet combinations
for all possible numbers of jets up to $N_{\max, Jet}$.
E.g. for `Nmax=3` jets, it will return
```
    [ [],                          # 0 jets
      [],                          # 1 jets
      [[0, 1]],                    # 2 jets
      [[0, 1], [0, 2], [1, 2]]     # 3 jets = Nmax jets
    ]
```
The idea is to create a look-up table at the beginning, instead of running the
`getJetCombinationIndices` function for each event individually.

<dl>
  <dt>Parameters:</dt>
  <dd><dl>
    <dt>Nmax_Jet : int</dt>
    <dd>
      <p>Maximum number of jets in one event accross all events in the tree.</p>
    </dd>
    </dl>
  </dd>
  <dt>Returns:</dt>
  <dd><dl>
    <dt>IndicesList : list</dt>
    <dd>
      <p>List of lists of possible index pairs.</p>
    </dd>
    </dl>
  </dd>
</dl>
