# Documentation for the `getMCInfo` file

Below you can find the documentation for the file `src/tools/getMCInfo.py`.

This module contains the function `getMCInfo`. A detailed example of the functioning, step by step, can be found in the jupyter notebook [here](../src/notebooks/makeNtuplesPAIReDjointMC.ipynb).

The module requires the following python packages:
* `numpy`
* `vector`
* `awkward`
* `tools.helpers`


## Function documentation

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
    <dt>label_BB : bool (array)</dt>
    <dd>
      <p>True label for BB PAIReD jets. See labeling scheme <a href="../notes/PAIReD_ROOT_file_content.md">here</a>.</p>
    </dd>
    <dt>label_CC : bool (array)</dt>
    <dd>
      <p>True label for CC PAIReD jets. See labeling scheme <a href="../notes/PAIReD_ROOT_file_content.md">here</a>.</p>
    </dd>
    <dt>label_LL : bool (array)</dt>
    <dd>
      <p>True label for LL PAIReD jets. See labeling scheme <a href="../notes/PAIReD_ROOT_file_content.md">here</a>.</p>
    </dd>
    <dt>label_ll : bool (array)</dt>
    <dd>
      <p>True label for ll PAIReD jets. See labeling scheme <a href="../notes/PAIReD_ROOT_file_content.md">here</a>.</p>
    </dd>
    <dt>label_cl : bool (array)</dt>
    <dd>
      <p>True label for cl PAIReD jets. See labeling scheme <a href="../notes/PAIReD_ROOT_file_content.md">here</a>.</p>
    </dd>
    <dt>label_cc : bool (array)</dt>
    <dd>
      <p>True label for cc PAIReD jets. See labeling scheme <a href="../notes/PAIReD_ROOT_file_content.md">here</a>.</p>
    </dd>
    <dt>label_bl : bool (array)</dt>
    <dd>
      <p>True label for bl PAIReD jets. See labeling scheme <a href="../notes/PAIReD_ROOT_file_content.md">here</a>.</p>
    </dd>
    <dt>label_bb : bool (array)</dt>
    <dd>
      <p>True label for bb PAIReD jets. See labeling scheme <a href="../notes/PAIReD_ROOT_file_content.md">here</a>.</p>
    </dd>
    </dl>
  </dd>
</dl>


