# Documentation for the `processEvents` file

Below you can find the documentation for the file `src/tools/processEvents.py`.

This module contains the function `processEvents`. A detailed example of the functioning, step by step, can be found in the jupyter notebook [here](../src/notebooks/makeNtuplesPAIReDjointMC.ipynb).

The module requires the following python packages:
* `numpy`
* `vector`
* `awkward`
* `tools.branchnames`
* `tools.helpers`
* `tools.getMCInfo`
* `tools.PAIReD_geometries`


## Function documentation

### `processEvents` <a name="processEvents">  </a>

**`processEvents(Events, physics_process=0, PAIReD_geometry="Ellipse")`**

Process the events to produce PAIReD jets with all relevant quantities.

<dl>
  <dt>Parameters:</dt>
  <dd><dl>
    <dt>Events : event tree object</dt>
    <dd>
      <p>Object containing all information of multiple simulated events.</p>
    </dd>
    </dl>
    <dt>physics_process : int, optional (default is 0)</dt>
    <dd>
      <p>Integer indicating the physics process in the file. <code>0</code> means unknown/not specified.</p>
    </dd>
    <dt>PAIReD_geometry : str, optional (default is "Ellipse")</dt>
    <dd>
      <p>String indicating the geometry of the PAIReD jet. Implemented choices: "Ellipse", "AK4" and "Dumbbell".</p>
    </dd>
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

