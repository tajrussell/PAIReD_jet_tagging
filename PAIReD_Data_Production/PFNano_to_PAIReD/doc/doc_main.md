# Documentation for the `processFileToPAIReD` file

Below you can find the documentation for the file `src/processFileToPAIReD.py`.

This module contains the following function which is used for the MC data conversion done in `PFNano_to_PAIReD`:
* [`makeNtuplesPAIReDjointMC()`](#makeNtuplesPAIReDjointMC): main function of the repository `PFNano_to_PAIReD` 

A detailed example of the functioning, step by step, can be found in the jupyter notebook [here](../src/notebooks/makeNtuplesPAIReDjointMC.ipynb).

The module requires the following python packages:
* `numpy`
* `awkward`
* `uproot`
* `vector`
* `time`


## Function documentation

### `makeNtuplesPAIReDjointMC` [main function] <a name="makeNtuplesPAIReDjointMC">  </a>

**`makeNtuplesPAIReDjointMC(inputFilePath, outputFilePath, batchsize = "100 MB", N_update = 10, N_events_to_process = -1, physics_process = 0, PAIReD_geometry = "Ellipse")`**

Performs the data conversion from [PFNano](https://github.com/cms-jet/PFNano) output files (extended NanoAODs) to a data format that the PAIReD tagger network can read in (ROOT file). Remark: this function is designed to deal with CMS MC simulations and **not** real data.

<dl>
  <dt>Parameters:</dt>
  <dd><dl>
    <dt>inputFilePath : str</dt>
    <dd>
      <p>Path to the input file (NanoAOD processed by PFNano, has to contain <b>all</b> particle flow candidates (PFCands)).</p>
    </dd>
    <dt>outputFilePath : str</dt>
    <dd>
      <p>Path and name of the output ROOT file.</p>
    </dd>
    <dt>batchsize : int or str, optional (default is "100 MB")</dt>
    <dd>
      <p>Size of the batches of events that are handled simultaneously. If <code>int</code> = number of events, if <code>str</code> = data size, e.g., <code>"100 MB"</code></p>
    </dd>
    <dt>N_update : int, optional (default is 10)</dt>
    <dd>
      <p>Number of batches to process before printing an update of the form:
        <code>process event batch 10 of 15</code></p>
    </dd>
    <dt>N_events_to_process : int, optional (default is -1)</dt>
    <dd>
      <p>Number of events to process. <code>-1</code> means all events are processed.</p>
    </dd>
    <dt>physics_process : int, optional (default is 0)</dt>
    <dd>
      <p>Integer indicating the physics process in the file. <code>0</code> means unknown/not specified.</p>
    </dd>
    <dt>PAIReD_geometry : str, optional (default is "Ellipse")</dt>
    <dd>
      <p>String indicating the geometry of the PAIReD jet. Implemented choices: "Ellipse", "AK4" and "Dumbbell".</p>
    </dd>
    </dl>
  </dd>
  <dt>Returns:</dt>
  <dd><dl>
    <dt></dt>
    <dd>
      <p>It creates an output ROOT file with the transformed data for the PAIReD tagger.</p>
    </dd>
    </dl>
  </dd>
</dl>

#### Examples
Here is an example that can be run in the `src` directory from the command line:
```
>>> python processFileToPAIReD.py ../data/example_PFNano_mcRun3_EE_allPF_noBTV.root ../data/example_PAIReD_mcRun3_EE.root

*****************************************************************
  Start makeNtuplesPAIReDjointMC()
*****************************************************************

** Open input file: ../data/example_PFNano_mcRun3_EE_allPF_noBTV.root
** Create output file: ../data/example_PAIReD_mcRun3_EE.root
** Create tree "tree" entitled "PAIReD jets for training" in output file
** Total number of events: 568
** Process all events
** Batch size: 100 MB
** Physics process: 0
** PAIReD geometry: Ellipse

** Start processing data
process event batch 0

** Written 1375 PAIReD jets to the output file in 4.23 seconds
** Computing time per PAIReD jet: 3.073 milliseconds
** Computing time per event: 7.440 milliseconds

*****************************************************************
  End makeNtuplesPAIReDjointMC()
*****************************************************************
```

#### Usage in the command line
You can see all the options for the usage in the command line by calling the `--help` option (`-h`):

```bash
>>> python makeNtuplesPAIReDjointMC.py -h
usage: processFileToPAIReD.py [-h] [--batchsize BATCHSIZE] [--N-update N_UPDATE] [-n NEVENTS] [-p PHYSICSPROCESS] [-g GEOMETRY] inputFilePath outputFilePath

Process MC PFNano files to the PAIReD data format for training.

positional arguments:
  inputFilePath         Path to the ROOT input file (NanoAOD processed by PFNano, has to contain particle flow candidates (PFcands))
  outputFilePath        Path and name of the ROOT output file

options:
  -h, --help            show this help message and exit
  --batchsize BATCHSIZE
                        Size of the batches of events that are handled simultaneously (default is '100 MB'). Can be an integer (=number of events) or a string (specifying the data size, e.g. '100 MB')
  --N-update N_UPDATE   Number of batches processed before printing an update in the terminal: e.g., 'process event batch 10 of 15' (default is 10)
  -n NEVENTS, --nevents NEVENTS
                        Number of events to be processed
  -p PHYSICSPROCESS, --physicsprocess PHYSICSPROCESS
                        Integer indicating the physics process. Default is 0.
  -g GEOMETRY, --geometry GEOMETRY
                        String indicating the geometry of the PAIReD jet. Default is Ellipse.
```
