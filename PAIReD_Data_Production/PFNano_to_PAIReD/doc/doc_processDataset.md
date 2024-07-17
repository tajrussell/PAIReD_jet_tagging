# Documentation for the `processDatasetToPAIReD` file

Below you can find the documentation for the file `src/processDatasetToPAIReD.py`.

This module contains the function `processDataset`.

The module requires the following python packages:
* `dbs`
* `argparse`
* `multiprocessing`
* `functools`
* `os`


## Function documentation

### `processDataset` <a name="processEvents">  </a>

**`processDataset(dataset, outputpath, storageprefix = "root://grid-cms-xrootd.physik.rwth-aachen.de/", batchsize = "100 MB", first_file = 0, N_files = -1, N_update = 10, N_events_to_process = -1, N_threads = 1, physics_process = 0, PAIReD_geometry = "Ellipse")`**

Process an entire data set to PAIReD jets.

<dl>
  <dt>Parameters:</dt>
  <dd><dl>
    <dt>dataset : str</dt>
    <dd>
      <p>Name of the dataset as shown in CMSDAS. Note that it has to be a NanoAOD set including PFCands.</p>
    </dd>
    <dt>outputpath : str</dt>
    <dd>
      <p>Path to the directory where the produced PAIReD jet files should be stored. If the path does not exist yet, it will be created.</p>
    </dd>
    <dt>storageprefix : str, optional (default is root://grid-cms-xrootd.physik.rwth-aachen.de/)</dt>
    <dd>
      <p>Storage prefix for the given data set. You can check on CMSDAS under <i>Sites</i> where the data set is stored.</p>
    </dd>
    <dt>batchsize : str, optional (default is 100MB)</dt>
    <dd>
      <p>Batchsize of the data amount to be processed by each worker simultaneously. Can be given with a unit, e.g. MB, or as an integer in which case it represents the number of events.</p>
    </dd>
    <dt>first_file : int, optional (default is 0)</dt>
    <dd>
      <p>First file of the list of files of the data set to be processed. Start at <code>0</code> if you want all files to be processed.</p>
    </dd>
    <dt>N_files : int, optional (default is -1)</dt>
    <dd>
      <p>Number of files to be processed. <code>-1</code> means up to the last file of the list.</p>
    </dd>
    <dt>N_update : int, optional (default is 10)</dt>
    <dd>
      <p>Number of batches to be processed by each worker before it prints an update on the progress.</p>
    </dd>
    <dt>N_event_to_process : int, optional (default is -1)</dt>
    <dd>
      <p>Number of events to be processed in each file. <code>-1</code> means all events are processed.</p>
    </dd>
    <dt>N_threads : int, optional (default is 1)</dt>
    <dd>
      <p>Number of threads/workers used simultaneously for processing.</p>
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
  <dt>Returns: nothing</dt>
</dl>

#### Examples
Here is an example on how the output might look like if you run the following command in the `src` directory:
```
>>> python processDatasetToPAIReD.py \
>>> /WminusH_Hto2C_WtoLNu_M-125_TuneCP5_13p6TeV_powheg-pythia8/jaschulz-Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2_BTV_Run3_2022_Comm_PFNANOAODv12-490aef546cce8088aae92a89578a4232/USER \
>>> ../data/test -u 5 -n 2 --start-file 0 -t 1 --physicsprocess 23 -g AK4

#################################################################
  Start processDataset()
#################################################################

** Process the following data set:
    /WminusH_Hto2C_WtoLNu_M-125_TuneCP5_13p6TeV_powheg-pythia8/jaschulz-Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2_BTV_Run3_2022_Comm_PFNANOAODv12-490aef546cce8088aae92a89578a4232/USER
** Storage prefix: root://grid-cms-xrootd.physik.rwth-aachen.de/
** Number of files in the data set: 98
** Number of files to be processed: 2
** Start at file: 0

Settings for processFileToPAIReD:
** Batch size: 100 MB
** N_update: 5
** Number of events to process: -1
** PAIReD geometry: AK4
** Create output directory if it does not exist yet: ../data/test
** Use only one thread

** Start processing the data set



*****************************************************************
  Start makeNtuplesPAIReDjointMC()
*****************************************************************

** Open input file: root://grid-cms-xrootd.physik.rwth-aachen.de//store/user/jaschulz/data/PFNano_Run3/mc_summer22EE_PFNANOAODv12/WminusH_Hto2C_WtoLNu_M-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2_BTV_Run3_2022_Comm_PFNANOAODv12/240326_092356/0000/PFNANOAODSIM_allPF_allBTV_1-1.root
** Create output file: ../data/test/PAIReD_1-1.root
** Create tree "tree" entitled "PAIReD jets for training" in output file
** Total number of events: 33705
** Process all events
** Batch size: 100 MB
** Physics process: 23
** PAIReD geometry: AK4

** Start processing data
process event batch 0
process event batch 5

** Written 68727 PAIReD jets to the output file in 65.39 seconds
** Computing time per PAIReD jet: 0.951 milliseconds
** Computing time per event: 1.940 milliseconds

*****************************************************************
  End makeNtuplesPAIReDjointMC()
*****************************************************************


*****************************************************************
  Start makeNtuplesPAIReDjointMC()
*****************************************************************

** Open input file: root://grid-cms-xrootd.physik.rwth-aachen.de//store/user/jaschulz/data/PFNano_Run3/mc_summer22EE_PFNANOAODv12/WminusH_Hto2C_WtoLNu_M-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2_BTV_Run3_2022_Comm_PFNANOAODv12/240326_092356/0000/PFNANOAODSIM_allPF_allBTV_1.root
** Create output file: ../data/test/PAIReD_1.root
** Create tree "tree" entitled "PAIReD jets for training" in output file
** Total number of events: 21507
** Process all events
** Batch size: 100 MB
** Physics process: 23
** PAIReD geometry: AK4

** Start processing data
process event batch 0

** Written 43899 PAIReD jets to the output file in 40.67 seconds
** Computing time per PAIReD jet: 0.927 milliseconds
** Computing time per event: 1.891 milliseconds

*****************************************************************
  End makeNtuplesPAIReDjointMC()
*****************************************************************



** Done processing


#################################################################
  End processDataset()
#################################################################
```

#### Usage in the command line
You can see all the options for the usage in the command line by calling the `--help` option (`-h`):

```bash
>>> python processDatasetToPAIReD.py -h
usage: processDatasetToPAIReD.py [-h] [-p STORAGEPREFIX] [-b BATCHSIZE] [-u N_UPDATE] [-n N_FILES]
                                 [--start-file START_FILE] [-e N_EVENTS] [-t N_THREADS]
                                 [--physicsprocess PHYSICSPROCESS] [-g GEOMETRY]
                                 dataset outputpath

Process MC PFNano data set to the PAIReD data format for training.

positional arguments:
  dataset               Name of the data set (NanoAOD processed by PFNano, has to contain particle flow
                        candidates (PFcands))
  outputpath            Path for the ROOT output files

options:
  -h, --help            show this help message and exit
  -p STORAGEPREFIX, --storageprefix STORAGEPREFIX
                        Prefix of the storgae site the data set is stored on (default for T2_DE_RWTH site)
  -b BATCHSIZE, --batchsize BATCHSIZE
                        Size of the batches of events that are handled simultaneously (default is '100 MB').
                        Can be an integer (=number of events) or a string (specifying the data size, e.g.
                        '100 MB')
  -u N_UPDATE, --N-update N_UPDATE
                        Number of batches processed before printing an update in the terminal: e.g.,
                        'process event batch 10 of 15' (default is 10)
  -n N_FILES, --N-files N_FILES
                        Number of files to be processed. -1 means all (default is -1)
  --start-file START_FILE
                        Index of first file to be processed (default is 0)
  -e N_EVENTS, --N-events N_EVENTS
                        Number of events to be processed per file. -1 means all (default is -1)
  -t N_THREADS, --N-threads N_THREADS
                        Number of threads (default is 1)
  --physicsprocess PHYSICSPROCESS
                        Integer indicating the physics process. Default is 0.
  -g GEOMETRY, --geometry GEOMETRY
                        String indicating the geometry of the PAIReD jet. Default is Ellipse.
```
