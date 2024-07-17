# Production of the data sets
In order to create a tagger for PAIReD jets, data sets of such PAIReD jets are required first of all. As this type of jet is new in CMS, PAIReD jets are not yet available in standardized file formats such as MiniAOD or NanoAOD. Accordingly, we have to generate the data sets ourselves. We do this in two steps:
1. We create a NanoAOD data set from a data set in MiniAOD format.
2. From these NanoAOD files we create PAIReD files with ntuples of PAIReD jets.
Both steps are described below. 

**Note:** Technically, it is unnecessary to do this in two steps. A single step from MiniAOD directly to PAIReD ntuples is possible, but was not implemented here due to an unfavorable choice of approach at the beginning of the project. This may be implemented in the future.

## 1. From MiniAOD to PFNanoAOD
We use the BTV Nano production workflow to produce NanoAOD data sets. This is necessary because the PF candidates are required as tagger input and these are not available in normal NanoAODs. You can find the corresponding [repo here](https://github.com/cms-btv-pog/btvnano-prod/tree/NanoAODv12_22Sep2023). Follow the instructions given there to set up CMSSW and crab on `lxplus`.

**Important:** Before sending jobs to produce data sets, a setting must be changed in the CMSSW code. The [following line](https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/custom_btv_cff.py#L658) in the `custom_btv_cff.py` file
```python
btvNano_addallPF_switch = cms.untracked.bool(False)
```
must be changed to  
```python
btvNano_addallPF_switch = cms.untracked.bool(True)
```
The reason for this is that in PAIReD jets there are potentially more PF candidates in the PAIReD jet than just those of the AK4 jet. Therefore, all PF candidates must be stored.
 
After the setup, a crab job for the production of NanoAOD data sets can be specified by a `yml` card following the instructions. It is important that the correct Python config file is specified in the `yml` card, an overview of this can be found under [How to Update](https://github.com/cms-btv-pog/btvnano-prod/tree/NanoAODv12_22Sep2023?tab=readme-ov-file#how-to-update) with the commands for generation.

The job is created using 
```bash
python3 crabby.py -c path/to/your_card.yml --make
```
and submitted using 
```bash
python3 crabby.py -c path/to/your_card.yml --submit
```

**Note:** I personally still used `lxplus7` which ended working recently. I have not tested the BTV-Nano workflow on the newer linux versions yet.

## 2. From PFNanoAOD to PAIReD
Once you have created the data sets in the NanoAOD file, the next step is to produce PAIReD jets from them. To do this, we use a specially developed framework `PFNano_to_PAIRED`. The following explains how to install `conda`, set up the required environment and download the framework. The use of the framework is then explained afterwards.

### Installing the framework
Go to the directory on your machine where you want to install the framework. Make sure that it is a reasonable place to install it to. Then, clone this repository to your current directory by running:
```bash
git clone https://github.com/JanGerritSchulz/PAIReD_jet_tagging.git
```

### Installing conda
If you don't have `Miniconda` installed already on your machine, you can install it by following these steps:

First, download `Miniconda` with this command:
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```
Next, run the installation
```bash
bash Miniconda3-latest-Linux-x86_64.sh
```
and follow the instructions. Choose a reasonable destination for the installation and make sure to choose `yes` at this point:
```
> Do you wish the installer to initialize Miniconda3
> by running conda init? [yes|no]
```
Optionally, you can disable the auto activation of the base environment at login:
```bash
conda config --set auto_activate_base false
```
Verify the installation is successful by running `conda info` and checking if the paths are pointing to your Miniconda installation.

If you cannot run the `conda` command, check if the conda path has been correctly set up in your `.bashrc`/`.zshrc` file. You may need to log out and log in again for the changes to take effect.

### Set up the conda environment
Once you have `Miniconda` installed, you can set up the `conda` environment `nano2paired` for the data processing to PAIReD files. For doing this, you need to go to main directory of the `PAIReD_jet_tagging` repository. If you just downloaded the repo, you can most likely get there with:
```bash
cd PAIReD_jet_tagging
```
When you are there, we can create and set up the environment by running:
```bash
conda env create -f environments/nano2paired.yml
```
After this is done, you can activate the environment with:
```bash
conda activate nano2paired
```

### Check if the setup works
In order to check whether everything worked fine, we will run a quick PAIReD production task. For this, we go to the `src` directory:
```bash
cd PAIReD_Data_Production/PFNano_to_PAIReD/src/
```
Make sure that your `nano2paired` environment is activated and run:
```bash
python processFileToPAIReD.py "../data/example_PFNano_mcRun3_EE_allPF_noBTV.root" "../data/example_PAIReD_mcRun3_EE.root"
```
If everything worked as expected, you should see the code running for a couple of seconds and get `End makeNtuplesPAIReDjoint()` as a final printed statement. If this is the case, you can take a closer look at the output and you will notice that you just produced some first PAIReD jets. Congratulations!

### Produce entire data sets of PAIReD jets
Now, we can start producing PAIReD jets on larger scale by processing an entire NanoAOD data set with a single command. For this, we use the [`src/processDatasetToPAIReD.py`](../../PAIReD_Data_Production/PFNano_to_PAIReD/src/processDatasetToPAIReD.py) script (see also the [documentation](../../PAIReD_Data_Production/PFNano_to_PAIReD/doc/doc_processDataset.md)). Make sure that your `nano2paired` environment is activated and you have a valid `voms-proxy` set up:
```bash
conda activate nano2paired
voms-proxy-init --voms cms -valid 192:00
```
In order to start the data processing, we simply run the python script with the required arguments:
```bash
python processDatasetToPAIReD.py {DATASET} {OUTPUTPATH} \
-p {STORAGEPREFIX} -t {NTHREADS} --physicsprocess {PHYS} \
-g {GEOMETRY}
```
As you can see there is a number of required and optional arguments that we can specify. Three of them are **absolutely necessary** for the code to work:
- `DATASET`: The name of the data set as it is shown on CMSDAS. This could for example be something like `/WminusH_Hto2C_WtoLNu_M-125_TuneCP5_13p6TeV_powheg-pythia8/jaschulz-Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2_BTV_Run3_2022_Comm_PFNANOAODv12-490aef546cce8088aae92a89578a4232/USER`.
- `OUTPUTPATH`: The path to the directory where the produced PAIReD files should be stored. If it does not exist, it will be created for you. Could for example be `../data/test`.
- `STORAGEPREFIX`: This one is a bit more tricky. It is the prefix for the storage site the NanoAOD set is stored on. The set default for this belongs here to the T2 site of RWTH Aachen. In case you produced your NanoAOD set on any other site, you will need to specify the path prefix here. As an example: The storage prefix for T2 RWTH is `root://grid-cms-xrootd.physik.rwth-aachen.de/`.

In addition to those required arguments, there are a bunch of optional ones to further specify the processing parameters or the PAIReD jet itself. Those include the batchsize (`-b`), the number of files to be processed (`-n`) or the number of threads to use in parallel (`-t`). The shape/geometry of the PAIReD jets can be specified with `-g` which accepts the implemented geometries (`Ellipse`, `AK4`, `Dumbbell`). Finally, a free variable `--physicsprocess` can be specified. The chosen value will be added an additional variable to each PAIReD jet object. It is meant to be used for keeping track of the different event types used in training, so that one can weight the jets according to the event class later on. There are a couple more arguments explained and summarized in the [documentation](../../PAIReD_Data_Production/PFNano_to_PAIReD/doc/doc_processDataset.md). There you can also find more information on the arguments mentioned here.

Here is a concrete example on how to run it:
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


### Subdividing in training/validation and test sets
The procedure described above will give you data sets of PAIReD jets. However for training a PAIReD tagger we will have to devide this into two parts:
- Training/validation: One set that is used in the training stage for optimizing the network parameters and evaluating the performance over the course of training.
- Test: One set for objectively testing the tagger performance after training with a previously untouched set of PAIReD jets.

This subdivision is done by hand here in the following manner. The files produced in the previous step are named `PAIReD_*.root` where `*` is a unique index for each file. In our convention, we will simply rename some files to `test_PAIReD_*.root` instead to mark them as belonging to the test data set.

We can easily do that by using the [`src/utils/addToTestset.sh`](../../PAIReD_Data_Production/PFNano_to_PAIReD/src/utils/addToTestset.sh) helpers script:
```bash
bash utils/addToTestset.sh {PAIReDPATH} {i1} {i2} {i3} {...}
```
This will rename the `PAIReD_*.root` files with indices `i1`, `i2` and so on, to `test_PAIReD_*.root` files. E.g. `PAIReD_i1.root` to `test_PAIReD_i1.root`.

This action can be reversed using the script [`src/utils/rmFromTestset.sh`](../../PAIReD_Data_Production/PFNano_to_PAIReD/src/utils/rmFromTestset.sh) which works in an analogous way.

File can even be removed entirely from both training/validation by setting them *inactive* with the script [`src/utils/addToInactive.sh`](../../PAIReD_Data_Production/PFNano_to_PAIReD/src/utils/addToInactive.sh).

### Counting the produced PAIReD jets
Finally, one can count how many labeled PAIReD jets one in a given directory `PAIReDPATH` by running the [`src/utils/countPAIReD.py`](../../PAIReD_Data_Production/PFNano_to_PAIReD/src/utils/countPAIReD.py) script:
```bash
python utils/countPAIReD.py {PAIReDPATH}
```