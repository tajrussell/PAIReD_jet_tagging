# Production of the data sets
In order to create a tagger for PAIReD jets, data sets of such PAIReD jets are required first of all. As this type of jet is new in CMS, PAIReD jets are not yet available in standardized file formats such as MiniAOD or NanoAOD. Accordingly, we have to generate the data sets ourselves. We do this in two steps:
1. We create a NanoAOD data set from a data set in MiniAOD format.
2. From these NanoAOD files we create PAIReD files with ntuples of PAIReD jets.
Both steps are described below. 

**Note:** Technically, it is unnecessary to do this in two steps. A single step from MiniAOD directly to PAIReD ntuples is possible, but was not implemented here due to an unfavorable choice of approach at the beginning of the project. This may be implemented in the future.

## 1. From MiniAOD to PFNanoAOD
We use the BTV Nano production workflow to produce NanoAOD data sets. This is necessary because the PF candidates are required as tagger input and these are not available in normal NanoAODs. You can find the corresponding [repo here](https://github.com/cms-btv-pog/btvnano-prod). Follow the instructions given there to set up CMSSW and crab on `lxplus`.

**Important:** Before sending jobs to produce data sets, a setting must be changed in the CMSSW code. The [following line](https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/custom_btv_cff.py#L658) in the `custom_btv_cff.py` file
```python
btvNano_addallPF_switch = cms.untracked.bool(False)
```
must be changed to  
```python
btvNano_addallPF_switch = cms.untracked.bool(True)
```
The reason for this is that in PAIReD jets there are potentially more PF candidates in the PAIReD jet than just those of the AK4 jet. Therefore, all PF candidates must be stored.
 
After the setup, a crab job for the production of NanoAOD data sets can be specified by a `yml` card following the instructions. It is important that the correct Python config file is specified in the `yml` card, an overview of this can be found under [How to Update](https://github.com/cms-btv-pog/btvnano-prod?tab=readme-ov-file#how-to-update) with the commands for generation.

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



### Subdividing in training/validation and test sets
