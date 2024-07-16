## Production of the data sets
In order to create a tagger for PAIReD jets, data sets of such PAIReD jets are required first of all. As this type of jet is new in CMS, PAIReD jets are not yet available in standardized file formats such as MiniAOD or NanoAOD. Accordingly, we have to generate the data sets ourselves. We do this in two steps:
1. We create a NanoAOD data set from a data set in MiniAOD format.
2. From these NanoAOD files we create PAIReD files with ntuples of PAIReD jets.
Both steps are described below. 

**Note:** Technically, it is unnecessary to do this in two steps. A single step from MiniAOD directly to PAIReD ntuples is possible, but was not implemented here due to an unfavorable choice of approach at the beginning of the project. This may be implemented in the future.

### 1. From MiniAOD to PFNanoAOD
We use the BTV Nano production workflow to produce NanoAOD data sets. This is necessary because the PF candidates are required as tagger input and these are not available in normal NanoAODs. You can find the corresponding [repo here](https://github.com/cms-btv-pog/btvnano-prod). Follow the instructions given there to set up CMSSW and crab on `lxplus`.

**Important:** Before sending jobs to produce data sets, a setting must be changed in the CMSSW code. This is this line in the `custom_btv_cff.py` file:
```python
btvNano_addallPF_switch = cms.untracked.bool(False)
```
This must be changed to 
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

**Note:** I personally still used `lxplus7` which was ended recently. I have not tested the BTV-Nano workflow on the newer linux versions yet.

### 2. From PFNanoAOD to PAIReD

### Subdividing in Training+Validation and Test
