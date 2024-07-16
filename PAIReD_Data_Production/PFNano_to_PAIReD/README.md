# PFNano_to_PAIReD
 Code to convert the output files of PFNano to a file format compatible with the PAIReD tagger network (see the [table of tree content](notes/PAIReD_ROOT_file_content.md) for the ouput ROOT file). See also the [table of data sets](notes/dataset_table.md) used and produced for trainings.

 This branch is currently under development and will contain in the end a working framework for producing an updated version of PAIReD jet construction. Compared to the [original PAIReD paper](https://doi.org/10.48550/arXiv.2311.11011), this version also contains finer background labels `ll`, `cl`, `cc`, `bl` and `bb` for PAIReD jets in addition to the usual `LL`, `CC` and `BB` labels. Furthermore several PAIReD jet geometries (Ellipse, AK4, Dumbbell) are implemented.

## Download
 In order to be able to use this module, you need to have a conda environment with `uproot`, `awkward` and `vector` installed.

### Installing the conda environment
 First, you will have to install Miniconda. See further instructions [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html).

 Then, create a new conda environment, e.g., named `nano2paired`:
 ```
conda create --name nano2paired python=3.10
 ```
 Activate the environment
 ```
conda activate nano2paired
 ```
 and install `uproot`, `vector` and `xrootd`:
 ```
conda install -c conda-forge uproot vector xrootd
 ```

### Cloning the Git repository
 Clone this repository to lxplus, e.g., via ssh:
 ```
git clone git@github.com:JanGerritSchulz/PFNano_to_PAIReD.git
 ```


## Documentation
 The documentation for the relevant files can be found here:
 * [Documentation for processFileToPAIReD](doc/doc_main.md)
 * [Documentation for tools/processEvents](doc/doc_processEvents.md)
 * [Documentation for tools/isInPAIReD/*](doc/doc_isInPAIReD.md)
 * [Documentation for tools/getMCInfo](doc/doc_getMCInfo.md)
 * [Documentation for tools/helpers](doc/doc_helpers.md)

Detailed examples of the functioning of the main function can be found, step by step, in the jupyter notebook [here](src/notebooks/makeNtuplesPAIReDjointMC.ipynb).

## Usage
 To perform a data conversion of a single file, one can call the <code>processFileToPAIReD</code>  python script in the <code>src</code> directory:
 ```
 python processFileToPAIReD.py "../data/example_PFNano_mcRun3_EE_allPF_noBTV.root" "../data/example_PAIReD_mcRun3_EE.root"
 ```
The first argument is the input file, the second is the output file. More on input arguments can be found in the [documentation](doc/doc_main.md).

## How to use PFNano
 The code of this repository is designed to convert [PFNano](https://github.com/cms-jet/PFNano) files to ROOT files that can be read by the PAIReD tagger. In order to get a meaningful output of all the jet pairs in the given input file, the file's content has to meet some criteria. Most importantly, it has to contain **all particle flow candidates** of an event (not only those inside, e.g., AK4 jets). This can be achieved if in the PFNano framework the function `addPFCands()` is called with the following arguments:
 ```python
 addPFCands(process,
            runOnMC = True,
            allPF = True,
            onlyAK4 = False,
            onlyAK8 = False)
 ```
 For example, the function `PFnano_customizeData_allPF(process)` could be used for this.

 Further information on how to produce PFNano files from MiniAODs in this [note](notes/produce_PFNano_files.md).

## Notes
* The particles are **not** cut in the following manner:
   * $|p_z| \leq 10000$
   * $|\eta| \leq 5$
   * $p_T > 0$
