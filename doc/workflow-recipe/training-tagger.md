# Training a PAIReD tagger

Here, we want to train a neural network on the classification of PAIReD jets and their mass regression. For this, we need to have dataset prepared. So, if you don't have a dataset yet, go back to the data production step described [here](./production-datasets.md). If you have your dataset ready, then let's get started on the neural network training!

#### Content
<!-- TOC -->

- [Set up conda environment](#set-up-the-conda-environment)
- [Training a network](#training-a-network)
    - [The network file](#the-network-file)
    - [The data configuration file](#the-data-configuration-file)
    - [Running the training](#running-the-training)
- [Testing a network](#testing-a-network)
- [Exporting a network to onnx](#exporting-a-network-to-onnx)

<!-- /TOC -->

## Set up the conda environment
For the training and inference, we will need a different conda environment than in the data production. If you don't have `conda` installed yet, go back to [Installing `conda`](./production-datasets.md#installing-conda). If you have conda installed, proceed with the following steps to set up the `weaver` environment.

Create a new environment:
```bash
conda create -n weaver python=3.10
```
Activate the environment:
```bash
conda activate weaver
```
Install `PyTorch` by following instructions for your OS/CUDA version on [the PyTorch website](https://pytorch.org/get-started). This could look for example like this:
```bash
pip install torch
```
Then, continue by installing the packages `weaver-core`, `numba` and `onnx`:
```bash
pip install weaver-core numba onnx
```
Once you have the environment installed, you are good to go for training your networks.

## Training a network
For training a neural network, you will need to set up and specify two important files. One is the data configuration file. It is a `yaml` file that specifies the inputs and outputs for your network as well as data selection and weighting of jets. The other one is the network file. This one is a `python` file that must contain two mandatory functions: `get_model()` and `get_loss()`. The former returns the network model and therefore specifies the architecture used. The latter defines the loss function which will be minimized in training.

### The network file
An example for a valid network file is [`PAIReD_ParT_sv_hybrid.py`](../../PAIReD_Tagger_Training/networks/PAIReD_ParT_sv_hybrid.py). In fact, if you want to train a hybrid PAIReD tagger that does classification and mass regression at the same time, you won't have to change anything and can simply use this network configuration. The used network architecture is the ParticleTransformer with PF candidates and SVs as inputs. More details on the implementation are given [here](../../PAIReD_Tagger_Training/notes/hybrid-tagger.md). If you want to do just classification without mass regression, you can use [`PAIReD_ParT_sv_classifier.py`](../../PAIReD_Tagger_Training/networks/PAIReD_ParT_sv_classifier.py).

### The data configuration file
An example for a valid data configuration file is [`PAIReD_hybrid_sv.train.yaml`](../../PAIReD_Tagger_Training/dataconfigs/PAIReD_hybrid_sv.train.yaml). If you check it out, you will find several sections where you can specify different things:
- `selection`: here you can define selection cuts that the jets must pass in order to be used.
- `preprocess`: here you can specify if your input data should be normalized before passing it to the network.
- `new_variables`: in this section you can define new variables based on the branches that exist in the ROOT data files.
- `inputs`: here you specify the input branches given to the network. Those come typically in subcategories. For ParT, you have for example the PF candidate features, the four-momenta of the PF candidates, the SV features and the four-momenta of the SVs.
- `labels`: here you first set the type of training (`hybrid` for a hybrid tagger, `simple` for a pure classifier). Then you set the target labels under `value` and if needed, the regression target under `value_custom`.
- `observers`: those are branches to be saved in the output ROOT files when running the inference/prediction.
- `weights`: here you can specify how your jets should be weighted. This should always be done in unbalanced datasets.

**Note:** If you use the given [`PAIReD_ParT_sv_hybrid.py`](../../PAIReD_Tagger_Training/networks/PAIReD_ParT_sv_hybrid.py), the first two `labels` **must** be the Higgs labels CC and BB. This is because the loss function treats the first two labels differently and uses the corresponding jets for the mass regression while the others are ignored for this.

### Running the training
For running the training, you actually just have to run a single python script in the terminal. Since it is a rather lengthy command, I would recommend writing it down in short bash script and running this script instead (see e.g. [training bash script `train.sh`](../../PAIReD_Tagger_Training/trainings/example_training_PAIReD_hybrid/train.sh)). The command to run could look like this:
```bash
python PAIReD_Tagger_Training/weaver-core-hybrid/weaver/train.py \
--train-mode "hybrid" \
--data-train \
"WpHcc:Path/to/data/WpHcc/PAIReD_*.root" \
"WmHcc:Path/to/data/WmHcc/PAIReD_*.root" \
"ZHcc:Path/to/data/ZHcc/PAIReD_*.root" \
"WpHbb:Path/to/data/WpHbb/PAIReD_*.root" \
"WmHbb:Path/to/data/WmHbb/PAIReD_*.root" \
"ZHbb:Path/to/data/ZHbb/PAIReD_*.root" \
"DY:Path/to/data/DY/PAIReD_*.root" \
--data-config PAIReD_Tagger_Training/dataconfig/PAIReD_hybrid_sv.train.yaml \
--network-config PAIReD_Tagger_Training/networks/PAIReD_ParT_sv_hybrid.py \
--batch-size 512 \
--start-lr 1e-4 \
--num-epochs 50 \
--optimizer ranger \
--gpus 0 \
--samples-per-epoch 3000000 \
--samples-per-epoch-val 300000 \
--num-workers 2 \
--fetch-step 0.02 \
--model-prefix Path/to/models/ \
--log-file Path/to/logfile.log \
--tensorboard "_NetworkName" \
--tensorboard-dir PAIReD_Tagger_Training/trainings \
--loss-option "factor_reg" 0.1 \
--loss-option "factor_err" 0.1 \
--loss-option "baseline_reg" 0.0 \
--loss-option "baseline_err" 0.0
```
Note that your `weaver` environment has to be activated for this. The most important given settings mean the following:
- `train-mode`: set `hybrid` for hybrid tagger and `cls` for a pure classifier.
- `data-train`: files with the training dataset. You can give multiple directories with different processes as shown here.
- `data-config`: path to your data configuration file.
- `network-config`: path to your network configuration file.
- `batch-size`: batchsize.
- `start-lr`: learning rate at the beginning of training (can change over time depending on the setting of `--lr-scheduler`).
- `num-epochs`: number of epochs to run the training over.
- `optimizer`: the optimizing algorithm to minimize the loss and optimize the network with (typically `ranger` is a fine choice).
- `gpus`: GPUs of the given machine that should be used. Note `0` means GPU with index 0 is taken. If you want no GPU (so run on CPU only) set it to `""`.
- `num-workers`: number of workers/threads for parallelization.
- `model-prefix`: path to save the pytorch models to.
- `log-file`: path for the log file.
- `loss-option`: set the parameters for the loss function, e.g. scaling factors or baselines. See [hybrid tagger note](../../PAIReD_Tagger_Training/notes/hybrid-tagger.md#loss-function) for more details.


## Testing a network
In order to get the predictions for jets from the test dataset, we have to run the inference task. This works very similar to the training task by running a command like in the [test bash script `test.sh`](../../PAIReD_Tagger_Training/trainings/example_training_PAIReD_hybrid/test.sh) or the command here:
```bash
python PAIReD_Tagger_Training/weaver-core-hybrid/weaver/train.py \
--predict \
--train-mode "hybrid" \
--data-test \
"WpHcc:Path/to/data/WpHcc/test_PAIReD_*.root" \
"WmHcc:Path/to/data/WmHcc/test_PAIReD_*.root" \
"ZHcc:Path/to/data/ZHcc/test_PAIReD_*.root" \
"WpHbb:Path/to/data/WpHbb/test_PAIReD_*.root" \
"WmHbb:Path/to/data/WmHbb/test_PAIReD_*.root" \
"ZHbb:Path/to/data/ZHbb/test_PAIReD_*.root" \
"DY:Path/to/data/DY/test_PAIReD_*.root" \
--data-config PAIReD_Tagger_Training/dataconfig/PAIReD_hybrid_sv.test.yaml \
--network-config PAIReD_Tagger_Training/networks/PAIReD_ParT_sv_hybrid.py \
--batch-size 512 \
--gpus 0 \
--num-workers 1 \
--model-prefix path/to/models/_best_epoch_state.pt \
--predict-output path/to/predictions/prediction.root
```
Note that:
- we added the `--predict` flag for running the prediction.
- we changed the data file names from `PAIReD_*.root` to `test_PAIReD_*.root` in order to use the test dataset.
- we also changed `--data-train` to `--data-test`.
- we exchanged the data-config file. The reason is that in the config file for training, we had some selection criteria for the jets, but we want to run the inference on all jets.
- the `--model-prefix` now points to the exact model file to use in the inference.

## Exporting a network to onnx
For exporting the trained model in the onnx format, we can run the command:
```bash
python PAIReD_Tagger_Training/weaver-core-hybrid/weaver/train.py \
--data-config PAIReD_Tagger_Training/dataconfig/PAIReD_hybrid_sv.test.yaml \
--network-config PAIReD_Tagger_Training/networks/PAIReD_ParT_sv_hybrid.py \
--model-prefix path/to/models/_best_epoch_state.pt \
--export-onnx path/to/export/model.onnx
```

Note that the dynamic axes for onnx seem to work only with `PyTorch` version below `2`.