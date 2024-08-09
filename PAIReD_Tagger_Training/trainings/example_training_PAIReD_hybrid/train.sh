#!/bin/bash
conda activate weaver

# name for the network (used in Tensorboard output)
NETNAME=myExampleNetwork
# hyperparameters
BATCHSIZE=512
LEARNINGRATE=5e-4
NEPOCHS=50
OPTIMIZER=ranger
NWORKERS=2
# path to the PAIReD_Tagger_Training directory
BASEDIR=../..
# path to the data set
DATADIR=/path/to/data

python $BASEDIR/weaver-core-hybrid/weaver/train.py \
--data-train \
"WpHcc:$DATA/WpHcc/PAIReD_*.root" \
"WmHcc:$DATA/WmHcc/PAIReD_*.root" \
"ZHcc:$DATA/ZHcc/PAIReD_*.root" \
"WpHbb:$DATA/WpHbb/PAIReD_*.root" \
"WmHbb:$DATA/WmHbb/PAIReD_*.root" \
"ZHbb:$DATA/ZHbb/PAIReD_*.root" \
"DY:$DATA/DY/PAIReD_*.root" \
--data-config $BASEDIR/dataconfig/PAIReD_hybrid_sv.train.yaml \
--network-config $BASEDIR/networks/PAIReD_ParT_sv_hybrid.py \
--batch-size $BATCHSIZE \
--start-lr $LEARNINGRATE \
--num-epochs $NEPOCHS \
--optimizer $OPTIMIZER \
--gpus 0 \
--samples-per-epoch 3000000 \
--samples-per-epoch-val 300000 \
--num-workers $NWORKERS --fetch-step 0.02 \
--model-prefix ./models/ \
--train-mode "hybrid" \
--log-file ./logfile.log \
--tensorboard "_$NETNAME" \
--tensorboard-dir $BASEDIR/trainings \
--loss-option "factor_reg" 0.1 \
--loss-option "factor_err" 0.1 \
--loss-option "baseline_reg" 0.0 \
--loss-option "baseline_err" 0.0
