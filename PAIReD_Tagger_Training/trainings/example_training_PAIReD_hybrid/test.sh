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
--predict \
--data-test \
"WpHcc:$DATA/WpHcc/test_PAIReD_*.root" \
"WmHcc:$DATA/WmHcc/test_PAIReD_*.root" \
"ZHcc:$DATA/ZHcc/test_PAIReD_*.root" \
"WpHbb:$DATA/WpHbb/test_PAIReD_*.root" \
"WmHbb:$DATA/WmHbb/test_PAIReD_*.root" \
"ZHbb:$DATA/ZHbb/test_PAIReD_*.root" \
"DY:$DATA/DY/test_PAIReD_*.root" \
--data-config $BASEDIR/dataconfig/PAIReD_hybrid_sv.test.yaml \
--network-config $BASEDIR/networks/PAIReD_ParT_sv_hybrid.py \
--batch-size $BATCHSIZE \
--gpus 0 \
--num-workers $NWORKERS \
--train-mode "hybrid" \
--model-prefix ./models/_best_epoch_state.pt \
--predict-output ./predictions/prediction.root \
--loss-option "factor_reg" 0.1 \
--loss-option "factor_err" 0.1 \
--loss-option "baseline_reg" 0.0 \
--loss-option "baseline_err" 0.0