#!/bin/bash
conda activate weaver

# path to the PAIReD_Tagger_Training directory
BASEDIR=../..
# path to the data set
DATADIR=/path/to/data

python $BASEDIR/weaver-core-hybrid/weaver/train.py \
--data-config $BASEDIR/dataconfig/PAIReD_hybrid_sv.test.yaml \
--network-config $BASEDIR/networks/PAIReD_ParT_sv_hybrid.py \
--model-prefix ./models/_best_epoch_state.pt \
--export-onnx ./model.onnx