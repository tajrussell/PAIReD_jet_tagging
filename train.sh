python PAIReD_Tagger_Training/weaver-core-hybrid/weaver/train.py --train-mode hybrid \
--data-train \
"DY:/HEP/export/home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/data/DY/PAIReD_*11.root" \
"XtoHH:/HEP/export/home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/data/XtoHH/*/PAIReD_*11.root" \
#"ttbar:/HEP/export/home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/data/TT/PAIReD_*.root" \
--data-test \
"DY:/HEP/export/home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/data/test/DY/PAIReD_*.root" \
"XtoHH:/HEP/export/home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/data/test/XtoHH/*/PAIReD_*.root" \
#"ttbar:/HEP/export/home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/data/test/TT/PAIReD_*.root" \
--train-val-split 0.88 \
--data-config /users/trussel1/PAIReD_jet_tagging/PAIReD_Tagger_Training/dataconfigs/PAIReD_hybrid_sv.train.yaml \
--network-config /users/trussel1/PAIReD_jet_tagging/PAIReD_Tagger_Training/networks/PAIReD_ParT_sv_hybrid.py \
--demo \
--batch-size 256 \
--start-lr 1e-4 \
--num-epochs 50 \
--optimizer ranger \
--gpus 0 \
--samples-per-epoch 300000 \
--samples-per-epoch-val 30000 \
--num-workers 11 \
--fetch-step 0.002 \
--model-prefix PAIReD_Tagger_Training/models \
--log-file PAIReD_Tagger_Training/log/training_0.log \
--tensorboard "_NetworkName" \
--tensorboard-dir PAIReD_Tagger_Training/trainings \
--loss-option "factor_reg" 0.1 \
--loss-option "factor_err" 0.1 \
--loss-option "baseline_reg" 0.0 \
--loss-option "baseline_err" 0.0