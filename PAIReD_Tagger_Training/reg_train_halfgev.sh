weaver --data-train \
"ZHvar:/HEP/export/home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/data/clustered_1GeV/ZHvar/haddedPAIReD_*.root" \
--data-test \
"ZHvar:/HEP/export/home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/data/clustered_1GeV/test/ZHvar/haddedPAIReD_*.root" \
--regression-mode \
--train-val-split 0.88 \
--data-config /users/trussel1/weaver-core/dataconfigs/PAIReD_regression_halfGeV.train.yaml \
--network-config /users/trussel1/weaver-core/networks/PAIReD_ParT_sv_regression.py \
--batch-size 256 \
--start-lr 5e-4 \
--num-epochs 80 \
--optimizer ranger \
--gpus 0,1,2,3,4,5 \
--samples-per-epoch 3000000 \
--samples-per-epoch-val 300000 \
--num-workers 11 \
--fetch-step 0.04 \
--model-prefix PAIReD_Tagger_Training/models/reg/halfgev_lowdr/models \
--log logs/training_0.log \
--tensorboard "_NetworkName" \
--no-remake-weights \
#--tensorboard-dir PAIReD_Tagger_Training/trainings \
#"ZHvar:/HEP/export/home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/data/clustered_1GeV/ZHvar/*/PAIReD_*.root" \
#"ZHvar:/HEP/export/home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/data/clustered_1GeV/test/ZHvar/*/PAIReD_*.root" \
#"DY:/HEP/export/home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/data/DY/PAIReD_*.root" \
#"DY:/HEP/export/home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/data/test/DY/PAIReD_*.root" \
