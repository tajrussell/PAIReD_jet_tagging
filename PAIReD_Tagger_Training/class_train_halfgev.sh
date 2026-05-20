weaver --data-train \
"DY:/HEP/export/home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/data/clustered_1GeV/DY/PAIReD_*.root" \
"TT:/HEP/export/home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/data/clustered_1GeV/TT/PAIReD_*.root" \
"ZHvar:/HEP/export/home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/data/clustered_1GeV/ZHvar/haddedPAIReD_*.root" \
--data-test \
"DY:/HEP/export/home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/data/clustered_1GeV/test/DY/PAIReD_*.root" \
"TT:/HEP/export/home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/data/clustered_1GeV/test/TT/PAIReD_*.root" \
"ZHvar:/HEP/export/home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/data/clustered_1GeV/test/ZHvar/haddedPAIReD_*.root" \
--train-val-split 0.88 \
--data-config /users/trussel1/weaver-core/dataconfigs/PAIReD_classification_halfGeV.train.yaml \
--network-config /users/trussel1/weaver-core/networks/PAIReD_ParT_sv_classifier.py \
--batch-size 512 \
--start-lr 1e-4 \
--num-epochs 80 \
--optimizer ranger \
--gpus 0,1,2,3 \
--samples-per-epoch 3000000 \
--samples-per-epoch-val 300000 \
--num-workers 11  \
--fetch-step 0.01 \
--model-prefix PAIReD_Tagger_Training/models/cls/halfgev_mjj/models \
--log logs/training_1.log \
--tensorboard "_NetworkName" \
--no-remake-weights \
#--tensorboard-dir PAIReD_Tagger_Training/trainings \
#"DY:/HEP/export/home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/data/DY/PAIReD_*.root" \
#"DY:/HEP/export/home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/data/test/DY/PAIReD_*.root" \
