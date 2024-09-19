#!/bin/bash
outdir="/home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/data/TT"
temp_input="PFNano.root"
temp_output="PAIReD.root"
export X509_USER_PROXY=/isilon/data/users/trussel1/x509up_user.pem
here=${PWD}

mkdir -p execute/job_$3
cd execute/job_$3

cp -r /home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/src .
ulimit -s unlimited
cd src

cp /isilon/hadoop/store/group/common/bruxhcc/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v2_BTV_Run3_2023_Comm_MINIAODv4/240916_150242/0000/$1 $temp_input
source /home/trussel1/miniconda3/etc/profile.d/conda.sh
conda activate nano2paired

python processFileToPAIReD.py $temp_input ${temp_output} --physicsprocess 66

cp $temp_output $outdir/$2

cd $here/execute/
rm -r job_$3
