#!/bin/bash
outdir="/home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/data/clustered_1GeV/XtoHH/MX-500-1000"
outdir_test="/home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/data/clustered_1GeV/test/XtoHH/MX-500-1000"
temp_input="PFNano.root"
temp_output="PAIReD.root"
temp_output_test="PAIReD_test.root"
export X509_USER_PROXY=/isilon/data/users/trussel1/x509up_user.pem
here=${PWD}

mkdir -p execute/job_X51CH_$3
cd execute/job_X51CH_$3

cp /home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/src/processFileToPAIReD.py .
cp -r /home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/src/tools .
ulimit -s unlimited
cd src

cp /home/trussel1/bruxhcc/XtoHH_MX-500to1000/$1 $temp_input
source /home/trussel1/miniconda3/etc/profile.d/conda.sh
conda activate nano2paired

python processFileToPAIReD.py $temp_input ${temp_output} --physicsprocess 25 -g Clustered_HighPT

cp $temp_output $outdir/$2
cp $temp_output_test $outdir_test/$2

cd $here/execute/
rm -r job_X51CH_$3
