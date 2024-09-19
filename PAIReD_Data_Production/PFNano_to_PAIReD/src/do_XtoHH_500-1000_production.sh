#!/bin/bash
outdir="/home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/data/XtoHH/MX-500-1000"
temp_input="PFNano.root"
temp_output="PAIReD.root"
export X509_USER_PROXY=/isilon/data/users/trussel1/x509up_user.pem
here=${PWD}

mkdir -p execute/job_$3
cd execute/job_$3

cp -r /home/trussel1/PAIReD_jet_tagging/PAIReD_Data_Production/PFNano_to_PAIReD/src .
ulimit -s unlimited
cd src

cp /home/trussel1/bruxhcc/XtoHH_MX-500to1000/$1 $temp_input
source /home/trussel1/miniconda3/etc/profile.d/conda.sh
conda activate nano2paired

python processFileToPAIReD.py $temp_input ${temp_output} --physicsprocess 2525

cp $temp_output $outdir/$2

cd $here/execute/
rm -r job_$3
