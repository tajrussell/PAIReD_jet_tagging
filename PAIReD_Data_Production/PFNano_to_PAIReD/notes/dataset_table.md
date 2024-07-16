# Tables of produced data sets / files

## First training with $m_H=125~\text{GeV}$ samples only
Since no variable mass samples were available in the beginning, the first studies are carried out with fixed Higgs mass samples. Here is an overview of the used data sets. The index number in the "Event type" column is the corresponding `MC_physics_process` value saved with the PAIReD jets in the training files. It helps keeping track of the data sets in training, especially regarding weights.

<table class="tg">
<thead>
  <tr>
    <th class="tg-0lax"></th>
    <th class="tg-dvid">Event type</th>
    <th class="tg-dvid">Data set name</th>
    <th class="tg-dvid">Parent data set<br>(MiniAOD)</th>
    <th class="tg-dvid">PFNano data set<br>(NanoAOD)</th>
    <th class="tg-dvid">Training files<br>(on RWTH Physics cluster)</th>
    <th class="tg-dvid">Number of PAIReD jets<br>(test set)<br>(ll:cc:bb:cl:bl:bb_ttbar)</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td class="tg-7ivu" rowspan="8">Hcc</td>
    <td class="tg-0pky">ttH<br><code>-</code></td>
    <td class="tg-0pky">ttHto2C_M-125_<br>TuneCP5_13p6TeV_powheg-pythia8 (124X)</td>
    <td class="tg-7od5"> <a href="https://cmsweb.cern.ch/das/request?view=list&amp;limit=50&amp;instance=prod%2Fglobal&amp;input=dataset%3D/ttHto2C_M-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v2/MINIAODSIM" target="_blank">DAS</a><br>(4.2M events)</td>
    <td class="tg-90e1"><a href="https://cmsweb.cern.ch/das/request?input=%2FttHto2C_M-125_TuneCP5_13p6TeV_powheg-pythia8%2Fjaschulz-Run3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v2_converted_to_PFNano_allPF_noBTV_NANOAODSIM-2b128d9a7ac74617343e5f44c7eb0fdf%2FUSER&instance=prod%2Fphys03" target="_blank">DAS</a><br>(4.2M events)</td>
    <td class="tg-90e1"></td>
    <td class="tg-90e1"></td>
  </tr>
  <tr>
    <td class="tg-0pky">ggZH<br><code>-</code></td>
    <td class="tg-0pky">ggZH_Hto2C_Zto2L_M-125_<br>TuneCP5_13p6TeV_powheg-pythia8 (124X)</td>
    <td class="tg-7od5"><a href="https://cmsweb.cern.ch/das/request?view=list&amp;limit=50&amp;instance=prod%2Fglobal&amp;input=dataset%3D/ggZH_Hto2C_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v2/MINIAODSIM" target="_blank">DAS</a><br>(3M events)</td>
    <td class="tg-90e1"><a href="https://cmsweb.cern.ch/das/request?input=%2FggZH_Hto2C_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8%2Fjaschulz-Run3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v2_converted_to_PFNano_allPF_noBTV_NANOAODSIM-2b128d9a7ac74617343e5f44c7eb0fdf%2FUSER&instance=prod%2Fphys03" target="_blank">DAS</a><br>(3M events)</td>
    <td class="tg-90e1"></td>
    <td class="tg-90e1"></td>
  </tr>
  <tr>
    <td class="tg-0lax" rowspan="4">WH<br><code>2/1</code></td>
    <td class="tg-0lax">WminusH_Hto2C_WtoLNu_M-125_<br>TuneCP5_13p6TeV_powheg-pythia8 (124X)</td>
    <td class="tg-yofg"><a href="https://cmsweb.cern.ch/das/request?view=list&amp;limit=50&amp;instance=prod%2Fglobal&amp;input=dataset%3D/WminusH_Hto2C_WtoLNu_M-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v2/MINIAODSIM" target="_blank">DAS</a><br>(3M events)</td>
    <td class="tg-ltxa"><a href="https://cmsweb.cern.ch/das/request?input=%2FWminusH_Hto2C_WtoLNu_M-125_TuneCP5_13p6TeV_powheg-pythia8%2Fjaschulz-Run3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v2_converted_to_PFNano_allPF_noBTV_NANOAODSIM-2b128d9a7ac74617343e5f44c7eb0fdf%2FUSER&instance=prod%2Fphys03" target="_blank">DAS</a><br>(3M events)</td>
    <td class="tg-ltxa"></td>
    <td class="tg-ltxa"></td>
  </tr>
  <tr>
    <td class="tg-0lax">WminusH_Hto2C_WtoLNu_M-125_<br>TuneCP5_13p6TeV_powheg-pythia8 (130X)</td>
    <td class="tg-yofg"><a href="https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset%3D%2FWminusH_Hto2C_WtoLNu_M-125_TuneCP5_13p6TeV_powheg-pythia8%2FRun3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2%2FMINIAODSIM" target="_blank">DAS</a><br>(3M events)</td>
    <td class="tg-ltxa"><a href="https://cmsweb.cern.ch/das/request?input=%2FWminusH_Hto2C_WtoLNu_M-125_TuneCP5_13p6TeV_powheg-pythia8%2Fjaschulz-Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2_BTV_Run3_2022_Comm_PFNANOAODv12-490aef546cce8088aae92a89578a4232%2FUSER&instance=prod%2Fphys03" target="_blank">DAS</a><br>(3M events)</td>
    <td class="tg-ltxa"><code>/net/scratch_cms3a/jgschulz/data/PAIReD/Hcc/WH/WminusH_130X</code></td>
    <td class="tg-ltxa"><b>5,942,568 (8%)</b><br>(585,039 : 2,496,761 : 0 : 2,860,768 : 0 : 0)</td>
  </tr>
  <tr>
    <td class="tg-0pky">WplusH_Hto2C_WtoLNu_M-125_<br>TuneCP5_13p6TeV_powheg-pythia8 (124X)</td>
    <td class="tg-7od5"><a href="https://cmsweb.cern.ch/das/request?view=list&amp;limit=50&amp;instance=prod%2Fglobal&amp;input=dataset%3D/WplusH_Hto2C_WtoLNu_M-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v2/MINIAODSIM" target="_blank">DAS</a><br>(2.9M events)</td>
    <td class="tg-90e1"><a href="https://cmsweb.cern.ch/das/request?input=%2FWplusH_Hto2C_WtoLNu_M-125_TuneCP5_13p6TeV_powheg-pythia8%2Fjaschulz-Run3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v2_converted_to_PFNano_allPF_noBTV_NANOAODSIM-2b128d9a7ac74617343e5f44c7eb0fdf%2FUSER&instance=prod%2Fphys03" target="_blank">DAS</a><br>(2.9M events)</td>
    <td class="tg-ltxa"></td>
    <td class="tg-ltxa"></td>
  </tr>
  <tr>
    <td class="tg-0pky">WplusH_Hto2C_WtoLNu_M-125_<br>TuneCP5_13p6TeV_powheg-pythia8 (130X)</td>
    <td class="tg-7od5"><a href="https://cmsweb.cern.ch/das/request?input=dataset%3D%2FWplusH_Hto2C_WtoLNu_M-125_TuneCP5_13p6TeV_powheg-pythia8%2FRun3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2%2FMINIAODSIM&instance=prod/global" target="_blank">DAS</a><br>(2.9M events)</td>
    <td class="tg-90e1"><a href="https://cmsweb.cern.ch/das/request?input=%2FWplusH_Hto2C_WtoLNu_M-125_TuneCP5_13p6TeV_powheg-pythia8%2Fjaschulz-Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2_BTV_Run3_2022_Comm_PFNANOAODv12-490aef546cce8088aae92a89578a4232%2FUSER&instance=prod%2Fphys03" target="_blank">DAS</a><br>(2.8M events)</td>
    <td class="tg-ltxa"><code>/net/scratch_cms3a/jgschulz/data/PAIReD/Hcc/WH/WplusH_130X</code></td>
    <td class="tg-ltxa"><b>5,118,382 (8%)</b><br>(577,530 : 2,103,581 : 0 : 2,437,271 : 0 : 0)</td>
  </tr>
  <tr>
    <td class="tg-0lax" rowspan="2">ZH<br><code>3</code></td>
    <td class="tg-0lax">ZH_Hto2C_Zto2L_M-125_<br>TuneCP5_13p6TeV_powheg-pythia8 (124X)</td>
    <td class="tg-yofg"><a href="https://cmsweb.cern.ch/das/request?view=list&amp;limit=50&amp;instance=prod%2Fglobal&amp;input=dataset%3D/ZH_Hto2C_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v2/MINIAODSIM" target="_blank">DAS</a><br>(3M events)</td>
    <td class="tg-ltxa"><a href="https://cmsweb.cern.ch/das/request?input=%2FZH_Hto2C_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8%2Fjaschulz-Run3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v2_converted_to_PFNano_allPF_noBTV_NANOAODSIM-2b128d9a7ac74617343e5f44c7eb0fdf%2FUSER&instance=prod%2Fphys03" target="_blank">DAS</a><br>(3M events)</td>
    <td class="tg-ltxa"></td>
    <td class="tg-90e1"></td>
  </tr>
  <tr><td class="tg-0lax">ZH_Hto2C_Zto2L_M-125_<br>TuneCP5_13p6TeV_powheg-pythia8 (130X)</td>
    <td class="tg-yofg"><a href="https://cmsweb.cern.ch/das/request?input=dataset%3D%2FZH_Hto2C_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8%2FRun3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2%2FMINIAODSIM&instance=prod/phys03" target="_blank">DAS</a><br>(2.9M events)</td>
    <td class="tg-ltxa"><a href="https://cmsweb.cern.ch/das/request?input=%2FZH_Hto2C_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8%2Fjaschulz-Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2_BTV_Run3_2022_Comm_PFNANOAODv12-490aef546cce8088aae92a89578a4232%2FUSER&instance=prod%2Fphys03" target="_blank">DAS</a><br>(1.6M events)</td>
    <td class="tg-ltxa"><code>/net/scratch_cms3a/jgschulz/data/PAIReD/Hcc/ZH/ZH_130X</code></td>
    <td class="tg-ltxa"><b>3,782,066 (9%)</b><br>(443,565 : 1,440,612 : 0 : 1,897,889 : 0 : 0)</td>
  </tr>

  
  <tr>
    <td class="tg-7ivu" rowspan="8">Hbb</td>
    <td class="tg-0pky">ttH<br><code>-</code></td>
    <td class="tg-0pky">ttHto2B_M-125_<br>TuneCP5_13p6TeV_powheg-pythia8 (124X)</td>
    <td class="tg-7od5"> <a href="https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset%3D%2FttHto2B_M-125_TuneCP5_13p6TeV_powheg-pythia8%2FRun3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v2%2FMINIAODSIM" target="_blank">DAS</a><br>(4.1M events)</td>
    <td class="tg-90e1"></td>
    <td class="tg-90e1"></td>
    <td class="tg-90e1"></td>
  </tr>
  <tr>
    <td class="tg-0pky">ggZH<br><code>-</code></td>
    <td class="tg-0pky">ggZH_Hto2B_Zto2L_M-125_<br>TuneCP5_13p6TeV_powheg-pythia8 (124X)</td>
    <td class="tg-7od5"><a href="https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset%3D%2FggZH_Hto2B_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8%2FRun3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v2%2FMINIAODSIM" target="_blank">DAS</a><br>(2.1M events)</td>
    <td class="tg-90e1"></td>
    <td class="tg-90e1"></td>
    <td class="tg-90e1"></td>
  </tr>
  <tr>
    <td class="tg-0lax" rowspan="4">WH<br><code>12/11</code></td>
    <td class="tg-0lax">WminusH_Hto2B_WtoLNu_M-125_<br>TuneCP5_13p6TeV_powheg-pythia8 (124X)</td>
    <td class="tg-yofg"><a href="https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset%3D%2FWminusH_Hto2B_WtoLNu_M-125_TuneCP5_13p6TeV_powheg-pythia8%2FRun3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v2%2FMINIAODSIM" target="_blank">DAS</a><br>(2M events)</td>
    <td class="tg-ltxa"><a href="https://cmsweb.cern.ch/das/request?input=%2FWminusH_Hto2B_WtoLNu_M-125_TuneCP5_13p6TeV_powheg-pythia8%2Fjaschulz-Run3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v2_converted_to_PFNano_allPF_noBTV_NANOAODSIM-2b128d9a7ac74617343e5f44c7eb0fdf%2FUSER&instance=prod%2Fphys03" target="_blank">DAS</a><br>(2M events)</td>
    <td class="tg-ltxa"></td>
    <td class="tg-ltxa"></td>
  </tr>
  <tr>
    <td class="tg-0lax">WminusH_Hto2B_WtoLNu_M-125_<br>TuneCP5_13p6TeV_powheg-pythia8 (130X)</td>
    <td class="tg-yofg"><a href="https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset%3D%2FWminusH_Hto2B_WtoLNu_M-125_TuneCP5_13p6TeV_powheg-pythia8%2FRun3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2%2FMINIAODSIM" target="_blank">DAS</a><br>(2M events)</td>
    <td class="tg-ltxa"><a href="https://cmsweb.cern.ch/das/request?input=%2FWminusH_Hto2B_WtoLNu_M-125_TuneCP5_13p6TeV_powheg-pythia8%2Fjaschulz-Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2_BTV_Run3_2022_Comm_PFNANOAODv12-490aef546cce8088aae92a89578a4232%2FUSER&instance=prod%2Fphys03" target="_blank">DAS</a><br>(1.9M events)</td>
    <td class="tg-ltxa"><code>/net/scratch_cms3a/jgschulz/data/PAIReD/Hbb/WH/WminusH_130X</code></td>
    <td class="tg-ltxa"><b>3,353,826 (13%)</b><br>(391,267 : 0 : 1,317,760 : 0 : 1,644,799 : 0)</td>
  </tr>
  <tr>
    <td class="tg-0pky">WplusH_Hto2B_WtoLNu_M-125_<br>TuneCP5_13p6TeV_powheg-pythia8 (124X)</td>
    <td class="tg-7od5"><a href="https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset%3D%2FWplusH_Hto2B_WtoLNu_M-125_TuneCP5_13p6TeV_powheg-pythia8%2FRun3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v2%2FMINIAODSIM" target="_blank">DAS</a><br>(2M events)</td>
    <td class="tg-90e1"><a href="https://cmsweb.cern.ch/das/request?input=%2FWplusH_Hto2B_WtoLNu_M-125_TuneCP5_13p6TeV_powheg-pythia8%2Fjaschulz-Run3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v2_converted_to_PFNano_allPF_noBTV_NANOAODSIM-2b128d9a7ac74617343e5f44c7eb0fdf%2FUSER&instance=prod%2Fphys03" target="_blank">DAS</a><br>(2M events)</td>
    <td class="tg-ltxa"></td>
    <td class="tg-ltxa"></td>
  </tr>
  <tr>
    <td class="tg-0pky">WplusH_Hto2B_WtoLNu_M-125_<br>TuneCP5_13p6TeV_powheg-pythia8 (130X)</td>
    <td class="tg-7od5"><a href="https://cmsweb.cern.ch/das/request?input=dataset%3D%2FWplusH_Hto2B_WtoLNu_M-125_TuneCP5_13p6TeV_powheg-pythia8%2FRun3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2%2FMINIAODSIM&instance=prod/global" target="_blank">DAS</a><br>(2M events)</td>
    <td class="tg-90e1"><a href="https://cmsweb.cern.ch/das/request?input=%2FWplusH_Hto2B_WtoLNu_M-125_TuneCP5_13p6TeV_powheg-pythia8%2Fjaschulz-Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2_BTV_Run3_2022_Comm_PFNANOAODv12-490aef546cce8088aae92a89578a4232%2FUSER&instance=prod%2Fphys03" target="_blank">DAS</a><br>(2M events)</td>
    <td class="tg-ltxa"><code>/net/scratch_cms3a/jgschulz/data/PAIReD/Hbb/WH/WplusH_130X</code></td>
    <td class="tg-ltxa"><b>3,429,063 (12%)</b><br>(451,013 : 0 : 1,323,027 : 0 : 1,655,023 : 0)</td>
  </tr>
  <tr>
    <td class="tg-0lax" rowspan="2">ZH<br><code>13</code></td>
    <td class="tg-0lax">ZH_Hto2B_Zto2L_M-125_<br>TuneCP5_13p6TeV_powheg-pythia8 (124X)</td>
    <td class="tg-yofg"><a href="https://cmsweb.cern.ch/das/request?input=dataset%3D%2FZH_Hto2B_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8%2FRun3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v1%2FMINIAODSIM&instance=prod/phys03" target="_blank">DAS</a><br>(2M events)<br>(Note: this one is <code>v1-v1</code> instead of <code>v1-v2</code>)</td>
    <td class="tg-ltxa"><a href="https://cmsweb.cern.ch/das/request?input=%2FZH_Hto2B_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8%2Fjaschulz-Run3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v1_BTV_Run3_2022_Comm_PFNANOAODv12-2ed60bf267d29da512b15586b8ac351a%2FUSER&instance=prod%2Fphys03" target="_blank">DAS</a><br>(2M events)</td>
    <td class="tg-ltxa"><code>/net/scratch_cms3a/jgschulz/data/PAIReD/Hbb/ZH/ZH_124X</code></td>
    <td class="tg-ltxa"><b>4,071,471 (13%)</b><br>(552,200 : 0 : 1,455,570 : 0 : 2,063,701 : 0)</td>
  </tr>
  <tr>
    <td class="tg-0lax">ZH_Hto2B_Zto2L_M-125_<br>TuneCP5_13p6TeV_powheg-pythia8 (130X)</td>
    <td class="tg-yofg"><a href="https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fphys03&input=dataset%3D%2FZH_Hto2B_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8%2FRun3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2%2FMINIAODSIM" target="_blank">DAS</a><br>(-M events)</td>
    <td class="tg-ltxa"></td>
    <td class="tg-ltxa"></td>
    <td class="tg-ltxa"></td>
  </tr>

  
  <tr>
    <td class="tg-7ivu" rowspan="3">QCD</td>
    <td class="tg-0lax" rowspan="3">4 Jets<br><code>21</code></td>
    <td class="tg-0lax">QCD-4Jets_HT-70to100_<br>TuneCP5_13p6TeV_madgraphMLM-pythia8 (124X)</td>
    <td class="tg-yofg"><a href="https://cmsweb.cern.ch/das/request?input=dataset%3D%2FQCD-4Jets_HT-70to100_TuneCP5_13p6TeV_madgraphMLM-pythia8%2FRun3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v2%2FMINIAODSIM&instance=prod/global" target="_blank">DAS</a><br>(66M events)</td>
    <td class="tg-ltxa"><a href="https://cmsweb.cern.ch/das/request?input=%2FQCD-4Jets_HT-70to100_TuneCP5_13p6TeV_madgraphMLM-pythia8%2Fjaschulz-Run3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v2_converted_to_PFNano_allPF_noBTV_NANOAODSIM-2b128d9a7ac74617343e5f44c7eb0fdf%2FUSER&instance=prod%2Fphys03" target="_blank">DAS</a><br>(10M events)</td>
    <td class="tg-ltxa"></td>
    <td class="tg-ltxa"></td>
  </tr>
  <tr>
    <td class="tg-0lax">QCD-4Jets_HT-200to400_<br>TuneCP5_13p6TeV_madgraphMLM-pythia8 (124X)</td>
    <td class="tg-yofg"><a href="https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset%3D%2FQCD-4Jets_HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8%2FRun3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v2%2FMINIAODSIM" target="_blank">DAS</a><br>(72M events)</td>
    <td class="tg-ltxa"><a href="https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fphys03&input=%2FQCD-4Jets_HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8%2Fjaschulz-Run3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v2_converted_to_PFNano_allPF_noBTV_NANOAODSIM-2b128d9a7ac74617343e5f44c7eb0fdf%2FUSER" target="_blank">DAS</a><br>(9M events)</td>
    <td class="tg-ltxa"></td>
    <td class="tg-ltxa"></td>
  </tr>
  <tr>
    <td class="tg-0lax">QCD-4Jets_HT-200to400_<br>TuneCP5_13p6TeV_madgraphMLM-pythia8 (130X)</td>
    <td class="tg-yofg"><a href="https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset%3D%2FQCD-4Jets_HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8%2FRun3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2%2FMINIAODSIM" target="_blank">DAS</a><br>(70M events)</td>
    <td class="tg-ltxa"><a href="https://cmsweb.cern.ch/das/request?input=%2FQCD-4Jets_HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8%2Fjaschulz-Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2_BTV_Run3_2022_Comm_PFNANOAODv12-490aef546cce8088aae92a89578a4232%2FUSER&instance=prod%2Fphys03" target="_blank">DAS</a><br>(5M events)</td>
    <td class="tg-ltxa"><code>/net/scratch_cms3a/jgschulz/data/PAIReD/QCD/4jets/HT-200to400_130X</code></td>
    <td class="tg-ltxa"><b>14,802,719 (2%)</b><br>(11,566,361 : 1,263,058 : 451,751 : 1,024,962 : 496,587 : 0)</td>
  </tr>

  <tr>
    <td class="tg-7ivu" rowspan="2">DY</td>
    <td class="tg-0lax" rowspan="2">DY + Jets<br><code>23</code></td>
    <td class="tg-0lax">DYto2L-4Jets_MLL-50_<br>TuneCP5_13p6TeV_madgraphMLM-pythia8 (130X)</td>
    <td class="tg-yofg"><a href="https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fphys03&input=dataset%3D%2FDYto2L-4Jets_MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8%2FRun3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2%2FMINIAODSIM">DAS</a><br>(240M events)</td>
    <td class="tg-ltxa"><a href="https://cmsweb.cern.ch/das/request?input=%2FDYto2L-4Jets_MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8%2Fjaschulz-Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2_BTV_Run3_2022_Comm_PFNANOAODv12-490aef546cce8088aae92a89578a4232%2FUSER&instance=prod%2Fphys03" target="_blank">DAS</a><br>(10M events)</td>
    <td class="tg-ltxa"></td>
    <td class="tg-ltxa"></td>
  </tr>
  <tr>
    <td class="tg-0lax">DYto2L-4Jets_MLL-50_2J_<br>TuneCP5_13p6TeV_madgraphMLM-pythia8 (130X)</td>
    <td class="tg-yofg"><a href="https://cmsweb.cern.ch/das/request?input=dataset%3D%2FDYto2L-4Jets_MLL-50_2J_TuneCP5_13p6TeV_madgraphMLM-pythia8%2FRun3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2%2FMINIAODSIM&instance=prod/phys03">DAS</a><br>(10M? events)</td>
    <td class="tg-ltxa"><a href="https://cmsweb.cern.ch/das/request?input=%2FDYto2L-4Jets_MLL-50_2J_TuneCP5_13p6TeV_madgraphMLM-pythia8%2Fjaschulz-Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2_BTV_Run3_2022_Comm_PFNANOAODv12-490aef546cce8088aae92a89578a4232%2FUSER&instance=prod%2Fphys03" target="_blank">DAS</a><br>(10M events)</td>
    <td class="tg-ltxa"><code>/net/scratch_cms3a/jgschulz/data/PAIReD/DY/DYto2L_2Jets_130X</code></td>
    <td class="tg-ltxa"><b>7,051,509 (5%)</b><br>(5,891,173 : 212,918 : 65,606 : 594,065 : 287,747 : 0)</td>
  </tr>

  <tr>
    <td class="tg-7ivu" rowspan="1">TT</td>
    <td class="tg-0lax" rowspan="1">TT<br><code>24</code></td>
    <td class="tg-0lax">TTtoLNu2Q_<br>TuneCP5_13p6TeV_powheg-pythia8 (130X)</td>
    <td class="tg-yofg"><a href="https://cmsweb.cern.ch/das/request?input=dataset%3D%2FTTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8%2FRun3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2%2FMINIAODSIM&instance=prod/global">DAS</a><br>(268M events)</td>
    <td class="tg-ltxa"><a href="https://cmsweb.cern.ch/das/request?input=%2FTTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8%2Fjaschulz-Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2_BTV_Run3_2022_Comm_PFNANOAODv12-490aef546cce8088aae92a89578a4232%2FUSER&instance=prod%2Fphys03" target="_blank">DAS</a><br>(20M events)</td>
    <td class="tg-ltxa"><code>/net/scratch_cms3a/jgschulz/data/PAIReD/TT/TTtoLNu2Q_130X</code></td>
    <td class="tg-ltxa"><b>14,604,714 (0%)</b><br>(1,735,944 : 147,655 : 0 : 920,182 : 7,938,798 : 3,862,135)</td>
  </tr>

  <tr>
    <td class="tg-7ivu" rowspan="4">V</td>
    <td class="tg-0lax" rowspan="4">W + Jets<br><code>22</code></td>
    <td class="tg-0lax">WtoLNu-4Jets_<br>TuneCP5_13p6TeV_madgraphMLM-pythia8 (130X)</td>
    <td class="tg-yofg"><a href="https://cmsweb.cern.ch/das/request?input=dataset%3D%2FWtoLNu-4Jets_TuneCP5_13p6TeV_madgraphMLM-pythia8%2FRun3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2%2FMINIAODSIM&instance=prod/global">DAS</a><br>(342M events)</td>
    <td class="tg-ltxa"><a href="https://cmsweb.cern.ch/das/request?input=%2FWtoLNu-4Jets_TuneCP5_13p6TeV_madgraphMLM-pythia8%2Fjaschulz-Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2_BTV_Run3_2022_Comm_PFNANOAODv12-490aef546cce8088aae92a89578a4232%2FUSER&instance=prod%2Fphys03" target="_blank">DAS</a><br>(10M events)</td>
    <td class="tg-ltxa"><code>/net/scratch_cms3a/jgschulz/data/PAIReD/V/WtoLNu_4Jets_130X</code></td>
    <td class="tg-ltxa"><b>1,397,908 (0%)</b><br>(1,191,427 : 41,799 : 10,929 : 140,092 : 13,661 : 0)</td>
  </tr>
  <tr>
    <td class="tg-0lax">WtoLNu-4Jets_2J_<br>TuneCP5_13p6TeV_madgraphMLM-pythia8 (130X)</td>
    <td class="tg-yofg"><a href="https://cmsweb.cern.ch/das/request?input=dataset%3D%2FWtoLNu-4Jets_2J_TuneCP5_13p6TeV_madgraphMLM-pythia8%2FRun3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2%2FMINIAODSIM&instance=prod/global">DAS</a><br>(36M events)</td>
    <td class="tg-ltxa"><a href="https://cmsweb.cern.ch/das/request?input=%2FWtoLNu-4Jets_2J_TuneCP5_13p6TeV_madgraphMLM-pythia8%2Fjaschulz-Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2_BTV_Run3_2022_Comm_PFNANOAODv12-490aef546cce8088aae92a89578a4232%2FUSER&instance=prod%2Fphys03">DAS</a><br>(20M events)</td>
    <td class="tg-ltxa"><code>/net/scratch_cms3a/jgschulz/data/PAIReD/V/WtoLNu_2Jets_130X</code></td>
    <td class="tg-ltxa"><b>10,675,837 (0%)</b><br>(8,562,203 : 319,143 : 57,533 : 1,629,107 : 107,851 : 0)</td>
  </tr>
  <tr>
    <td class="tg-0lax">WtoLNu-4Jets_3J_<br>TuneCP5_13p6TeV_madgraphMLM-pythia8 (130X)</td>
    <td class="tg-yofg"><a href="https://cmsweb.cern.ch/das/request?input=dataset%3D%2FWtoLNu-4Jets_3J_TuneCP5_13p6TeV_madgraphMLM-pythia8%2FRun3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2%2FMINIAODSIM&instance=prod/global">DAS</a><br>(28M events)</td>
    <td class="tg-ltxa"><a href="https://cmsweb.cern.ch/das/request?input=%2FWtoLNu-4Jets_3J_TuneCP5_13p6TeV_madgraphMLM-pythia8%2Fjaschulz-Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2_BTV_Run3_2022_Comm_PFNANOAODv12-490aef546cce8088aae92a89578a4232%2FUSER&instance=prod%2Fphys03">DAS</a><br>(10M events)</td>
    <td class="tg-ltxa"><code>/net/scratch_cms3a/jgschulz/data/PAIReD/V/WtoLNu_3Jets_130X</code></td>
    <td class="tg-ltxa"></td>
  </tr>
  <tr>
    <td class="tg-0lax">WtoLNu-4Jets_4J_<br>TuneCP5_13p6TeV_madgraphMLM-pythia8 (130X)</td>
    <td class="tg-yofg"><a href="https://cmsweb.cern.ch/das/request?input=dataset%3D%2FWtoLNu-4Jets_4J_TuneCP5_13p6TeV_madgraphMLM-pythia8%2FRun3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2%2FMINIAODSIM&instance=prod/global">DAS</a><br>(5M events)</td>
    <td class="tg-ltxa"></td>
    <td class="tg-ltxa"></td>
    <td class="tg-ltxa"></td>
  </tr>
</tbody>
</table>
