#!/bin/bash

#make log and err directories
mkdir -p log
mkdir -p err

#compile the code (if needed probable a good idea)
#g++ MakeSingleParticleMassPlots_5bins.C $(root-config --cflags --libs) -Wall -O2 -o MassPlots5bins.exe
#g++ temp1_5bins.C $(root-config --cflags --libs) -Wall -O2 -o temp1
#g++ temp_5bins.C $(root-config --cflags --libs) -Wall -O2 -o temp_5bins
g++ temp_5bins_differentialBinning.C $(root-config --cflags --libs) -Wall -O2 -o temp_diffBin

#submit the jobs with different arguments and redirect output 

./temp_diffBin 1 0.1 SB_out_20Nov_bin1.root newCode_newData/sidebands_01_newData_allbins.root  > log/outputb1.log 2> err/errorb1.err &
./temp_diffBin 2 0.1 SB_out_20Nov_bin2.root newCode_newData/sidebands_01_newData_allbins.root  > log/outputb2.log 2> err/errorb2.err &
./temp_diffBin 3 0.1 SB_out_20Nov_bin3.root newCode_newData/sidebands_01_newData_allbins.root  > log/outputb3.log 2> err/errorb3.err &
./temp_diffBin 4 0.1 SB_out_20Nov_bin4.root newCode_newData/sidebands_01_newData_allbins.root  > log/outputb4.log 2> err/errorb4.err &
./temp_diffBin 5 0.1 SB_out_20Nov_bin5.root newCode_newData/sidebands_01_newData_allbins.root  > log/outputb5.log 2> err/errorb5.err &


: '
./temp_5bins 0.1 XXXSB_out_19Nov_01.root newCode_newData/sidebands_01_newData.root > log/output01.log 2> err/error01.err &
./temp_5bins 0.01 XXXSB_out_19Nov_001.root newCode_newData/sidebands_001_newData.root > log/output001.log 2> err/error001.err &
./temp_5bins 0.05 XXXSB_out_19Nov_005.root newCode_newData/sidebands_005_newData.root > log/output005.log 2> err/error005.err &
./temp_5bins 0.2 XXXSB_out_19Nov_02.root newCode_newData/sidebands_02_newData.root > log/output02.log 2> err/error02.err &
./temp_5bins 0.4 XXXSB_out_19Nov_04.root newCode_newData/sidebands_04_newData.root > log/output04.log 2> err/error04.err &
'

: '
#./temp1 0.1 newCode_newData/sidebands_01_newData_allbins.root TH2F_08Nov_yesskip.root TF1_outputs_12Nov_01_allbins.root > log/output01.log 2> err/error01.err &
./temp1 0.1 newCode_newData/sidebands_01_newData.root TH2F_08Nov_yesskip.root TF1_newData_01.root > log/output01.log 2> err/error01.err &
./temp1 0.01 newCode_newData/sidebands_001_newData.root TH2F_12Nov_001.root TF1_outputs_12Nov_001.root > log/output001.log 2> err/error001.err &
./temp1 0.2 newCode_newData/sidebands_02_newData.root TH2F_12Nov_02.root TF1_outputs_12Nov_02.root > log/output02.log 2> err/error02.err &
./temp1 0.4 newCode_newData/sidebands_04_newData.root TH2F_12Nov_04.root TF1_outputs_12Nov_04.root > log/output04.log 2> err/error04.err &
./temp1 0.05 newCode_newData/sidebands_005_newData.root TH2F_12Nov_005.root TF1_outputs_12Nov_005.root > log/output005.log 2> err/error005.err &
'

#after MassPlots you need to run updateCreateTemplates.C to get output TF1 

: '
./MassPlots5bins.exe 0.4 TH2F_12Nov_04.root > log/output04.log 2> err/error04.err &  
./MassPlots5bins.exe 0.2 TH2F_12Nov_02.root > log/output02.log 2> err/error02.err &  
./MassPlots5bins.exe 0.05 TH2F_12Nov_005.root > log/output005.log 2> err/error005.err &  
./MassPlots5bins.exe 0.01 TH2F_12Nov_001.root > log/output001.log 2> err/error001.err &  
'
