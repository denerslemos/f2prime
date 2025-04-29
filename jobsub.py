#!/usr/bin/env python

import os.path

codehome = "/afs/cern.ch/work/d/ddesouza/UIC/SPRACE/CMSSW_13_0_5/src/f2prime"
cmssource = "/afs/cern.ch/work/d/ddesouza/UIC/SPRACE/CMSSW_13_0_5/src/"
outfolder = "/eos/user/d/ddesouza/f2prime"

os.system("mkdir -p "+str(outfolder))

''' p going direction '''

#os.system("mkdir -p "+str(outfolder)+"/HM250/pgoing/ && python3 HTCondor_submit.py -i DATA_SAMPLES/HM250/pgoing/listoffiles_pPb_DATA_HM250_pgoing -v DATA_SAMPLES/HM250/pgoing/V0_HM250_pgoing -o "+str(outfolder)+"/HM250/pgoing/f2prime_HM250_pgoing -n 3 -s f2primeHM250_pgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/HM185/pgoing/PD1/ && python3 HTCondor_submit.py -i DATA_SAMPLES/HM185/pgoing/listoffiles_pPb_DATA_HM185_pgoing_PD1 -v DATA_SAMPLES/HM185/pgoing/V0_HM185_PD1_pgoing -o "+str(outfolder)+"/HM185/pgoing/PD1/f2prime_HM185_PD1_pgoing -n 2 -s f2primeHM185_PD1_pgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/HM185/pgoing/PD2/ && python3 HTCondor_submit.py -i DATA_SAMPLES/HM185/pgoing/listoffiles_pPb_DATA_HM185_pgoing_PD2 -v DATA_SAMPLES/HM185/pgoing/V0_HM185_PD2_pgoing -o "+str(outfolder)+"/HM185/pgoing/PD2/f2prime_HM185_PD2_pgoing -n 2 -s f2primeHM185_PD2_pgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/HM185/pgoing/PD3/ && python3 HTCondor_submit.py -i DATA_SAMPLES/HM185/pgoing/listoffiles_pPb_DATA_HM185_pgoing_PD3 -v DATA_SAMPLES/HM185/pgoing/V0_HM185_PD3_pgoing -o "+str(outfolder)+"/HM185/pgoing/PD3/f2prime_HM185_PD3_pgoing -n 2 -s f2primeHM185_PD3_pgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/HM185/pgoing/PD4/ && python3 HTCondor_submit.py -i DATA_SAMPLES/HM185/pgoing/listoffiles_pPb_DATA_HM185_pgoing_PD4 -v DATA_SAMPLES/HM185/pgoing/V0_HM185_PD4_pgoing -o "+str(outfolder)+"/HM185/pgoing/PD4/f2prime_HM185_PD4_pgoing -n 2 -s f2primeHM185_PD4_pgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/HM185/pgoing/PD5/ && python3 HTCondor_submit.py -i DATA_SAMPLES/HM185/pgoing/listoffiles_pPb_DATA_HM185_pgoing_PD5 -v DATA_SAMPLES/HM185/pgoing/V0_HM185_PD5_pgoing -o "+str(outfolder)+"/HM185/pgoing/PD5/f2prime_HM185_PD5_pgoing -n 2 -s f2primeHM185_PD5_pgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/HM185/pgoing/PD6/ && python3 HTCondor_submit.py -i DATA_SAMPLES/HM185/pgoing/listoffiles_pPb_DATA_HM185_pgoing_PD6 -v DATA_SAMPLES/HM185/pgoing/V0_HM185_PD6_pgoing -o "+str(outfolder)+"/HM185/pgoing/PD6/f2prime_HM185_PD6_pgoing -n 2 -s f2primeHM185_PD6_pgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/MB/pgoing/PD1/ && python3 HTCondor_submit.py -i DATA_SAMPLES/MB/pgoing/listoffiles_pPb_DATA_MB_pgoing_PD1 -v DATA_SAMPLES/MB/pgoing/V0_MB_PD1_pgoing -o "+str(outfolder)+"/MB/pgoing/PD1/f2prime_MB_PD1_pgoing -n 0 -s f2primeMB_PD1_pgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/MB/pgoing/PD2/ && python3 HTCondor_submit.py -i DATA_SAMPLES/MB/pgoing/listoffiles_pPb_DATA_MB_pgoing_PD2 -v DATA_SAMPLES/MB/pgoing/V0_MB_PD2_pgoing -o "+str(outfolder)+"/MB/pgoing/PD2/f2prime_MB_PD2_pgoing -n 0 -s f2primeMB_PD2_pgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/MB/pgoing/PD3/ && python3 HTCondor_submit.py -i DATA_SAMPLES/MB/pgoing/listoffiles_pPb_DATA_MB_pgoing_PD3 -v DATA_SAMPLES/MB/pgoing/V0_MB_PD3_pgoing -o "+str(outfolder)+"/MB/pgoing/PD3/f2prime_MB_PD3_pgoing -n 0 -s f2primeMB_PD3_pgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/MB/pgoing/PD4/ && python3 HTCondor_submit.py -i DATA_SAMPLES/MB/pgoing/listoffiles_pPb_DATA_MB_pgoing_PD4 -v DATA_SAMPLES/MB/pgoing/V0_MB_PD4_pgoing -o "+str(outfolder)+"/MB/pgoing/PD4/f2prime_MB_PD4_pgoing -n 0 -s f2primeMB_PD4_pgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/MB/pgoing/PD5/ && python3 HTCondor_submit.py -i DATA_SAMPLES/MB/pgoing/listoffiles_pPb_DATA_MB_pgoing_PD5 -v DATA_SAMPLES/MB/pgoing/V0_MB_PD5_pgoing -o "+str(outfolder)+"/MB/pgoing/PD5/f2prime_MB_PD5_pgoing -n 0 -s f2primeMB_PD5_pgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/MB/pgoing/PD6/ && python3 HTCondor_submit.py -i DATA_SAMPLES/MB/pgoing/listoffiles_pPb_DATA_MB_pgoing_PD6 -v DATA_SAMPLES/MB/pgoing/V0_MB_PD6_pgoing -o "+str(outfolder)+"/MB/pgoing/PD6/f2prime_MB_PD6_pgoing -n 0 -s f2primeMB_PD6_pgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/MB/pgoing/PD7/ && python3 HTCondor_submit.py -i DATA_SAMPLES/MB/pgoing/listoffiles_pPb_DATA_MB_pgoing_PD7 -v DATA_SAMPLES/MB/pgoing/V0_MB_PD7_pgoing -o "+str(outfolder)+"/MB/pgoing/PD7/f2prime_MB_PD7_pgoing -n 0 -s f2primeMB_PD7_pgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/MB/pgoing/PD8/ && python3 HTCondor_submit.py -i DATA_SAMPLES/MB/pgoing/listoffiles_pPb_DATA_MB_pgoing_PD8 -v DATA_SAMPLES/MB/pgoing/V0_MB_PD8_pgoing -o "+str(outfolder)+"/MB/pgoing/PD8/f2prime_MB_PD8_pgoing -n 0 -s f2primeMB_PD8_pgoing -c "+str(cmssource)+" -p "+str(codehome))

'''Pb going direction'''

#os.system("mkdir -p "+str(outfolder)+"/HM250/Pbgoing/ && python3 HTCondor_submit.py -i DATA_SAMPLES/HM250/Pbgoing/listoffiles_pPb_DATA_HM250_Pbgoing -v DATA_SAMPLES/HM250/Pbgoing/V0_HM250_Pbgoing -o "+str(outfolder)+"/HM250/Pbgoing/f2prime_HM250_Pbgoing -n 3 -s f2primeHM250_Pbgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/HM185/Pbgoing/PD1/ && python3 HTCondor_submit.py -i DATA_SAMPLES/HM185/Pbgoing/listoffiles_pPb_DATA_HM185_Pbgoing_PD1 -v DATA_SAMPLES/HM185/Pbgoing/V0_HM185_PD1_Pbgoing -o "+str(outfolder)+"/HM185/Pbgoing/PD1/f2prime_HM185_PD1_Pbgoing -n 2 -s f2primeHM185_PD1_Pbgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/HM185/Pbgoing/PD2/ && python3 HTCondor_submit.py -i DATA_SAMPLES/HM185/Pbgoing/listoffiles_pPb_DATA_HM185_Pbgoing_PD2 -v DATA_SAMPLES/HM185/Pbgoing/V0_HM185_PD2_Pbgoing -o "+str(outfolder)+"/HM185/Pbgoing/PD2/f2prime_HM185_PD2_Pbgoing -n 2 -s f2primeHM185_PD2_Pbgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/HM185/Pbgoing/PD3/ && python3 HTCondor_submit.py -i DATA_SAMPLES/HM185/Pbgoing/listoffiles_pPb_DATA_HM185_Pbgoing_PD3 -v DATA_SAMPLES/HM185/Pbgoing/V0_HM185_PD3_Pbgoing -o "+str(outfolder)+"/HM185/Pbgoing/PD3/f2prime_HM185_PD3_Pbgoing -n 2 -s f2primeHM185_PD3_Pbgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/HM185/Pbgoing/PD4/ && python3 HTCondor_submit.py -i DATA_SAMPLES/HM185/Pbgoing/listoffiles_pPb_DATA_HM185_Pbgoing_PD4 -v DATA_SAMPLES/HM185/Pbgoing/V0_HM185_PD4_Pbgoing -o "+str(outfolder)+"/HM185/Pbgoing/PD4/f2prime_HM185_PD4_Pbgoing -n 2 -s f2primeHM185_PD4_Pbgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/HM185/Pbgoing/PD5/ && python3 HTCondor_submit.py -i DATA_SAMPLES/HM185/Pbgoing/listoffiles_pPb_DATA_HM185_Pbgoing_PD5 -v DATA_SAMPLES/HM185/Pbgoing/V0_HM185_PD5_Pbgoing -o "+str(outfolder)+"/HM185/Pbgoing/PD5/f2prime_HM185_PD5_Pbgoing -n 2 -s f2primeHM185_PD5_Pbgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/HM185/Pbgoing/PD6/ && python3 HTCondor_submit.py -i DATA_SAMPLES/HM185/Pbgoing/listoffiles_pPb_DATA_HM185_Pbgoing_PD6 -v DATA_SAMPLES/HM185/Pbgoing/V0_HM185_PD6_Pbgoing -o "+str(outfolder)+"/HM185/Pbgoing/PD6/f2prime_HM185_PD6_Pbgoing -n 2 -s f2primeHM185_PD6_Pbgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/MB/Pbgoing/PD1/ && python3 HTCondor_submit.py -i DATA_SAMPLES/MB/Pbgoing/listoffiles_pPb_DATA_MB_Pbgoing_PD1 -v DATA_SAMPLES/MB/Pbgoing/V0_MB_PD1_Pbgoing -o "+str(outfolder)+"/MB/Pbgoing/PD1/f2prime_MB_PD1_Pbgoing -n 0 -s f2primeMB_PD1_Pbgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/MB/Pbgoing/PD2/ && python3 HTCondor_submit.py -i DATA_SAMPLES/MB/Pbgoing/listoffiles_pPb_DATA_MB_Pbgoing_PD2 -v DATA_SAMPLES/MB/Pbgoing/V0_MB_PD2_Pbgoing -o "+str(outfolder)+"/MB/Pbgoing/PD2/f2prime_MB_PD2_Pbgoing -n 0 -s f2primeMB_PD2_Pbgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/MB/Pbgoing/PD3/ && python3 HTCondor_submit.py -i DATA_SAMPLES/MB/Pbgoing/listoffiles_pPb_DATA_MB_Pbgoing_PD3 -v DATA_SAMPLES/MB/Pbgoing/V0_MB_PD3_Pbgoing -o "+str(outfolder)+"/MB/Pbgoing/PD3/f2prime_MB_PD3_Pbgoing -n 0 -s f2primeMB_PD3_Pbgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/MB/Pbgoing/PD4/ && python3 HTCondor_submit.py -i DATA_SAMPLES/MB/Pbgoing/listoffiles_pPb_DATA_MB_Pbgoing_PD4 -v DATA_SAMPLES/MB/Pbgoing/V0_MB_PD4_Pbgoing -o "+str(outfolder)+"/MB/Pbgoing/PD4/f2prime_MB_PD4_Pbgoing -n 0 -s f2primeMB_PD4_Pbgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/MB/Pbgoing/PD5/ && python3 HTCondor_submit.py -i DATA_SAMPLES/MB/Pbgoing/listoffiles_pPb_DATA_MB_Pbgoing_PD5 -v DATA_SAMPLES/MB/Pbgoing/V0_MB_PD5_Pbgoing -o "+str(outfolder)+"/MB/Pbgoing/PD5/f2prime_MB_PD5_Pbgoing -n 0 -s f2primeMB_PD5_Pbgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/MB/Pbgoing/PD6/ && python3 HTCondor_submit.py -i DATA_SAMPLES/MB/Pbgoing/listoffiles_pPb_DATA_MB_Pbgoing_PD6 -v DATA_SAMPLES/MB/Pbgoing/V0_MB_PD6_Pbgoing -o "+str(outfolder)+"/MB/Pbgoing/PD6/f2prime_MB_PD6_Pbgoing -n 0 -s f2primeMB_PD6_Pbgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/MB/Pbgoing/PD7/ && python3 HTCondor_submit.py -i DATA_SAMPLES/MB/Pbgoing/listoffiles_pPb_DATA_MB_Pbgoing_PD7 -v DATA_SAMPLES/MB/Pbgoing/V0_MB_PD7_Pbgoing -o "+str(outfolder)+"/MB/Pbgoing/PD7/f2prime_MB_PD7_Pbgoing -n 0 -s f2primeMB_PD7_Pbgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/MB/Pbgoing/PD8/ && python3 HTCondor_submit.py -i DATA_SAMPLES/MB/Pbgoing/listoffiles_pPb_DATA_MB_Pbgoing_PD8 -v DATA_SAMPLES/MB/Pbgoing/V0_MB_PD8_Pbgoing -o "+str(outfolder)+"/MB/Pbgoing/PD8/f2prime_MB_PD8_Pbgoing -n 0 -s f2primeMB_PD8_Pbgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/MB/Pbgoing/PD9/ && python3 HTCondor_submit.py -i DATA_SAMPLES/MB/Pbgoing/listoffiles_pPb_DATA_MB_Pbgoing_PD9 -v DATA_SAMPLES/MB/Pbgoing/V0_MB_PD9_Pbgoing -o "+str(outfolder)+"/MB/Pbgoing/PD9/f2prime_MB_PD9_Pbgoing -n 0 -s f2primeMB_PD9_Pbgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/MB/Pbgoing/PD10/ && python3 HTCondor_submit.py -i DATA_SAMPLES/MB/Pbgoing/listoffiles_pPb_DATA_MB_Pbgoing_PD10 -v DATA_SAMPLES/MB/Pbgoing/V0_MB_PD10_Pbgoing -o "+str(outfolder)+"/MB/Pbgoing/PD10/f2prime_MB_PD10_Pbgoing -n 0 -s f2primeMB_PD10_Pbgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/MB/Pbgoing/PD11/ && python3 HTCondor_submit.py -i DATA_SAMPLES/MB/Pbgoing/listoffiles_pPb_DATA_MB_Pbgoing_PD11 -v DATA_SAMPLES/MB/Pbgoing/V0_MB_PD11_Pbgoing -o "+str(outfolder)+"/MB/Pbgoing/PD11/f2prime_MB_PD11_Pbgoing -n 0 -s f2primeMB_PD11_Pbgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/MB/Pbgoing/PD12/ && python3 HTCondor_submit.py -i DATA_SAMPLES/MB/Pbgoing/listoffiles_pPb_DATA_MB_Pbgoing_PD12 -v DATA_SAMPLES/MB/Pbgoing/V0_MB_PD12_Pbgoing -o "+str(outfolder)+"/MB/Pbgoing/PD12/f2prime_MB_PD12_Pbgoing -n 0 -s f2primeMB_PD12_Pbgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/MB/Pbgoing/PD13/ && python3 HTCondor_submit.py -i DATA_SAMPLES/MB/Pbgoing/listoffiles_pPb_DATA_MB_Pbgoing_PD13 -v DATA_SAMPLES/MB/Pbgoing/V0_MB_PD13_Pbgoing -o "+str(outfolder)+"/MB/Pbgoing/PD13/f2prime_MB_PD13_Pbgoing -n 0 -s f2primeMB_PD13_Pbgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/MB/Pbgoing/PD14/ && python3 HTCondor_submit.py -i DATA_SAMPLES/MB/Pbgoing/listoffiles_pPb_DATA_MB_Pbgoing_PD14 -v DATA_SAMPLES/MB/Pbgoing/V0_MB_PD14_Pbgoing -o "+str(outfolder)+"/MB/Pbgoing/PD14/f2prime_MB_PD14_Pbgoing -n 0 -s f2primeMB_PD14_Pbgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/MB/Pbgoing/PD15/ && python3 HTCondor_submit.py -i DATA_SAMPLES/MB/Pbgoing/listoffiles_pPb_DATA_MB_Pbgoing_PD15 -v DATA_SAMPLES/MB/Pbgoing/V0_MB_PD15_Pbgoing -o "+str(outfolder)+"/MB/Pbgoing/PD15/f2prime_MB_PD15_Pbgoing -n 0 -s f2primeMB_PD15_Pbgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/MB/Pbgoing/PD16/ && python3 HTCondor_submit.py -i DATA_SAMPLES/MB/Pbgoing/listoffiles_pPb_DATA_MB_Pbgoing_PD16 -v DATA_SAMPLES/MB/Pbgoing/V0_MB_PD16_Pbgoing -o "+str(outfolder)+"/MB/Pbgoing/PD16/f2prime_MB_PD16_Pbgoing -n 0 -s f2primeMB_PD16_Pbgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/MB/Pbgoing/PD17/ && python3 HTCondor_submit.py -i DATA_SAMPLES/MB/Pbgoing/listoffiles_pPb_DATA_MB_Pbgoing_PD17 -v DATA_SAMPLES/MB/Pbgoing/V0_MB_PD17_Pbgoing -o "+str(outfolder)+"/MB/Pbgoing/PD17/f2prime_MB_PD17_Pbgoing -n 0 -s f2primeMB_PD17_Pbgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/MB/Pbgoing/PD18/ && python3 HTCondor_submit.py -i DATA_SAMPLES/MB/Pbgoing/listoffiles_pPb_DATA_MB_Pbgoing_PD18 -v DATA_SAMPLES/MB/Pbgoing/V0_MB_PD18_Pbgoing -o "+str(outfolder)+"/MB/Pbgoing/PD18/f2prime_MB_PD18_Pbgoing -n 0 -s f2primeMB_PD18_Pbgoing -c "+str(cmssource)+" -p "+str(codehome))

#os.system("mkdir -p "+str(outfolder)+"/MB/Pbgoing/PD19/ && python3 HTCondor_submit.py -i DATA_SAMPLES/MB/Pbgoing/listoffiles_pPb_DATA_MB_Pbgoing_PD19 -v DATA_SAMPLES/MB/Pbgoing/V0_MB_PD19_Pbgoing -o "+str(outfolder)+"/MB/Pbgoing/PD19/f2prime_MB_PD19_Pbgoing -n 0 -s f2primeMB_PD19_Pbgoing -c "+str(cmssource)+" -p "+str(codehome))

os.system("mkdir -p "+str(outfolder)+"/MB/Pbgoing/PD20/ && python3 HTCondor_submit.py -i DATA_SAMPLES/MB/Pbgoing/listoffiles_pPb_DATA_MB_Pbgoing_PD20 -v DATA_SAMPLES/MB/Pbgoing/V0_MB_PD20_Pbgoing -o "+str(outfolder)+"/MB/Pbgoing/PD20/f2prime_MB_PD20_Pbgoing -n 0 -s f2primeMB_PD20_Pbgoing -c "+str(cmssource)+" -p "+str(codehome))
