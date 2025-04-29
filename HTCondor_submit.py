#!/usr/bin/env python

'''Add important stuff'''
import os
import os.path
import optparse
import subprocess

''' Inputs for the skim code '''
usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-i', '--infiles', dest='infiles', help='input list of files', default='', type='string')
parser.add_option('-v', '--infilesv0', dest='infilesv0', help='V0 input list of files', default='', type='string')
parser.add_option('-o', '--outfiles', dest='outfiles', help='output files', default='', type='string')
parser.add_option('-n', '--multiplicity', dest='multiplicity', help='0 for no cut, 1 for MB [0,185], 2 for HM 1 to 6 [185, 250] and 3 for HM 7 [250, inf]', default='0', type='int')
parser.add_option('-s', '--subfiles', dest='subfiles', help='HTCondor submission file', default='', type='string')
parser.add_option('-c', '--cmssw', dest='cmsswdir', help='CMSSW directory', default='$CMSSW_BASE/src', type='string')
parser.add_option('-p', '--pwd', dest='pwddir', help='Local directory', default='$(pwd)', type='string')

(opt, args) = parser.parse_args()
inFiles = opt.infiles
inFilesV0 = opt.infilesv0
outFiles = opt.outfiles
mult = opt.multiplicity
subFiles = opt.subfiles
cmsswDir = opt.cmsswdir
pwdDir = opt.pwddir


''' Read list of files '''
listOfFiles = open(inFiles+'.txt', 'r')
Lines = listOfFiles.readlines()
print ("Number of files/jobs: "+str(len(Lines)))

''' Make .sh file automatically '''
script_content = """#!/bin/bash

echo "Setup CMSSW (ROOT version)"
cd """+cmsswDir+"""
eval `scramv1 runtime -sh`
cd """+pwdDir+"""
mkdir -p cond
echo "Submit skim jobs at "
echo PWD: $PWD

./K0Star $1 $2 $3 $4
"""

# Save to a shell script file
with open("sub_skim.sh", "w") as f:
    f.write(script_content)

os.chmod("sub_skim.sh", 0o755)

''' Start the write submission file '''
fsubfile = open(subFiles+".sub", "w")
command_lines = '''universe   = vanilla
getenv     = True
executable = sub_skim.sh
+JobFlavour           = "tomorrow"
requirements = ((OpSysAndVer =?= "AlmaLinux9") && (CERNEnvironment =?= "qa"))
RequestCpus = 2
transfer_input_files  = voms_proxy.txt
environment = "X509_USER_PROXY=voms_proxy.txt"
'''

''' Loop over files '''
i=0
for line in Lines:
    outtempfiles = open(inFiles+"_part"+str(i)+".txt", "w")
    outtempfiles.write(line)
    outtempfiles.close()
    temp = '''
log        = cond/'''+subFiles+'''_part_'''+str(i)+'''.log
output     = cond/'''+subFiles+'''_part_'''+str(i)+'''.out
error      = cond/'''+subFiles+'''_part_'''+str(i)+'''.err
arguments = '''+inFiles+'''_part'''+str(i)+'''.txt '''+inFilesV0+'''.txt '''+outFiles+'''_'''+str(i)+'''.root '''+str(mult)+'''
queue
'''
    command_lines += temp
    i=i+1

fsubfile.write(command_lines)
fsubfile.close()
subprocess.call(["condor_submit", subFiles+".sub"])
