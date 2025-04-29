# F2Prime(1525) -> K0s + K0s skim code using HTCondor

From Dener De Souza Lemos

## Intructions

Setup CMSSW (just for root versioning)
```
export SCRAM_ARCH=el9_amd64_gcc12
cmsrel CMSSW_13_0_5
cd CMSSW_13_0_5/src
cmsenv
```
Inside of the src folder, download the code using
```
git clone git@github.com:denerslemos/f2prime.git
cd f2prime
mkdir cond
```
Once this steps are done you can compile the code with
```
g++ -O2 f2prime.C `root-config --libs` `root-config --cflags` -o f2prime
```
This will create the executable: ```f2prime``` 

You will need your VOMS certificate, do it using
```
voms-proxy-init -rfc -voms cms --out voms_proxy.txt --hours 200
```
that creates a certificate file valid for 200 hours: voms_proxy.txt

After that, you will submit jobs by using [HTCondor_submit.py](https://github.com/denerslemos/f2prime/blob/main/HTCondor_submit.py):
```
python3 HTCondor_submit.py -i input_text_file -v v0_input_file -o output_name_file -n X -s Y -c Z -p Q
```
 - input_text_file: is the text file (use it without the .txt extension) with inputs and can be found in the folder DATA_SAMPLES each .root input will be a job
 - v0_input_file: is the K0s text file (use it without the .txt extension) with inputs that are also at DATA_SAMPLES. For each PD the code will match all V0 .root files with one input_text_file
 - output_name_file: output file name (use it without the .root extension), it will automatically include a counter for each input. You can use paths to save on EOS.
 - X: 0 or 1 for no multiplicity cut for MB sample, 2 for HM185 [185,250] and 3 for HM250 [250,inf]
 - Y: name for the submission files, I have used HTcondor_sub_ + some information from the sample, PD, MB, ... + pgoing or Pbgoing.
 - Z: CMSSW folder you give ```cmsenv```
 - Q: PWD of the folder you are at (or where you wanna run your code)
It will automatically include a counter for each input

You can also edit the ```jobsub.py``` to select what you wanna submit. For that, you must edit the first few lines for your own repositories and uncomment/comment the jobs you wanna submit. Please DO NOT submit all jobs at the same time.
