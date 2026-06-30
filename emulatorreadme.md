Install:

```
cmsrel CMSSW_16_1_0
cd CMSSW_16_1_0/src
cmsenv
git cms-merge-topic pkflana:new_data_format
scram b -j8
```

Before you run the emulator at all, you must decide which version of alignment corrections you want to use. The corrections used will be those in the folder:\
L1Trigger/CSCTriggerPrimitives/data/GEMCSC/AlignmentCorrection\
Currently, I have Konstantin's 0 placeholders in the folder AlignmentCorrection_0 and the current best alignment in AlignmentCorrection_non0, if you wish to use one of these, first rename
the desired folder to AlignmentCorrection. Any other corrections you make/want to test need to be put into a folder with that name.

There's currently two versions of running the emulator, running the "bare" emulator, which here means the emulator in the official cmssw + the reader + the emulated LCT reports GEM hit roll
information. This requires the file you run over to have the RAW datatier. An example command (assuming you're running from src) is:
```
cmsRun L1Trigger/CSCTriggerPrimitives/test/runCSCTriggerPrimitiveProducer_cfg.py mc=True inputFiles="file:/eos/user/p/pflanaga/singlemu2_50_eta_restricted/1/step3_1.root"
```
For the output of this, slopeex is the new 6-bit slope

The second version is running the version of the emulator that matches LCTs to emtftracks and then to reco muons, helping to single out high quality muons. This version requires both the RAW
and RECO datatiers. An example command is:
```
cmsRun GEMCSCTriggerTest/CSCSlopeFinder/test/run_unpack_CSCGEM_matcher.py mc=True inputFiles="file:/eos/user/p/pflanaga/singlemu2_50_eta_restricted/1/step3_1.root"
```
This version can either match data LCTs to data reco muons and also match data LCTs to emulated LCTs, or it can match emulated LCTs directly to data reco muons. This is controlled with 
the option "useEmulatorLCTs" where True means matching emulated directly to data. This allows you to test the difference made by making changes to how LCTs are emulated more directly.

In both versions, you must make sure you have the correct global tag and set the correct option for mc/emulator. Most of the other options are fairly self-explanatory and I have barely had
any cause to use them.