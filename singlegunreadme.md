The scripts I used to create the muon gun samples are included. I ran this in CMSSW_15_1_0, you'll have to edit all eos paths in both files. I two arguments in the bash file are folder number\
(for parallelization) and config used. The configs are located in Configuration/Generator/python, you can edit configs there or make new ones

```
#!/bin/bash
cd /eos/user/p/pflanaga/cmssw_clean/CMSSW_15_1_0/src/
cmsenv
cd /eos/user/p/pflanaga/singlemu2_50_eta_restricted
mkdir $1
cd $1

cmsDriver.py $2  -s GEN,SIM -n 500 --conditions 150X_mcRun3_2025_realistic_v2 \
--beamspot DBrealistic --datatier GEN-SIM --eventcontent FEVTDEBUG --geometry DB:Extended --era Run3_2025 \
--no_exec --fileout file:step1_$1.root --python_filename step1_cfg_$1.py --customise_commands process.RandomNumberGeneratorService.generator.initialSeed="int($1)"

cmsRun step1_cfg_$1.py

cmsDriver.py step2  -s DIGI:pdigi_valid,L1,DIGI2RAW,HLT:@relval2022 --conditions 150X_mcRun3_2025_realistic_v2 \
--datatier GEN-SIM-DIGI-RAW -n 500 --eventcontent FEVTDEBUGHLT --geometry DB:Extended --era Run3_2025 \
--no_exec --filein  file:step1_$1.root  --fileout file:step2_$1.root --python_filename step2_cfg_$1.py --customise_commands process.RandomNumberGeneratorService.generator.initialSeed="int($1)"

cmsRun step2_cfg_$1.py

rm step1_$1.root

cmsDriver.py step3  -s RAW2DIGI,L1Reco,RECO,RECOSIM,PAT,NANO,VALIDATION:@standardValidationNoHLT+@miniAODValidation,DQM:@standardDQMFakeHLT+@miniAODDQM+@nanoAODDQM \
--conditions 150X_mcRun3_2025_realistic_v2 --datatier GEN-SIM-RECO,MINIAODSIM,NANOAODSIM,DQMIO -n 500 \
--eventcontent FEVTDEBUGHLT,MINIAODSIM,NANOEDMAODSIM,DQM --geometry DB:Extended --era Run3_2025 --no_exec \
--filein  file:step2_$1.root  --fileout file:step3_$1.root --python_filename step3_cfg_$1.py --customise_commands process.RandomNumberGeneratorService.generator.initialSeed="int($1)"

cmsRun step3_cfg_$1.py

rm step2_$1.root
rm step3_$1_in*.root
rm *.py
```

and for submitting to condor:

```
import htcondor
import os
schedd = htcondor.Schedd()
col = htcondor.Collector()
credd = htcondor.Credd()
credd.add_user_cred(htcondor.CredTypes.Kerberos, None)

time = "workday"
configs = ["SingleMuPt2_50_Eta1p5_2p5_cfi","SingleMuPt2_50_Eta-1p5_-2p5_cfi"]
before = 0
after = 0
arguments = []
for i in range(1,1001):
        si = str(i)
        if os.path.exists("/eos/user/p/pflanaga/singlemu2_50_eta_restricted/"+si):
            continue
        if i<501:
            config = configs[1]
            before += 1
        else:
            config = configs[0]
            after += 1
        arguments.append({"arguments": si+' '+config})
    
job = htcondor.Submit({
    "executable": "simulation_reco.sh",
    "arguments": "$(arguments)",
    "output": "output/$(ClusterId).$(ProcId).out",
    "error": "error/$(ClusterId).$(ProcId).err",
    "log": "log/$(ClusterId).$(ProcId).log",
    "universe": "vanilla",
    "Requirements": '(OpSysAndVer =?= "AlmaLinux9")',
    "+MaxRuntime": "28800",
    "max_retries": "1"
})
job['MY.SendCredential'] = "true"
print(before,after)
submit_result = schedd.submit(job,itemdata = iter(arguments))
```

I have not yet looked into saving less to save space, I will do this later and update.

