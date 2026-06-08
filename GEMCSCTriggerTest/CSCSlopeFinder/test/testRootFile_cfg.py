import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.Eras.Era_Run3_cff import Run3
from Configuration.Eras.Era_Run2_2018_cff import Run2_2018
options = VarParsing('analysis')
options.register("unpack", True, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Set to True when you want to unpack the CSC DAQ data.")
options.register("selectCSCs", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Set to True when you want to (un)select certain CSCs.")
options.register("maskedChambers", "", VarParsing.multiplicity.list, VarParsing.varType.string,
                 "Chambers you want to explicitly mask.")
options.register("selectedChambers", "", VarParsing.multiplicity.list, VarParsing.varType.string,
                 "Chambers you want to explicitly mask.")
options.register("unpackGEM", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Set to True when you want to unpack the GEM DAQ data.")
options.register("l1", True, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Set to True when you want to re-emulate the CSC trigger primitives.")
options.register("l1GEM", True, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Set to True when you want to re-emulate the GEM trigger primitives.")
options.register("mc", True, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Set to True when running on MC.")
options.register("dqm", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Set to True when you want to run the CSC DQM")
options.register("dqmGEM", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Set to True when you want to run the GEM DQM")
options.register("useEmtfGEM", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Set to True when you want to use GEM clusters from the EMTF in the DQM")
options.register("useB904ME11", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Set to True when using B904 ME1/1 data.")
options.register("useB904ME21", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Set to True when using B904 ME2/1 data (also works for ME3/1 and ME4/1).")
options.register("useB904ME234s2", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Set to True when using B904 ME4/2 data (also works for MEX/2 and ME1/3).")
options.register('useB904ME11PositiveEndcap',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
                 "Set to True when using data from ME1/1 set as Positive Endcap chamber in B904.")
options.register('useB904ME11NegativeEndcap',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
                 "Set to True when using data from ME1/1 set as Negative Endcap chamber in B904.")
options.register('useB904GE11Short',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
                 "Set to True when using data from GE1/1 Short super chamber in B904.")
options.register('useB904GE11Long',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
                 "Set to True when using data from GE1/1 Long super chamber in B904.")
options.register("run3", True, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Set to True when using Run-3 data.")
options.register("runCCLUTOTMB", True, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Set to True when using the CCLUT OTMB algorithm.")
options.register("runCCLUTTMB", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Set to True when using the CCLUT TMB algorithm.")
options.register("runME11ILT", True, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Set to True when running the GEM-CSC integrated local trigger algorithm in ME1/1.")
options.register("runME21ILT", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Set to True when running the GEM-CSC integrated local trigger algorithm in ME2/1.")
options.register("saveEdmOutput", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Set to True if you want to keep the EDM ROOT after unpacking and re-emulating.")
options.register("preTriggerAnalysis", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Set to True if you want to print out more details about CLCTs and LCTs in the offline CSC DQM module.")
options.register("dropNonMuonCollections", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Option to drop most non-muon collections generally considered unnecessary for GEM/CSC analysis")
options.register("dqmOutputFile", "step_DQM.root", VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "Name of the DQM output file. Default: step_DQM.root")
options.parseArguments()

process_era = Run3
if not options.run3:
      process_era = Run2_2018

process = cms.Process("L1CSCTPG", process_era)
process.load("Configuration/StandardSequences/GeometryRecoDB_cff")
process.load("Configuration/StandardSequences/MagneticField_cff")
process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.DQMSaverAtRunEnd_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load("EventFilter.CSCRawToDigi.cscUnpacker_cfi")
process.load('EventFilter.GEMRawToDigi.muonGEMDigis_cfi')
process.load('EventFilter.L1TRawToDigi.emtfStage2Digis_cfi')
process.load("L1Trigger.CSCTriggerPrimitives.cscTriggerPrimitiveDigis_cfi")
process.load("CalibMuon.CSCCalibration.CSCL1TPLookupTableEP_cff")
process.load('L1Trigger.L1TGEM.simGEMDigis_cff')
process.load("DQM.L1TMonitor.L1TdeCSCTPG_cfi")
process.load("DQM.L1TMonitor.L1TdeGEMTPG_cfi")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoMuon.TrackingTools.MuonServiceProxy_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')

process.maxEvents = cms.untracked.PSet(
      input = cms.untracked.int32(options.maxEvents)
)

# Set single-threaded mode correctly
process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(1),
    wantSummary = cms.untracked.bool(True),
    TryToContinue = cms.untracked.vstring('ProductNotFound')
)

# Source: your ROOT file
process.source = cms.Source(
      "PoolSource",
      fileNames = cms.untracked.vstring(options.inputFiles),
      inputCommands = cms.untracked.vstring(
            'keep *',
            'drop CSCDetIdCSCShowerDigiMuonDigiCollection_simCscTriggerPrimitiveDigis_*_*',
            'drop LCTDebugobjects_*_*_*',
            'drop floatBXVector_gtStage2Digis_CICADAScore_RECO' #This was causing a fail for some reason in 2024F/G? 
      )
)

if options.unpackGEM:
      process.source.labelRawDataLikeMC = cms.untracked.bool(False)

# Dummy output (optional)
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("testOutput.root")
)

# Keep the process alive with an empty path
process.p = cms.Path()  # No modules to run

# Add output step if you want to write something (optional)
process.end = cms.EndPath(process.out)
