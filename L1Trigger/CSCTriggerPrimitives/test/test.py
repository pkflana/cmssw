import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
        "root://cms-xrd-global.cern.ch//store/data/Run2026A/Muon0/RAW-RECO/ZMu-PromptReco-v1/000/401/699/00000/47f8adfa-a0ca-4bcc-8692-7be0aea9f23d.root"
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

process.test = cms.EDAnalyzer("EventContentAnalyzer")

process.p = cms.Path(process.test)
