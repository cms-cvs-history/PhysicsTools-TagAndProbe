import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("PhysicsTools.TagAndProbe.Electron_TagProbeEDMAnalysis_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:demo.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.p = cms.Path(process.demo)
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.demo.TagProbeType = 1
process.demo.FitFileName = 'test_electroneff_GsfToIso.root'


