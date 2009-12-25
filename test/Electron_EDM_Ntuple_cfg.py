import FWCore.ParameterSet.Config as cms

process = cms.Process("ElectronEff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("CondCore.DBCommon.CondDBSetup_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = cms.string('IDEAL_V11::All')
#process.GlobalTag.globaltag = cms.string('CRAFT_30X::All')
#process.GlobalTag.globaltag = "STARTUP_V7::All"
process.GlobalTag.globaltag = cms.string('MC_31X_V3::All')
process.load("PhysicsTools.TagAndProbe.tag_probe_electron_cfi")


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/mc/Summer09/Zee/GEN-SIM-RECO/MC_31X_V3_SD_Ele15-v1/0003/C62D4337-F1AB-DE11-A96B-0018F3D0970E.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000




process.TPEdm = cms.EDProducer("TagProbeEDMNtuple",
    isMC = cms.untracked.bool(False),
    allProbeTruthMatchMapTags = cms.untracked.VInputTag(cms.InputTag("SuperClustersMatch"), cms.InputTag("GsfElectronsMatch"), cms.InputTag("IsolationMatch"), cms.InputTag("IdMatch")),
    checkExactOverlap = cms.untracked.bool(False),
    triggerDelRMatch = cms.untracked.double(0.3),
    triggerDelPtRelMatch = cms.untracked.double(0.3),
    allProbeCandTags = cms.untracked.VInputTag(cms.InputTag("theSuperClusters"), cms.InputTag("theGsfElectrons"), cms.InputTag("theIsolation"), cms.InputTag("theId")),
    # Truth Matching tags
    passProbeTruthMatchMapTags = cms.untracked.VInputTag(cms.InputTag("GsfElectronsMatch"), cms.InputTag("IsolationMatch"), cms.InputTag("IdMatch"), cms.InputTag("HLTMatch")),
    # Tag & Probe Electron Candidate Collections
    tagCandTags = cms.untracked.VInputTag(cms.InputTag("theHLT"), cms.InputTag("theHLT"), cms.InputTag("theHLT"), cms.InputTag("theHLT")),
    # Tag & Probe Muon Association Map 
    tagProbeMapTags = cms.untracked.VInputTag(cms.InputTag("tpMapSuperClusters"), cms.InputTag("tpMapGsfElectrons"), cms.InputTag("tpMapIsolation"), cms.InputTag("tpMapId")),
    # Type of tag-probe candidates, use "Muon" or "Electron"
    # For the moment this only affects the kind of particle
    # used for storing MC truth information.
    tagProbeType = cms.untracked.string('Electron'),
    # Truth Map Tags
    tagTruthMatchMapTags = cms.untracked.VInputTag(cms.InputTag("HLTMatch"), cms.InputTag("HLTMatch"), cms.InputTag("HLTMatch"), cms.InputTag("HLTMatch")),
    # Store some generic information about the event
    # in case we want it
    mcParticles = cms.untracked.vint32(23, 11, 22),
    trackTags = cms.untracked.VInputTag(cms.InputTag("generalTracks")),
    # Pass Probe Electron Candidate Collections
    passProbeCandTags = cms.untracked.VInputTag(cms.InputTag("theGsfElectrons"), cms.InputTag("theIsolation"), cms.InputTag("theId"), cms.InputTag("theHLT")),
    verticesTag = cms.untracked.InputTag("offlinePrimaryVertices"),
    mcParents = cms.untracked.vint32(0, 0, 0),
    BestProbeCriteria = cms.untracked.vstring("HighestProbePt","HighestProbePt","HighestProbePt","HighestProbePt")
)


# Only keep events where a Tag-Probe pair was found (useful for background)
process.TPFilter = cms.EDFilter("TagProbeEDMFilter")

process.p1 = cms.Path( process.lepton_cands*process.TPEdm*process.TPFilter )

process.outpath = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("demo.root"),
    outputCommands = cms.untracked.vstring(
	  "drop *",
	  "keep *_TPEdm_*_*" 
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p1')
    )
)
    
process.the_end = cms.EndPath( process.outpath )
