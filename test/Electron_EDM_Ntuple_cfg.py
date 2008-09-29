import FWCore.ParameterSet.Config as cms

process = cms.Process("ElectronEff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")

#process.load("Configuration.StandardSequences.MagneticField_cff")
#process.GlobalTag.globaltag = cms.string('IDEAL_V6::All')

process.load("CondCore.DBCommon.CondDBSetup_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.load("PhysicsTools.TagAndProbe.tag_probe_electron_cfi")


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 20



process.TPEdm = cms.EDFilter("TagProbeEDMNtuple",
    allProbeTruthMatchMapTags = cms.untracked.VInputTag(cms.InputTag("SuperClustersMatch"), cms.InputTag("GsfElectronsMatch"), cms.InputTag("IsolationMatch"), cms.InputTag("IdMatch")),
    checkExactOverlap = cms.untracked.bool(False),
    triggerDelRMatch = cms.untracked.double(0.3),
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
    mcParents = cms.untracked.vint32(0, 0, 0)
)

process.outpath = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_TPEdm_*_*'),
    fileName = cms.untracked.string('test_EDM_ntuple.root')
)

process.p1 = cms.Path(process.lepton_cands+process.TPEdm)
process.the_end = cms.EndPath(process.outpath)
process.PoolSource.fileNames = ['/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/04419036-F385-DD11-B3A7-001617C3B6E8.root',
        '/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/0A28F869-F285-DD11-AF3C-001617DBD5B2.root',
        '/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/162C4B5E-F585-DD11-872A-001617C3B64C.root',
        '/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/205E6CE3-F485-DD11-9D53-001617C3B76A.root',
        '/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/562BAFA1-F585-DD11-B931-001617DBD224.root',
        '/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/565AFE10-EF85-DD11-8353-000423D6B42C.root',
        '/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/5C66302A-F185-DD11-81D3-000423D98834.root',
        '/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/66F60641-F685-DD11-A493-000423D987FC.root',
        '/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/6E6A6E2D-F485-DD11-B707-001617DBD472.root',
        '/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/70191C8D-F485-DD11-8280-001617E30D06.root',
        '/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/9A204B65-F385-DD11-9CF1-000423D98B6C.root',
        '/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/A6E8BAB0-F085-DD11-9AB1-000423D986C4.root',
        '/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/B4D8FA14-F485-DD11-A41D-001617C3B76A.root',
        '/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/B6CA9FAB-F185-DD11-B66B-001617E30D0A.root',
        '/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/D8DAE3BF-F385-DD11-9C8E-001617C3B65A.root',
        '/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/EA02E08F-F285-DD11-8AF3-000423D9870C.root',
        '/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0001/6E9B44E2-0487-DD11-BFA7-001617C3B78C.root']

