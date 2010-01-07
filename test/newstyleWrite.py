import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_3_4_1/RelValJpsiMM/GEN-SIM-RECO/STARTUP3X_V14-v1/0004/DC4D9E50-82ED-DE11-BF03-003048D375AA.root',
        '/store/relval/CMSSW_3_4_1/RelValJpsiMM/GEN-SIM-RECO/STARTUP3X_V14-v1/0004/88F9B6DC-85ED-DE11-A295-001D09F2514F.root',
        '/store/relval/CMSSW_3_4_1/RelValJpsiMM/GEN-SIM-RECO/STARTUP3X_V14-v1/0004/66B9D232-85ED-DE11-B5EC-0019B9F7312C.root',
        '/store/relval/CMSSW_3_4_1/RelValJpsiMM/GEN-SIM-RECO/STARTUP3X_V14-v1/0004/509FE9BB-83ED-DE11-A25F-001D09F253C0.root',
        '/store/relval/CMSSW_3_4_1/RelValJpsiMM/GEN-SIM-RECO/STARTUP3X_V14-v1/0004/4294CBAC-B5ED-DE11-B7E8-003048D2BE08.root',
    )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )    

process.tagMuons = cms.EDFilter("MuonRefSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("isGlobalMuon && pt > 3 && abs(eta) < 2.4"), 
)

process.trkProbes = cms.EDFilter("MuonRefSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("isTrackerMuon && pt > 3 && abs(eta) < 2.4"), 
)
process.trkPassingGlb = cms.EDFilter("MuonRefSelector",
    src = cms.InputTag("muons"),
    cut = cms.string(process.trkProbes.cut.value() + " && isGlobalMuon"), 
)


process.allTagsAndProbes = cms.Sequence(
    process.tagMuons +
    process.trkProbes +
    process.trkPassingGlb
)

process.tagProbes = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagMuons@+ trkProbes@-"), # charge coniugate states are implied
    cut   = cms.string("2.5 < mass < 3.8"),
)

process.muMcMatch = cms.EDFilter("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(13),
    src = cms.InputTag("muons"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles")
)

process.MakeHisto = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tagProbes"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
        eta = cms.string("eta()"),
        pt  = cms.string("pt()"),
    ),
    flags = cms.PSet(
        passing = cms.InputTag("trkPassingGlb"),
    ),
    isMC = cms.bool(True),
    tagMatches = cms.InputTag("muMcMatch"),
    probeMatches  = cms.InputTag("muMcMatch"),
    motherPdgId = cms.int32(443),
    makeMCUnbiasTree = cms.bool(False),
    #checkMotherInUnbiasEff = cms.bool(True),
    #allProbes     = cms.InputTag("trkProbes"),
)

process.tagAndProbe = cms.Path( 
    process.allTagsAndProbes *
    process.tagProbes * 
    process.muMcMatch *
    process.MakeHisto
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("testNewWrite.root"))
