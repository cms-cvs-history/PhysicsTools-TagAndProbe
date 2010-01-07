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

process.tagProbes = cms.EDProducer("TagProbeProducer",
    TagCollection   = cms.InputTag("tagMuons"),
    ProbeCollection = cms.InputTag("trkProbes"),
    PassingProbeCollection = cms.InputTag("trkPassingGlb"),
    MassMinCut = cms.untracked.double(2.8),
    MassMaxCut = cms.untracked.double(3.5),
)

process.muMcMatch = cms.EDFilter("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(13),
    src = cms.InputTag("muons"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles")
)

process.TPEdm = cms.EDProducer("TagProbeEDMNtuple",
    tagProbeType = cms.untracked.string( "Muon" ),
    genParticlesTag = cms.untracked.InputTag("genParticles"),
    mcParticles = cms.untracked.vint32(13),
    mcParents   = cms.untracked.vint32( 0),
    
    ## Tag & Probe Muon Candidate Collections
    ## all this will be filled once for each T&P measurement
    tagCandTags = cms.untracked.VInputTag(cms.InputTag("tagMuons")),
    passProbeCandTags = cms.untracked.VInputTag(cms.InputTag("trkPassingGlb")),
    allProbeCandTags  = cms.untracked.VInputTag(cms.InputTag("trkProbes")),
    tagTruthMatchMapTags       = cms.untracked.VInputTag(cms.InputTag("muMcMatch")),
    passProbeTruthMatchMapTags = cms.untracked.VInputTag(cms.InputTag("muMcMatch")),
    allProbeTruthMatchMapTags  = cms.untracked.VInputTag(cms.InputTag("muMcMatch")),
    tagProbeMapTags = cms.untracked.VInputTag(cms.InputTag("tagProbes")),
    BestProbeCriteria = cms.untracked.vstring("OneProbe"),
    BestProbeInvMass  = cms.untracked.vdouble(3.09),

    ## Ignore the following, we don't care
    triggerEventTag = cms.untracked.InputTag( "hltTriggerSummaryAOD" ),
    hltL1Tag        = cms.untracked.InputTag( "hltSingleMuIsoLevel1Seed" ),
    hltTag          = cms.untracked.InputTag( "hltSingleMuIsoL3IsoFiltered" ),
    triggerDelRMatch     = cms.untracked.double( 0.15 ),
    triggerDelPtRelMatch = cms.untracked.double( 0.15 )
)

process.MakeHisto = cms.EDAnalyzer("TagProbeEDMAnalysis",      
      # AFAIK these are the only meaningful parameters for now that we just make the Histo
      # --------------------------------------------
      Mode = cms.untracked.string("Write"),
      FitFileName = cms.untracked.string( "histo_output.root" ),

      TagProbeType = cms.untracked.int32(0),

      MCTruthParentId = cms.untracked.int32(443),
      Weight = cms.untracked.double(1.0),

      CalculateEffSideBand = cms.untracked.bool( True ), ## Calculate and store effs using SB
      CalculateEffFitter   = cms.untracked.bool( True ), ## Calculate and store effs from Roofit
      CalculateEffTruth    = cms.untracked.bool( True ), ## Calculate and store true effs
      UnbinnedFit          = cms.untracked.bool( True ),
      Do2DFit              = cms.untracked.bool( True ),

      NumBinsMass         = cms.untracked.int32( 20 ),
      MassLow             = cms.untracked.double( 2.8 ),
      MassHigh            = cms.untracked.double( 3.5 ),

      NameVar1             = cms.untracked.string( "pt" ),
      Var1BinBoundaries   = cms.untracked.vdouble( 3, 4.5, 6, 10, 20 ),
      NameVar2             = cms.untracked.string( "eta" ),
      Var2BinBoundaries   = cms.untracked.vdouble( -2.4,-1.3,-0.8,0.8,1.3,2.4),

      # All the following is useless now
      # --------------------------------------------
      GaussLineShape = cms.untracked.PSet(
        GaussMean        = cms.untracked.vdouble( 3.09, 2.9,  3.1 ),
        GaussSigma       = cms.untracked.vdouble( 0.03, 0.01, 0.05 )
      ),

      CMSBkgLineShape = cms.untracked.PSet(
        CMSBkgAlpha           = cms.untracked.vdouble( 0 ),
        CMSBkgBeta            = cms.untracked.vdouble( 0 ),
        CMSBkgPeak            = cms.untracked.vdouble( 3.09 ),
        CMSBkgGamma           = cms.untracked.vdouble( 0, -10, 10 )
      ),

      Efficiency        = cms.untracked.vdouble( 0.99,0.5,1.0 ),    
      NumSignal         = cms.untracked.vdouble( 4000.0,-10.0,30000.0 ),    
      NumBkgPass        = cms.untracked.vdouble( 4000.0,0.0,  10000.0 ),    
      NumBkgFail        = cms.untracked.vdouble( 1000.0,-10.0, 7000.0 ),    

      SBSPeak            = cms.untracked.double( 3.1 ),   ## Mass peak
      SBSStanDev         = cms.untracked.double( 2 )        ## SD from peak for subtraction
)

process.tagAndProbe = cms.Path( 
    process.allTagsAndProbes *
    process.tagProbes * 
    process.muMcMatch *
    process.TPEdm * 
    process.MakeHisto
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("testWrite.root"))
