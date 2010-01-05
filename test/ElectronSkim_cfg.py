import FWCore.ParameterSet.Config as cms
process = cms.Process("SKIM")
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load("CondCore.DBCommon.CondDBSetup_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")


# configure HLT
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0 AND((40 OR 41) AND NOT (36 OR 37 OR 38 OR 39))')

# require physics declared
process.physDecl = cms.EDFilter("PhysDecl",
                                applyfilter = cms.untracked.bool(True)
                                )

# require scraping filter
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
                                    applyfilter = cms.untracked.bool(True),
                                    debugOn = cms.untracked.bool(False),
                                    numtrack = cms.untracked.uint32(10),
                                    thresh = cms.untracked.double(0.2)
                                    )


# require primary vertex
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
    vertexCollection = cms.InputTag('offlinePrimaryVertices'),
    minimumNumberOfTracks = cms.uint32(1) ,
    maxAbsZ = cms.double(20), 
    maxd0 = cms.double(2) 
)




#############   Set the number of events #############
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5000)
)


#############   Define the source file ###############

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
     '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/F80561F5-A2E6-DE11-B342-000423D996C8.root',
        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/F4767383-9FE6-DE11-8099-000423D985B0.root',
        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/EE646B4A-99E6-DE11-AF39-001D09F2AF96.root',
        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/E826F9B2-AAE6-DE11-94AB-003048D3750A.root',
        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/E466542D-A8E6-DE11-89BB-000423D9853C.root',
        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/DEF9A66D-ABE6-DE11-9F1D-000423D33970.root',
        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/DCCA559F-9CE6-DE11-97E4-000423D98804.root',
        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/DC5C34E5-B3E6-DE11-9526-003048D2BE08.root',
        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/D67048AB-A3E6-DE11-9AB0-000423D94990.root',
        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/BC780CBB-A5E6-DE11-AF91-001D09F2932B.root',
        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/B2C169DD-A0E6-DE11-80FE-001D09F2525D.root',
        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/AA8AF135-A0E6-DE11-9C48-001D09F2426D.root',
        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/A81DC948-A2E6-DE11-AF25-000423D94524.root',
        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/A459C2E4-A7E6-DE11-925D-001D09F252E9.root',
        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/A286ECF8-A2E6-DE11-9C21-001617C3B778.root',
        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/A256FE95-A1E6-DE11-9F5D-000423D98B6C.root',
        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/98271EC4-AAE6-DE11-A34B-003048D37514.root',
        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/8EEB3E47-A2E6-DE11-8D31-000423D951D4.root',
        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/86A0E2DE-A0E6-DE11-945E-001D09F2462D.root',
        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/80E7934F-9DE6-DE11-B5F3-001D09F25041.root',
        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/808DAB04-9EE6-DE11-A45D-0019B9F704D6.root',
        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/6AA2CC3F-A7E6-DE11-883A-001D09F295FB.root',
        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/62A1C6E9-A7E6-DE11-B136-0019B9F730D2.root',
        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/36FB7AF0-9BE6-DE11-BED4-001D09F28755.root',
        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/32F48C17-A5E6-DE11-ABCC-001617DBD224.root',
        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/282EFE2C-A8E6-DE11-AFFD-000423D990CC.root',
        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/24F54541-A7E6-DE11-A90F-000423D99AA2.root',
        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/0ECDEC97-A1E6-DE11-BB90-000423D986C4.root',
        '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/124/009/0E502584-9FE6-DE11-B954-000423D9989E.root')
)



#  SuperClusters  ################
process.load("RecoEcal.EgammaClusterProducers.hybridSuperClusters_cfi")
process.load("RecoEcal.EgammaClusterProducers.multi5x5SuperClustersWithPreshower_cfi")


process.EBSuperClusters = cms.EDFilter("SuperClusterSelector",
    src = cms.InputTag("hybridSuperClusters"),
    cut = cms.string('abs( eta ) < 1.4442')
)

process.EESuperClusters = cms.EDFilter("SuperClusterSelector",
    src = cms.InputTag("multi5x5SuperClustersWithPreshower"),
    cut = cms.string('abs( eta ) > 1.560 & abs( eta ) < 2.5')
)

process.allSuperClusters = cms.EDFilter("SuperClusterMerger",
   src = cms.VInputTag(cms.InputTag("EBSuperClusters"), cms.InputTag("EESuperClusters"))
)


process.sc_sequence = cms.Sequence( ( (process.hybridSuperClusters*process.EBSuperClusters) +
                                      (process.multi5x5SuperClustersWithPreshower*
                                       process.EESuperClusters))* process.allSuperClusters
                                    )


## process.theGsfElectrons = cms.EDFilter("GsfElectronSelector",
##     src = cms.InputTag("gsfElectrons"),
##     cut = cms.string('((abs( eta ) < 1.4442) || (abs( eta ) > 1.560 & abs( eta ) < 3.0)) &&  et  > 15.0')
## )
## process.skimPath = cms.Path(process.HLTElectrons*process.theGsfElectrons)



process.skimPath = cms.Path( process.hltLevel1GTSeed *process.physDecl *
                             process.scrapingVeto * process.primaryVertexFilter *
                             process.sc_sequence )





#############  output module if just want to skim by HLT path ##############
process.out = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('skimPath')), 
    fileName = cms.untracked.string('SkimSuperClusters.root')
)
process.p = cms.EndPath(process.out)




#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

