import FWCore.ParameterSet.Config as cms 
 
from PhysicsTools.HepMCCandAlgos.genParticles_cfi import * 
from Geometry.CMSCommonData.cmsIdealGeometryXML_cff import * 
from Geometry.CaloEventSetup.CaloGeometry_cff import * 
from Configuration.EventContent.EventContent_cff import * 
 
process = cms.Process("ElectronEff") 
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("CondCore.DBCommon.CondDBSetup_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('MC_31X_V3::All')

readFiles=[]


readFiles.extend( [
       '/store/mc/Summer09/Zee/GEN-SIM-RECO/MC_31X_V3-v1/0028/DAA8926E-CE8B-DE11-9654-0030487DF78A.root',
       '/store/mc/Summer09/Zee/GEN-SIM-RECO/MC_31X_V3-v1/0027/E6394B5A-338B-DE11-A710-00151764221C.root',
       '/store/mc/Summer09/Zee/GEN-SIM-RECO/MC_31X_V3-v1/0027/B00D1565-438B-DE11-A443-001E682F8528.root',
       '/store/mc/Summer09/Zee/GEN-SIM-RECO/MC_31X_V3-v1/0027/481B38F7-338B-DE11-94AB-00304853CC4A.root',
       '/store/mc/Summer09/Zee/GEN-SIM-RECO/MC_31X_V3-v1/0024/40D768B7-1D89-DE11-9248-001517642260.root',
       '/store/mc/Summer09/Zee/GEN-SIM-RECO/MC_31X_V3-v1/0022/F6C3053C-8C88-DE11-B969-00144F0D84D8.root',
       '/store/mc/Summer09/Zee/GEN-SIM-RECO/MC_31X_V3-v1/0022/F4536933-EA87-DE11-9A0E-000423D2F35E.root',
       '/store/mc/Summer09/Zee/GEN-SIM-RECO/MC_31X_V3-v1/0022/F4502B77-BA87-DE11-94C6-003048984614.root',
       '/store/mc/Summer09/Zee/GEN-SIM-RECO/MC_31X_V3-v1/0022/F2D492C1-8687-DE11-B0C5-00144F9E7CDE.root',
       '/store/mc/Summer09/Zee/GEN-SIM-RECO/MC_31X_V3-v1/0022/F01A0DFC-8387-DE11-B1C5-001E6837DFEA.root',
       '/store/mc/Summer09/Zee/GEN-SIM-RECO/MC_31X_V3-v1/0022/EEDB882D-0289-DE11-9283-001517641EB0.root',
       '/store/mc/Summer09/Zee/GEN-SIM-RECO/MC_31X_V3-v1/0022/EE23B82A-DF87-DE11-B617-00E081402C6F.root',
       '/store/mc/Summer09/Zee/GEN-SIM-RECO/MC_31X_V3-v1/0022/E8E30907-4288-DE11-A856-0030487CD620.root',
       '/store/mc/Summer09/Zee/GEN-SIM-RECO/MC_31X_V3-v1/0022/E854CB2D-8587-DE11-934C-00096BB5BECE.root',
       '/store/mc/Summer09/Zee/GEN-SIM-RECO/MC_31X_V3-v1/0022/E69C476B-0189-DE11-9725-001517641EB0.root',
       '/store/mc/Summer09/Zee/GEN-SIM-RECO/MC_31X_V3-v1/0022/D854EB0A-4288-DE11-8185-001A921D5A7D.root',
       '/store/mc/Summer09/Zee/GEN-SIM-RECO/MC_31X_V3-v1/0022/D4FE7EF2-078D-DE11-A17E-001E0B469778.root',
       '/store/mc/Summer09/Zee/GEN-SIM-RECO/MC_31X_V3-v1/0022/D4B2359B-0289-DE11-B902-000423CA664C.root'])

files=["12095AEE-D8AB-DE11-8964-00304867905A.root",
    "30516DA9-F9AB-DE11-9562-003048678FFA.root",
    "663EA598-F5AB-DE11-BB78-001A92810AE6.root",
    "A0582EBF-E7AB-DE11-9733-0018F3D0969A.root",
    "AE047ABD-D4AB-DE11-B943-0018F3D096EC.root",
    "C62D4337-F1AB-DE11-A96B-0018F3D0970E.root"    
]
FL=['/store/mc/Summer09/Zee/GEN-SIM-RECO/MC_31X_V3_SD_Ele15-v1/0003/%s'%i for i in files]



process.source = cms.Source("PoolSource",
     fileNames = cms.untracked.vstring(readFiles)
)

process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(5000) 
) 
 
process.MessageLogger.destinations = ['cout', 'cerr'] 
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
 


#  SuperClusters  ################

 
process.HybridSuperClusters = cms.EDProducer("ConcreteEcalCandidateProducer",
    src = cms.InputTag("correctedHybridSuperClusters"),
    particleType = cms.string('gamma')
)
process.EBSuperClusters = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("HybridSuperClusters"),
    cut = cms.string('et>20.0 & abs( eta ) < 1.4442')
)



process.EndcapSuperClusters = cms.EDProducer("ConcreteEcalCandidateProducer",
    src = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower"),
    particleType = cms.string('gamma')
)
process.EESuperClusters = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("EndcapSuperClusters"),
    cut = cms.string('et>20.0 & abs( eta ) > 1.560 & abs( eta ) < 2.5')
)

process.allSuperClusters = cms.EDFilter("CandViewMerger",
   src = cms.VInputTag(cms.InputTag("EBSuperClusters"), cms.InputTag("EESuperClusters"))
)


process.sc_sequence = cms.Sequence( (process.HybridSuperClusters *
                                     process.EBSuperClusters +
                                     process.EndcapSuperClusters *
                                     process.EESuperClusters) * process.allSuperClusters)




# Cut-based Robust electron ID  ######
from RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi import *
import RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi
process.eid = RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi.eidCutBasedExt.clone()
process.eidTight = RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi.eidCutBasedExt.clone()
process.eidTight.electronQuality = cms.string('tight')


#  electron isolation  ################
#
process.load("RecoEgamma.EgammaIsolationAlgos.eleIsoDepositTk_cff")
process.load("RecoEgamma.EgammaIsolationAlgos.eleIsoDepositEcalFromHits_cff")
process.load("RecoEgamma.EgammaIsolationAlgos.eleIsoDepositHcalFromTowers_cff")
process.load("RecoEgamma.EgammaIsolationAlgos.eleIsoFromDepsModules_cff")



process.mceff = cms.EDAnalyzer("ElectronMCEfficiency",
    SuperClusters = cms.InputTag('allSuperClusters'),                        
    electronIDSourceLoose = cms.InputTag('eid'),
    electronIDSourceTight = cms.InputTag('eidTight'),                          
    tkIsoTag = cms.InputTag('eleIsoFromDepsTk'),
    ecalIsoTag = cms.InputTag('eleIsoFromDepsEcalFromHits'),
    hcalIsoTag = cms.InputTag('eleIsoFromDepsHcalFromTowers'),
    hltFilter  =  cms.untracked.InputTag( "hltL1NonIsoHLTNonIsoSingleElectronEt15LTITrackIsolFilter","","HLT"),
    MCTruthParentId = cms.untracked.vint32(23,22),                           
    ElectronPtCut = cms.untracked.double(20.0),
    deltaR = cms.untracked.double(0.3),
    HistOutFile = cms.untracked.string('ntuple_ZeeSEP.root')
)

process.p = cms.Path(
    #process.sc_sequence +
    #
    process.HybridSuperClusters +
    process.EBSuperClusters +
    process.EndcapSuperClusters +
    process.EESuperClusters + process.allSuperClusters+
    #
    process.eleIsoDepositTk + process.eleIsoFromDepsTk +
    process.eleIsoDepositEcalFromHits +
    process.eleIsoFromDepsEcalFromHits +
    process.eleIsoDepositHcalFromTowers +  
    process.eleIsoFromDepsHcalFromTowers +
    process.eid + process.eidTight + process.mceff
    )
