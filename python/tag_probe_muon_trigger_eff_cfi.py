import FWCore.ParameterSet.Config as cms

muonCands = cms.EDFilter(
    "MuonRefSelector",
    src = cms.InputTag("muons"),
    cut = cms.string('isGlobalMuon > 0 &  pt > 10.0 & abs(eta)<2.1')
 )

## idCands = cms.EDFilter("MuonSelectorWmunu",
## 	src = cms.InputTag("muonCands")
## )

## tagCands = cms.EDProducer("eTriggerCandidateCollection",
##     InputProducer = cms.InputTag('muonCands'),
##     hltTag = cms.untracked.InputTag("HLT_Mu9","","HLT"),
##     triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","","HLT")
##     )


tagCands = cms.EDProducer("TrgMatchedMuonRefProducer",
	ProbeCollection = cms.InputTag("muonCands"),
	triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","","HLT"),
	muonFilterTags = cms.untracked.VInputTag(
	cms.InputTag("HLT_Mu9","","HLT")
 #       cms.InputTag("hltSingleMuPrescale3L3PreFiltered","","HLT")
	),
	usePtMatching = cms.untracked.bool(False)
)

# 3. Make a collection of opposite sign global muons, the probes
probeCands = cms.EDFilter(
     "MuonRefSelector",
    src = cms.InputTag("muons"),
    cut = cms.string('isGlobalMuon > 0 & pt > 10.0 & abs(eta)<2.1 ')
)


passProbeCands  = cms.EDProducer("TrgMatchedMuonRefProducer",
	ProbeCollection = cms.InputTag("probeCands"),
	triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","","HLT"),
	muonFilterTags = cms.untracked.VInputTag(
	cms.InputTag("HLT_Mu9","","HLT")
 #       cms.InputTag("hltSingleMuPrescale3L3PreFiltered","","HLT")
	),
	usePtMatching = cms.untracked.bool(False)
)
## # 4. Passing probes are those that fired the trigger
## passProbeCands = cms.EDProducer("eTriggerCandidateCollection",
##     InputProducer = cms.InputTag('probeCands'),
##     hltTag = cms.untracked.InputTag("HLT_Mu9","","HLT"),
##     triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","","HLT")
##     )

## cms.EDProducer(
##     "TrgMatchedMuonRefProducer",
##     ProbeCollection = cms.InputTag("probeCands"),
##     triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD"),
##     muonFilterTags = cms.untracked.VInputTag(
##         cms.InputTag("hltSingleMuNoIsoL3PreFiltered9","","HLT")
##   #      cms.InputTag("hltSingleMuPrescale3L3PreFiltered","","HLT")
## 	),
##     usePtMatching = cms.untracked.bool(True)
## )

# Make the tag probe association map
muonTagProbeMap = cms.EDProducer("TagProbeProducer",
    MassMaxCut = cms.untracked.double(120.0),
    TagCollection = cms.InputTag("tagCands"),
    MassMinCut = cms.untracked.double(50.0),
    ProbeCollection = cms.InputTag("probeCands")
)

# find generator particles matching by DeltaR
tagMuonMatch = cms.EDFilter("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(13),
    src = cms.InputTag("tagCands"),
    distMin = cms.double(0.15),
    matched = cms.InputTag("genParticles")
)

allProbeMuonMatch = cms.EDFilter("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(13),
    src = cms.InputTag("probeCands"),
    distMin = cms.double(0.15),
    matched = cms.InputTag("genParticles")
)

passProbeMuonMatch = cms.EDFilter("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(13),
    src = cms.InputTag("passProbeCands"),
    distMin = cms.double(0.15),
    matched = cms.InputTag("genParticles")
)

muon_cands = cms.Sequence(muonCands+tagCands+probeCands+passProbeCands+muonTagProbeMap+tagMuonMatch+allProbeMuonMatch+passProbeMuonMatch)
