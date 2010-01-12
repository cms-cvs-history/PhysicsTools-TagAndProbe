import FWCore.ParameterSet.Config as cms

numbers=range(1,432)
numbers.remove(66)
numbers.remove(206)
numbers.remove(371)
numbers.remove(381)
numbers.remove(408)
numbers.remove(108)
numbers.remove(159)
numbers.remove(210)
numbers.remove(315)

source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    ['dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/resilient/makouski/minBiasFiltered/SkimSuperClusters_%i.root'%i for i in numbers]
   )
)


