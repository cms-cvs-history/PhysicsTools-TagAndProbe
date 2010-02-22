import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )    

process.load("PhysicsTools.TagAndProbe.TagProbeFitTreeAnalyzer_cfi")
process.TagProbeFitTreeAnalyzer.OutputFileName = cms.untracked.string("testNewAnalyzer.root")
process.TagProbeFitTreeAnalyzer.Variables.mass = ["Tag-Probe Mass", "2.5", "3.8", "GeV/c^{2}"]
process.TagProbeFitTreeAnalyzer.Variables.pt   = ["Probe p_{T}",      "3",  "20", "GeV/c"    ]
process.TagProbeFitTreeAnalyzer.Variables.eta  = ["Probe #eta",    "-2.4", "2.4", ""         ]
process.TagProbeFitTreeAnalyzer.PDFs = cms.untracked.PSet(
    gaussPlusExpo = cms.untracked.vstring(
        "Gaussian::signal(mass, mean[3.09,2.9,3.2], sigma[0.03,0.01,0.05])",
        "Exponential::backgroundPass(mass, slope[0,-10,10])",
        "Exponential::backgroundFail(mass, slope)",
        "efficiency[0.9,0,1]",
        "signalFractionInPassing[0.9]"
    )
)
process.TagProbeFitTreeAnalyzer.Efficiencies = cms.untracked.PSet(
    pt = cms.untracked.PSet(
        EfficiencyCategoryAndState = cms.untracked.vstring("passing","pass"),
        UnbinnedVariables = cms.untracked.vstring("mass"),
        BinnedVariables = cms.untracked.PSet(
            pt = cms.untracked.vdouble(3.0, 4.5, 6.0, 10.0, 20.0),
        ),
        BinToPDFmap = cms.untracked.vstring("gaussPlusExpo"), # single model for all bins
    ),
    pt_mcTrue = cms.untracked.PSet(
        EfficiencyCategoryAndState = cms.untracked.vstring("passing","pass"),
        UnbinnedVariables = cms.untracked.vstring("mass"),
        BinnedVariables = cms.untracked.PSet(
            mcTrue = cms.untracked.vstring("true"),
            pt = cms.untracked.vdouble(3.0, 4.5, 6.0, 10.0, 20.0),
        ),
        # no BinToPDFmap: no fit
    ),
    eta = cms.untracked.PSet(
        EfficiencyCategoryAndState = cms.untracked.vstring("passing","pass"),
        UnbinnedVariables = cms.untracked.vstring("mass"),
        BinnedVariables = cms.untracked.PSet(
            eta = cms.untracked.vdouble(-2.4,-1.3,-0.8,0.8,1.3,2.4),
        ),
        BinToPDFmap = cms.untracked.vstring("gaussPlusExpo"), # single model for all bins
    ),
    eta_mcTrue = cms.untracked.PSet(
        EfficiencyCategoryAndState = cms.untracked.vstring("passing","pass"),
        UnbinnedVariables = cms.untracked.vstring("mass"),
        BinnedVariables = cms.untracked.PSet(
            mcTrue = cms.untracked.vstring("true"),
            eta = cms.untracked.vdouble(-2.4,-1.3,-0.8,0.8,1.3,2.4),
        ),
        # no BinToPDFmap: no fit
    ),
    pt_eta = cms.untracked.PSet(
        EfficiencyCategoryAndState = cms.untracked.vstring("passing","pass"),
        UnbinnedVariables = cms.untracked.vstring("mass"),
        BinnedVariables = cms.untracked.PSet(
            pt = cms.untracked.vdouble(3.0, 4.5, 6.0, 10.0, 20.0),
            eta = cms.untracked.vdouble(-2.4,-1.3,-0.8,0.8,1.3,2.4),
        ),
        BinToPDFmap = cms.untracked.vstring("gaussPlusExpo"), # single model for all bins
    ),
    pt_eta_mcTrue = cms.untracked.PSet(
        EfficiencyCategoryAndState = cms.untracked.vstring("passing","pass"),
        UnbinnedVariables = cms.untracked.vstring("mass"),
        BinnedVariables = cms.untracked.PSet(
            mcTrue = cms.untracked.vstring("true"),
            pt = cms.untracked.vdouble(3.0, 4.5, 6.0, 10.0, 20.0),
            eta = cms.untracked.vdouble(-2.4,-1.3,-0.8,0.8,1.3,2.4),
        ),
        # no BinToPDFmap: no fit
    )
)


process.fitness = cms.Path(
    process.TagProbeFitTreeAnalyzer
)

