import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer("TagProbeEDMAnalysis",
    # Efficiency/Fitting variables                
    CalculateEffSideBand = cms.untracked.bool(True),
    CalculateEffFitter = cms.untracked.bool(True), ## effs from Roofit
    UnbinnedFit = cms.untracked.bool(True),
    CalculateEffTruth = cms.untracked.bool(True), ## true effs
    Do2DFit = cms.untracked.bool(True),
    # if want to float signal shape parameters in each bin: VERY SLOW                  
    floatAllShapeParameters  = cms.untracked.bool(False),                  
    # Type of Efficiency : 0 => SC-->GsfElectron
    # 1 ==> GsfElectron-->isolation
    # 2 ==> isolation-->id
    # 3 ==> id-->HLT
    TagProbeType = cms.untracked.int32(0),
    # if want to use reconstructed values for truth-matched candidates                   
    useRecoVarsForTruthMatchedCands = cms.untracked.bool(True),
    # Root file to eff histograms to
    FitFileName = cms.untracked.string('demo.root'),                   
    # Variables for sideband subtraction
    SBSPeak = cms.untracked.double(91.1876),
    SBSStanDev = cms.untracked.double(5.0), ## SD from peak for subtraction                  
    # Variable Specifications for SB subtractions and Roofit                      
    # Choose binned or unbinned fitter ...
    # Note that the unbinned fit will not fit weighted data,
    # if you wish to use weights, use the binned fit.
    # Variables and binning for the eff hists
    # Valid variables names for the eff binning are:
    # "pt","p","px","py","pz","e","et","eta" and "phi"
    # This way of declaring the bin will overide any other
    # If omitted the defaults are var1 = pt and var2 = eta
    NameVar1 = cms.untracked.string('pt'),
    NumBinsVar1 = cms.untracked.int32(4),
    Var1Low = cms.untracked.double(20.0),                   
    Var1High = cms.untracked.double(100.0),
    NameVar2 = cms.untracked.string('eta'),                      
    NumBinsVar2 = cms.untracked.int32(5),
    Var2Low = cms.untracked.double(-2.5),
    Var2High = cms.untracked.double(2.5),              
    # There is also an option to read the variables in 
    # via a file. This allows for much greater binning flexability
    #
    # Fitter variables - for the Roofit fitter
    # If you want the variable to float in the fit fill
    # three array elements {default, range_low, range_high}
    # If the variable should be fixed, fill one element {value}                      
    #      
    # Background variables
    ## Background variables
    CMSBkgLineShape = cms.untracked.PSet(
      CMSBkgAlpha           = cms.untracked.vdouble( 62., 50, 70 ),
      CMSBkgBeta            = cms.untracked.vdouble( 0.001, 0.0, 0.1 ),
      CMSBkgPeak            = cms.untracked.vdouble( 91.1876 ),
      CMSBkgGamma           = cms.untracked.vdouble( 0.05,0.0, 0.1 )
    ),
    NumBkgPass = cms.untracked.vdouble(1000.0, 0.0, 1000000.0),
    NumBkgFail = cms.untracked.vdouble(10.0, 0.0, 1000000.0),                      
    # Signal variables
    ## Fitter variables - for the Roofit fitter
    ## If you want the variable to float in the fit fill
    ## three array elements {default, range_low, range_high}
    ## If the variable should be fixed, fill one element {value}
    ## Signal variables
    ZLineShape = cms.untracked.PSet(
     ZMean        = cms.untracked.vdouble( 91.1876, 90, 92),
     ZWidth       = cms.untracked.vdouble( 2.495,   1.5, 10.0 ),
     ZSigma       = cms.untracked.vdouble( 0.75,    0.1, 5.0 ),
     ZWidthL      = cms.untracked.vdouble( 15.0,    1,   30.0 ),
     ZWidthR      = cms.untracked.vdouble( 4.0,     1, 15.0 ),
     ZBifurGaussFrac    = cms.untracked.vdouble( 0.10, 0.0, 1.0 )
    ),
    ## Other possible signal line shapes, uncomment to use                                 
    ##  CBLineShape = cms.untracked.PSet(
    ##	CBMean      = cms.untracked.vdouble( 9.4603,9.0,10.0 ),
    ##	CBSigma     = cms.untracked.vdouble( 0.75,0.01,5.0 ),
    ##	CBAlpha     = cms.untracked.vdouble( 8.0,0.0,30.0 ),
    ##	CBN         = cms.untracked.vdouble( 1.0,0.0,30.0 )
    ## ),
    ## GaussLineShape = cms.untracked.PSet(
    ##	GaussMean     = cms.untracked.vdouble( 9.4603,9.0,10.0 ),
    ##	GaussSigma    = cms.untracked.vdouble( 0.75,0.01,5.0 )
    ## ),                      
    NumSignal = cms.untracked.vdouble(4000.0, 0.0, 1000000.0),                        
    # Mass window for fitting
    # untracked int32 NumBinsMass         = 60
    # untracked double MassLow            = 60.0
    # untracked double MassHigh           = 120.0
    NumBinsMass = cms.untracked.int32(60),
    MassLow = cms.untracked.double(60.0),
    MassHigh = cms.untracked.double(120.0),                      
    # Efficiency variables
    Efficiency = cms.untracked.vdouble(0.9, 0.0, 1.0),
    # Make some plots of tree variables ...
    quantities = cms.untracked.vstring('TPmass'),
    # Binning for the above plots 
    XBins = cms.untracked.vuint32(60),
    XMin = cms.untracked.vdouble(60.0),                  
    XMax = cms.untracked.vdouble(120.0),
    logY = cms.untracked.vuint32(1),
    # file names and condition for the above plots                 
    outputFileNames = cms.untracked.vstring('Zmass_pass.eps'),                           
    conditions = cms.untracked.vstring('TPppass==1')
)


