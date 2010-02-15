import FWCore.ParameterSet.Config as cms
process = cms.Process("Fit")

# Add your own files here
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )    

process.RunFit = cms.EDAnalyzer("TagProbeEDMAnalysis",  
      ## Efficiency/Fitting variables
      CalculateEffSideBand = cms.untracked.bool( True ), ## Calculate and store effs using SB
      CalculateEffFitter   = cms.untracked.bool( True ), ## Calculate and store effs from Roofit
      CalculateEffTruth    = cms.untracked.bool( True ), ## Calculate and store true effs

      ## Set mode to read from files ...
      Mode = cms.untracked.string("ReadNew"),

      ReadFromFiles     = cms.untracked.vstring("testNewWrite.root"),
      ReadFromDirectory = cms.untracked.string("GsfToIso"),
      FitFileName       = cms.untracked.string("fit_output.root"), ##<< actually not used

      UnbinnedFit  = cms.untracked.bool( True ),
      Do2DFit      = cms.untracked.bool( True ),
      #Do2DFit      = cms.untracked.bool( False ),
      HasWeights   = cms.untracked.bool( False ), ## <<< Specify there are no weights

      ## Mass window for fitting
      NumBinsMass         = cms.untracked.int32( 20 ), # matters only for the plots, the fit is unbinned
      MassLow             = cms.untracked.double( 60.0 ),
      MassHigh            = cms.untracked.double( 120.0 ),
    
      ## Variable Specifications for SB subtractions and Roofit
      NameVar1             = cms.untracked.string( "pt" ),
      Var1BinBoundaries   = cms.untracked.vdouble( 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120),
      NameVar2             = cms.untracked.string( "eta" ),
      Var2BinBoundaries   = cms.untracked.vdouble( -2.4,-1.3,-0.8,0.8,1.3,2.4),

      PassingProbeName = cms.untracked.string("passing"), # <<< You have to specify the column that says if a probe passed or not

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

      ## Background variables
      CMSBkgLineShape = cms.untracked.PSet(
        CMSBkgAlpha           = cms.untracked.vdouble( 62., 50, 70 ),
        CMSBkgBeta            = cms.untracked.vdouble( 0.001, 0.0, 0.1 ),
        CMSBkgPeak            = cms.untracked.vdouble( 91.1876 ),
        CMSBkgGamma           = cms.untracked.vdouble( 0.05,0.0, 0.1 )
      ),

      ## Efficiency variables
      Efficiency        = cms.untracked.vdouble( 0.90,0.0,1.0 ),    
      NumSignal         = cms.untracked.vdouble( 4000.0,-10.0,1000000.0 ),    
      NumBkgPass        = cms.untracked.vdouble( 4000.0,-10.0,1000000.0 ),    
      NumBkgFail        = cms.untracked.vdouble( 1000.0,-10.0,1000000.0  ),    

      ## Variables for sideband subtraction
      SBSPeak            = cms.untracked.double( 91.18 ),   ## Mass peak
      SBSStanDev         = cms.untracked.double(  5.0 ),      ## SD from peak for subtraction

      # All the following is useless now that we just read and fit
      # but it's required...
      # --------------------------------------------
      MCTruthParentId = cms.untracked.int32(22),
)

process.p = cms.Path(process.RunFit)
process.TFileService = cms.Service("TFileService", fileName = cms.string("testNewFit.root"))
