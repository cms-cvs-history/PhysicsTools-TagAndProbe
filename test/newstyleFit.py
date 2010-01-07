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
      ReadFromDirectory = cms.untracked.string("MakeHisto"),
      FitFileName       = cms.untracked.string("fit_output.root"), ##<< actually not used

      UnbinnedFit  = cms.untracked.bool( True ),
      Do2DFit      = cms.untracked.bool( True ),
      #Do2DFit      = cms.untracked.bool( False ),
      HasWeights   = cms.untracked.bool( False ), ## <<< Specify there are no weights

      ## Mass window for fitting
      NumBinsMass         = cms.untracked.int32( 20 ), # matters only for the plots, the fit is unbinned
      MassLow             = cms.untracked.double( 2.8 ),
      MassHigh            = cms.untracked.double( 3.5 ),
    
      ## Variable Specifications for SB subtractions and Roofit
      NameVar1             = cms.untracked.string( "pt" ),
      Var1BinBoundaries   = cms.untracked.vdouble( 3, 4.5, 6, 10, 20),
      NameVar2             = cms.untracked.string( "eta" ),
      Var2BinBoundaries   = cms.untracked.vdouble( -2.4,-1.3,-0.8,0.8,1.3,2.4),

      PassingProbeName = cms.untracked.string("passing"), # <<< You have to specify the column that says if a probe passed or not

      GaussLineShape = cms.untracked.PSet(
	GaussMean        = cms.untracked.vdouble( 3.09, 2.9,  3.2  ),
	GaussSigma       = cms.untracked.vdouble( 0.03, 0.01, 0.05 ),
      ),

      ## Background variables
      CMSBkgLineShape = cms.untracked.PSet(
	CMSBkgAlpha           = cms.untracked.vdouble( 0 ), # fix these two to zero
	CMSBkgBeta            = cms.untracked.vdouble( 0 ), # so it's just an exp
	CMSBkgPeak            = cms.untracked.vdouble( 3.09 ),
	CMSBkgGamma           = cms.untracked.vdouble( 0, -10, 10 )
      ),

      ## Efficiency variables
      Efficiency        = cms.untracked.vdouble( 0.90,0.0,1.0 ),    
      NumSignal         = cms.untracked.vdouble( 4000.0,-10.0,30000.0 ),    
      NumBkgPass        = cms.untracked.vdouble( 4000.0,-10.0,10000.0 ),    
      NumBkgFail        = cms.untracked.vdouble( 1000.0,-10.0,7000.0  ),    

      ## Variables for sideband subtraction
      SBSPeak            = cms.untracked.double( 3.09 ),   ## Mass peak
      SBSStanDev         = cms.untracked.double(  .03 ),      ## SD from peak for subtraction

      # All the following is useless now that we just read and fit
      # but it's required...
      # --------------------------------------------
      MCTruthParentId = cms.untracked.int32(443),
)

process.p = cms.Path(process.RunFit)
process.TFileService = cms.Service("TFileService", fileName = cms.string("testNewFit.root"))
