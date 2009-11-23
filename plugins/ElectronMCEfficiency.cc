#ifndef PhysicsTools_TagAndProbe_ElectronMCEfficiency_h
#define PhysicsTools_TagAndProbe_ElectronMCEfficiency_h

#include <string>
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/ClusterShape.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Math/interface/deltaR.h" // reco::deltaR

#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TFile.h"
#include "TTree.h"

/* ******************************************************** */
/**  Store Monte Carlo generated  information to compute efficiency */
/* ******************************************************** */


class ElectronMCEfficiency : public edm::EDAnalyzer
{
 public:
  explicit ElectronMCEfficiency (const edm::ParameterSet&);
  ~ElectronMCEfficiency();

 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void boolResults(  const edm::Event&, const edm::EventSetup&, 
			     const reco::Candidate&);  
  virtual bool CheckAcceptance( const edm::Event&, 
			       const reco::Candidate&);
  const edm::ValueMap<double>& getValueMap(const edm::Event&, edm::InputTag&);

  bool isFromResDecay(int);

  bool overlap( const reco::Candidate &, 
			      const reco::Candidate &) const;

  bool cutDecision ( int, int, double, double, double, double, double,
		     double, double);

  bool CheckTriggerMatch( edm::Handle<trigger::TriggerEvent>,double, double);
  bool CheckSuperClusterMatch( const edm::Event&, const reco::Candidate&);


  // ----------member data ---------------------------

  TFile*      hOutputFile ;
  TTree *     myTree;
  std::string fOutputFileName ;


  edm::InputTag  tkIsoTag_;
  edm::InputTag  ecalIsoTag_;
  edm::InputTag  hcalIsoTag_;
  edm::InputTag  elecIDSrcLoose_;
  edm::InputTag  elecIDSrcTight_;
  edm::InputTag  triggerSummaryLabel_;
  edm::InputTag  hltFilter_;
  edm::InputTag  SuperClusters_;

  std::vector<int> truthParentId_;
  double ElectronPtCut_;
  double dRWindow_;

  float mRes;
  float Res_px;
  float Res_py;
  float Res_pz;
  float Res_E;
  float Res_Pt;
  float Res_Et;
  float Res_Eta;   
  float Res_Phi;
  float Res_Vx;
  float Res_Vy;
  float Res_Vz;
     
  float probepx;
  float probepy;
  float probepz;
  float probeE;
  float probePt;
  float probeEt;
  float probeEta;    
  float probePhi;
  int probeCharge;
  float probeVx;
  float probeVy;
  float probeVz;

  float probe_trackiso;
  float probe_ecaliso;
  float probe_hcaliso;
  float probe_SCEta;
  float probe_SCE;
  float probe_SCEt;
  float probe_deltaEta;
  float probe_deltaPhi;
  float probe_sigmaEtaEta;

  bool isInAcceptance;
  bool isSuperCluster;
  bool isGsfElectron;
  bool isIsolated;
  bool isRobustLoose;
  bool isRobustTight;
  bool isPassWoff;
  bool isPassZoff;
  bool isTriggered;
  bool isPassWon;
  bool isPassZon;
};
#endif






ElectronMCEfficiency::ElectronMCEfficiency(const edm::ParameterSet& iConfig)
{
  // electron id
  elecIDSrcLoose_     = iConfig.getParameter<edm::InputTag>("electronIDSourceLoose");
  elecIDSrcTight_     = iConfig.getParameter<edm::InputTag>("electronIDSourceTight");
  tkIsoTag_           = iConfig.getParameter<edm::InputTag>("tkIsoTag");
  ecalIsoTag_         = iConfig.getParameter<edm::InputTag>("ecalIsoTag");
  hcalIsoTag_         = iConfig.getParameter<edm::InputTag>("hcalIsoTag");
  SuperClusters_      = iConfig.getParameter<edm::InputTag>("SuperClusters");


  // trigger 
  const edm::InputTag dSummaryObj( "hltTriggerSummaryAOD","","HLT" );
  triggerSummaryLabel_ = 
    iConfig.getUntrackedParameter<edm::InputTag>("triggerSummaryLabel",  dSummaryObj );

  const edm::InputTag 
    dHLT("hltL1NonIsoHLTNonIsoSingleElectronEt15LTITrackIsolFilter","","HLT");
  hltFilter_ = iConfig.getUntrackedParameter<edm::InputTag>("hltFilter",dHLT);


  // MC truth parent Id
  std::vector<int>       dEmptyIntVec;
  truthParentId_  = iConfig.getUntrackedParameter< std::vector<int> >("MCTruthParentId", 
								      dEmptyIntVec);

  // pT cut
  ElectronPtCut_ = iConfig.getUntrackedParameter< double >("ElectronPtCut", 20.0);


  // deltaR matching window
  dRWindow_ = iConfig.getUntrackedParameter< double >("deltaR", 0.3);

  // get name of output file with histogramms
  fOutputFileName = iConfig.getUntrackedParameter<std::string>("HistOutFile",
							       "demo.root"); 
}





ElectronMCEfficiency::~ElectronMCEfficiency()
{
  // Clean up
}





  void ElectronMCEfficiency::beginJob( const edm::EventSetup& iSetup)
  {
    // Open output ROOT file and TTree
    hOutputFile   = new TFile( fOutputFileName.c_str(), "RECREATE" ) ; 
    myTree = new TTree("probe","Probe MC Truth Tree");

    myTree->Branch("mRes",        &mRes,        "mRes/F");
    myTree->Branch("Res_px",      &Res_px,      "Res_px/F");
    myTree->Branch("Res_py",      &Res_py,      "Res_py/F");
    myTree->Branch("Res_pz",      &Res_pz,      "Res_pz/F");
    myTree->Branch("Res_E",       &Res_E,       "Res_E/F");
    myTree->Branch("Res_Pt",      &Res_Pt,      "Res_Pt/F");
    myTree->Branch("Res_Et",      &Res_Et,      "Res_Et/F");
    myTree->Branch("Res_Eta",     &Res_Eta,     "Res_Eta/F");    
    myTree->Branch("Res_Phi",     &Res_Phi,     "Res_Phi/F");
    myTree->Branch("Res_Vx",      &Res_Vx,      "Res_Vx/F");
    myTree->Branch("Res_Vy",      &Res_Vy,      "Res_Vy/F");
    myTree->Branch("Res_Vz",      &Res_Vz,      "Res_Vz/F");

    myTree->Branch("probepx",     &probepx,     "probepx/F");
    myTree->Branch("probepy",     &probepy,     "probepy/F");
    myTree->Branch("probepz",     &probepz,     "probepz/F");
    myTree->Branch("probeE",      &probeE,      "probeE/F");
    myTree->Branch("probePt",     &probePt,     "probePt/F");
    myTree->Branch("probeEt",     &probeEt,     "probeEt/F");    
    myTree->Branch("probeEta",    &probeEta,    "probeEta/F");    
    myTree->Branch("probePhi",    &probePhi,    "probePhi/F");
    myTree->Branch("probeCharge", &probeCharge, "probeCharge/I");
    myTree->Branch("probeVx",     &probeVx,     "probeVx/F");
    myTree->Branch("probeVy",     &probeVy,     "probeVy/F");
    myTree->Branch("probeVz",     &probeVz,     "probeVz/F");

    myTree->Branch("probe_trackiso",    &probe_trackiso,    "probe_trackiso/F");
    myTree->Branch("probe_ecaliso",     &probe_ecaliso,     "probe_ecaliso/F");
    myTree->Branch("probe_hcaliso",     &probe_hcaliso,     "probe_hcaliso/F");
    myTree->Branch("probe_SCEta",       &probe_SCEta,       "probe_SCEta/F");
    myTree->Branch("probe_SCE",         &probe_SCE,         "probe_SCE/F");
    myTree->Branch("probe_SCEt",        &probe_SCEt,        "probe_SCEt/F");
    myTree->Branch("probe_deltaEta",    &probe_deltaEta,    "probe_deltaEta/F");
    myTree->Branch("probe_deltaPhi",    &probe_deltaPhi,    "probe_deltaPhi/F");
    myTree->Branch("probe_sigmaiEtaiEta", &probe_sigmaEtaEta, "probe_sigmaEtaEta/F");

    myTree->Branch("isInAcceptance",     &isInAcceptance, "isInAcceptance/O");
    myTree->Branch("isSuperCluster",     &isSuperCluster, "isSuperCluster/O");
    myTree->Branch("isGsfElectron",     &isGsfElectron, "isGsfElectron/O");
    myTree->Branch("isIsolated",        &isIsolated,     "isIsolated/O");
    myTree->Branch("isRobustLoose",     &isRobustLoose,  "isRobustLoose/O");
    myTree->Branch("isRobustTight",     &isRobustTight,  "isRobustTight/O");
    myTree->Branch("isPassWoff",        &isPassWoff,     "isPassWoff/O");
    myTree->Branch("isPassZoff",        &isPassZoff,     "isPassZoff/O");
    myTree->Branch("isTriggered",       &isTriggered,    "isTriggered/O");
    myTree->Branch("isPassWon",         &isPassWon,      "isPassWon/O");
    myTree->Branch("isPassZon",         &isPassZon,      "isPassZon/O");
  }




  void ElectronMCEfficiency::endJob()
  {
    hOutputFile->SetCompressionLevel(2);
     hOutputFile->cd();
     myTree->Write();

    delete myTree;
    hOutputFile->Close();
    delete hOutputFile;
  }







void ElectronMCEfficiency::analyze(const edm::Event& iEvent, 
				    const edm::EventSetup& iSetup) {

  mRes               = -1.;
  Res_px             = -99999.;
  Res_py             = -99999.;
  Res_pz             = -99999.;
  Res_E              = -1.;
  Res_Pt             = -1.;
  Res_Et             = -1.;
  Res_Eta            = -10.;
  Res_Phi            = -10.;
  Res_Vx             = -10.;
  Res_Vy             = -10.;
  Res_Vz             = -10.;

  probepx            = -99999.;
  probepy            = -99999.;
  probepz            = -99999.;
  probeE             = -1.;
  probePt            = -1.;
  probeEt            = -1.;
  probeEta           = -10.;
  probePhi           = -10.;
  probeCharge        = -10;
  probeVx            = -10.;
  probeVy            = -10.;
  probeVz            = -10.;   

  probe_trackiso = 99999.;
  probe_ecaliso  = 99999.;
  probe_hcaliso  =  99999.;
  probe_SCEta        = -10.;
  probe_SCE          = -1.;
  probe_SCEt         = -1.;
  probe_deltaEta     = -1.;
  probe_deltaPhi     = -1.;
  probe_sigmaEtaEta  = -1.;

  isInAcceptance     = false;
  isSuperCluster     = false;
  isGsfElectron      = false;
  isIsolated         = false;
  isRobustLoose      = false;
  isRobustTight      = false;
  isPassWoff         = false;
  isPassZoff         = false;
  isTriggered        = false;
  isPassWon          = false;
  isPassZon          = false;

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel("genParticles", genParticles);


  // ------------ trigger objects 
  edm::Handle<trigger::TriggerEvent> triggerObj;
  iEvent.getByLabel(triggerSummaryLabel_,triggerObj); 
  if(!triggerObj.isValid()) { 
    edm::LogInfo("TriggerEvent") << " objects not found"; 
  }
  
  size_t nZ = genParticles->size();
  if( nZ < 1 ) return;
  const reco::Candidate *Res(0);

  for(size_t i = 0; i < nZ; ++ i) {

    Res = &((*genParticles)[i]);
    size_t ndau = 0;
    if(!(Res==0)) ndau = Res->numberOfDaughters();
    if( ndau<1 || !isFromResDecay(Res->pdgId()) ) continue;


    for(size_t j = 0; j < ndau; ++ j) {
      const reco::Candidate *d = Res->daughter( j );
      if( d==0 ) continue;
      if( ! (abs(d->pdgId()) == 11) ) continue;

      isInAcceptance = CheckAcceptance(iEvent, *d);
      isSuperCluster = CheckSuperClusterMatch( iEvent, *d);

      mRes    = Res->mass();
      Res_px  = Res->px();
      Res_py  = Res->py();
      Res_pz  = Res->pz();
      Res_E   = Res->energy();
      Res_Pt  = Res->pt();
      Res_Et  = Res->et();
      Res_Eta = Res->eta();   
      Res_Phi = Res->phi();
      Res_Vx  = Res->vx();
      Res_Vy  = Res->vy();
      Res_Vz  = Res->vz();
     
      probepx     = d->px();
      probepy     = d->py();
      probepz     = d->pz();
      probeE      = d->energy();
      probePt     = d->pt();
      probeEt     = d->et();
      probeEta    = d->eta();    
      probePhi    = d->phi();
      probeCharge = d->charge();
      probeVx     = d->vx();
      probeVy     = d->vy();
      probeVz     = d->vz();

      boolResults( iEvent, iSetup, *d);

      // check trigger matching
      isTriggered = CheckTriggerMatch( triggerObj, probeEta, probePhi );

      isPassWoff    = isPassWoff && isInAcceptance;
      isPassZoff    = isPassZoff && isInAcceptance;
      isPassWon     = isPassWoff && isTriggered;
      isPassZon     = isPassZoff && isTriggered;

      // fill the tree for each candidate
      myTree->Fill();
    }
  }

}
// ---------------------------------








void ElectronMCEfficiency::boolResults( const edm::Event& iEvent, 
					const edm::EventSetup& iSetup, 
					const reco::Candidate& ele) {

  // --------- Read the electron collection in the event
  edm::Handle<edm::View<reco::GsfElectron> > electrons;
  iEvent.getByLabel("gsfElectrons", electrons);

  // ---------- electron id
  edm::Handle<edm::ValueMap<float> >  idhandleLoose;
  edm::Handle<edm::ValueMap<float> >  idhandleTight;
  iEvent.getByLabel( elecIDSrcLoose_, idhandleLoose );
  iEvent.getByLabel( elecIDSrcTight_, idhandleTight );

    
  // ---------- electron isolation
  const edm::ValueMap<double>& tkIsoMap = getValueMap(iEvent, tkIsoTag_);
  const edm::ValueMap<double>& ecalIsoMap = getValueMap(iEvent, ecalIsoTag_);
  const edm::ValueMap<double>& hcalIsoMap = getValueMap(iEvent, hcalIsoTag_);


  // ------------- cluster shape
  const edm::InputTag ebRecHit("reducedEcalRecHitsEB");
  const edm::InputTag eeRecHit("reducedEcalRecHitsEE");
  EcalClusterLazyTools lazyTools( iEvent, iSetup, ebRecHit, eeRecHit); 
   

  // loop over electron collection
  int index =0;
  for(edm::View<reco::GsfElectron>::const_iterator  
	elec = electrons->begin(); elec != electrons->end();++elec) {

    edm::RefToBase<reco::GsfElectron> electronRef(electrons, index);
    const reco::GsfElectron *x = 
      dynamic_cast<const reco::GsfElectron *>(electronRef.get());
      

    // ------------ isolation & id variables -----------------
    float electronIdLoose = (*idhandleLoose)[electronRef];
    float electronIdTight = (*idhandleTight)[electronRef];
    probe_trackiso = tkIsoMap[electronRef];
    probe_ecaliso  = ecalIsoMap[electronRef];
    probe_hcaliso  = hcalIsoMap[electronRef];


    // ---------- dPhi, dEta, sigmaEtaEta -------------
    reco::SuperClusterRef sc = elec->superCluster();
    std::vector<float> vCov2 = lazyTools.localCovariances(*(sc->seed()));
    probe_SCEta = sc.get()->eta();
    probe_SCE   = sc.get()->energy();
    probe_SCEt  = probe_SCE / cosh(probe_SCEta);
    probe_deltaEta  = elec->deltaEtaSuperClusterTrackAtVtx();
    probe_deltaPhi  = elec->deltaPhiSuperClusterTrackAtVtx();
    probe_sigmaEtaEta  = sqrt(vCov2[0]);


    int eClass = 2;
    if(electronRef->isEB()) eClass = 0;
    else if(electronRef->isEE()) eClass = 1;

    bool passW = cutDecision ( 0, eClass, probe_SCEt, probe_deltaEta, 
			       probe_deltaPhi, probe_sigmaEtaEta, probe_trackiso,
			       probe_ecaliso, probe_hcaliso);

    bool passZ = cutDecision ( 1, eClass, probe_SCEt, probe_deltaEta, 
			       probe_deltaPhi, probe_sigmaEtaEta, probe_trackiso,
			       probe_ecaliso, probe_hcaliso);

    if( overlap( *x, ele) ) {
      isGsfElectron = true;
      if( probe_trackiso / probe_SCEt < 0.2) isIsolated = true;
      if( electronIdLoose > 0.0 ) isRobustLoose = true;
      if( electronIdTight > 0.0 ) isRobustTight = true;
      isPassWoff = passW;
      isPassZoff = passZ;
      break;
    }
    else {
      probe_trackiso = 99999.;
      probe_ecaliso  = 99999.;
      probe_hcaliso  =  99999.;
      probe_SCEta        = -10.;
      probe_SCE          = -1.;
      probe_SCEt         = -1.;
      probe_deltaEta     = -1.;
      probe_deltaPhi     = -1.;
      probe_sigmaEtaEta  = -1.;

    }		

    ++index;
  }

}





// --------------- apply analysis cuts here --------------------
// classification: 0 == barrel, 1 == endcaps
// WZ:             0 == W, 1 == Z

bool ElectronMCEfficiency::cutDecision ( int WZ, int classification, 
					 double pt, double deta, 
					 double dphi, double sietaeta, double tkiso,
					 double ecaliso, double hcaliso) {
  
  double ptCut_=20.0, deltaEtaCut_=1.0, deltaPhiCut_=1.0, 
    sigmaEtaEtaCut_=1.0, tkIsoCut_=10000.0, ecalIsoCut_=10000.0, 
    hcalIsoCut_=10000.0;


  // for W selection
  if( WZ == 0 ) {
    ptCut_ = 30.0;
    if( classification ==0 ) {  // barrel
      deltaEtaCut_      = 0.0040;
      deltaPhiCut_      = 0.025;
      sigmaEtaEtaCut_   = 0.0099;
      tkIsoCut_         = 2.2;
      ecalIsoCut_       = 4.2;
      hcalIsoCut_       = 2.0;
    }
    else if( classification == 1 ) {
      deltaEtaCut_     = 0.0066;
      deltaPhiCut_     = 0.020;
      sigmaEtaEtaCut_  = 0.028;
      tkIsoCut_        = 1.1;
      ecalIsoCut_      = 3.4;
      hcalIsoCut_      = 1.3;
    }
  }

  // for Z selection
  if( WZ == 1 ) {
    ptCut_ = 20.0;
    if( classification ==0 ) {  // barrel
      deltaEtaCut_      = 0.0071;
      deltaPhiCut_      = 0.5;
      sigmaEtaEtaCut_   = 0.010;
      tkIsoCut_         = 7.2;
      ecalIsoCut_       = 5.7;
      hcalIsoCut_       = 8.1;
    }
    else if( classification == 1 ) {
      deltaEtaCut_     = 0.0066;
      deltaPhiCut_     = 0.5;
      sigmaEtaEtaCut_  = 0.028;
      tkIsoCut_        = 5.1;
      ecalIsoCut_      = 5.0;
      hcalIsoCut_      = 3.4;
    }
  }


  bool decision = (pt > ptCut_) && 
    (fabs(deta) < deltaEtaCut_) && 
    (fabs(dphi) < deltaPhiCut_) &&
    (sietaeta < sigmaEtaEtaCut_) && (tkiso < tkIsoCut_) && 
    (ecaliso < ecalIsoCut_) && (hcaliso < hcalIsoCut_);

  return decision;
}






////////// check the MC truth Id of the parent ///////////////////

bool ElectronMCEfficiency::isFromResDecay(int pdgId) {
  
  bool result = false;  
  for(int j=0; j< (int) truthParentId_.size(); j++) {
    if(pdgId == truthParentId_[j]) {
      result = true;
      break;
    }
  }

  return result;
}





////////// Apply event selection cuts ///////////////////

bool ElectronMCEfficiency::CheckAcceptance( const edm::Event& iEvent, 
					   const reco::Candidate& ele ) {

  bool result = true;

  float eEta = ele.eta();
  float ePt  = ele.pt();


  // electron pT cut
  if( ePt < ElectronPtCut_ ) result = false;

  // electron acceptance
  if( !((fabs(eEta)<1.4442) || 
	(fabs(eEta)>1.560 && fabs(eEta)<2.5)) ) result = false;

  return result;
}






//  ************* Utility: check overlap **********************
bool ElectronMCEfficiency::overlap( const reco::Candidate & e, 
				    const reco::Candidate & c) const {

  if( fabs(deltaR(e,c) ) < dRWindow_ ) return true;
    
  return false;
}
//--------------------------------------------




//little labour saving function to get the reference to the ValueMap
const edm::ValueMap<double>& ElectronMCEfficiency::getValueMap(const edm::Event& ev, 
						       edm::InputTag& inTag)
{
  edm::Handle<edm::ValueMap<double> > handle;
  ev.getByLabel(inTag,handle);
  return *(handle.product());
}






// ------------- perform trigger matching ----------------------

bool ElectronMCEfficiency::CheckTriggerMatch( edm::Handle<trigger::TriggerEvent> 
					      triggerObj, double eta, double phi) {

  bool result = false;
  const trigger::TriggerObjectCollection & toc(triggerObj->getObjects());
  const int trigindex = triggerObj->filterIndex( hltFilter_ );
  if ( trigindex >= triggerObj->sizeFilters() ) return false; 

  const trigger::Keys & l1k = triggerObj->filterKeys( trigindex );
  if ( l1k.size() <= 0 ) return false; 


  for (trigger::Keys::const_iterator ki = l1k.begin(); ki !=l1k.end(); ++ki ) {
    
    if (reco::deltaR( eta, phi, toc[*ki].eta(),toc[*ki].phi()) < dRWindow_) 
      result = true;      /////// ---------- trigger match                 
  }
  
  return result;
}




////////// Does the MC generated electron match a super cluster ///////////////////

bool ElectronMCEfficiency::CheckSuperClusterMatch( const edm::Event& iEvent, 
						   const reco::Candidate& ele ) {

  bool result = false;

  edm::Handle<edm::View<reco::Candidate> > SuperClusters;
  iEvent.getByLabel( SuperClusters_, SuperClusters);

  for(edm::View<reco::Candidate>::const_iterator  
	sc = SuperClusters->begin(); sc != SuperClusters->end();++sc) {

    if( overlap( *sc, ele) ) result = true; 
  }
  return result; 
}



//define this as a plug-in
DEFINE_FWK_MODULE( ElectronMCEfficiency );

