#ifndef PhysicsTools_TagAndProbe_ElectronEfficiency_h
#define PhysicsTools_TagAndProbe_ElectronEfficiency_h

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/ClusterShape.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"

#include "DataFormats/Math/interface/deltaR.h" // reco::deltaR
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/JetExtendedAssociation.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"

//Hcal Noise Objects
#include "RecoMET/METAlgorithms/interface/HcalNoiseRBXArray.h"
#include "DataFormats/METReco/interface/HcalNoiseHPD.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"


#include "TFile.h"
#include "TTree.h"
#include "TNamed.h"
#include <vector>
#include <string>
#include <map>


using namespace edm;
using namespace reco;
using namespace std;



class ElectronEfficiency : public edm::EDAnalyzer
{
 public:
  explicit ElectronEfficiency (const edm::ParameterSet&);
  ~ElectronEfficiency();

 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  
  void buildTree();
  void clearTreeVectors(void);


  // ----------member data ---------------------------

  TFile*      hOutputFile ;
  TTree *     mTree;

  //---- configurable parameters --------  
  std::string fOutputFileName ;

  double mJetPtMin;
  std::string mJetsName;
  std::string mJetsIDName;
  std::string mJetExtender;
  std::string mMetName;
  std::string mMetNoHFName;
  edm::InputTag  mSuperClusters;
  edm::InputTag mHcalNoiseTag;

  int mRunNo,mEvtNo,mLumi,mBunch, mNSC, mNGsf, mNPV, mLooseHcalNoise, mTightHcalNoise;
  double mMET, mMETnoHF, mSumET, mSumETnoHF;

  std::vector<int>    *mNtrkVtx,*mNtrkCalo,*mN90,*mN90Hits,*mPVntracks;
  std::vector<double> *mSCeta,*mSCphi,*mSCe,*mSCet, *mSCx, *mSCy, 
    *mSCz, *mSCtheta, *mSCpx, *mSCpy, *mSCpz, *mSCpt, *mSCE9overE25,
    *mSCSigmaiEtaiEta, *mSCSigmaiEtaiPhi, *mSCSigmaiPhiiPhi;

  std::vector<int>    *mGSF_charge, *mGSF_numberOfBrems, *mGSF_Classification;
  std::vector<bool>   *mGSF_isEB, *mGSF_isEE, *mGSF_isGap;
  std::vector<double> *mGSF_px, *mGSF_py, *mGSF_pz, *mGSF_e,
    *mGSF_pt, *mGSF_et, *mGSF_eta, *mGSF_phi, *mGSF_vx, *mGSF_vy,
    *mGSF_vz, *mGSF_deltaEta, *mGSF_deltaPhi, *mGSF_EoverP,
    *mGSF_sigmaEtaEta, *mGSF_sigmaiEtaiEta, *mGSF_e1x5, *mGSF_e2x5Max, 
    *mGSF_e5x5, *mGSF_H1overE, *mGSF_H2overE, *mGSF_HoverE, 
    *mGSF_dr03TkSumPt, *mGSF_dr03EcalRecHitSumEt, *mGSF_dr03HcalDepth1TowerSumEt,
    *mGSF_dr03HcalDepth2TowerSumEt, *mGSF_dr03HcalTowerSumEt,
    *mGSF_dr04TkSumPt, *mGSF_dr04EcalRecHitSumEt, *mGSF_dr04HcalDepth1TowerSumEt,
    *mGSF_dr04HcalDepth2TowerSumEt, *mGSF_dr04HcalTowerSumEt, *mGSF_fbrem;

  std::vector<double> *mE,*mPt,*mEta,*mEtaD,*mPhi,*mY,*mEmf;
  std::vector<double> *mTrkCaloPt,*mTrkCaloEta,*mTrkCaloPhi;
  std::vector<double> *mTrkVtxPt,*mTrkVtxEta,*mTrkVtxPhi;
  std::vector<double> *mfHPD,*mfRBX,*mEtaMoment,*mPhiMoment;
  std::vector<double> *mPVx,*mPVy,*mPVz,*mPVchi2,*mPVndof;
  std::vector<double> *mfHcalNoise;
 
  typedef math::PtEtaPhiELorentzVectorF LorentzVector;
};
#endif





ElectronEfficiency::ElectronEfficiency(const edm::ParameterSet& iConfig)
{
  mSuperClusters          = iConfig.getParameter<edm::InputTag>            ("SuperClusters");
  mJetsName               = iConfig.getParameter<std::string>              ("jets");
  mJetsIDName             = iConfig.getParameter<std::string>              ("jetsID");
  mMetName                = iConfig.getParameter<std::string>              ("met");
  mMetNoHFName            = iConfig.getParameter<std::string>              ("metNoHF");
  mJetExtender            = iConfig.getParameter<std::string>              ("jetExtender");
  mHcalNoiseTag           = iConfig.getParameter<edm::InputTag>            ("hcalNoiseTag");
  mJetPtMin               = iConfig.getParameter<double>                   ("minJetPt");

  // get name of output file with histogramms
  fOutputFileName = iConfig.getUntrackedParameter<std::string>("HistOutFile","demo.root");

}





ElectronEfficiency::~ElectronEfficiency()
{
  // Clean up
}





void ElectronEfficiency::beginJob( const edm::EventSetup& iSetup)
{
  // Open output ROOT file and TTree
  hOutputFile   = new TFile( fOutputFileName.c_str(), "RECREATE" ) ; 
  mTree = new TTree("probe","Probe Data Tree");
  buildTree();
}




void ElectronEfficiency::endJob()
{
  hOutputFile->SetCompressionLevel(2);
  hOutputFile->cd();
  mTree->Write();
  
  delete mTree;
  hOutputFile->Close();
  delete hOutputFile;
}




void ElectronEfficiency::analyze(const edm::Event& event, 
				 const edm::EventSetup& iSetup){
  

  mRunNo = 0;
  mEvtNo = 0;
  mLumi  = 0;
  mBunch = 0;
  mNSC   = 0;
  mNGsf  = 0; 
  mNPV   = 0; 
  mLooseHcalNoise = 0;
  mTightHcalNoise = 0;
  
  clearTreeVectors();
  mRunNo = event.id().run();
  mEvtNo = event.id().event();
  mLumi  = event.luminosityBlock();
  mBunch = event.bunchCrossing();

  ////////////Vertices//////////////
  edm::Handle<reco::VertexCollection> recVtxs;
  event.getByLabel("offlinePrimaryVertices",recVtxs);
  for(unsigned int ind=0;ind<recVtxs->size();ind++) 
    {
      if (!((*recVtxs)[ind].isFake())) 
        {
	  mNPV += 1;
	  mPVx->push_back( (*recVtxs)[ind].x() );
	  mPVy->push_back( (*recVtxs)[ind].y() );
	  mPVz->push_back( (*recVtxs)[ind].z() );
        }
    }




  //-------- get SC collection ----------
  edm::Handle<reco::SuperClusterCollection> SuperClusters;
  event.getByLabel( mSuperClusters, SuperClusters);
  int counter = 0;

  for(reco::SuperClusterCollection::const_iterator  
	sc = SuperClusters->begin(); sc != SuperClusters->end();++sc) {

    counter++;

    mSCeta->push_back(sc->eta());
    mSCphi->push_back(sc->phi());
    mSCe  ->push_back(sc->energy());
    mSCet ->push_back(sc->energy()/cosh(sc->eta()));   

    mSCx->push_back( sc->x() );
    mSCy->push_back( sc->y() );
    mSCz->push_back( sc->z() );
    mSCtheta->push_back( sc->position().Theta() );
    mSCpx->push_back( sc->energy() * sin(sc->position().Theta()) * cos(sc->phi()) );
    mSCpy->push_back( sc->energy() * sin(sc->position().Theta()) * sin(sc->phi()) );
    mSCpz->push_back( sc->energy() * cos(sc->position().Theta()) );
    mSCpt->push_back( sc->energy()/cosh(sc->eta()) );

  
    const edm::InputTag ebRecHit("reducedEcalRecHitsEB");
    const edm::InputTag eeRecHit("reducedEcalRecHitsEE");
    
    EcalClusterLazyTools lazyTools( event, iSetup, ebRecHit, eeRecHit); 
    std::vector<float> vCov = lazyTools.localCovariances(*(sc->seed()));
    mSCE9overE25->push_back( lazyTools.e3x3(*(sc->seed()))/
      lazyTools.e5x5(*(sc->seed())) );
    mSCSigmaiEtaiEta->push_back(  sqrt(vCov[0]) );
    mSCSigmaiEtaiPhi->push_back(  sqrt(vCov[1]) );
    mSCSigmaiPhiiPhi->push_back(  sqrt(vCov[2]) );
  }  

  mNSC = counter;




  // --------- Read the electron collection in the event
  edm::Handle<edm::View<reco::GsfElectron> > electrons;
  event.getByLabel("gsfElectrons", electrons);
  counter = 0;

  for(edm::View<reco::GsfElectron>::const_iterator  
	elec = electrons->begin(); elec != electrons->end();++elec) {

    counter++;

    mGSF_px->push_back(     elec->px() );
    mGSF_py ->push_back(    elec->py() );
    mGSF_pz->push_back(     elec->pz() );
    mGSF_e->push_back(      elec->energy() );
    mGSF_pt->push_back(     elec->pt() );
    mGSF_et->push_back(     elec->et() );
    mGSF_eta->push_back(    elec->eta() );    
    mGSF_phi->push_back(    elec->phi() );
    mGSF_charge->push_back( elec->charge() );
    mGSF_vx->push_back(     elec->vx() );
    mGSF_vy->push_back(     elec->vy() );
    mGSF_vz->push_back(     elec->vz() );

    mGSF_deltaEta->push_back(     elec->deltaEtaSuperClusterTrackAtVtx() );
    mGSF_deltaPhi->push_back(     elec->deltaPhiSuperClusterTrackAtVtx() );
    mGSF_EoverP->push_back(       elec->eSuperClusterOverP() );
    mGSF_sigmaEtaEta->push_back(  elec->sigmaEtaEta() ); 
    mGSF_sigmaiEtaiEta->push_back(elec->sigmaIetaIeta() ); 
    mGSF_e1x5->push_back(         elec->e1x5() ); 
    mGSF_e2x5Max->push_back(      elec->e2x5Max() ); 
    mGSF_e5x5->push_back(         elec->e5x5() ); 
    mGSF_H1overE->push_back(      elec->hcalDepth1OverEcal()  );
    mGSF_H2overE->push_back(      elec->hcalDepth2OverEcal() ); 
    mGSF_HoverE->push_back(       elec->hcalOverEcal() ); 


    // dR=0.3 isolations
    mGSF_dr03TkSumPt->push_back(  elec->dr03TkSumPt() );
    mGSF_dr03EcalRecHitSumEt->push_back(  elec->dr03EcalRecHitSumEt() );
    mGSF_dr03HcalDepth1TowerSumEt->push_back(  elec->dr03HcalDepth1TowerSumEt() );
    mGSF_dr03HcalDepth2TowerSumEt->push_back(  elec->dr03HcalDepth2TowerSumEt() );
    mGSF_dr03HcalTowerSumEt->push_back(  elec->dr03HcalTowerSumEt() );

    // dR=0.4 isolations
    mGSF_dr04TkSumPt->push_back(  elec->dr04TkSumPt() );
    mGSF_dr04EcalRecHitSumEt->push_back(  elec->dr04EcalRecHitSumEt()  );
    mGSF_dr04HcalDepth1TowerSumEt->push_back(  elec->dr04HcalDepth1TowerSumEt() );
    mGSF_dr04HcalDepth2TowerSumEt->push_back(  elec->dr04HcalDepth2TowerSumEt() );
    mGSF_dr04HcalTowerSumEt->push_back(  elec->dr04HcalTowerSumEt() );
 

    // electron in barrel or endcap ?
    mGSF_isEB->push_back(  elec->isEB() );
    mGSF_isEE->push_back(  elec->isEE() );
    mGSF_isGap->push_back(  elec->isGap() );
    mGSF_fbrem->push_back(  elec->fbrem() );
    mGSF_numberOfBrems->push_back(   elec->numberOfBrems() );
    mGSF_Classification->push_back(   elec->classification() );
  }

  mNGsf = counter;


  ///////// HcalNoiseCollection //////////////                 
  edm::Handle<HcalNoiseRBXCollection> rbxColl;
  event.getByLabel(mHcalNoiseTag,rbxColl);
  if(!rbxColl.isValid()) {
    throw edm::Exception(edm::errors::ProductNotFound)
      << " could not find HcalNoiseRBXCollection named " << "hcalnoise" << ".\n";
    return;
  }

  // loop over the RBXs 
  std::map<CaloTowerDetId, double> hcalNoise;
  for(HcalNoiseRBXCollection::const_iterator rit=rbxColl->begin(); 
      rit!=rbxColl->end(); ++rit) {
    HcalNoiseRBX rbx    = (*rit);
    const std::vector<HcalNoiseHPD> hpds = rbx.HPDs();
    for(unsigned ihpd=0; ihpd<hpds.size(); ihpd++){
      const edm::RefVector<CaloTowerCollection> noiseCTowers = hpds[ihpd].caloTowers();
      for(unsigned int itow=0; itow<noiseCTowers.size(); itow++){
        hcalNoise.insert( std::pair<CaloTowerDetId, double>( noiseCTowers[itow]->id(), 
							     noiseCTowers[itow]->hadEnergy()) );
      }
    }
  }


  ////////////// Jets //////
  edm::Handle<CaloJetCollection> jets;
  event.getByLabel(mJetsName,jets);
  edm::Handle<JetExtendedAssociation::Container> jetExtender;
  event.getByLabel(mJetExtender,jetExtender);
  edm::Handle<edm::ValueMap<reco::JetID> > jetsID;
  event.getByLabel(mJetsIDName,jetsID);
  if ((*jets).size() < 1) return;
  for(unsigned int ind=0;ind<(*jets).size();ind++) 
    {
      if ((*jets)[ind].pt() < mJetPtMin) continue;
      LorentzVector TrkCaloP4 = JetExtendedAssociation::tracksAtCaloP4(*jetExtender,
								       (*jets)[ind]);
      LorentzVector TrkVtxP4  = JetExtendedAssociation::tracksAtVertexP4(*jetExtender, 
									 (*jets)[ind]);
      edm::RefToBase<reco::Jet> jetRef(edm::Ref<reco::CaloJetCollection>(jets,ind));
      mPt        ->push_back((*jets)[ind].pt());
      mEta       ->push_back((*jets)[ind].eta());
      mEtaD      ->push_back((*jets)[ind].detectorP4().eta());
      mY         ->push_back((*jets)[ind].y()); 
      mPhi       ->push_back((*jets)[ind].phi());
      mE         ->push_back((*jets)[ind].energy());
      mN90       ->push_back((*jets)[ind].n90());
      mfHPD      ->push_back((*jetsID)[jetRef].fHPD);
      mfRBX      ->push_back((*jetsID)[jetRef].fRBX); 
      mEmf       ->push_back((*jets)[ind].emEnergyFraction()); 	
      mEtaMoment ->push_back((*jets)[ind].etaetaMoment()); 
      mPhiMoment ->push_back((*jets)[ind].phiphiMoment());
      mNtrkVtx   ->push_back(JetExtendedAssociation::tracksAtVertexNumber(*jetExtender,
									  (*jets)[ind]));
      mNtrkCalo  ->push_back(JetExtendedAssociation::tracksAtCaloNumber(*jetExtender,
									(*jets)[ind])); 
      mTrkCaloPt ->push_back(TrkCaloP4.pt());
      mTrkCaloEta->push_back(TrkCaloP4.eta());
      mTrkCaloPhi->push_back(TrkCaloP4.phi());
      mTrkVtxPt  ->push_back(TrkVtxP4.pt());
      mTrkVtxEta ->push_back(TrkVtxP4.eta());
      mTrkVtxPhi ->push_back(TrkVtxP4.phi());
      
      double jetEneNoise=0.0;
      std::vector< CaloTowerPtr >jTowers = (*jets)[ind].getCaloConstituents ();
      for(unsigned int itow=0; itow<jTowers.size(); itow++) {
	std::map<CaloTowerDetId, double>::iterator thisTow = 
	  hcalNoise.find(jTowers[itow]->id());
	if( thisTow != hcalNoise.end() ) jetEneNoise += jTowers[itow]->energy();
      }
      mfHcalNoise->push_back( jetEneNoise/(*jets)[ind].energy() );
    } 


  ////////////// MET //////
  edm::Handle<CaloMETCollection> met;
  event.getByLabel(mMetName,met);
  if (met->size() == 0)
    {
      mMET   = -1;
      mSumET = -1;
    }
  else
    {
      mMET   = (*met)[0].et();
      mSumET = (*met)[0].sumEt();
    }

  edm::Handle<CaloMETCollection> metNoHF;
  event.getByLabel(mMetNoHFName,metNoHF);
  if (metNoHF->size() == 0)
    {
      mMETnoHF   = -1;
      mSumETnoHF = -1;
    }
  else
    {
      mMETnoHF   = (*metNoHF)[0].et();
      mSumETnoHF = (*metNoHF)[0].sumEt();
    }



  mTree->Fill();
  
}
// ---------------------------------

void ElectronEfficiency::clearTreeVectors(void){
  mSCeta->clear();
  mSCphi->clear();
  mSCe->clear();
  mSCet->clear();
  mSCx->clear(); 
  mSCy->clear();
  mSCz->clear(); 
  mSCtheta->clear(); 
  mSCpx->clear(); 
  mSCpy->clear(); 
  mSCpz->clear(); 
  mSCpt->clear(); 
  mSCE9overE25->clear();
  mSCSigmaiEtaiEta->clear();
  mSCSigmaiEtaiPhi->clear(); 
  mSCSigmaiPhiiPhi->clear();

  mGSF_px->clear();
  mGSF_py ->clear();
  mGSF_pz->clear();
  mGSF_e->clear();
  mGSF_pt->clear();
  mGSF_et->clear();
  mGSF_eta->clear();    
  mGSF_phi->clear();
  mGSF_charge->clear();
  mGSF_vx->clear();
  mGSF_vy->clear();
  mGSF_vz->clear();
  mGSF_deltaEta->clear();
  mGSF_deltaPhi->clear();
  mGSF_EoverP->clear();
  mGSF_sigmaEtaEta->clear(); 
  mGSF_sigmaiEtaiEta->clear(); 
  mGSF_e1x5->clear(); 
  mGSF_e2x5Max->clear(); 
  mGSF_e5x5->clear(); 
  mGSF_H1overE->clear();
  mGSF_H2overE->clear(); 
  mGSF_HoverE->clear(); 
  mGSF_dr03TkSumPt->clear();
  mGSF_dr03EcalRecHitSumEt->clear();
  mGSF_dr03HcalDepth1TowerSumEt->clear();
  mGSF_dr03HcalDepth2TowerSumEt->clear();
  mGSF_dr03HcalTowerSumEt->clear();
  mGSF_dr04TkSumPt->clear();
  mGSF_dr04EcalRecHitSumEt->clear();
  mGSF_dr04HcalDepth1TowerSumEt->clear();
  mGSF_dr04HcalDepth2TowerSumEt->clear();
  mGSF_dr04HcalTowerSumEt->clear();
  mGSF_isEB->clear();
  mGSF_isEE->clear();
  mGSF_isGap->clear();
  mGSF_fbrem->clear();
  mGSF_numberOfBrems->clear();
  mGSF_Classification->clear();

  mPt        ->clear();
  mEta       ->clear();
  mEtaD      ->clear();
  mY         ->clear();
  mPhi       ->clear();
  mE         ->clear();
  mEmf       ->clear();
  mEtaMoment ->clear();
  mPhiMoment ->clear();
  mNtrkVtx   ->clear();
  mNtrkCalo  ->clear();
  mTrkCaloPt ->clear();
  mTrkCaloEta->clear();
  mTrkCaloPhi->clear();
  mTrkVtxPt  ->clear();
  mTrkVtxEta ->clear();
  mTrkVtxPhi ->clear();
  mN90       ->clear(); 
  mfHPD      ->clear();
  mfRBX      ->clear();
  mfHcalNoise->clear();
  mPVx       ->clear();
  mPVy       ->clear();
  mPVz       ->clear();
}



void ElectronEfficiency::buildTree(){

  mPt         = new std::vector<double>();
  mEta        = new std::vector<double>();
  mEtaD       = new std::vector<double>();
  mY          = new std::vector<double>();
  mPhi        = new std::vector<double>();
  mE          = new std::vector<double>();
  mEmf        = new std::vector<double>();
  mEtaMoment  = new std::vector<double>();
  mPhiMoment  = new std::vector<double>();
  mNtrkVtx    = new std::vector<int>   ();
  mNtrkCalo   = new std::vector<int>   ();
  mTrkCaloPt  = new std::vector<double>();
  mTrkCaloEta = new std::vector<double>();
  mTrkCaloPhi = new std::vector<double>();
  mTrkVtxPt   = new std::vector<double>();
  mTrkVtxEta  = new std::vector<double>();
  mTrkVtxPhi  = new std::vector<double>();
  mN90        = new std::vector<int>   ();
  mfHPD       = new std::vector<double>();
  mfRBX       = new std::vector<double>();
  mfHcalNoise = new std::vector<double>();

  mPVx        = new std::vector<double>();
  mPVy        = new std::vector<double>();
  mPVz        = new std::vector<double>();

  mSCeta = new std::vector<double>();
  mSCphi = new std::vector<double>();
  mSCe = new std::vector<double>();
  mSCet = new std::vector<double>();
  mSCx = new std::vector<double>(); 
  mSCy = new std::vector<double>();
  mSCz = new std::vector<double>(); 
  mSCtheta = new std::vector<double>(); 
  mSCpx = new std::vector<double>(); 
  mSCpy = new std::vector<double>(); 
  mSCpz = new std::vector<double>(); 
  mSCpt = new std::vector<double>(); 
  mSCE9overE25 = new std::vector<double>();
  mSCSigmaiEtaiEta = new std::vector<double>();
  mSCSigmaiEtaiPhi = new std::vector<double>(); 
  mSCSigmaiPhiiPhi = new std::vector<double>();


  mGSF_px = new std::vector<double>();
  mGSF_py  = new std::vector<double>();
  mGSF_pz = new std::vector<double>();
  mGSF_e = new std::vector<double>();
  mGSF_pt = new std::vector<double>();
  mGSF_et = new std::vector<double>();
  mGSF_eta = new std::vector<double>();    
  mGSF_phi = new std::vector<double>();
  mGSF_charge = new std::vector<int>();
  mGSF_vx = new std::vector<double>();
  mGSF_vy = new std::vector<double>();
  mGSF_vz = new std::vector<double>();
  mGSF_deltaEta = new std::vector<double>();
  mGSF_deltaPhi = new std::vector<double>();
  mGSF_EoverP = new std::vector<double>();
  mGSF_sigmaEtaEta = new std::vector<double>(); 
  mGSF_sigmaiEtaiEta = new std::vector<double>(); 
  mGSF_e1x5 = new std::vector<double>(); 
  mGSF_e2x5Max = new std::vector<double>(); 
  mGSF_e5x5 = new std::vector<double>(); 
  mGSF_H1overE = new std::vector<double>();
  mGSF_H2overE = new std::vector<double>(); 
  mGSF_HoverE = new std::vector<double>(); 
  mGSF_dr03TkSumPt = new std::vector<double>();
  mGSF_dr03EcalRecHitSumEt = new std::vector<double>();
  mGSF_dr03HcalDepth1TowerSumEt = new std::vector<double>();
  mGSF_dr03HcalDepth2TowerSumEt = new std::vector<double>();
  mGSF_dr03HcalTowerSumEt = new std::vector<double>();
  mGSF_dr04TkSumPt = new std::vector<double>();
  mGSF_dr04EcalRecHitSumEt = new std::vector<double>();
  mGSF_dr04HcalDepth1TowerSumEt = new std::vector<double>();
  mGSF_dr04HcalDepth2TowerSumEt = new std::vector<double>();
  mGSF_dr04HcalTowerSumEt = new std::vector<double>();
  mGSF_isEB = new std::vector<bool>();
  mGSF_isEE = new std::vector<bool>();
  mGSF_isGap = new std::vector<bool>();
  mGSF_fbrem = new std::vector<double>();
  mGSF_numberOfBrems = new std::vector<int>();
  mGSF_Classification = new std::vector<int>();


  mTree->Branch("event_evtNo"      ,&mEvtNo               ,"mEvtNo/I");
  mTree->Branch("event_runNo"      ,&mRunNo               ,"mRunNo/I");
  mTree->Branch("event_lumi"       ,&mLumi                ,"mLumi/I");
  mTree->Branch("event_bunch"      ,&mBunch               ,"mBunch/I"); 

  mTree->Branch("SC_eta"           ,"vector<double>"      ,&mSCeta);
  mTree->Branch("SC_phi"           ,"vector<double>"      ,&mSCphi);
  mTree->Branch("SC_e"             ,"vector<double>"      ,&mSCe);
  mTree->Branch("SC_et"            ,"vector<double>"       ,&mSCet);
  mTree->Branch("SC_x"             ,"vector<double>"      ,&mSCx); 
  mTree->Branch("SC_y"             ,"vector<double>"      ,&mSCy);
  mTree->Branch("SC_z"             ,"vector<double>"      ,&mSCz); 
  mTree->Branch("SC_theta"         ,"vector<double>"      ,&mSCtheta); 
  mTree->Branch("SC_px"            ,"vector<double>"      ,&mSCpx); 
  mTree->Branch("SC_py"            ,"vector<double>"      ,&mSCpy); 
  mTree->Branch("SC_pz"            ,"vector<double>"      ,&mSCpz); 
  mTree->Branch("SC_pt"            ,"vector<double>"      ,&mSCpt); 
  mTree->Branch("SC_e9e25"         ,"vector<double>"      ,&mSCE9overE25);
  mTree->Branch("SC_sigmaiEtaiEta" ,"vector<double>"      ,&mSCSigmaiEtaiEta);
  mTree->Branch("SC_sigmaiEtaiPhi" ,"vector<double>"      ,&mSCSigmaiEtaiPhi); 
  mTree->Branch("SC_sigmaiPhiiPhi" ,"vector<double>"      ,&mSCSigmaiPhiiPhi);
  mTree->Branch("nSC"              ,&mNSC                 ,"mNSC/I");


  mTree->Branch("Gsf_px"           ,"vector<double>"      ,&mGSF_px);
  mTree->Branch("Gsf_py"           ,"vector<double>"      ,&mGSF_py);
  mTree->Branch("Gsf_pz"           ,"vector<double>"      ,&mGSF_pz);
  mTree->Branch("Gsf_e"            ,"vector<double>"      ,&mGSF_e);
  mTree->Branch("Gsf_pt"           ,"vector<double>"      ,&mGSF_pt);
  mTree->Branch("Gsf_et"           ,"vector<double>"      ,&mGSF_et);
  mTree->Branch("Gsf_eta"          ,"vector<double>"      ,&mGSF_eta);    
  mTree->Branch("Gsf_phi"          ,"vector<double>"      ,&mGSF_phi);
  mTree->Branch("Gsf_charge"       ,"vector<int>"         ,&mGSF_charge);
  mTree->Branch("Gsf_vx"           ,"vector<double>"      ,&mGSF_vx);
  mTree->Branch("Gsf_vy"           ,"vector<double>"      ,&mGSF_vy);
  mTree->Branch("Gsf_vz"           ,"vector<double>"      ,&mGSF_vz);
  mTree->Branch("Gsf_deltaEta"     ,"vector<double>"      ,&mGSF_deltaEta);
  mTree->Branch("Gsf_deltaPhi"     ,"vector<double>"      ,&mGSF_deltaPhi);
  mTree->Branch("Gsf_EoverP"       ,"vector<double>"      ,&mGSF_EoverP);
  mTree->Branch("Gsf_sigmaEtaEta"  ,"vector<double>"      ,&mGSF_sigmaEtaEta); 
  mTree->Branch("Gsf_sigmaiEtaiEta","vector<double>"      ,&mGSF_sigmaiEtaiEta); 
  mTree->Branch("Gsf_e1x5"         ,"vector<double>"      ,&mGSF_e1x5); 
  mTree->Branch("Gsf_e2x5Max"      ,"vector<double>"      ,&mGSF_e2x5Max); 
  mTree->Branch("Gsf_e5x5"         ,"vector<double>"      ,&mGSF_e5x5); 
  mTree->Branch("Gsf_H1overE"      ,"vector<double>"      ,&mGSF_H1overE);
  mTree->Branch("Gsf_H2overE"      ,"vector<double>"      ,&mGSF_H2overE); 
  mTree->Branch("Gsf_HoverE"       ,"vector<double>"      ,&mGSF_HoverE); 
  mTree->Branch("Gsf_dr03TrackIso" ,"vector<double>"      ,&mGSF_dr03TkSumPt);
  mTree->Branch("Gsf_dr03EcalIso"  ,"vector<double>"      ,&mGSF_dr03EcalRecHitSumEt);
  mTree->Branch("Gsf_dr03HcalIso1" ,"vector<double>"      ,&mGSF_dr03HcalDepth1TowerSumEt);
  mTree->Branch("Gsf_dr03HcalIso2" ,"vector<double>"      ,&mGSF_dr03HcalDepth2TowerSumEt);
  mTree->Branch("Gsf_dr03HcalIso"  ,"vector<double>"      ,&mGSF_dr03HcalTowerSumEt);
  mTree->Branch("Gsf_dr04TrackIso" ,"vector<double>"      ,&mGSF_dr04TkSumPt);
  mTree->Branch("Gsf_dr04EcalIso"  ,"vector<double>"      ,&mGSF_dr04EcalRecHitSumEt);
  mTree->Branch("Gsf_dr04HcalIso1" ,"vector<double>"      ,&mGSF_dr04HcalDepth1TowerSumEt);
  mTree->Branch("Gsf_dr04HcalIso2" ,"vector<double>"      ,&mGSF_dr04HcalDepth2TowerSumEt);
  mTree->Branch("Gsf_dr04HcalIso"  ,"vector<double>"      ,&mGSF_dr04HcalTowerSumEt);
  mTree->Branch("Gsf_isEB"         ,"vector<bool>"        ,&mGSF_isEB);
  mTree->Branch("Gsf_isEE"         ,"vector<bool>"        ,&mGSF_isEE);
  mTree->Branch("Gsf_isGap"        ,"vector<bool>"        ,&mGSF_isGap);
  mTree->Branch("Gsf_fbrem"        ,"vector<double>"      ,&mGSF_fbrem);
  mTree->Branch("Gsf_numberOfBrems","vector<int>"         ,&mGSF_numberOfBrems);
  mTree->Branch("Gsf_classification","vector<int>"        ,&mGSF_Classification);
  mTree->Branch("nGsf"              ,&mNGsf               ,"mNGsf/I");


  mTree->Branch("jet_pt"         ,"vector<double>"      ,&mPt);
  mTree->Branch("jet_eta"        ,"vector<double>"      ,&mEta);
  mTree->Branch("jet_etaDetector","vector<double>"      ,&mEtaD);
  mTree->Branch("jet_y"          ,"vector<double>"      ,&mY);
  mTree->Branch("jet_phi"        ,"vector<double>"      ,&mPhi);
  mTree->Branch("jet_e"          ,"vector<double>"      ,&mE);
  mTree->Branch("jet_emf"        ,"vector<double>"      ,&mEmf);
  mTree->Branch("jet_etaMoment"  ,"vector<double>"      ,&mEtaMoment);
  mTree->Branch("jet_phiMoment"  ,"vector<double>"      ,&mPhiMoment);
  mTree->Branch("jet_nTrkVtx"    ,"vector<int>"        ,&mNtrkVtx);
  mTree->Branch("jet_nTrkCalo"   ,"vector<int>"        ,&mNtrkCalo);
  mTree->Branch("jet_TrkCaloPt"  ,"vector<double>"      ,&mTrkCaloPt);
  mTree->Branch("jet_TrkCaloEta" ,"vector<double>"      ,&mTrkCaloEta);
  mTree->Branch("jet_TrkCaloPhi" ,"vector<double>"      ,&mTrkCaloPhi);
  mTree->Branch("jet_TrkVtxPt"   ,"vector<double>"      ,&mTrkVtxPt);
  mTree->Branch("jet_TrkVtxEta"  ,"vector<double>"      ,&mTrkVtxEta);
  mTree->Branch("jet_TrkVtxPhi"  ,"vector<double>"      ,&mTrkVtxPhi);
  mTree->Branch("jet_n90"        ,"vector<int>"        ,&mN90);
  mTree->Branch("jet_fHPD"       ,"vector<double>"      ,&mfHPD);
  mTree->Branch("jet_fRBX"       ,"vector<double>"      ,&mfRBX);  
  mTree->Branch("jet_fHcalNoise" ,"vector<double>"      ,&mfHcalNoise);

  mTree->Branch("event_PVx"        ,"vector<double>"      ,&mPVx);
  mTree->Branch("event_PVy"        ,"vector<double>"      ,&mPVy);
  mTree->Branch("event_PVz"        ,"vector<double>"      ,&mPVz);
  mTree->Branch("event_nPV"                ,&mNPV                 ,"mNPV/I");

  mTree->Branch("event_met"        ,&mMET                ,"mMET/F");
  mTree->Branch("event_sumet"      ,&mSumET              ,"mSumET/F");
  mTree->Branch("event_metNoHF"    ,&mMETnoHF            ,"mMETnoHF/F");
  mTree->Branch("event_sumetNoHF"  ,&mSumETnoHF          ,"mSumETnoHF/F");


}


//define this as a plug-in
DEFINE_FWK_MODULE( ElectronEfficiency );

