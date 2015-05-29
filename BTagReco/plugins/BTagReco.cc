// -*- C++ -*-
//
// Package:    BTagReco/BTagReco
// Class:      BTagReco
// 
/**\class BTagReco BTagReco.cc BTagReco/BTagReco/plugins/BTagReco.cc

 Description:
 This class aims to map miniAOD infos into flat tree for stand-alone analyses (using only root)

 Implementation:
 The implementation follows some rules. 
 If you modify part of the code, you are kindly invited to be coherent with the style it is written.
 For instance, pay attention to  
  the way you write comments
  the indentation
  the use of methods and functions in the main code 

 Usefull reference:
 https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD
 Variables are stored in a tree described by CTree, implemented in TreeVariables.h
*/
//
// Original Author:  Francesco Romeo
//         Created:  Mon, 19 Jan 2015 08:54:10 GMT
//
//
/////
//   Headers
/////
//System and event
#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TStopwatch.h"
//Gen info
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//Trigger
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
//Vertex
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
//Muon
#include "DataFormats/PatCandidates/interface/Muon.h"
//Electron
#include "DataFormats/PatCandidates/interface/Electron.h"
//Tau
#include "DataFormats/PatCandidates/interface/Tau.h"
//Photon
#include "DataFormats/PatCandidates/interface/Photon.h"
//Jet and Met
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
//Packed
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
//Track builder infos
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexUpdator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexTrackCompatibilityEstimator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexSmoother.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
//Math
#include "DataFormats/Math/interface/deltaR.h"
//Track extrapolator
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalTrajectoryExtrapolatorToLine.h"
//Store info
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
//User defined classes
#include "BTagRunII/BTagReco/interface/TreeVariables.h"
#include "BTagRunII/BTagReco/interface/ObjEvtFunctions.h"
#include "BTagRunII/BTagReco/interface/JetVariables.h"
/////
//   Namespace
/////
using namespace reco;
using namespace edm;
using namespace std;
/////
//   Class declaration
/////
class BTagReco : public edm::EDAnalyzer {
 public:
 explicit BTagReco(const edm::ParameterSet&);
 ~BTagReco();
 static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
 private:
 //Default methods
 virtual void beginJob() override;
 virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
 virtual void endJob() override;
 //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
 //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
 //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
 //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
 //Collections of objects
 edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
 edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
 edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
 edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
 edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
 edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
 edm::EDGetTokenT<pat::MuonCollection> muonToken_;
 edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
 edm::EDGetTokenT<pat::TauCollection> tauToken_;
 edm::EDGetTokenT<pat::PhotonCollection> photonToken_;
 edm::EDGetTokenT<pat::JetCollection> jetToken_;
 edm::EDGetTokenT<pat::JetCollection> fatjetToken_;
 edm::EDGetTokenT<pat::METCollection> metToken_;
 edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
 edm::EDGetTokenT<pat::PackedCandidateCollection> lostTracksToken_;
 //Values for begin job
 int evt_totnum;
 //Values for the whole analysis
 const double mindr_p3;
 const double mindr_p5;
 const bool   first_jet_highest_btag;
 const bool   first_jet_lowest_btag;
 //Watch time and cpu for the analysis
 TStopwatch* stopwatch;
 //Tree
 CTree *tree;
 const edm::Service<TFileService> fs;
};
/////
//   Constructors and destructor
/////
BTagReco::BTagReco(const edm::ParameterSet& iConfig):
 //Collections of objects
 prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
 packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
 triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
 triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
 triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
 vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
 muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
 electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
 tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
 photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
 jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
 fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
 metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
 pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
 lostTracksToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("lostTracks"))),
 //Values for the whole analysis
 mindr_p3(iConfig.getUntrackedParameter<double>("mindr_p3")),
 mindr_p5(iConfig.getUntrackedParameter<double>("mindr_p5")),
 first_jet_highest_btag(iConfig.getParameter<bool>("first_jet_highest_btag")),
 first_jet_lowest_btag(iConfig.getParameter<bool>("first_jet_lowest_btag")),
 //TTree 
 tree(new CTree(fs->make<TTree>("tree", "tree")))   
{
 //Now do what ever initialization is needed
 tree->make_branches();
 stopwatch = new TStopwatch();
}

BTagReco::~BTagReco()
{
 // do anything here that needs to be done at desctruction time
 // (e.g. close files, deallocate resources etc.)
 delete stopwatch;
}
/////
//   Member functions
/////
// ------------ method called for each event  ------------
void BTagReco::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
 using namespace edm;
 evt_totnum++; 
 /////
 //   Handle iEvent.getByToken
 /////
 Handle<edm::View<reco::GenParticle> > pruned;
 iEvent.getByToken(prunedGenToken_,pruned);
 Handle<edm::View<pat::PackedGenParticle> > packed;
 iEvent.getByToken(packedGenToken_,packed);
 edm::Handle<edm::TriggerResults> triggerBits;
 iEvent.getByToken(triggerBits_, triggerBits);
 edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
 iEvent.getByToken(triggerObjects_, triggerObjects);
 edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
 iEvent.getByToken(triggerPrescales_, triggerPrescales);
 edm::Handle<reco::VertexCollection> vertices;
 iEvent.getByToken(vtxToken_, vertices);
 edm::Handle<pat::MuonCollection> muons;
 iEvent.getByToken(muonToken_, muons);
 edm::Handle<pat::ElectronCollection> electrons;
 iEvent.getByToken(electronToken_, electrons);
 edm::Handle<pat::TauCollection> taus;
 iEvent.getByToken(tauToken_, taus);
 edm::Handle<pat::PhotonCollection> photons;
 iEvent.getByToken(photonToken_, photons);
 edm::Handle<pat::JetCollection> jets;
 iEvent.getByToken(jetToken_, jets);
 PFJetIDSelectionFunctor pfLooseJetID(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE),
                         pfTightJetID(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT);
 pat::strbitset passLooseCuts(pfLooseJetID.getBitTemplate()),
                passTightCuts(pfTightJetID.getBitTemplate());
 edm::Handle<pat::JetCollection> fatjets;
 iEvent.getByToken(fatjetToken_, fatjets);
 edm::Handle<pat::METCollection> mets;
 iEvent.getByToken(metToken_, mets);
 edm::Handle<pat::PackedCandidateCollection> pfs;
 iEvent.getByToken(pfToken_, pfs);
 edm::Handle<pat::PackedCandidateCollection> lostrks;
 iEvent.getByToken(lostTracksToken_, lostrks);
 edm::ESHandle<TransientTrackBuilder> ttrkbuilder;
 iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrkbuilder);
 KalmanVertexFitter vtxFitter(true);
 /////
 //   Primary vertex
 /////
 if(vertices->empty()) return; // skip the event if no PV found
 const reco::Vertex &PV = vertices->front();
 //tree->loop_initialize();
 /////
 //   Initial skim
 /////
 vector<int>  candpfcpvindices;
 vector<const reco::Candidate*> looseleps; 
 vector<const reco::Candidate*> tightleps; 
 //Muons
 for(const pat::Muon &mu : *muons){
  if(!is_loose_muon(mu)) continue;
  //if(!is_cand_inpv(mindr_p3, mu, *pfs, candpfcpvindices)) continue; 
  looseleps.push_back((const reco::Candidate*)&mu);
  if(!is_tight_muon(mu,PV)) continue;
  tightleps.push_back((const reco::Candidate*)&mu);
 }
 //Electrons
 for(const pat::Electron &ele : *electrons){
  if(!is_loose_electron(ele,PV)) continue;
  //if(!is_cand_inpv(mindr_p3, ele, *pfs, candpfcpvindices)) continue;
  bool matchelemu = false;
  for(uint gl=0; gl<looseleps.size(); gl++) if(deltaR(looseleps[gl]->p4(),ele.p4())<mindr_p3) matchelemu = true;
  if(matchelemu) continue;
  looseleps.push_back((const reco::Candidate*)&ele);
  if(!is_tight_electron(ele,PV)) continue;
  tightleps.push_back((const reco::Candidate*)&ele);
 }
 int lep_numl = looseleps.size();
 //tree->lep_numl = lep_numl;
 int lep_numt = tightleps.size();
 //tree->lep_numt = lep_numt;
 if(lep_numl>1) tree->lep_dichprodl = looseleps[0]->charge()*looseleps[1]->charge();
 if(lep_numt>1) tree->lep_dichprodt = tightleps[0]->charge()*tightleps[1]->charge();
 //Jets
 int jet_pos = 0;
 int jet_num = 0;
 vector<pair<double,int> > jet_csv_pos;
 for(const pat::Jet &j : *jets){
  if(!is_good_jet(j)){jet_pos++; continue;}
  bool jetmatchedlepts = false;
  for(uint gl=0; gl<tightleps.size(); gl++) if(deltaR(tightleps[gl]->p4(),j.p4())<mindr_p5) jetmatchedlepts = true;
  if(jetmatchedlepts){jet_pos++; continue;}
  double csvcurrjet = j.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags");
  jet_csv_pos.push_back(make_pair(csvcurrjet,jet_pos));
  jet_pos++;
  jet_num++;
 }
 //tree->jet_num = jet_num; 
 //We want at least one jet per event (For our study, we want indeed to consider one jet per event)
 if(jet_num==0) return; 
 //first_jet_highest_btag means that the first jet is the one with the highest value of the b tag discriminator (for signal)
 //if(first_jet_highest_btag) sort(jet_csv_pos.rbegin(), jet_csv_pos.rend());
 //first_jet_lowest_btag means that the first jet is the one with the lowest value of the b tag discriminator (for bkg)
 //if(first_jet_lowest_btag)  sort(jet_csv_pos.begin(), jet_csv_pos.end());
 //Now we take the first jet and we consider only one jet per each event
 const pat::Jet & j = (*jets)[jet_csv_pos[0].second]; 
 //Gen level association: 5 for b jets, !5 for non b jets
 if(first_jet_highest_btag && abs(j.partonFlavour())!=5) return; //The first jet must be a b at parton level (for signal) 
 if(first_jet_lowest_btag  && (abs(j.partonFlavour()) ==5 || abs(j.partonFlavour()) ==4)) return; //The first jet must not be a b at parton level (for bkg) 
 tree->loop_initialize();
 tree->lep_numl      = lep_numl;
 tree->lep_numt      = lep_numt;
 tree->jet_num       = jet_num;
 tree->partonFlavour = j.partonFlavour();
 /////
 //   Take relevant info only for the first jet 
 /////
 //kinematic
 tree->jet_pt  = j.pt();
 tree->jet_eta = j.eta();
 tree->jet_phi = j.phi();
 tree->jet_en  = j.energy();
 //B tag prop
 double jet_csv_double = j.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags");
 tree->jet_csv = jet_csv_double;
 //Jet trks
 int jet_ndaus_int     = 0;
 int jet_chtrks_int    = 0;
 int jet_chtrkspv_int  = 0;
 int jet_chtrksnpv_int = 0;
 vector<Track> jetchtrks;
 vector<Track> jetchtrksnpv; 
 int jet_chtrkspvtt_int  = 0;
 int jet_chtrksnpvtt_int = 0;
 get_jettrks(j, PV, *ttrkbuilder, jet_ndaus_int, jet_chtrks_int, jet_chtrkspv_int, jet_chtrksnpv_int, jetchtrks, jetchtrksnpv, jet_chtrkspvtt_int, jet_chtrksnpvtt_int); 
 tree->jet_ndaus       = jet_ndaus_int;
 tree->jet_chtrks      = jet_chtrks_int;
 tree->jet_chtrkspv    = jet_chtrkspv_int;
 tree->jet_chtrksnpv   = jet_chtrksnpv_int;
 tree->jet_chtrkspvtt  = jet_chtrkspvtt_int;
 tree->jet_chtrksnpvtt = jet_chtrksnpvtt_int;
 /////Vertex compatibility 
 //chi2
 double jet_chi2tot_double  = -999;
 double jet_chi2ndf_double  = -1;
 double jet_chi2pval_double = -1;
 get_chi2info(jetchtrks, *ttrkbuilder, jet_chi2tot_double, jet_chi2ndf_double, jet_chi2pval_double);
 tree->jet_chi2tot  = jet_chi2tot_double;
 tree->jet_chi2ndf  = jet_chi2ndf_double;
 tree->jet_chi2pval = jet_chi2pval_double;
 //Two trk info
 double jet_num2v_double     = 0;
 double jet_numno2v_double   = 0;
 double jet_dca3d2t_double   = 0;
 double jet_dca3dno2t_double = 0;
 double jet_dca2d2t_double   = 0;
 double jet_dca2dno2t_double = 0;
 get_2trksinfo(jetchtrks, *ttrkbuilder, jet_num2v_double, jet_numno2v_double, jet_dca3d2t_double, jet_dca3dno2t_double, jet_dca2d2t_double, jet_dca2dno2t_double); 
 tree->jet_num2v       = jet_num2v_double;
 tree->jet_numno2v     = jet_numno2v_double;
 tree->jet_num2vno2v   = jet_num2v_double+jet_numno2v_double;
 ////Get the track IP info
 //Variables
 int track_3D_pos = 0;
 int track_2D_pos = 0;
 int track_1D_pos = 0;
 int track_num    = 0; 
 vector<pair<double,int> > trk_3D_IP;
 vector<pair<double, int>> trk_2D_IP;
 vector<pair<double, int>> trk_1D_IP;
 //Take the reco direction
 //GlobalVector gv3D(j.px(),j.py(),j.pz());
 //GlobalVector gv2D(j.px(),j.py(),0);
 //GlobalVector gv1D(0,0,j.pz());
 //Take the gen direction (from parton)
 const reco::GenParticle * jet_genparton = j.genParton(); 
 GlobalVector gv3D(jet_genparton->px(),jet_genparton->py(),jet_genparton->pz());
 GlobalVector gv2D(jet_genparton->px(),jet_genparton->py(),0);
 GlobalVector gv1D(0,0,jet_genparton->pz());
 //Take the gen direction (from jet)
 //const reco::GenJet * jet_genjet = j.genJet();
 //GlobalVector gv3D(jet_genjet->px(),jet_genjet->py(),jet_genjet->pz());
 //GlobalVector gv2D(jet_genjet->px(),jet_genjet->py(),0);
 //GlobalVector gv1D(0,0,jet_genjet->pz());


 double trk_IP3D_val  = -999;
 double trk_IP3D_sig  = -999;
 double trk_IP2D_val  = -999;
 double trk_IP2D_sig  = -999;
 double trk_sIP3D_val = -999;
 double trk_sIP3D_sig = -999;
 double trk_sIP2D_val = -999;
 double trk_sIP2D_sig = -999;
 double trk_IP1D_val  = -999;
 double trk_IP1D_sig  = -999;
 double trk_sIP1D_val = -999;
 double trk_sIP1D_sig = -999;
 double trk_IP3D_err  = -999;
 double trk_sIP3D_err = -999;
 double trk_IP2D_err  = -999;
 double trk_sIP2D_err = -999;
 double trk_IP1D_err  = -999;
 double trk_sIP1D_err = -999;
 //Get track 3D info
 for(uint i=0;i<jetchtrks.size();i++)
 {
  Track trk            = jetchtrks[i];
  TransientTrack ttrk  = ttrkbuilder->build(&trk);
  IPToolsValues3D(ttrk,PV,gv3D,trk_IP3D_val,trk_IP3D_sig,trk_sIP3D_val,trk_sIP3D_sig,trk_IP3D_err,trk_sIP3D_err);
  trk_3D_IP.push_back(make_pair(trk_IP3D_val,track_3D_pos)); 
  track_3D_pos++;
  track_num++;
 }
 if(track_num==0) return;
 sort(trk_3D_IP.rbegin(), trk_3D_IP.rend());
 for(uint k=0;k<jetchtrks.size();k++)
 {
  int track3D           = trk_3D_IP[k].second;
  Track trk3D           = jetchtrks[track3D];
  TransientTrack ttrk3D = ttrkbuilder->build(&trk3D);
  IPToolsValues3D(ttrk3D,PV,gv3D,trk_IP3D_val,trk_IP3D_sig,trk_sIP3D_val,trk_sIP3D_sig,trk_IP3D_err,trk_sIP3D_err);
  tree->trk_IP3D_val[k]  = trk_IP3D_val;
  tree->trk_IP3D_sig[k]  = trk_IP3D_sig;
  tree->trk_IP3D_err[k]  = trk_IP3D_err;
  tree->trk_sIP3D_val[k] = trk_sIP3D_val;
  tree->trk_sIP3D_sig[k] = trk_sIP3D_sig;
  tree->trk_sIP3D_err[k] = trk_sIP3D_err;
 } 
 //Get track 2D info 
 for(uint i=0;i<jetchtrks.size();i++)
 {
  Track trk2D            = jetchtrks[i];
  TransientTrack ttrk2D  = ttrkbuilder->build(&trk2D);
  IPToolsValues2D(ttrk2D,PV,gv2D,trk_IP2D_val,trk_IP2D_sig,trk_sIP2D_val,trk_sIP2D_sig,trk_IP2D_err,trk_sIP2D_err);
  trk_2D_IP.push_back(make_pair(trk_IP2D_val,track_2D_pos));
  track_2D_pos++;
 }
 sort(trk_2D_IP.rbegin(),trk_2D_IP.rend());
 for(uint k=0;k<jetchtrks.size();k++)
 {
  int track2D           = trk_2D_IP[k].second;
  Track trk2D           = jetchtrks[track2D];
  TransientTrack ttrk2D = ttrkbuilder->build(&trk2D);
  IPToolsValues2D(ttrk2D,PV,gv2D,trk_IP2D_val,trk_IP2D_sig,trk_sIP2D_val,trk_sIP2D_sig, trk_IP2D_err, trk_sIP2D_err);
  tree->trk_IP2D_val[k]  = trk_IP2D_val;
  tree->trk_IP2D_sig[k]  = trk_IP2D_sig;
  tree->trk_IP2D_err[k]  = trk_IP2D_err;
  tree->trk_sIP2D_val[k] = trk_sIP2D_val;
  tree->trk_sIP2D_sig[k] = trk_sIP2D_sig;
  tree->trk_sIP2D_err[k] = trk_sIP2D_err;
 }
 //Get track 1D info 
 for(uint i=0;i<jetchtrks.size();i++)
 {
  Track trk1D            = jetchtrks[i];
  TransientTrack ttrk1D  = ttrkbuilder->build(&trk1D);
  IPToolsValues1D(ttrk1D,PV,gv1D,trk_IP1D_val,trk_IP1D_sig,trk_sIP1D_val,trk_sIP1D_sig, trk_IP1D_err, trk_sIP1D_err);
  trk_1D_IP.push_back(make_pair(trk_IP1D_val,track_1D_pos));
  track_1D_pos++;
 }
 sort(trk_1D_IP.rbegin(),trk_1D_IP.rend());
 for(uint k=0;k<jetchtrks.size();k++)
 {
  int track1D           = trk_1D_IP[k].second;
  Track trk1D           = jetchtrks[track1D];
  TransientTrack ttrk1D = ttrkbuilder->build(&trk1D);
  IPToolsValues1D(ttrk1D,PV,gv1D,trk_IP1D_val,trk_IP1D_sig,trk_sIP1D_val,trk_sIP1D_sig, trk_IP1D_err, trk_sIP1D_err);
  tree->trk_IP1D_val[k]  = trk_IP1D_val;
  tree->trk_IP1D_sig[k]  = trk_IP1D_sig;
  tree->trk_IP1D_err[k]  = trk_IP1D_err;
  tree->trk_sIP1D_val[k] = trk_sIP1D_val;
  tree->trk_sIP1D_sig[k] = trk_sIP1D_sig;
  tree->trk_sIP1D_err[k] = trk_sIP1D_err;
 }
 /////
 //   Fill tree
 /////
 tree->tree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void 
BTagReco::beginJob()
{
 stopwatch->Start();
 evt_totnum = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
BTagReco::endJob() 
{
 stopwatch->Stop(); 
 cout<<endl;
 cout<<"Rapid job summary "<<endl;
 cout<<evt_totnum<<" events analysed in "<<stopwatch->RealTime()<<" seconds"<<endl;
 cout<<endl;
}
// ------------ method called when starting to processes a run  ------------
/*
void 
BTagReco::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
BTagReco::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
BTagReco::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
BTagReco::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BTagReco::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(BTagReco);
/*
 /////
 //   New variables defined refiting the primary vertex
 /////
 //Refiting the primary vertex with all the tracks     
 vector<Track> pvtrks;
 vector<TransientTrack> pvttrks;
 for(uint i=0; i<pfs->size(); i++){
  const pat::PackedCandidate & c = (*pfs)[i];
  if(c.charge()!=0 && c.numberOfHits()>0 && c.fromPV()==pat::PackedCandidate::PVUsedInFit){
   Track trk = Track(c.pseudoTrack());
   pvtrks.push_back(trk);
   TransientTrack ttrk = ttrkbuilder->build(&trk);
   pvttrks.push_back(ttrk);
  }
 }
 for(uint i=0; i<lostrks->size(); i++){
  const pat::PackedCandidate & c = (*lostrks)[i];
  if(c.charge()!=0 && c.numberOfHits()>0 && c.fromPV()==pat::PackedCandidate::PVUsedInFit){
   Track trk = Track(c.pseudoTrack());
   pvtrks.push_back(trk);
   TransientTrack ttrk = ttrkbuilder->build(&trk);
   pvttrks.push_back(ttrk);
  }
 }
 double Spvx=0;
 double Spvy=0;
 double Spvz=0;
 Spvx = PV.position().x();
 Spvy = PV.position().y();
 Spvz = PV.position().z();
 TransientVertex pv = vtxFitter.vertex(pvttrks);
 double Rpvx=0;
 double Rpvy=0;
 double Rpvz=0;
 Rpvx = pv.position().x();
 Rpvy = pv.position().y();
 Rpvz = pv.position().z();
 tree->PVxpvx = Spvx-Rpvx;
 tree->PVypvy = Spvy-Rpvy;
 tree->PVzpvz = Spvz-Rpvz;
 //Refiting the primary vertex with all the tracks, but the jet tracks     
 vector<Track> pvtrks_wt_jettrk;
 vector<TransientTrack> pvttrks_wt_jettrk;
 for(uint t=0; t<pvtrks.size(); t++){
  for(uint jt = 0; jt<jetchtrks.size(); jt++){
   if(pvtrks[t].pt() !=  jetchtrks[jt].pt()){
    pvtrks_wt_jettrk.push_back(pvtrks[t]);
    TransientTrack ttrk = ttrkbuilder->build(&pvtrks[t]);
    pvttrks_wt_jettrk.push_back(ttrk);
   }
  }
 }
 if(pvttrks_wt_jettrk.size()>2){
  TransientVertex pvnb = vtxFitter.vertex(pvttrks_wt_jettrk);
  if(pvnb.isValid()){
   double RpvX=0;
   double RpvY=0;
   double RpvZ=0;
   RpvX = pvnb.position().x();
   RpvY = pvnb.position().y();
   RpvZ = pvnb.position().z();
   tree->PVxnob = Spvx-RpvX;
   tree->PVynob = Spvy-RpvY;
   tree->PVznob = Spvz-RpvZ;
   tree->pvxnob = Rpvx-RpvX;
   tree->pvynob = Rpvy-RpvY;
   tree->pvznob = Rpvz-RpvZ;
  }
 }
*/
