//Jet and Met
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
//Math
#include <algorithm>
using namespace reco;
using namespace edm;
using namespace std;
/////
//   Lepton ID 
/////
//Ask candidate to come from PV
bool is_cand_inpv(const double mindR, const reco::Candidate& cand, const pat::PackedCandidateCollection& pcc, vector<int>& goodcandpfcindices){
 bool iscandinpv = false;
 double mindr = mindR;
 double goodcandpfcindex = -1;
 for(uint i=0; i<pcc.size(); i++){
  const pat::PackedCandidate cpf = pcc[i];
  if(deltaR(cand.p4(),cpf.p4())<mindr && //dR is the standard geometrical way to associate 
    (fabs(cand.pt()-cpf.pt())/cand.pt())<0.05 && //Check in pT, because ele,tau are usually faked by jets (many trks) and dR may not be enough
    cpf.charge()!=0 && cpf.numberOfHits()>0 && //Leptons are charged and built from tracks, also to be consistent with PV tracks  
    cpf.fromPV()==pat::PackedCandidate::PVUsedInFit //Coming from PV
    ){
   bool wasnotagoodcandindex = false; //Check that cpf was not already used for association
   for(uint vgi=0; vgi<goodcandpfcindices.size(); vgi++) if(goodcandpfcindices[vgi]==int(i)) wasnotagoodcandindex = true;
   if(!wasnotagoodcandindex){
    iscandinpv = true;
    mindr = deltaR(cand.p4(),cpf.p4());
    goodcandpfcindex = int(i);
   }
  }
 }
 if(goodcandpfcindex!=-1) goodcandpfcindices.push_back(goodcandpfcindex);
 return iscandinpv;
}
//Get transient track of a reco candidate
TransientTrack get_ttrk(const reco::Candidate* cand, const pat::PackedCandidateCollection& pcc, const TransientTrackBuilder& ttrkbuilder){
 Track trk;
 double mindr = 0.3;
 for(uint i=0; i<pcc.size(); i++){
  const pat::PackedCandidate cpf = pcc[i];
  if(deltaR(cand->p4(),cpf.p4())<mindr && (fabs(cand->pt()-cpf.pt())/cand->pt())<0.05){
   mindr = deltaR(cand->p4(),cpf.p4());
   trk = Track(cpf.pseudoTrack());
  }
 } 
 //cout<<"trk check "<<cand->pt()<<" "<<trk.pt()<<endl;
 TransientTrack ttrk = ttrkbuilder.build(&trk);
 //Wrong trk.charge()???
 //cout<<"check charge "<<cand->charge()<<" "<<trk.charge()<<" "<<ttrk.charge()<<endl;
 return ttrk;
}
/////
//   Relative isolation with delta-beta correction
/////
double rel_iso_dbc(const pat::Muon& lepton){
 return ( (lepton.chargedHadronIso() + 
           std::max(0.0, lepton.neutralHadronIso() + lepton.photonIso() - 0.5*lepton.puChargedHadronIso()))/lepton.pt() );
}
double rel_iso_dbc(const pat::Electron& lepton){
 return ( (lepton.chargedHadronIso() + 
           std::max(0.0, lepton.neutralHadronIso() + lepton.photonIso() - 0.5*lepton.puChargedHadronIso()))/lepton.pt() );
}
/////
//   Muon ID
/////
//https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopMUO
//lepton + jets / single-top
bool is_tight_muon(const pat::Muon& mu, const reco::Vertex& vtx){
 if(mu.track().isNull() || mu.globalTrack().isNull() || mu.muonBestTrack().isNull() || mu.innerTrack().isNull()) return false;
 return(
  mu.pt()>30 &&
  TMath::Abs(mu.eta()) < 2.1 &&
  mu.isPFMuon() &&
  mu.isGlobalMuon() &&
  mu.normChi2() < 10 &&
  mu.track()->hitPattern().trackerLayersWithMeasurement() > 5 &&
  mu.globalTrack()->hitPattern().numberOfValidMuonHits() > 0 &&
  mu.muonBestTrack()->dxy(vtx.position()) < 0.2 &&
  mu.muonBestTrack()->dz(vtx.position())< 0.5 &&
  mu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0 &&
  mu.numberOfMatchedStations() > 1 &&
  rel_iso_dbc(mu) < 0.12
 );
}
//https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopMUO
//dilepton
bool is_loose_muon(const pat::Muon& mu){
 return(
  mu.pt()>20 &&
  mu.isPFMuon() &&
  TMath::Abs(mu.eta()) < 2.4 &&
  (mu.isGlobalMuon() || mu.isTrackerMuon()) &&
  rel_iso_dbc(mu) < 0.2
 );
}
/////
//   Ele ID
/////
bool is_tight_electron(const pat::Electron& ele, const reco::Vertex& vtx){
 if(ele.gsfTrack().isNull()) return false;
 const float ae = TMath::Abs(ele.eta());
 return(
  ele.gsfTrack().isNonnull() &&
  ele.pt() > 30 && ae < 2.5 && !(ae>1.4442 && ae<1.5660) &&
  TMath::Abs(ele.gsfTrack()->dxy(vtx.position())) < 0.02 &&
  TMath::Abs(ele.gsfTrack()->dz(vtx.position())) < 0.2   &&
  ele.passConversionVeto() &&
  //throws error:
  //pat::Electron: the ID mvaTrigV0 can't be found in this pat::Electron.
  //The available IDs are: 'eidLoose' 'eidRobustHighEnergy' 'eidRobustLoose' 'eidRobustTight' 'eidTight'
  //ele.electronID("mvaTrigV0") > 0.5 &&
  ele.electronID("eidLoose") > 0.5 &&
  //ele.gsfTrack()->trackerExpectedHitsInner().numberOfHits() <= 0 &&
  rel_iso_dbc(ele) < 0.1
 );
}
bool is_loose_electron(const pat::Electron& ele, const reco::Vertex& vtx){
 if(ele.gsfTrack().isNull()) return false;
 const float ae = TMath::Abs(ele.eta());
 return(
  ele.gsfTrack().isNonnull() &&
  ele.pt() > 20 && ae < 2.5 &&
  TMath::Abs(ele.gsfTrack()->dxy(vtx.position())) < 0.04 &&
  TMath::Abs(ele.gsfTrack()->dz(vtx.position())) < 0.2   &&
  ele.passConversionVeto() &&
  //throws error:
  //pat::Electron: the ID mvaTrigV0 can't be found in this pat::Electron.
  //The available IDs are: 'eidLoose' 'eidRobustHighEnergy' 'eidRobustLoose' 'eidRobustTight' 'eidTight'
  //ele.electronID("mvaTrigV0") > 0.5 &&
  ele.electronID("eidLoose") > 0.5 &&
  //ele.gsfTrack()->trackerExpectedHitsInner().numberOfHits() <= 0 &&
  rel_iso_dbc(ele) < 0.15
);
}
/////
//   Jet id
/////
//Namespace to identify the pileup nature of a jet 
//working points from RecoJets/JetProducers/python/PileupJetIDCutParams_cfi.py
//Eta Categories  0-2.5 2.5-2.75 2.75-3.0 3.0-5.0
namespace pu_mva{
 float full_chs_loose[4][4] = {
  //full_5x_chs_wp Loose id
  {-0.98,-0.95,-0.94,-0.94}, //pt 0-10
  {-0.98,-0.95,-0.94,-0.94}, //pt 10-20
  {-0.89,-0.77,-0.69,-0.75}, //pt 20-30
  {-0.89,-0.77,-0.69,-0.57}, //pt 30-50
 };
 int eta_idx(float eta){//Note that in principle this is not needed as we always require jet to be within eta<2.4
  const float ae = TMath::Abs(eta);
  if(ae < 2.5){
   return 0;
  }else if(ae < 2.75){
   return 1;
  }else if(ae < 3.0){
   return 2;
  }else if(ae < 5.0) {
   return 3;
  }else{
   return 3;
  }
 }
 int pt_idx(float pt){
  if(pt < 10){
   return 0;
  }else if(pt < 20){
   return 1;
  }else if(pt < 30){
   return 2;
  }else if(pt < 50){
   return 3;
  }else{
   //edm::LogWarning("jet_pu_id") << "pt outside range " << pt;
   return -1;
  }
 }
 bool pass_id(const pat::Jet& x, float mva) {
  int vpt_idx  = pt_idx(x.pt());
  int veta_idx = eta_idx(x.eta());
  //hard jet
  if(vpt_idx==-1) return true;
  if(mva>full_chs_loose[vpt_idx][veta_idx]) return true;
  return false;
 }
}
bool is_good_jet(const pat::Jet &j){
 PFJetIDSelectionFunctor pfLooseJetID(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE),
                         pfTightJetID(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT);
 pat::strbitset passLooseCuts(pfLooseJetID.getBitTemplate()),
                passTightCuts(pfTightJetID.getBitTemplate());
 return( j.pt()>30 && fabs(j.eta())<2.5 &&
         pfLooseJetID(j,passLooseCuts) && pu_mva::pass_id(j, j.userFloat("pileupJetId:fullDiscriminant")) );
}
//Find jet associated to a b quark
void jettobquark(double& mindr, const reco::Candidate& cand, const pat::JetCollection& jets, vector<int>& goodjetindices, vector<const reco::Candidate*> leps){
 int goodjetindex = -1;
 int jet_pos = 0;
 for(const pat::Jet &j : jets){
  if(!is_good_jet(j)){jet_pos++; continue;}
  bool jetmatchedlepts = false;
  for(uint gl=0; gl<leps.size(); gl++) if(deltaR(leps[gl]->p4(),j.p4())<0.5) jetmatchedlepts = true;
  if(jetmatchedlepts){jet_pos++; continue;}
  if(deltaR(cand.p4(),j.p4())<mindr //dR is the standard geometrical way to associate 
     //&& (fabs(cand.pt()-j.pt())/cand.pt())<0.1 //Check in pT as well
    ){
   bool wasnotagoodcandindex = false; //Check that cpf was not already used for association
   for(uint vgi=0; vgi<goodjetindices.size(); vgi++) if(goodjetindices[vgi]==jet_pos) wasnotagoodcandindex = true;
   if(!wasnotagoodcandindex){
    mindr = deltaR(cand.p4(),j.p4());
    goodjetindex = jet_pos;
   }
  }
  jet_pos++;
 }
 goodjetindices.push_back(goodjetindex);
}
