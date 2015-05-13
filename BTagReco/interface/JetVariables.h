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
#include "TLorentzVector.h"
KalmanVertexFitter vtxFitter(true);
//Namespaces
using namespace reco;
using namespace edm;
using namespace std;
//Check that the track is a good track
bool is_goodtrk(Track trk,const reco::Vertex& vtx){
 bool isgoodtrk = false;
 if(trk.pt()>1 &&
   trk.hitPattern().numberOfValidHits()>=8 && 
   trk.hitPattern().numberOfValidPixelHits()>=2 &&
   trk.normalizedChi2()<5 &&
   std::abs(trk.dxy(vtx.position()))<0.2 &&
   std::abs(trk.dz(vtx.position()))<17    
   ) isgoodtrk = true;
 return isgoodtrk;
}
//Get transient tracks from track
TransientTrack get_ttrk(Track trk, const TransientTrackBuilder& ttrkbuilder){
 TransientTrack ttrk;
 ttrk = ttrkbuilder.build(&trk);
 return ttrk;
}
vector<TransientTrack> get_ttrks(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder){
 vector<TransientTrack> ttrks;
 for(uint tr=0; tr<trks.size(); tr++){
  TransientTrack ttrk = ttrkbuilder.build(&trks[tr]);
  ttrks.push_back(ttrk);
 }
 return ttrks;
}
//Get transient pv from tracks
TransientVertex get_tv(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder){
 TransientVertex tv;
 vector<TransientTrack> ttrks = get_ttrks(trks,ttrkbuilder);
 if(ttrks.size()>=2) tv = vtxFitter.vertex(ttrks);
 return tv;
}
TransientVertex get_ttv(vector<TransientTrack> ttrks, const TransientTrackBuilder& ttrkbuilder){
 TransientVertex tv;
 if(ttrks.size()>=2) tv = vtxFitter.vertex(ttrks);
 return tv;
}
//Get chi2 information from trk vertex
void get_chi2info(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, double& chi2tot, double& chi2ndf, double& chi2pval){
 if(trks.size()>=2){
  TransientVertex trks_tv = get_tv(trks, ttrkbuilder);
  if(trks_tv.isValid()){
   chi2tot  = trks_tv.totalChiSquared();
   chi2ndf  = trks_tv.degreesOfFreedom();
   chi2pval = TMath::Prob(trks_tv.totalChiSquared(),trks_tv.degreesOfFreedom());
  }
 }
}
//DCA between two trks
pair<double,double> dca2trks(Track tkA, Track tkB, const TransientTrackBuilder& ttrkbuilder){
 double dca3d2trks_sig = 0;
 double dca2d2trks_sig = 0;
 TransientTrack ttkA = get_ttrk(tkA, ttrkbuilder); 
 TransientTrack ttkB = get_ttrk(tkB, ttrkbuilder); 
 if(ttkA.impactPointTSCP().isValid() && ttkB.impactPointTSCP().isValid()){
  //Minimum distance
  FreeTrajectoryState state1 = ttkA.impactPointTSCP().theState();
  FreeTrajectoryState state2 = ttkB.impactPointTSCP().theState();
  TwoTrackMinimumDistance minDist;
  minDist.calculate(state1, state2);
  if(minDist.status()){
   //3D distance
   //const float dist3D = minDist.distance();
   std::pair<GlobalPoint,GlobalPoint> pcas = minDist.points();
   GlobalPoint pca1 = pcas.first;
   GlobalPoint pca2 = pcas.second;
   ROOT::Math::SVector<double, 3> distanceVector(pca1.x()-pca2.x(), pca1.y()-pca2.y(), pca1.z()-pca2.z());
   const float twoTkDis3D = ROOT::Math::Mag(distanceVector);
   //3D err distance
   float mass     = 0.139526;
   float massigma = mass*1e-6;
   float chi2 = 0.0f, ndf = 0.0f;
   KinematicParticleFactoryFromTransientTrack pFactory;
   RefCountedKinematicParticle tkAParticle = pFactory.particle(ttkA, mass, chi2, ndf, massigma);
   RefCountedKinematicParticle tkBParticle = pFactory.particle(ttkB, mass, chi2, ndf, massigma);
   float sig[6];
   sig[0] = tkAParticle->stateAtPoint(pca1).kinematicParametersError().matrix()(0,0);
   sig[1] = tkAParticle->stateAtPoint(pca1).kinematicParametersError().matrix()(0,1);
   sig[2] = tkAParticle->stateAtPoint(pca1).kinematicParametersError().matrix()(1,1);
   sig[3] = tkAParticle->stateAtPoint(pca1).kinematicParametersError().matrix()(0,2);
   sig[4] = tkAParticle->stateAtPoint(pca1).kinematicParametersError().matrix()(1,2);
   sig[5] = tkAParticle->stateAtPoint(pca1).kinematicParametersError().matrix()(2,2);
   ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > pca1Cov(sig, sig+6);
   sig[0] = tkBParticle->stateAtPoint(pca2).kinematicParametersError().matrix()(0,0);
   sig[1] = tkBParticle->stateAtPoint(pca2).kinematicParametersError().matrix()(0,1);
   sig[2] = tkBParticle->stateAtPoint(pca2).kinematicParametersError().matrix()(1,1);
   sig[3] = tkBParticle->stateAtPoint(pca2).kinematicParametersError().matrix()(0,2);
   sig[4] = tkBParticle->stateAtPoint(pca2).kinematicParametersError().matrix()(1,2);
   sig[5] = tkBParticle->stateAtPoint(pca2).kinematicParametersError().matrix()(2,2);
   ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > pca2Cov(sig, sig+6);
   ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > totCov = pca1Cov + pca2Cov;
   double twoTkDist3DEr = TMath::Sqrt(ROOT::Math::Similarity(totCov, distanceVector)) / twoTkDis3D;
   dca3d2trks_sig = twoTkDis3D/twoTkDist3DEr;
   //2D distance and err distance
   distanceVector(2) = 0.0;
   double twoTauDis2D    = ROOT::Math::Mag(distanceVector);
   double twoTauDist2DEr = TMath::Sqrt(ROOT::Math::Similarity(totCov, distanceVector)) / twoTauDis2D;
   dca2d2trks_sig = twoTauDis2D/twoTauDist2DEr;
  }//if(minDist.status())
 }//ttkA.impactPointTSCP().isValid() && ttkB.impactPointTSCP().isValid()
 return make_pair(dca3d2trks_sig,dca2d2trks_sig);
}
//Get Num of two-trk vertices
void get_2trksinfo(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, double& num2v, double& numno2v, double& dca3d2t, double& dca3dno2t, double& dca2d2t, double& dca2dno2t){
 for(uint t=0; t<trks.size(); t++){
  for(uint t2=t+1; t2<trks.size(); t2++){
   vector<Track> twotrks;
   twotrks.push_back(trks[t]);  
   twotrks.push_back(trks[t2]);    
   TransientVertex tv = get_tv(twotrks, ttrkbuilder);
   if(tv.isValid() && TMath::Prob(tv.totalChiSquared(),tv.degreesOfFreedom())>0.05){
    num2v++;
    pair<double,double> dca2trks3d2d = dca2trks(trks[t], trks[t2], ttrkbuilder);
    dca3d2t += dca2trks3d2d.first;
    dca2d2t += dca2trks3d2d.second;
   }else{
    numno2v++;
    pair<double,double> dca2trks3d2d = dca2trks(trks[t], trks[t2], ttrkbuilder);
    dca3dno2t += dca2trks3d2d.first;
    dca2dno2t += dca2trks3d2d.second;
   }
  }
 }
}
//Get jet tracks
void get_jettrks(const pat::Jet& jet, const reco::Vertex& vtx, const TransientTrackBuilder& ttrkbuilder, int& jet_ndaus, int& jet_chtrks, int& jet_chtrkspv, int& jet_chtrksnpv, vector<Track>& jetchtrks, vector<Track>& jetchtrksnpv, int& jet_chtrkspvtt, int& jet_chtrksnpvtt){
 //Access jet daughters
 vector<CandidatePtr> jdaus(jet.daughterPtrVector());
 sort(jdaus.begin(), jdaus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) {return p1->pt() > p2->pt();});
 jet_ndaus = jdaus.size(); 
 for(int jd=0; jd<jet_ndaus; jd++){
  const pat::PackedCandidate &jcand = dynamic_cast<const pat::PackedCandidate &>(*jdaus[jd]);
  if(deltaR(jcand.p4(),jet.p4())>0.4) continue;
  Track trk = Track(jcand.pseudoTrack());
  bool isgoodtrk = is_goodtrk(trk,vtx);
  //Minimal conditions for a track 
  //if(jcand.charge()!=0 && jcand.numberOfHits()>0)
  if(isgoodtrk && jcand.charge()!=0 && jcand.fromPV()>1)
  {
   jet_chtrks++; 
   jetchtrks.push_back(trk);
   //Other conditions on jet daughters
   //Using fromPV method
   if(jcand.fromPV()==pat::PackedCandidate::PVUsedInFit){
    jet_chtrkspv++;
   }else{
    jet_chtrksnpv++;
    jetchtrksnpv.push_back(trk);
   }
  }//Ch trks 
 }//Loop on jet daus 
}
//Get info for IP between a reco candidate and a tv
void get_ipinfo_llepjettv(TransientTrack ttrk, TransientVertex tv, GlobalVector gv,
 double& aip3Dsig, double& sip3Dsig, double& aip2Dsig, double& sip2Dsig,
 double& aip3Dsig_val, double& sip3Dsig_val, double& aip2Dsig_val, double& sip2Dsig_val,
 double& aip3Dsig_err, double& sip3Dsig_err, double& aip2Dsig_err, double& sip2Dsig_err
){
 pair<bool, Measurement1D> currIP = IPTools::absoluteImpactParameter3D(ttrk,tv);
 if(currIP.first){
  aip3Dsig     = currIP.second.value()/currIP.second.error();
  aip3Dsig_val = currIP.second.value();
  aip3Dsig_err = currIP.second.error();
 } 
 currIP = IPTools::signedImpactParameter3D(ttrk,gv,tv);
 if(currIP.first){
  sip3Dsig     = currIP.second.value()/currIP.second.error();
  sip3Dsig_val = currIP.second.value();
  sip3Dsig_err = currIP.second.error();
 }
 currIP = IPTools::absoluteTransverseImpactParameter(ttrk,tv);
 if(currIP.first){
  aip2Dsig     = currIP.second.value()/currIP.second.error();
  aip2Dsig_val = currIP.second.value();
  aip2Dsig_err = currIP.second.error();
 }
 currIP = IPTools::signedTransverseImpactParameter(ttrk,gv,tv);
 if(currIP.first){
  sip2Dsig     = currIP.second.value()/currIP.second.error();
  sip2Dsig_val = currIP.second.value();
  sip2Dsig_err = currIP.second.error();
 }
}
