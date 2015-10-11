//New class
#ifndef TREEVARIABLES
#define TREEVARIABLES
/////
//   Headers
/////
#include <TTree.h>
#include <string>
#include <map>
#include <cmath>
#include <iostream>
/////
//   Constants
/////
#define DEF_SIZE1D 100
#define DEF_VAL_INT -9999
#define DEF_VAL_FLOAT -9999.0f
#define DEF_VAL_DOUBLE -9999.0d
#define FLOAT_EPS 0.0000001f
#define DOUBLE_EPS 0.0000001d
/////
//   Functions
/////
#define INIT_1DARRAY(x,n,y) for(int i=0;i<n;i++) {x[i]=y;}
#define INIT_2DARRAY(x,n,m,y) for(int i=0;i<n;i++) { for(int j=0;j<m;j++) { x[i][j]=y; } }
inline bool is_undef(int x) { return x==DEF_VAL_INT; };
inline bool is_undef(float x) { return fabs(x-DEF_VAL_FLOAT) < FLOAT_EPS; };
inline bool is_undef(double x) { return fabs(x-DEF_VAL_DOUBLE) < DOUBLE_EPS; }
/////
//   Class declaration
/////
class CTree{
public:
 CTree(TTree* _tree) { tree = _tree; };
 TTree* tree;
 /////
 //   Helper functions for accessing branches
 /////
 template <typename T>
 T get_address(const std::string name){
  auto* br = tree->GetBranch(name.c_str());
  if(br==0){
   std::cerr << "ERROR: get_address CTree " << "branch " << name << " does not exist" << std::endl;
   throw std::exception();
  }
  auto* p = br->GetAddress();
  return reinterpret_cast<T>(p);
 }
 /////
 //   Declare variables
 /////
 double sumhcsvjets, prodhcsvjets;
 double duno;
 int ngenb;
 int lep_numl, lep_numt, lep_dichprodl, lep_dichprodt;
 int jet_num, jet_ngenb, jet_ngenbdr;
 double jet_pt, jet_eta, jet_phi, jet_en;
 double jet_csv;
 double partonFlavour;
 double partonFlavourdR;
 double deltaR_reco_gen,deltaR_reco_jet, deltaR_gen_jet;
 //Num trk info
 double jet_ndaus, jet_chtrks, jet_chtrkspv, jet_chtrksnpv, jet_chtrkspvtt, jet_chtrksnpvtt;
 //Chi2 info
 double jet_chi2tot, jet_chi2ndf, jet_chi2pval;
 //Two trk info
 double RV_kvf_ntrks,RV_nojet_ntrks;
 //Vertex info
 double RV_kvf_3D,RV_avf_bs_3D,RV_avf_3D,RV_kvf_bs_3D,RV_RVnojet_3D,RV_RVaddjtrcks_NPV_3D;
 double jet_num2v, jet_numno2v, jet_num2vno2v, jet_dca3d2t, jet_dca3dno2t, jet_dca3d2tno2t, jet_dca2d2t, jet_dca2dno2t, jet_dca2d2tno2t;
 //Primary Vertex Info
 double diff_PVx_RVx_kvf, diff_PVy_RVy_kvf, diff_PVz_RVz_kvf,diff_PVx_RVx_kvf_bs, diff_PVy_RVy_kvf_bs, diff_PVz_RVz_kvf_bs,diff_PVx_RVx_avf, diff_PVy_RVy_avf, diff_PVz_RVz_avf,diff_PVx_RVx_avf_bs, diff_PVy_RVy_avf_bs, diff_PVz_RVz_avf_bs,PVx, PVy, PVz, PV_nojet_x, PV_nojet_y, PV_nojet_z, diff_RV_RVnob_x, diff_RV_RVnob_y, diff_RV_RVnob_z, pv_ntrks, RVnojet_chi2_ndf, diff_pv_RVnojet_chi2_ndf, RV_nchi2_ndf, diff_RV_RVnojet_chi2_ndf;
 // Get the IP info of the Jet tracks
double trk_IP3D_val_AE[DEF_SIZE1D], trk_IP3D_val_TE[DEF_SIZE1D], trk_IP3D_sig_AE[DEF_SIZE1D], trk_IP3D_sig_TE[DEF_SIZE1D],track_IP3D_val_AE[DEF_SIZE1D], track_IP3D_val_TE[DEF_SIZE1D], track_IP3D_sig_AE[DEF_SIZE1D], track_IP3D_sig_TE[DEF_SIZE1D],track_pt[DEF_SIZE1D];
double PV_x,PV_y,PV_z,trk_sIP1D_val_x[DEF_SIZE1D], trk_sIP1D_val_y[DEF_SIZE1D] ,trk_sIP1D_val_z[DEF_SIZE1D], diff_RP_PV_x[DEF_SIZE1D], diff_RP_PV_y[DEF_SIZE1D], diff_RP_PV_z[DEF_SIZE1D], IPVec_PV_x[DEF_SIZE1D],  IPVec_PV_y[DEF_SIZE1D],  IPVec_PV_z[DEF_SIZE1D],diff_chi2_ndf_RV_addjettrks_RVnob;
//IP variables 
double trk_IP3D_val_RV_jetV[DEF_SIZE1D], trk_IP3D_sig_RV_jetV[DEF_SIZE1D], trk_IP2D_val_RV_jetV[DEF_SIZE1D], trk_IP2D_sig_RV_jetV[DEF_SIZE1D],trk_sIP3D_val_RV_jetV[DEF_SIZE1D], trk_sIP3D_sig_RV_jetV[DEF_SIZE1D], trk_sIP2D_val_RV_jetV[DEF_SIZE1D], trk_sIP2D_sig_RV_jetV[DEF_SIZE1D],trk_IP3D_err_RV_jetV[DEF_SIZE1D], trk_sIP3D_err_RV_jetV[DEF_SIZE1D], trk_IP2D_err_RV_jetV[DEF_SIZE1D], trk_sIP2D_err_RV_jetV[DEF_SIZE1D];
double RV_x,RV_y,RV_z,RVnojtrk_x,RVnojtrk_y,RVnojtrk_z;
double trk_IP3D_val[DEF_SIZE1D], trk_IP3D_sig[DEF_SIZE1D], trk_IP2D_val[DEF_SIZE1D], trk_IP2D_sig[DEF_SIZE1D], trk_IP1D_val[DEF_SIZE1D],trk_IP1D_sig[DEF_SIZE1D], trk_sIP3D_val[DEF_SIZE1D], trk_sIP3D_sig[DEF_SIZE1D], trk_sIP2D_val[DEF_SIZE1D], trk_sIP2D_sig[DEF_SIZE1D],trk_sIP1D_val[DEF_SIZE1D], trk_sIP1D_sig[DEF_SIZE1D], trk_IP3D_err[DEF_SIZE1D], trk_sIP3D_err[DEF_SIZE1D], trk_IP2D_err[DEF_SIZE1D], trk_sIP2D_err[DEF_SIZE1D],trk_IP1D_err[DEF_SIZE1D],trk_sIP1D_err[DEF_SIZE1D];
///// ////   Initialise
//Decay Length and Flight Distance
double FD,FD3D_val,FD3D_sig,FD2D_val,FD2D_sig,DL1D_val[DEF_SIZE1D],DL1D_sig[DEF_SIZE1D],sDL3D_val[DEF_SIZE1D],sDL3D_sig[DEF_SIZE1D],sDL2D_val[DEF_SIZE1D],sDL2D_sig[DEF_SIZE1D],sDL1D_val[DEF_SIZE1D],sDL1D_sig[DEF_SIZE1D],DL3D_val[DEF_SIZE1D],DL3D_sig[DEF_SIZE1D],DL2D_val[DEF_SIZE1D],DL2D_sig[DEF_SIZE1D]; 
///////
 void loop_initialize(void){
  duno        = 1.;
  ngenb           = DEF_VAL_INT;
  lep_numl        = DEF_VAL_INT;
  lep_numt        = DEF_VAL_INT;
  lep_dichprodl   = DEF_VAL_INT;
  lep_dichprodt   = DEF_VAL_INT;
  jet_num         = DEF_VAL_INT;
  jet_ngenb       = DEF_VAL_INT;
  jet_ngenbdr     = DEF_VAL_INT;
  jet_pt          = DEF_VAL_DOUBLE;
  jet_eta         = DEF_VAL_DOUBLE;
  jet_phi         = DEF_VAL_DOUBLE;
  jet_en          = DEF_VAL_DOUBLE;
  jet_csv         = DEF_VAL_DOUBLE;
  partonFlavour   = DEF_VAL_DOUBLE;
  partonFlavourdR = DEF_VAL_DOUBLE;
  deltaR_reco_gen = DEF_VAL_DOUBLE;
  deltaR_reco_jet = DEF_VAL_DOUBLE;
  deltaR_gen_jet  = DEF_VAL_DOUBLE;
  //Num trk info
  jet_ndaus       = DEF_VAL_DOUBLE;
  jet_chtrks      = DEF_VAL_DOUBLE;
  jet_chtrkspv    = DEF_VAL_DOUBLE;
  jet_chtrksnpv   = DEF_VAL_DOUBLE;
  jet_chtrkspvtt  = DEF_VAL_DOUBLE;
  jet_chtrksnpvtt = DEF_VAL_DOUBLE;
  //Chi2 info
  jet_chi2tot     = DEF_VAL_DOUBLE;
  jet_chi2ndf     = DEF_VAL_DOUBLE;
  jet_chi2pval    = DEF_VAL_DOUBLE;
  //Vertex info 
  RV_kvf_bs_3D           = DEF_VAL_DOUBLE;
  RV_avf_3D              = DEF_VAL_DOUBLE;
  RV_avf_bs_3D           = DEF_VAL_DOUBLE;
  RV_kvf_3D              = DEF_VAL_DOUBLE;
  RV_RVaddjtrcks_NPV_3D  = DEF_VAL_DOUBLE;
  RV_RVnojet_3D          = DEF_VAL_DOUBLE;
  //Two trk info
  jet_num2v       = DEF_VAL_DOUBLE;
  jet_numno2v     = DEF_VAL_DOUBLE;
  jet_num2vno2v   = DEF_VAL_DOUBLE;
  jet_dca3d2t     = DEF_VAL_DOUBLE;
  jet_dca3dno2t   = DEF_VAL_DOUBLE;
  jet_dca2d2t     = DEF_VAL_DOUBLE;
  jet_dca2dno2t   = DEF_VAL_DOUBLE;
  jet_dca3d2tno2t = DEF_VAL_DOUBLE;
  jet_dca2d2tno2t = DEF_VAL_DOUBLE;
  RV_kvf_ntrks    = DEF_VAL_DOUBLE;
  RV_nojet_ntrks  = DEF_VAL_DOUBLE;
  //Flight Distance
  FD3D_val        = DEF_VAL_DOUBLE;
  FD2D_val        = DEF_VAL_DOUBLE;
  FD3D_sig        = DEF_VAL_DOUBLE;
  FD2D_sig        = DEF_VAL_DOUBLE;
 // Primary Vertex Information 
  diff_PVx_RVx_kvf          = -9999.;
  diff_PVy_RVy_kvf          = -9999.;
  diff_PVz_RVz_kvf          = -9999.;
  diff_PVx_RVx_kvf_bs       = -9999.;
  diff_PVy_RVy_kvf_bs       = -9999.;
  diff_PVz_RVz_kvf_bs       = -9999.;
  diff_PVx_RVx_avf          = -9999.;
  diff_PVy_RVy_avf          = -9999.;
  diff_PVz_RVz_avf          = -9999.;
  diff_PVx_RVx_avf_bs       = -9999.;
  diff_PVy_RVy_avf_bs       = -9999.;
  diff_PVz_RVz_avf_bs       = -9999.;
  FD                        = DEF_VAL_DOUBLE;
  PVy                       = DEF_VAL_DOUBLE;
  PVx                       = DEF_VAL_DOUBLE;
  PVz                       = DEF_VAL_DOUBLE;
  PV_nojet_x                = DEF_VAL_DOUBLE;
  PV_nojet_y                = DEF_VAL_DOUBLE;
  PV_nojet_z                = DEF_VAL_DOUBLE;
  diff_RV_RVnob_x           = DEF_VAL_DOUBLE;
  diff_RV_RVnob_y           = DEF_VAL_DOUBLE;
  diff_RV_RVnob_z           = DEF_VAL_DOUBLE;
  pv_ntrks                  = DEF_VAL_DOUBLE;
  RVnojet_chi2_ndf          = DEF_VAL_DOUBLE;
  diff_pv_RVnojet_chi2_ndf  = DEF_VAL_DOUBLE;
  RV_nchi2_ndf              = DEF_VAL_DOUBLE;
  diff_RV_RVnojet_chi2_ndf  = DEF_VAL_DOUBLE;
  PV_x                      = DEF_VAL_DOUBLE;
  PV_y                      = DEF_VAL_DOUBLE;
  PV_z                      = DEF_VAL_DOUBLE;
  RV_x                      = DEF_VAL_DOUBLE;
  RV_y                      = DEF_VAL_DOUBLE;
  RV_z                      = DEF_VAL_DOUBLE;
  RVnojtrk_x                = DEF_VAL_DOUBLE;
  RVnojtrk_y                = DEF_VAL_DOUBLE;
  RVnojtrk_z                = DEF_VAL_DOUBLE;   
 // Get the IP info of the Jet Tracks
  diff_chi2_ndf_RV_addjettrks_RVnob = DEF_VAL_DOUBLE;
  INIT_1DARRAY(track_pt,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_IP3D_val,DEF_SIZE1D,DEF_VAL_DOUBLE); 
  INIT_1DARRAY(trk_IP3D_sig,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_IP2D_val,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_IP2D_sig,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_IP3D_val_RV_jetV,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_IP3D_sig_RV_jetV,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_IP2D_val_RV_jetV,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_IP2D_sig_RV_jetV,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_IP1D_val,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_IP1D_sig,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_sIP3D_val,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_sIP3D_sig,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_sIP3D_val_RV_jetV,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_sIP3D_sig_RV_jetV,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_sIP2D_val,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_sIP2D_sig,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_sIP2D_val_RV_jetV,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_sIP2D_sig_RV_jetV,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_sIP1D_val,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_sIP1D_sig,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_IP3D_err,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_sIP3D_err,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_IP2D_err,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_sIP2D_err,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_IP3D_err_RV_jetV,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_sIP3D_err_RV_jetV,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_IP2D_err_RV_jetV,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_sIP2D_err_RV_jetV,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_IP1D_err,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_sIP1D_err,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_sIP1D_val_x,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_sIP1D_val_y,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(trk_sIP1D_val_z,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(diff_RP_PV_x,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(diff_RP_PV_y,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(diff_RP_PV_z,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(IPVec_PV_x,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(IPVec_PV_y,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(IPVec_PV_z,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(DL1D_val,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(DL1D_sig,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(DL2D_val,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(DL2D_sig,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(DL3D_val,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(DL3D_sig,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(sDL3D_val,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(sDL3D_val,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(sDL3D_sig,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(sDL2D_val,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(sDL2D_sig,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(sDL1D_val,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(sDL1D_sig,DEF_SIZE1D,DEF_VAL_DOUBLE);
  } 
 /////
 //   Set branches
 /////
 void make_branches(void){
  tree->Branch("duno", &duno, "duno/D");
  tree->Branch("sumhcsvjets", &sumhcsvjets, "sumhcsvjets/D");
  tree->Branch("prodhcsvjets", &prodhcsvjets, "prodhcsvjets/D");
  tree->Branch("ngenb", &ngenb, "ngenb/I");
  tree->Branch("lep_numl", &lep_numl, "lep_numl/I");
  tree->Branch("lep_numt", &lep_numt, "lep_numt/I");
  tree->Branch("lep_dichprodl", &lep_dichprodl, "lep_dichprodl/I");
  tree->Branch("lep_dichprodt", &lep_dichprodt, "lep_dichprodt/I");
  tree->Branch("jet_num", &jet_num, "jet_num/I");
  tree->Branch("jet_ngenb", &jet_ngenb, "jet_ngenb/I");
  tree->Branch("jet_ngenbdr", &jet_ngenbdr, "jet_ngenbdr/I");
  tree->Branch("jet_pt", &jet_pt, "jet_pt/D");
  tree->Branch("jet_eta", &jet_eta, "jet_eta/D");
  tree->Branch("jet_phi", &jet_phi, "jet_phi/D");
  tree->Branch("jet_en", &jet_en, "jet_en/D");
  tree->Branch("jet_csv", &jet_csv, "jet_csv/D");
  tree->Branch("partonFlavour", &partonFlavour, "partonFlavour/D");
  tree->Branch("partonFlavourdR", &partonFlavourdR, "partonFlavourdR/D");
  tree->Branch("deltaR_reco_gen", &deltaR_reco_gen, "deltaR_reco_gen/D");
  tree->Branch("deltaR_reco_jet", &deltaR_reco_jet, "deltaR_reco_jet/D");
  tree->Branch("deltaR_gen_jet", &deltaR_gen_jet, "deltaR_gen_jet/D");
  //Num trk info
  tree->Branch("jet_ndaus", &jet_ndaus, "jet_ndaus/D");
  tree->Branch("jet_chtrks", &jet_chtrks, "jet_chtrks/D");
  tree->Branch("jet_chtrkspv", &jet_chtrkspv, "jet_chtrkspv/D");
  tree->Branch("jet_chtrksnpv", &jet_chtrksnpv, "jet_chtrksnpv/D");  
  tree->Branch("jet_chtrkspvtt", &jet_chtrkspvtt, "jet_chtrkspvtt/D");
  tree->Branch("jet_chtrksnpvtt", &jet_chtrksnpvtt, "jet_chtrksnpvtt/D");
  //Chi2 info
  tree->Branch("jet_chi2tot", &jet_chi2tot, "jet_chi2tot/D");
  tree->Branch("jet_chi2ndf", &jet_chi2ndf, "jet_chi2ndf/D");
  tree->Branch("jet_chi2pval", &jet_chi2pval, "jet_chi2pval/D");
  //Two trk info
  tree->Branch("jet_num2v", &jet_num2v, "jet_num2v/D");
  tree->Branch("jet_numno2v", &jet_numno2v, "jet_numno2v/D");
  tree->Branch("jet_num2vno2v", &jet_num2vno2v, "jet_num2vno2v/D");
  tree->Branch("jet_dca3d2t", &jet_dca3d2t, "jet_dca3d2t/D");
  tree->Branch("jet_dca3dno2t", &jet_dca3dno2t, "jet_dca3dno2t/D");
  tree->Branch("jet_dca2d2t", &jet_dca2d2t, "jet_dca2d2t/D");
  tree->Branch("jet_dca2dno2t", &jet_dca2dno2t, "jet_dca2dno2t/D");
  tree->Branch("jet_dca3d2tno2t", &jet_dca3d2tno2t, "jet_dca3d2tno2t/D");
  tree->Branch("jet_dca2d2tno2t", &jet_dca2d2tno2t, "jet_dca2d2tno2t/D");
  //Primary Vertex information
  tree->Branch("RV_kvf_bs_3D",&RV_kvf_bs_3D,"RV_kvf_bs_3D/D");
  tree->Branch("RV_kvf_3D",&RV_kvf_3D,"RV_kvf_3D/D");
  tree->Branch("RV_avf_3D",&RV_avf_3D,"RV_avf_3D/D");
  tree->Branch("RV_avf_bs_3D",&RV_avf_bs_3D,"RV_avf_bs_3D/D");
  tree->Branch("RV_RVnojet_3D",&RV_RVnojet_3D,"RV_RVnojet_3D/D");
  tree->Branch("RV_RVaddjtrcks_NPV_3D",&RV_RVaddjtrcks_NPV_3D,"RV_RVaddjtrcks_NPV_3D/D");
  tree->Branch("diff_PVx_RVx_kvf",&diff_PVx_RVx_kvf,"diff_PVx_RVx_kvf/D");
  tree->Branch("diff_PVy_RVy_kvf",&diff_PVy_RVy_kvf,"diff_PVy_RVy_kvf/D");
  tree->Branch("diff_PVz_RVz_kvf",&diff_PVz_RVz_kvf,"diff_PVz_RVz_kvf/D");
  tree->Branch("diff_PVx_RVx_kvf_bs",&diff_PVx_RVx_kvf_bs,"diff_PVx_RVx_kvf_bs/D");
  tree->Branch("diff_PVy_RVy_kvf_bs",&diff_PVy_RVy_kvf_bs,"diff_PVy_RVy_kvf_bs/D");
  tree->Branch("diff_PVz_RVz_kvf_bs",&diff_PVz_RVz_kvf_bs,"diff_PVz_RVz_kvf_bs/D");
  tree->Branch("diff_PVx_RVx_avf",&diff_PVx_RVx_avf,"diff_PVx_RVx_kvf/D");
  tree->Branch("diff_PVy_RVy_avf",&diff_PVy_RVy_avf,"diff_PVy_RVy_kvf/D");
  tree->Branch("diff_PVz_RVz_avf",&diff_PVz_RVz_avf,"diff_PVz_RVz_kvf/D");
  tree->Branch("diff_PVx_RVx_avf_bs",&diff_PVx_RVx_avf_bs,"diff_PVx_RVx_avf_bs/D");
  tree->Branch("diff_PVy_RVy_avf_bs",&diff_PVy_RVy_avf_bs,"diff_PVy_RVy_avf_bs/D");
  tree->Branch("diff_PVz_RVz_avf_bs",&diff_PVz_RVz_avf_bs,"diff_PVz_RVz_avf_bs/D");
  tree->Branch("PVx", &PVx ,"PVx/D");
  tree->Branch("PVy", &PVy ,"PVy/D");
  tree->Branch("PVz", &PVz ,"PVz/D");
  tree->Branch("PV_nojet_x", &PV_nojet_x , "PV_nojet_x/D");
  tree->Branch("PV_nojet_y", &PV_nojet_y , "PV_nojet_y/D");
  tree->Branch("PV_nojet_z", &PV_nojet_z , "PV_nojet_z/D");
  tree->Branch("diff_RV_RVnob_x", &diff_RV_RVnob_x ,"diff_RV_RVnob_x/D");
  tree->Branch("diff_RV_RVnob_y", &diff_RV_RVnob_y ,"diff_RV_RVnob_y/D");
  tree->Branch("diff_RV_RVnob_z", &diff_RV_RVnob_z ,"diff_RV_RVnob_z/D");
  tree->Branch("pv_ntrks", &pv_ntrks , "pv_ntrks/D");
  tree->Branch("RVnojet_chi2_ndf", &RVnojet_chi2_ndf , "RVnojet_chi2_ndf/D");
  tree->Branch("diff_pv_RVnojet_chi2_ndf", &diff_pv_RVnojet_chi2_ndf , "diff_pv_RVnojet_chi2_ndf/D");
  tree->Branch("RV_nchi2_ndf", &RV_nchi2_ndf , "RV_nchi2_ndf/D");
  tree->Branch("diff_RV_RVnojet_chi2_ndf", &diff_RV_RVnojet_chi2_ndf , "diff_RV_RVnojet_chi2_ndf/D");
  tree->Branch("diff_chi2_ndf_RV_addjettrks_RVnob", &diff_chi2_ndf_RV_addjettrks_RVnob, "diff_chi2_ndf_RV_addjettrks_RVnob/D");
  tree->Branch("RV_x",&RV_x,"RV_x/D");
  tree->Branch("RV_y",&RV_y,"RV_y/D");
  tree->Branch("RV_z",&RV_z,"RV_z/D");
  tree->Branch("RVnojtrk_x",&RVnojtrk_x,"RVnojtrk_x/D");
  tree->Branch("RVnojtrk_y",&RVnojtrk_y,"RVnojtrk_y/D");
  tree->Branch("RVnojtrk_z",&RVnojtrk_z,"RVnojtrk_z/D");
  tree->Branch("RV_kvf_ntrks",&RV_kvf_ntrks,"RV_kvf_ntrks/D");
  tree->Branch("RV_nojet_ntrks",&RV_nojet_ntrks,"RV_nojet_ntrks/D");
  //Get the IP info of the tracks
  tree->Branch("trk_IP3D_val", &trk_IP3D_val, "trk_IP3D_val[25]/D");
  tree->Branch("trk_IP3D_sig", &trk_IP3D_sig, "trk_IP3D_sig[25]/D");
  tree->Branch("trk_IP3D_val_RV_jetV", &trk_IP3D_val_RV_jetV, "trk_IP3D_val_RV_jetV[25]/D");
  tree->Branch("trk_IP3D_sig_RV_jetV", &trk_IP3D_sig_RV_jetV, "trk_IP3D_sig[25]_RV_jetV/D");
  tree->Branch("trk_IP2D_val", &trk_IP2D_val, "trk_IP2D_val[25]/D");
  tree->Branch("trk_IP2D_sig", &trk_IP2D_sig, "trk_IP2D_sig[25]/D");
  tree->Branch("trk_IP2D_val_RV_jetV", &trk_IP2D_val_RV_jetV, "trk_IP2D_val_RV_jetV[25]/D");
  tree->Branch("trk_IP2D_sig_RV_jetV", &trk_IP2D_sig_RV_jetV, "trk_IP2D_sig_RV_jetV[25]/D");
  tree->Branch("trk_IP1D_val", &trk_IP1D_val, "trk_IP1D_val[25]/D");
  tree->Branch("trk_IP1D_sig", &trk_IP1D_sig, "trk_IP1D_sig[25]/D");
  tree->Branch("trk_sIP3D_val", &trk_sIP3D_val, "trk_sIP3D_val[25]/D");
  tree->Branch("trk_sIP3D_sig", &trk_sIP3D_sig, "trk_sIP3D_sig[25]/D");
  tree->Branch("trk_sIP3D_val_RV_jetV", &trk_sIP3D_val_RV_jetV, "trk_sIP3D_val_RV_jetV[25]/D");
  tree->Branch("trk_sIP3D_sig_RV_jetV", &trk_sIP3D_sig_RV_jetV, "trk_sIP3D_sig_RV_jetV[25]/D");
  tree->Branch("trk_sIP2D_val", &trk_sIP2D_val, "trk_sIP2D_val[25]/D");
  tree->Branch("trk_sIP2D_sig", &trk_sIP2D_sig, "trk_sIP2D_sig[25]/D");
  tree->Branch("trk_sIP2D_val_RV_jetV", &trk_sIP2D_val_RV_jetV, "trk_sIP2D_val_RV_jetV[25]/D");
  tree->Branch("trk_sIP2D_sig_RV_jetV", &trk_sIP2D_sig_RV_jetV, "trk_sIP2D_sig_RV_jetV[25]/D");
  tree->Branch("trk_sIP1D_val", &trk_sIP1D_val, "trk_sIP1D_val[25]/D");
  tree->Branch("trk_sIP1D_sig", &trk_sIP1D_sig, "trk_sIP1D_sig[25]/D");
  tree->Branch("trk_IP3D_err", &trk_IP3D_err, "trk_IP3D_err[25]/D");
  tree->Branch("trk_sIP3D_err", &trk_sIP3D_err, "trk_sIP3D_err[25]/D");
  tree->Branch("trk_IP3D_err_RV_jetV", &trk_IP3D_err_RV_jetV, "trk_IP3D_err_RV_jetV[25]/D");
  tree->Branch("trk_sIP3D_err_RV_jetV", &trk_sIP3D_err_RV_jetV, "trk_sIP3D_err_RV_jetV[25]/D");
  tree->Branch("trk_IP2D_err", &trk_IP2D_err, "trk_IP2D_err[25]/D"); 
  tree->Branch("trk_sIP2D_err", &trk_sIP2D_err, "trk_sIP2D_err[25]/D");
  tree->Branch("trk_IP2D_err_RV_jetV", &trk_IP2D_err_RV_jetV, "trk_IP2D_err_RV_jetV[25]/D"); 
  tree->Branch("trk_sIP2D_err_RV_jetV", &trk_sIP2D_err_RV_jetV, "trk_sIP2D_err_RV_jetV[25]/D");
  tree->Branch("track_pt", &track_pt, "track_pt[25]/D");
  tree->Branch("trk_IP1D_err", &trk_IP1D_err, "trk_IP1D_err[25]/D");
  tree->Branch("trk_sIP1D_err", &trk_sIP1D_err, "trk_sIP1D_err[25]/D");
  tree->Branch("PV_x", &PV_x, "PV_x/D"); 
  tree->Branch("PV_y", &PV_x, "PV_y/D");
  tree->Branch("PV_z", &PV_z, "PV_z/D");
  tree->Branch("trk_sIP1D_val_x", &trk_sIP1D_val_x,"trk_sIP1D_val_x[25]/D");
  tree->Branch("trk_sIP1D_val_y", &trk_sIP1D_val_y,"trk_sIP1D_val_y[25]/D");
  tree->Branch("trk_sIP1D_val_z", &trk_sIP1D_val_z,"trk_sIP1D_val_z[25]/D");
  tree->Branch("diff_RP_PV_x", &diff_RP_PV_x,"diff_RP_PV_x[25]/D");
  tree->Branch("diff_RP_PV_y", &diff_RP_PV_y,"diff_RP_PV_y[25]/D");
  tree->Branch("diff_RP_PV_z", &diff_RP_PV_z,"diff_RP_PV_z[25]/D");
  tree->Branch("IPVec_PV_x", &IPVec_PV_x,"IPVec_PV_x[25]/D");
  tree->Branch("IPVec_PV_y", &IPVec_PV_y,"IPVec_PV_y[25]/D");
  tree->Branch("IPVec_PV_z", &IPVec_PV_z,"IPVec_PV_z[25]/D");
  //ee->Branch("FD3D_val", &FD3D_val,"FD3D_val[25]/D");
  tree->Branch("FD3D_val", &FD3D_val,"FD3D_val/D");
  tree->Branch("FD3D_sig", &FD3D_sig,"FD3D_sig/D");
  tree->Branch("FD2D_val", &FD2D_val,"FD2D_val/D");
  tree->Branch("FD2D_sig", &FD2D_sig,"FD2D_sig/D");
  tree->Branch("sDL1D_val", &sDL1D_val,"sDL1D_val[25]/D");
  tree->Branch("sDL1D_sig", &sDL1D_sig,"sDL1D_sig[25]/D");
  tree->Branch("sDL3D_val", &sDL3D_val,"sDL3D_val[25]/D");
  tree->Branch("sDL3D_sig", &sDL3D_sig,"sDL3D_sig[25]/D");
  tree->Branch("sDL2D_val", &sDL2D_val,"sDL2D_val[25]/D");
  tree->Branch("sDL2D_sig", &sDL2D_sig,"sDL2D_sig[25]/D");
  tree->Branch("DL1D_val", &DL1D_val,"DL1D_val[25]/D");
  tree->Branch("DL3D_sig", &DL3D_sig,"DL3D_sig[25]/D");
  tree->Branch("DL2D_sig", &DL2D_sig,"DL2D_sig[25]/D");
  tree->Branch("DL1D_sig", &DL1D_sig,"DL1D_sig[25]/D");
  tree->Branch("DL3D_val", &DL3D_val,"DL3D_val[25]/D");
  tree->Branch("DL2D_val", &DL2D_val,"DL2D_val[25]/D");
  tree->Branch("DL1D_val", &DL1D_val,"DL1D_val[25]/D");
 }
 /////
 //   Set branch address
 /////
 //Connects the branches of an existing TTree to variables used when loading the file
 void set_branch_addresses(void){
  tree->SetBranchAddress("duno", &duno);
  tree->SetBranchAddress("sumhcsvjets", &sumhcsvjets);
  tree->SetBranchAddress("prodhcsvjets", &prodhcsvjets);
  tree->SetBranchAddress("ngenb", &ngenb);
  tree->SetBranchAddress("lep_numl", &lep_numl);
  tree->SetBranchAddress("lep_numt", &lep_numt);
  tree->SetBranchAddress("lep_dichprodl", &lep_dichprodl);
  tree->SetBranchAddress("lep_dichprodt", &lep_dichprodt);
  tree->SetBranchAddress("jet_num", &jet_num);
  tree->SetBranchAddress("jet_ngenb", &jet_ngenb);
  tree->SetBranchAddress("jet_ngenbdr", &jet_ngenbdr);
  tree->SetBranchAddress("jet_pt", &jet_pt);
  tree->SetBranchAddress("jet_eta", &jet_eta);
  tree->SetBranchAddress("jet_phi", &jet_phi);
  tree->SetBranchAddress("jet_en", &jet_en);
  tree->SetBranchAddress("jet_csv", &jet_csv);
  tree->SetBranchAddress("partonFlavour", &partonFlavour);
  tree->SetBranchAddress("partonFlavourdR", &partonFlavourdR);
  tree->SetBranchAddress("deltaR_reco_gen", &deltaR_reco_gen);
  tree->SetBranchAddress("deltaR_reco_jet", &deltaR_reco_jet);
  tree->SetBranchAddress("deltaR_gen_jet", &deltaR_gen_jet);
  //Num trk
  tree->SetBranchAddress("jet_ndaus", &jet_ndaus);
  tree->SetBranchAddress("jet_chtrks", &jet_chtrks);
  tree->SetBranchAddress("jet_chtrkspv", &jet_chtrkspv);
  tree->SetBranchAddress("jet_chtrksnpv", &jet_chtrksnpv);
  tree->SetBranchAddress("jet_chtrkspvtt", &jet_chtrkspvtt);
  tree->SetBranchAddress("jet_chtrksnpvtt", &jet_chtrksnpvtt);
  //Chi2 info
  tree->SetBranchAddress("jet_chi2tot", &jet_chi2tot);
  tree->SetBranchAddress("jet_chi2ndf", &jet_chi2ndf);
  tree->SetBranchAddress("jet_chi2pval", &jet_chi2pval);
  //Two trk info
  tree->SetBranchAddress("jet_num2v", &jet_num2v);
  tree->SetBranchAddress("jet_numno2v", &jet_numno2v);
  tree->SetBranchAddress("jet_num2vno2v", &jet_num2vno2v);
  tree->SetBranchAddress("jet_dca3d2t", &jet_dca3d2t);
  tree->SetBranchAddress("jet_dca3dno2t", &jet_dca3dno2t);
  tree->SetBranchAddress("jet_dca2d2t", &jet_dca2d2t);
  tree->SetBranchAddress("jet_dca2dno2t", &jet_dca2dno2t);
  tree->SetBranchAddress("jet_dca3d2tno2t", &jet_dca3d2tno2t);
  tree->SetBranchAddress("jet_dca2d2tno2t", &jet_dca2d2tno2t);
  //Primary Vertex Info
  tree->SetBranchAddress("RV_kvf_bs_3D", &RV_kvf_bs_3D);
  tree->SetBranchAddress("RV_avf_3D", &RV_avf_3D);
  tree->SetBranchAddress("RV_avf_bs_3D", &RV_avf_bs_3D);
  tree->SetBranchAddress("RV_kvf_3D", &RV_kvf_3D);
  tree->SetBranchAddress("RV_RVaddjtrcks_NPV_3D",&RV_RVaddjtrcks_NPV_3D);
  tree->SetBranchAddress("RV_RVnojet_3D",&RV_RVnojet_3D);
  tree->SetBranchAddress("diff_PVx_RVx_kvf", &diff_PVx_RVx_kvf);
  tree->SetBranchAddress("diff_PVy_RVy_kvf", &diff_PVy_RVy_kvf);
  tree->SetBranchAddress("diff_PVz_RVz_kvf", &diff_PVz_RVz_kvf);
  tree->SetBranchAddress("diff_PVx_RVx_kvf_bs", &diff_PVx_RVx_kvf_bs);
  tree->SetBranchAddress("diff_PVy_RVy_kvf_bs", &diff_PVy_RVy_kvf_bs);
  tree->SetBranchAddress("diff_PVz_RVz_kvf_bs", &diff_PVz_RVz_kvf_bs);
  tree->SetBranchAddress("diff_PVx_RVx_kvf", &diff_PVx_RVx_kvf);
  tree->SetBranchAddress("diff_PVy_RVy_kvf", &diff_PVy_RVy_kvf);
  tree->SetBranchAddress("diff_PVz_RVz_kvf", &diff_PVz_RVz_kvf);
  tree->SetBranchAddress("diff_PVx_RVx_avf_bs", &diff_PVx_RVx_avf_bs);
  tree->SetBranchAddress("diff_PVy_RVy_avf_bs", &diff_PVy_RVy_avf_bs);
  tree->SetBranchAddress("diff_PVz_RVz_avf_bs", &diff_PVz_RVz_avf_bs);
  tree->SetBranchAddress("PVx", &PVx);
  tree->SetBranchAddress("PVy", &PVy);
  tree->SetBranchAddress("PVz", &PVz);
  tree->SetBranchAddress("PV_nojet_x", &PV_nojet_x);
  tree->SetBranchAddress("PV_nojet_y", &PV_nojet_y);
  tree->SetBranchAddress("PV_nojet_z", &PV_nojet_z);
  tree->SetBranchAddress("diff_RV_RVnob_x", &diff_RV_RVnob_x);
  tree->SetBranchAddress("diff_RV_RVnob_y", &diff_RV_RVnob_y);
  tree->SetBranchAddress("diff_RV_RVnob_z", &diff_RV_RVnob_z);
  tree->SetBranchAddress("pv_ntrks", &pv_ntrks);
  tree->SetBranchAddress("RVnojet_chi2_ndf", &RVnojet_chi2_ndf); 
  tree->SetBranchAddress("diff_pv_RVnojet_chi2_ndf", &diff_pv_RVnojet_chi2_ndf);
  tree->SetBranchAddress("RV_nchi2_ndf", &RV_nchi2_ndf);
  tree->SetBranchAddress("diff_RV_RVnojet_chi2_ndf", &diff_RV_RVnojet_chi2_ndf);
  tree->SetBranchAddress("diff_chi2_ndf_RV_addjettrks_RVnob", &diff_chi2_ndf_RV_addjettrks_RVnob); 
  //Get the IP info of the Jet Tracks
  tree->SetBranchAddress("trk_IP3D_val",&trk_IP3D_val);
  tree->SetBranchAddress("trk_IP3D_sig",&trk_IP3D_sig);
  tree->SetBranchAddress("trk_IP3D_val_RV_jetV",&trk_IP3D_val_RV_jetV);
  tree->SetBranchAddress("trk_IP3D_sig_RV_jetV",&trk_IP3D_sig_RV_jetV);
  tree->SetBranchAddress("trk_IP2D_val",&trk_IP2D_val);
  tree->SetBranchAddress("trk_IP2D_sig",&trk_IP2D_sig);
  tree->SetBranchAddress("trk_IP2D_val_RV_jetV",&trk_IP2D_val_RV_jetV);
  tree->SetBranchAddress("trk_IP2D_sig_RV_jetV",&trk_IP2D_sig_RV_jetV);
  tree->SetBranchAddress("trk_IP1D_val",&trk_IP1D_val);
  tree->SetBranchAddress("trk_IP1D_sig",&trk_IP1D_sig);
  tree->SetBranchAddress("trk_sIP3D_val",&trk_sIP3D_val);
  tree->SetBranchAddress("trk_sIP3D_sig",&trk_sIP3D_sig);
  tree->SetBranchAddress("trk_sIP2D_val",&trk_sIP2D_val);
  tree->SetBranchAddress("trk_sIP2D_sig",&trk_sIP2D_sig);
  tree->SetBranchAddress("trk_sIP3D_val_RV_jetV",&trk_sIP3D_val_RV_jetV);
  tree->SetBranchAddress("trk_sIP3D_sig_RV_jetV",&trk_sIP3D_sig_RV_jetV);
  tree->SetBranchAddress("trk_sIP2D_val_RV_jetV",&trk_sIP2D_val_RV_jetV);
  tree->SetBranchAddress("trk_sIP2D_sig_RV_jetV",&trk_sIP2D_sig_RV_jetV);
  tree->SetBranchAddress("trk_sIP1D_val",&trk_sIP1D_val);
  tree->SetBranchAddress("trk_sIP1D_sig",&trk_sIP1D_sig);
  tree->SetBranchAddress("trk_IP3D_err",&trk_IP3D_err);
  tree->SetBranchAddress("trk_sIP3D_err",&trk_sIP3D_err);
  tree->SetBranchAddress("trk_IP2D_err",&trk_IP2D_err);
  tree->SetBranchAddress("trk_sIP2D_err",&trk_sIP2D_err);
  tree->SetBranchAddress("trk_IP3D_err_RV_jetV",&trk_IP3D_err_RV_jetV);
  tree->SetBranchAddress("trk_sIP3D_err_RV_jetV",&trk_sIP3D_err_RV_jetV);
  tree->SetBranchAddress("trk_IP2D_err_RV_jetV",&trk_IP2D_err_RV_jetV);
  tree->SetBranchAddress("trk_sIP2D_err_RV_jetV",&trk_sIP2D_err_RV_jetV);
  tree->SetBranchAddress("trk_IP1D_err",&trk_IP1D_err);
  tree->SetBranchAddress("trk_sIP1D_err",&trk_sIP1D_err);
  tree->SetBranchAddress("PV_x", &PV_x);
  tree->SetBranchAddress("PV_y", &PV_y);
  tree->SetBranchAddress("PV_z", &PV_z);
  tree->SetBranchAddress("RV_x", &RV_x);
  tree->SetBranchAddress("RV_y", &RV_y);
  tree->SetBranchAddress("RV_z", &RV_z);
  tree->SetBranchAddress("RVnojtrk_x", &RVnojtrk_x);
  tree->SetBranchAddress("RVnojtrk_y", &RVnojtrk_y);
  tree->SetBranchAddress("RVnojtrk_z", &RVnojtrk_z);
  tree->SetBranchAddress("trk_sIP1D_val_x", &trk_sIP1D_val_x);
  tree->SetBranchAddress("trk_sIP1D_val_y", &trk_sIP1D_val_y);
  tree->SetBranchAddress("trk_sIP1D_val_z", &trk_sIP1D_val_z);
  tree->SetBranchAddress("diff_RP_PV_x", &diff_RP_PV_x);
  tree->SetBranchAddress("diff_RP_PV_y", &diff_RP_PV_y);
  tree->SetBranchAddress("diff_RP_PV_z", &diff_RP_PV_z);
  tree->SetBranchAddress("IPVec_PV_x", &IPVec_PV_x);
  tree->SetBranchAddress("IPVec_PV_y", &IPVec_PV_y);
  tree->SetBranchAddress("IPVec_PV_z", &IPVec_PV_z);
  tree->SetBranchAddress("track_pt", &track_pt);
  tree->SetBranchAddress("FD3D_val", &FD3D_val);
  tree->SetBranchAddress("FD3D_sig", &FD3D_sig);
  tree->SetBranchAddress("FD2D_val", &FD2D_val);
  tree->SetBranchAddress("FD2D_sig", &FD2D_sig);
  tree->SetBranchAddress("DL1D_val", &DL1D_val);
  tree->SetBranchAddress("DL1D_sig", &DL1D_sig);
  tree->SetBranchAddress("DL3D_val", &DL3D_val);
  tree->SetBranchAddress("DL3D_sig", &DL3D_sig);
  tree->SetBranchAddress("DL2D_val", &DL2D_val);
  tree->SetBranchAddress("DL2D_sig", &DL2D_sig);
  tree->SetBranchAddress("sDL3D_val", &sDL3D_val);
  tree->SetBranchAddress("sDL3D_sig", &sDL3D_sig);
  tree->SetBranchAddress("sDL2D_val", &sDL2D_val);
  tree->SetBranchAddress("sDL2D_sig", &sDL2D_sig);
  tree->SetBranchAddress("sDL1D_val", &sDL1D_val);
  tree->SetBranchAddress("sDL1D_sig", &sDL1D_sig);
  tree->SetBranchAddress("RV_kvf_ntrks",&RV_kvf_ntrks);
  tree->SetBranchAddress("RV_nojet_ntrks",&RV_nojet_ntrks);
 }
 };
#endif
