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
 //Num trk info
 double jet_ndaus, jet_chtrks, jet_chtrkspv, jet_chtrksnpv, jet_chtrkspvtt, jet_chtrksnpvtt;
 //Chi2 info
 double jet_chi2tot, jet_chi2ndf, jet_chi2pval;
 //Two trk info
 double jet_num2v, jet_numno2v, jet_num2vno2v, jet_dca3d2t, jet_dca3dno2t, jet_dca3d2tno2t, jet_dca2d2t, jet_dca2dno2t, jet_dca2d2tno2t;
 //Primary Vertex Info
 double PVxpvx, PVypvy, PVzpvz;
 /////
 ////   Initialise
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
 // Primary Vertex Information 
  PVxpvx = -9999.;
  PVypvy = -9999.;
  PVzpvz = -9999.;
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
  tree->Branch("PVxpvx",&PVxpvx,"PVxpvx/D");
  tree->Branch("PVypvy",&PVypvy,"PVypvy/D");
  tree->Branch("PVzpvz",&PVzpvz,"PVzpvz/D");
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
  tree->SetBranchAddress("PVxpvx", &PVxpvx);
  tree->SetBranchAddress("PVypvy", &PVypvy);
  tree->SetBranchAddress("PVzpvz", &PVzpvz);
   
 }
};
#endif
