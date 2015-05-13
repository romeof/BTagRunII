/**
This Macro
1. Calls trees, fix cuts and creates ntuples with cut data  

Need to specify
0. See Declare constants
*/
/////
//   To run: root -l SkimNtuple.cc+  
/////
/////
//   Prepare Root and Roofit
/////
#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"
#include "TLegend.h"
#include "TTree.h"
#include "TTreePlayer.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <iostream>
using namespace std;
////
//   Declare constants
/////
//Path - samples 
const string path     = "/afs/cern.ch/work/s/sshaheen/CMSSW_7_2_3/src/BTagRunII/BTagReco/rootfiles/signal/";
const char *samples[] = {"sgnl_nocsv_sub"};
const string specsel  = "";
//Selection
const string suffisso   = "goodevt_2mupt2010_check_abs_sgnl_nocsv"; 
//const string selection  = "mu_num==2 && ele_num==0 && lep_charge[0]==lep_charge[1] && lep_pt[0]>20 && lep_pt[1]>10";

const string selection = "abs(partonFlavour[0]) ==5 && abs(partonFlavour[0]) !=9999";
//const string selection = "abs(partonFlavour[0]) !=4 && abs(partonFlavour[0]) !=5 && abs(partonFlavour[0]) !=9999";
const string dotroot    = ".root"; 
/////
//   Declare functions 
/////
TFile* Call_TFile(string rootpla);
/////
//   Main function
/////
void SkimNtuple(){
 //For all the samples
 vector<string> rootplas(samples, samples + sizeof(samples)/sizeof(samples[0]));
 for(uint i=0; i<rootplas.size(); i++){
  //Call old ntupla
  TFile* f = Call_TFile(rootplas[i]); TTree* tree; f->GetObject("demo/tree",tree);
  //Create new ntupla
  string newfil  = specsel+rootplas[i]+"_"+suffisso+dotroot; 
  TFile *newfile = new TFile(newfil.c_str(),"recreate");
  newfile->cd();
  //Do cut
  string varCut  = selection;
  TTree* newtree = tree->CopyTree(varCut.c_str());
  //Save
  newtree->Write(); 
  newfile->Write();
  newfile->Close();
 }
}
/////
//   Call TFile to be read
/////
TFile* Call_TFile(string rootpla){
 string file_name = path+specsel+rootpla+dotroot;
 TFile* f = new TFile(file_name.c_str(),"update");
 return f;
}
//const string selection  = "l1pvval!=-9999 && nmus==2 && neles==0 && cands_pt[0]>25 && cands_pt[1]>25";
// && cands_pt[0]>25 && cands_pt[1]>25 && mus_pt[1]==cands_pt[1]";// && cands_pt[0]>20 && cands_pt[1]>10";
// ((nmus==2 && mus_pt[0]>20 && mus_pt[1]>20) || (neles==2 && eles_pt[0]>20 && eles_pt[1]>20) || (nmus+neles==2 && mus_pt[0]>20 && meles_pt[0]>20))";
//2mu
//const string suffisso   = "2mu"; 
//const string selection  = "l1pvval!=-9999 && nmus==2 && neles==0";
//sig 2L
//const string suffisso   = "sr_2lep"; 
//const string selection  = "nLepGood==2 && LepGood_pt[0]>20 && LepGood_pt[1]>20 && LepGood_charge[0]*LepGood_charge[1] > 0 && nJet25>=4 && (nBJetLoose25 >= 2 || nBJetMedium25 >= 1)";
//const string suffisso   = "sr_3lep"; 
//const string selection  = "nLepGood==3 && LepGood_pt[0]>20 && LepGood_pt[1]>20 && LepGood_pt[2]>20 && nJet25>=4 && (nBJetLoose25 >= 2 || nBJetMedium25 >= 1)";

