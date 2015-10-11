/**
This Macro   
1. Plots variables of signal and background overlaid

Need to specify
1. See Declare Constants
2. There is a difference between int first[arraycomp]; int second[arraycomp]; and double first[arraycomp]; double second[arraycomp];
*/
/////
//   To run: root -l SigBkgVarsOverlaid.cc+
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
#include "TGaxis.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <iostream>
using namespace std;
/////
//   Declare constants
/////
const string path     = "/afs/cern.ch/work/s/sshaheen/CMSSW_7_2_3/src/BTagRunII/BTagReco/macros_panga/";
const char *samples[] = { "sig_RVnew_","bak_RVnew_","charm_RVnew_"};
const string selection  = "";
const int numvar        = 100;
//declare variables to draw

//const char *varfirst[]        = {"diff_PVx_RVx_kvf","diff_PVy_RVy_kvf","diff_PVz_RVz_kvf","diff_PVx_RVx_kvf_bs","diff_PVy_RVy_kvf_bs","diff_PVz_RVz_kvf_bs","diff_PVx_RVx_avf","diff_PVy_RVy_avf","diff_PVz_RVz_avf","diff_PVx_RVx_avf_bs","diff_PVy_RVy_avf_bs","diff_PVz_RVz_avf_bs"};
//const char *varfirst[]        = {"PVz","PVz","PVz","PVz","PVz","PVz","PVz","PVz"};
//const char *varsecond[]       = {"RV_z","RV_kvf_bs_3D","RV_avf_3D","RV_avf_bs_3D","RV_z","RV_kvf_bs_3D","RV_avf_3D","RV_avf_bs_3D"};
const char *varfirst[]          = {"diff_chi2_ndf_RV_addjettrks_RVnob"};
const char *varsecond[]         = {"1","diff_chi2_ndf_RV_addjettrks_RVnob","1","PVx","PVy","PVz"};
//const char *vartitle[]        = {"(diff_PVx_RVx_kvf)/PVx","(diff_PVy_RVy_kvf)/PVy","(diff_PVz_RVz_kvf)/PVz","(diff_PVx_RVx_kvf_bs)/PVx","(diff_PVy_RVy_kvf_bs)/PVy","(diff_PVz_RVz_kvf_bs)/PVz","(diff_PVx_RVx_avf)/PVx","(diff_PVy_RVy_avf)/PVy","(diff_PVz_RVz_avf)/PVz","(diff_PVx_RVx_avf_bs)/PVx","(diff_PVy_RVy_avf_bs)/PVy","(diff_PVz_RVz_avf_bs)/PVz"};
const char *vartitle[]        = {"#chi^{2}/ndf_RVaddjtrk_RVnojet","#chi^{2}/ndf(RV-RVaddjtrck)","diff_RV_RVnob_z","(diff_RV_RVnob_x)/PV_x","(diff_RV_RVnob_y)/PV_y","(diff_RV_RVnob_z/PV_z)"};
const double inRange[numvar]  = {-5,-100,-100,-0.1,-0.1,-0.1,-0.1,-0.1,-0.02,-0.02,-0.02,-0.02};
const double endRange[numvar] = {7,100,100,0.1,0.1,0.1,0.1,0.1,0.1,0.02,0.02,0.02,};
const int    bin[numvar]      = {100,100,100,100,100,100,100,100,100,100,100,100};
double ylogscale[numvar]      = {1,1,1,1,1,1,1,1,1,1,1,1};
double setminimum[numvar]     = {0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01};
double setmaximum[numvar]     = {100,100,100,100,100,100,100,100,100,100,100,100};
double overflow[numvar]       = {5,500,1000,500,500,1000,500,500,1000,500,500,1000};
double overflow1[numvar]      = {5,500,1000,500,500,1000,500,500,1000,500,500,1000};
double overflowl[numvar]      = {-500,-500,-1000,-500,-500,-1000,-500,-500,-1000,-500,-500,-1000};
//0ther options
bool saveplots = true;

//Variables
/////
//   Declare functions 
/////
TFile* Call_TFile(string rootpla);
string convertToString(double number);
void setTDRStyle();
/////
//   Main function
/////
void SigBkgVarsOverlaid(){
 //Preliminarly
 setTDRStyle();
 //Loop over variables
 vector<string> varfirsts(varfirst, varfirst + sizeof(varfirst)/sizeof(varfirst[0])); 
 vector<string> varseconds(varsecond, varsecond + sizeof(varsecond)/sizeof(varsecond[0]));
 vector<string> vartitles(vartitle, vartitle + sizeof(vartitle)/sizeof(vartitle[0])); 
 const uint variables_size = varfirsts.size();
 cout<<"the size is"<<variables_size;
 for(uint vars=0; vars<variables_size; vars++)
 //for(uint vars=0; vars<1; vars++)
 {
  TCanvas* c1 = new TCanvas(vartitles[vars].c_str(),vartitles[vars].c_str(),200,200,800,600);
cout<<"the variable is"<<varfirsts[vars]<<endl; 
 if(ylogscale[vars]==1) c1->SetLogy();
  TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  //jet_csv //TLegend *leg = new TLegend(0.5, 0.7, 0.7, 0.9);
  leg->SetHeader("Samples");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.05);
  //leg->SetTextSize(0.035);
  //Do plots
  //Loop over samples
  vector<string> rootplas(samples, samples + sizeof(samples)/sizeof(samples[0]));
  const uint rootplas_size = rootplas.size();
  for(uint smp=0; smp<rootplas_size; smp++){
   cout<<"Rootpla "<<rootplas[smp]<<endl;
   //Declare histo
   TH1F *hist = new TH1F(rootplas[smp].c_str(),rootplas[smp].c_str(),bin[vars],inRange[vars],endRange[vars]); 
   hist->SetTitle("");
   hist->SetMarkerStyle(1);
   hist->GetXaxis()->SetTitle(vartitles[vars].c_str());
   hist->GetYaxis()->SetTitle("Percentage");
   if(smp==0){
    hist->SetMarkerColor(kRed);
    hist->SetLineColor(kRed);
   }else if(smp==1){
    hist->SetMarkerColor(kBlue);
    hist->SetLineColor(kBlue);
   }
   else{
    hist->SetMarkerColor(kGreen);
    hist->SetLineColor(kGreen);
    }
   //Make plot
   TFile* f = Call_TFile((rootplas[smp]).c_str()); TTree* tree; f->GetObject("tree",tree);
   //int arraysize=100;
   double first;
   double second;
   double ratio=0;
   //double first[arraycomp];
   int a=0;
   tree->SetBranchAddress(varfirsts[vars].c_str(),&first); 
   tree->SetBranchAddress(varseconds[vars].c_str(),&second); 
   for(int en=0; en<tree->GetEntries(); en++)
    
   //cout<<"the panga is"<<first[vars]<<endl;
   //for(int en=0; en<25; en++)
   {
    tree->GetEntry(en);
    if(first==-9999) continue;
    if(varseconds[vars]=="1") 
    {
    if(first>overflow1[vars]) first=overflow1[vars];
    if(first<overflowl[vars]) first=overflowl[vars];
     //hist->Fill(first[vars]);
    hist->Fill(first);  
    }
    //if(vars<=3){
    if(varseconds[vars]!="1"){
    if(vars<=2){
    ratio=first-second;
   //Overflow the last bin
    if(ratio>overflow[vars]) 
   {
   ratio=overflow[vars];
   }
   }else{
   ratio=first/second;
   //cout<<"It should not work "<<endl;
   if(ratio>overflow[vars]) ratio=overflow[vars];
   //cout<<"vars";
   }
   //ratio= numerator/denominator;
   
   //if(ratio>overflow[vars]) ratio=overflow[vars];
   //cout<<"the ratio is "<<ratio<<endl;
    hist->Fill(ratio); 
   } 
    }
  // }
   cout<<"the entries with -999 values are"<<a<<endl;
   cout<<"Entries of "<<rootplas[smp]<<" is "<<hist->Integral()<<endl;
   double scale = hist->Integral();
   hist->Scale(1/scale*100);
   //Draw plot
   gStyle->SetOptStat("mr");     
   gStyle->SetStatColor(kWhite);
   gStyle->SetStatX(0.9); //Starting position on X axis
   gStyle->SetStatW(0.3); //Horizontal size 
   //jet_csv
   //gStyle->SetStatX(0.7); //Starting position on X axis
   //gStyle->SetStatW(0.2); //Horizontal size 
   if(smp==0){
    gStyle->SetStatTextColor(kRed);
    gStyle->SetStatY(0.7); //Starting position on Y axis
    gStyle->SetStatFontSize(0.2); //Vertical Size
   }else if(smp==1){
    gStyle->SetStatTextColor(kBlue);
    gStyle->SetStatY(0.5); //Starting position on Y axis
    gStyle->SetStatFontSize(0.2); //Vertical Size
   } else{
    gStyle->SetStatTextColor(kGreen);
    gStyle->SetStatY(0.3); //Starting position on Y axis
    gStyle->SetStatFontSize(0.2); //Vertical Size
    }
   if(smp==0){
    leg->AddEntry(hist,"b jet","L"); 
    hist->SetMinimum(setminimum[vars]);
    hist->SetMaximum(setmaximum[vars]);
    hist->Draw(""); 
   }else if(smp==1){
    leg->AddEntry(hist,"light Flavor jets","L"); 
    hist->SetMinimum(setminimum[vars]);
    hist->SetMaximum(setmaximum[vars]);
    hist->Draw("sames");
   }else{
   leg->AddEntry(hist,"c jets","L");
    hist->SetMinimum(setminimum[vars]);
    hist->SetMaximum(setmaximum[vars]);
    hist->Draw("sames");
   }
  } 
  leg->Draw();
  string num = convertToString(double(vars));
  string namefile = "plots_RPVnojet/"+varfirsts[vars]+"_"+num+".pdf";
  if(saveplots) c1->SaveAs(namefile.c_str());
 }
}
/////
//   Call TFile to be read
/////
TFile* Call_TFile(string rootpla){
 string dotroot   = ".root";
 string file_name = path+rootpla+dotroot;
 cout<<"file_name "<<file_name<<endl;
 TFile* f = new TFile(file_name.c_str(),"update");
 return f;
}
/////
//   Convert number to string 
/////
string convertToString(double number){
 char buffer[5];
 sprintf(buffer,"%g",number);
 return string(buffer);
}
/////
//   Set setTDRStyle_modified (from link https://twiki.cern.ch/twiki/pub/CMS/TRK10001/setTDRStyle_modified.C)
/////
void setTDRStyle(){
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

  // For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);
  
  // For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);
  // For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  // For the histo:
  tdrStyle->SetHistFillColor(0);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  tdrStyle->SetErrorX(0.);
//  tdrStyle->SetErrorMarker(20);

  tdrStyle->SetMarkerStyle(20);

  //For the fit/function:
  tdrStyle->SetOptFit(0);
  tdrStyle->SetFitFormat("5.4g");
  //tdrStyle->SetFuncColor(1);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  //For the date:
  tdrStyle->SetOptDate(0);

  // For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetTextAlign(11);
  tdrStyle->SetStatBorderSize(0);;

  // Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.165);
  tdrStyle->SetPadLeftMargin(0.1);
  tdrStyle->SetPadRightMargin(0.05);

  // For the Global title:
  //  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.5);
  tdrStyle->SetTitleH(0.05); // Set the height of the title box
  //tdrStyle->SetTitleW(0); // Set the width of the title box
  tdrStyle->SetTitleX(0.15); // Set the position of the title box
  tdrStyle->SetTitleY(1.0); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  tdrStyle->SetTitleBorderSize(0);

  // For the axis titles:
  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.045, "X");
  tdrStyle->SetTitleSize(0.06, "YZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(1.5);
  //tdrStyle->SetTitleYOffset(1.0);
  tdrStyle->SetTitleOffset(0.75, "Y"); // Another way to set the Offset

  // For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

  // For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  // Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);
  tdrStyle->cd();
}
