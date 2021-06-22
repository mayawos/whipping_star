#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <cstring>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TMath.h"
#include "TSystem.h"
#include "TMatrixT.h"
#include "TMatrixD.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TError.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TDecompSVD.h"
#include "TText.h"
#include "TPaveStats.h"

#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBNfit.h"
#include "SBNfit3pN.h"
#include "SBNgenerate.h"
#include "prob.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;


void DrawDataMCAndSyst(TString cName, TH1D* MC1, TH1D* MC2, TH1D* data, std::string var, std::string title);
void plot_one(TMatrixD matrix, std::string tag);
bool verbose = true;

int main(int argc, char* argv[])
{
  std::string xml = "example.xml";
  std::string tag = "example";
  int setno = 1;
  bool sys = false;
  bool detsys = false;
  bool mcerr = false;
  bool addext = false;
  
  const struct option longopts[] =
    {
      {"xml", 		required_argument, 	0, 'x'},
      {"tag", 		required_argument, 	0, 't'},
      {"setno", 	required_argument, 	0, 'n'},
      {"sys",	        no_argument, 		0, 's'},
      {"detsys",	no_argument, 		0, 'd'},
      {"mcerr",	        no_argument, 		0, 'm'},
      {"addext",	no_argument, 		0, 'a'},
      {0,   		        no_argument,     	0,  0},
    };
  
  int iarg = 0;
  opterr=1;
  int index;
  
  while(iarg != -1)
    {
      iarg = getopt_long(argc,argv, ":x:t:n:dmas:", longopts, &index);
      
      switch(iarg)
	{
	case 'x':
	  xml = optarg;
	  break;
	case 't':
	  tag = optarg;
	  break;
	case 'n':
	  setno = atoi(optarg);
	  break;
	case 's':
	  sys = true;
	  break;
	case 'd':
	  detsys = true;
	  break;
	case 'm':
	  mcerr = true;
	  break;
	case 'a':
	  addext = true;
	  break;
	case '?':
	case 'h':
	  std::cout<<"Allowed arguments:"<<std::endl;
	  std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
	  return 0;
	}
    }
  
  gStyle->SetPalette(87);

  std::string detsyslabel = "";
  std::string mcerrlabel = "";
  std::string zerobinlabel = "";
  if( mcerr ) mcerrlabel = "_with_mc_err"; 
  if( addext ) zerobinlabel = "_with_zerobin_err"; 
  if( detsys ) detsyslabel = "_detsys"; 
  //Lçoad up our covariance matricies we calculated in example1 (we could also load up single variation ones)
  TString bf = Form("../bin/constrained_%s_mc_fakedata%d%s%s%s_BF.SBNspec.root",tag.c_str(),setno,mcerrlabel.c_str(),zerobinlabel.c_str(),detsyslabel.c_str());
  TString data = Form("/uboone/app/users/davidc/SBNFit/build/bin/fakedata1/0126/10bin/constrained_np_zp_numu_reco_e_H1_mc_fakedata%d_DATA_CV.SBNspec.root",setno);
  std::cout << "bf = " << bf << std::endl;
  //TString data = Form("../bin/%s_fakedata%d.SBNspec.root",tag.c_str(),setno);
  //TString run1 = Form("../examples/%s_run1%s%s%s.SBNspec.root",tag.c_str(),mcerrlabel.c_str(),zerobinlabel.c_str(),detsyslabel.c_str());
  //TString run2 = Form("../examples/%s_run2%s%s%s.SBNspec.root",tag.c_str(),mcerrlabel.c_str(),zerobinlabel.c_str(),detsyslabel.c_str());
  //TString run3 = Form("../examples/%s_run3%s%s%s.SBNspec.root",tag.c_str(),mcerrlabel.c_str(),zerobinlabel.c_str(),detsyslabel.c_str());
  //TString runall = Form("../examples/%s%s%s%s.SBNspec.root",tag.c_str(),mcerrlabel.c_str(),zerobinlabel.c_str(),detsyslabel.c_str());
  TFile * fbf = new TFile(bf,"read");
  TFile * fdata = new TFile(data,"read");
  //TFile * frun1 = new TFile(run1,"read");
  //TFile * frun2 = new TFile(run2,"read");
  //TFile * frun3 = new TFile(run3,"read");
  //TFile * frunall = new TFile(runall,"read");
  TH1D *histo_1eNp_lee_bf = (TH1D*)fbf->Get("nu_uBooNE_1eNp_lee");
  TH1D *histo_1e0p_lee_bf = (TH1D*)fbf->Get("nu_uBooNE_1e0p_lee");
  TH1D *histo_1eNp_bg_bf = (TH1D*)fbf->Get("nu_uBooNE_1eNp_bg");
  TH1D *histo_1e0p_bg_bf = (TH1D*)fbf->Get("nu_uBooNE_1e0p_bg");
  histo_1eNp_lee_bf->Add(histo_1eNp_bg_bf);
  histo_1e0p_lee_bf->Add(histo_1e0p_bg_bf);
  TH1D *histo_1eNp_data = (TH1D*)fdata->Get("nu_uBooNE_1eNp_data");
  TH1D *histo_1e0p_data = (TH1D*)fdata->Get("nu_uBooNE_1e0p_data");
  //TH1D *histo_1eNp_lee_run1 = (TH1D*)frun1->Get("nu_uBooNE_1eNp_lee");
  //TH1D *histo_1eNp_lee_run2 = (TH1D*)frun2->Get("nu_uBooNE_1eNp_lee");
  //TH1D *histo_1eNp_lee_run3 = (TH1D*)frun3->Get("nu_uBooNE_1eNp_lee");
  //TH1D *histo_1eNp_lee_runall = (TH1D*)frunall->Get("nu_uBooNE_1eNp_lee");
  //TH1D *histo_1eNp_lee_run123_sum = (TH1D*)histo_1eNp_lee_run1->Clone("nu_uBooNE_1eNp_lee_sum");
  //histo_1eNp_lee_run123_sum->Sumw2();
  //histo_1eNp_lee_run123_sum->Add(histo_1eNp_lee_run2);
  //histo_1eNp_lee_run123_sum->Add(histo_1eNp_lee_run3);
  /*
  TH1D *histo_1eNp_bg_run1 = (TH1D*)frun1->Get("nu_uBooNE_1eNp_bg");
  TH1D *histo_1eNp_bg_run2 = (TH1D*)frun2->Get("nu_uBooNE_1eNp_bg");
  TH1D *histo_1eNp_bg_run3 = (TH1D*)frun3->Get("nu_uBooNE_1eNp_bg");
  TH1D *histo_1eNp_bg_runall = (TH1D*)frunall->Get("nu_uBooNE_1eNp_bg");
  TH1D *histo_1eNp_bg_run123_sum = (TH1D*)histo_1eNp_bg_run1->Clone("nu_uBooNE_1eNp_bg_sum");
  std::cout << "run1 1eNp bg abs err = " << histo_1eNp_bg_run123_sum->GetBinError(2) << std::endl;
  std::cout << "run1 1eNp bg = " << histo_1eNp_bg_run123_sum->GetBinContent(2) << std::endl;
  std::cout << "run1 1eNp bg frac err = " << histo_1eNp_bg_run123_sum->GetBinError(2)/histo_1eNp_bg_run123_sum->GetBinContent(2) << std::endl;
  histo_1eNp_bg_run123_sum->Sumw2();
  histo_1eNp_bg_run123_sum->Add(histo_1eNp_bg_run2);
  std::cout << "run1+2 1eNp bg abs err = " << histo_1eNp_bg_run2->GetBinError(2) << ", " << histo_1eNp_bg_run123_sum->GetBinError(2) << std::endl;
  std::cout << "run1+2 1eNp bg = " << histo_1eNp_bg_run2->GetBinContent(2) << ", " << histo_1eNp_bg_run123_sum->GetBinContent(2) << std::endl;
  std::cout << "run1+2 1eNp bg frac err = " << histo_1eNp_bg_run2->GetBinError(2)/histo_1eNp_bg_run2->GetBinContent(2) << ", " << histo_1eNp_bg_run123_sum->GetBinError(2)/histo_1eNp_bg_run123_sum->GetBinContent(2) << std::endl;
  histo_1eNp_bg_run123_sum->Add(histo_1eNp_bg_run3);
  std::cout << "run1+2+3 1eNp bg abs err = " << histo_1eNp_bg_run3->GetBinError(2) << ", " << histo_1eNp_bg_run123_sum->GetBinError(2) << std::endl;
  std::cout << "run1+2+3 1eNp bg = " << histo_1eNp_bg_run3->GetBinContent(2) << ", " << histo_1eNp_bg_run123_sum->GetBinContent(2) << std::endl;
  std::cout << "run1+2+3 1eNp bg frac err = " << histo_1eNp_bg_run3->GetBinError(2)/histo_1eNp_bg_run3->GetBinContent(2) << ", " << histo_1eNp_bg_run123_sum->GetBinError(2)/histo_1eNp_bg_run123_sum->GetBinContent(2) << std::endl;
  std::cout << "run123 1eNp bg abs err = " << histo_1eNp_bg_runall->GetBinError(2) << ", " << histo_1eNp_bg_runall->GetBinContent(2) << std::endl;
  std::cout << "run123 1eNp bg frac err = " << histo_1eNp_bg_runall->GetBinError(2)/histo_1eNp_bg_runall->GetBinContent(2) << ", " << histo_1eNp_bg_runall->GetBinError(2)/histo_1eNp_bg_runall->GetBinContent(2) << std::endl;
  
  TH1D *histo_1e0p_lee_run1 = (TH1D*)frun1->Get("nu_uBooNE_1e0p_lee");
  TH1D *histo_1e0p_lee_run2 = (TH1D*)frun2->Get("nu_uBooNE_1e0p_lee");
  TH1D *histo_1e0p_lee_run3 = (TH1D*)frun3->Get("nu_uBooNE_1e0p_lee");
  TH1D *histo_1e0p_lee_runall = (TH1D*)frunall->Get("nu_uBooNE_1e0p_lee");
  TH1D *histo_1e0p_lee_run123_sum = (TH1D*)histo_1e0p_lee_run1->Clone("nu_uBooNE_1e0p_lee_sum");
  histo_1e0p_lee_run123_sum->Sumw2();
  histo_1e0p_lee_run123_sum->Add(histo_1e0p_lee_run2);
  histo_1e0p_lee_run123_sum->Add(histo_1e0p_lee_run3);
  
  TH1D *histo_1e0p_bg_run1 = (TH1D*)frun1->Get("nu_uBooNE_1e0p_bg");
  TH1D *histo_1e0p_bg_run2 = (TH1D*)frun2->Get("nu_uBooNE_1e0p_bg");
  TH1D *histo_1e0p_bg_run3 = (TH1D*)frun3->Get("nu_uBooNE_1e0p_bg");
  TH1D *histo_1e0p_bg_runall = (TH1D*)frunall->Get("nu_uBooNE_1e0p_bg");
  TH1D *histo_1e0p_bg_run123_sum = (TH1D*)histo_1e0p_bg_run1->Clone("nu_uBooNE_1e0p_bg_sum");
  histo_1e0p_bg_run123_sum->Sumw2();
  histo_1e0p_bg_run123_sum->Add(histo_1e0p_bg_run2);
  histo_1e0p_bg_run123_sum->Add(histo_1e0p_bg_run3);

  TH1D *h_numu_data1 = (TH1D*)histo_1e0p_bg_run1->Clone("dummy");
  h_numu_data1->Reset();
*/
  TString tag1 = tag + mcerrlabel + zerobinlabel + detsyslabel + "_BFplots_1eNp"; 
  DrawDataMCAndSyst(tag1, histo_1eNp_lee_bf, histo_1eNp_bg_bf, histo_1eNp_data, "Reconstructed Visible Energy [GeV]", "#nu_{e} 1eNp0#pi selection");
  TString tag2 = tag + mcerrlabel + zerobinlabel + detsyslabel + "_BFplots_1e0p"; 
  DrawDataMCAndSyst(tag2, histo_1e0p_lee_bf, histo_1e0p_bg_bf, histo_1e0p_data, "Reconstructed Visible Energy [GeV]", "#nu_{e} 1e0p0#pi selection");
  
  return 0;
  
  
}

void DrawDataMCAndSyst(TString cName, TH1D* MC1, TH1D* MC2, TH1D* data, std::string var, std::string title){
  
  TCanvas can(cName,cName,1200,800);
  
  TH1D* tmpMC1 = (TH1D*)MC1->Clone("tempCV");
  TH1D* tmpMC2 = (TH1D*)MC2->Clone("tempCVnolee");
  TH1D* sysErr1 = (TH1D*)MC1->Clone("tempSysCV");
  TH1D* sysErr2 = (TH1D*)MC2->Clone("tempSysCVnolee");
  TH1D* tmpData = (TH1D*)data->Clone("tempData");

  TH1D* tmpratio1 = (TH1D*)MC1->Clone("tempRatio1");
  TH1D* tmpMC1_staterr = (TH1D*)MC1->Clone("tempCV_staterr2");
  for( int bin = 1; bin < tmpMC1_staterr->GetNbinsX()+1; bin++){ tmpMC1_staterr->SetBinError(bin,sqrt(tmpMC1_staterr->GetBinContent(bin))); }
  tmpratio1->Reset();
  tmpratio1->Divide(tmpMC2,tmpMC1_staterr,1.0,1.0,"B");
  if(verbose) std::cout << "nbins data, mc = " << tmpData->Integral() << ", " << tmpMC1->Integral() << std::endl;

  TH1D* tmpratio2 = (TH1D*)MC1->Clone("tempRatio2");
  TH1D* tmpMC2_staterr = (TH1D*)MC2->Clone("tempCV_staterr2");
  for( int bin = 1; bin < tmpMC2_staterr->GetNbinsX()+1; bin++){ tmpMC2_staterr->SetBinError(bin,sqrt(tmpMC2_staterr->GetBinContent(bin))); }
  tmpratio2->Reset();
  tmpratio2->Divide(tmpMC1,tmpMC2_staterr,1.0,1.0,"B");
  if(verbose) std::cout << "nbins data, mc nolee = " << tmpData->Integral() << ", " << tmpMC2->Integral() << std::endl;

  TH1D* tmpsys = (TH1D*)MC2->Clone("tempSysRatio");
  TH1D* tmpsys2 = (TH1D*)MC1->Clone("tempSysRatio1");

  for( int bin = 1; bin < tmpsys->GetNbinsX()+1; bin++){
    tmpsys->SetBinContent(bin,1.0);
    tmpsys->SetBinError(bin,tmpMC2->GetBinError(bin)/tmpMC2->GetBinContent(bin));
    std::cout << "tempSys MC2 " << tmpsys->GetBinError(bin) << std::endl;
    //if(bin==1) std::cout << "tempSys MC2 " << tmpMC2->GetBinError(bin) << "/" << tmpMC2->GetBinContent(bin) << " = " << tmpMC2->GetBinError(bin)/tmpMC2->GetBinContent(bin) << std::endl;
    tmpsys2->SetBinContent(bin,1.0);
    tmpsys2->SetBinError(bin,tmpMC1->GetBinError(bin)/tmpMC1->GetBinContent(bin));
    if(verbose) std::cout << "tempSys MC1 " << tmpsys2->GetBinError(bin) << std::endl;
    //tmpratio->SetBinError(bin,0.01);
    sysErr1->SetBinContent(bin,tmpMC1->GetBinContent(bin));
    sysErr1->SetBinError(bin,tmpMC1->GetBinError(bin));
    sysErr2->SetBinContent(bin,tmpMC2->GetBinContent(bin));
    sysErr2->SetBinError(bin,tmpMC2->GetBinError(bin));
    //std::cout << "temp ratio = " << tmpratio->GetBinContent(bin) << std::endl; 
  }
  for( int bin = 1; bin < tmpMC1->GetNbinsX()+1; bin++){ tmpMC1->SetBinError(bin,sqrt(tmpMC1->GetBinContent(bin))); }
  for( int bin = 1; bin < tmpMC2->GetNbinsX()+1; bin++){ tmpMC2->SetBinError(bin,sqrt(tmpMC2->GetBinContent(bin))); }
  
  // Upper plot will be in pad1
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.0, 1.0, 1.0);
  pad1->SetBottomMargin(0.13); // Upper and lower plot are joined
  pad1->SetGridx(2);        // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad
  //TH1 *frame= gPad->DrawFrame(100., 0.5, 1100., 1.5);
  double max = 1.1;
  double max1 = 1.1;
  double maxMC1 = tmpMC1->GetMaximum();
  double maxMC2 = tmpMC2->GetMaximum();
  double maxData = tmpData->GetMaximum();
  if(maxMC1 > maxData) max1 *= maxMC1;
  if(maxData > maxMC1 ) max1 *= maxData;
  if(maxMC2 > max1) max *= maxMC2;
  else max = max1;

  //title = title+"; Reconstructed Visible Energy [GeV]; Events/0.1 GeV";
  std::cout << "title: " << title << std::endl;
  if(std::string(cName).find("1eNp") != std::string::npos ) max=19;
  else if(std::string(cName).find("1e0p") != std::string::npos ) max=14.5;

  sysErr1->SetTitleSize(0.05);
  sysErr1->SetTitleFont(62);
  sysErr1->SetTitle(title.c_str());
  sysErr1->SetStats(0);
  //sysErr1->GetYaxis()->SetTitleSize(20);
  //sysErr1->GetYaxis()->SetTitleOffset(0.0);
  ///sysErr1->GetXaxis()->SetTitleOffset(0.0);
  //sysErr1->GetYaxis()->SetTitleFont(42);
  //sysErr1->GetXaxis()->SetTitleSize(20);
  //sysErr1->GetXaxis()->SetTitleFont(42);
  //sysErr1->GetXaxis()->SetLabelSize(0);

  Int_t trans_blue1 = TColor::GetColorTransparent(kBlue, 0.7);
  Int_t trans_red1 = TColor::GetColorTransparent(kRed, 0.7);

  tmpMC1->SetLineColor(trans_blue1);
  tmpMC1->SetLineWidth(3);
  tmpMC1->SetLineStyle(1);
  tmpMC2->SetLineColor(trans_red1);
  tmpMC2->SetLineWidth(3);
  tmpMC2->SetLineStyle(1);
  tmpData->SetLineWidth(2);
  tmpData->SetMarkerStyle(20);
  tmpData->SetLineColor(kBlack);
  if(std::string(cName).find("nowgt") != std::string::npos && std::string(cName).find("nue") != std::string::npos ) tmpData->SetLineColor(kBlack);

  Int_t trans_blue = TColor::GetColorTransparent(kBlue-10, 0.3);
  Int_t trans_red = TColor::GetColorTransparent(kRed-10, 0.3);

  sysErr1->SetLineColor(kBlue-10);
  sysErr1->SetLineWidth(1);
  sysErr1->SetFillColor(trans_blue);
  sysErr1->SetFillStyle(1001);

  sysErr2->SetLineColor(kRed-10);
  sysErr2->SetLineWidth(1);
  sysErr2->SetFillColor(trans_red);
  sysErr2->SetFillStyle(1001);

  sysErr1->SetMaximum(1.*max);
  sysErr1->SetMinimum(0.);

  sysErr1->GetYaxis()->SetTitleFont(62);
  sysErr1->GetYaxis()->SetTitleSize(0.04);
  sysErr1->GetYaxis()->SetTitle("Events/0.1 GeV");
  sysErr1->GetXaxis()->SetTitle("Reconstructed Visible Energy [GeV]");
  sysErr1->GetXaxis()->SetTitleSize(0.04);
  sysErr1->Draw("E2");
  sysErr2->Draw("E2 same");
  tmpMC1->Draw("hist same"); 
  tmpMC2->Draw("hist same"); 
  sysErr1->Draw("E2 same");
  tmpData->Draw("PE same");
  
  if(std::string(cName).find("nowgt") != std::string::npos && std::string(cName).find("nue") != std::string::npos ) tmpData->Draw("hist same"); 

  double x1=0.0; 
  double y1=0.0; 
  double x2=0.0; 
  double y2=0.0; 
  if(title == "#nu_{#mu} Selection"){
    x1=0.7;
    x2=0.9;
    y1=0.75;
    y2=0.9;
  }else{
    x1=0.6;
    x2=0.9;
    y1=0.7;
    y2=0.9;
  } 
  TLegend  *legend = new TLegend(x1,y1,x2,y2); // we need different positions for the legend to not 
  // get the plot titles for the legend
  if(title == "#nu_{#mu} Selection"){
    if(std::string(cName).find("nowgt") != std::string::npos ){
      legend->AddEntry(tmpMC1,"#nu_{#mu} CC, post-constraint","l"); 
      legend->AddEntry(tmpMC2,"#nu_{#mu} CC, no Genie wgt","l"); 
      legend->AddEntry(tmpData,"fake data","l");
    }else{
      legend->AddEntry(tmpMC1,"old tune, unconstrained","l"); 
      legend->AddEntry(tmpMC2,"new tune, post-constraint","l");
    }
  }else{
    if(std::string(cName).find("nowgt") != std::string::npos ){
      legend->AddEntry(tmpMC1,"#nu_e, Genie Wgt","l"); 
      legend->AddEntry(tmpMC2,"#nu_e, no Genie wgt","l"); 
      legend->AddEntry(tmpData,"#nu_e, post-constraint","l");
    }else if(std::string(cName).find("before") != std::string::npos ){
      legend->AddEntry(tmpMC1,"Old Tune, unconstrained","l"); 
      legend->AddEntry(tmpMC2,"New Tune, unconstrained","l"); 
    }else if(std::string(cName).find("BF") != std::string::npos ){
      legend->AddEntry(tmpMC1,"LEE, #mu = 3.05","l"); 
      legend->AddEntry(tmpMC2,"#nu_{e} intrinsic+BG","l"); 
      legend->AddEntry(tmpData,"fakedata","pel"); 
    }else if(std::string(cName).find("trueSignal") != std::string::npos ){
      legend->AddEntry(tmpMC1,"LEE, #mu = 3.5","l"); 
      legend->AddEntry(tmpMC2,"#nu_{e} intrinsic+BG","l"); 
      legend->AddEntry(tmpData,"fakedata","pel"); 
    }else if(std::string(cName).find("unconstrained") != std::string::npos ){
      legend->SetHeader("Before Constraint","C");
      legend->AddEntry(tmpMC1,Form("Run by Run, #nu_{e}+bkg: %.2f",tmpMC1->Integral()),"le"); 
      legend->AddEntry(tmpMC2,Form("Run 3, #nu_{e}+bkg: %.2f",tmpMC2->Integral()),"le"); 
    }else if(std::string(cName).find("constrained") != std::string::npos ){
      legend->SetHeader("Post-constraint","C");
      legend->AddEntry(tmpMC1,Form("Run by Run, #nu_{e}+bkg: %.2f",tmpMC1->Integral()),"le"); 
      legend->AddEntry(tmpMC2,Form("Run 3, #nu_{e}+bkg: %.2f",tmpMC2->Integral()),"le"); 
    }
  } 
  legend->Draw("same");
/*
  // lower plot will be in pad
  can.cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.45);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.25);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad

  tmpsys2->SetTitle("");
  tmpsys2->GetYaxis()->SetTitle("Run 3/Run 1-3");
  tmpsys2->GetXaxis()->SetTitle(var.c_str());
  tmpsys2->SetStats(0);
  tmpsys->SetLineColor(kRed-10);
  tmpsys->SetLineWidth(1);
  tmpsys->SetFillColor(kRed-10);
  tmpsys->SetFillStyle(1001);
  tmpsys2->SetLineColor(kBlue-10);
  tmpsys2->SetLineWidth(1);
  tmpsys2->SetFillColor(kBlue-10);
  tmpsys2->SetFillStyle(1001);
  tmpsys2->GetYaxis()->SetTitleSize(25);
  tmpsys2->GetYaxis()->SetTitleFont(43);
  tmpsys2->GetXaxis()->SetTitleSize(25);
  tmpsys2->GetXaxis()->SetTitleOffset(3.5);
  tmpsys2->GetYaxis()->SetTitleOffset(1.1);
  tmpsys2->GetXaxis()->SetLabelSize(0.1);
  tmpsys2->GetYaxis()->SetLabelSize(0.08);
  tmpsys2->GetYaxis()->SetNdivisions(210, kTRUE);
  tmpsys2->GetXaxis()->SetTitleFont(43);
  std::cout << "@@@@@@ tmpratio2, tmpratio1 = " << tmpratio2->GetMaximum() << ", " << tmpratio1->GetMaximum() << ", " << tmpratio2->GetMinimum() << ", " << tmpratio1->GetMinimum() << std::endl;
  if(tmpratio2->GetMaximum() >= tmpratio1->GetMaximum() ) tmpsys2->SetMaximum(1.1*tmpratio2->GetMaximum()); 
  if(tmpratio1->GetMaximum() >= tmpratio2->GetMaximum() ) tmpsys2->SetMaximum(1.1*tmpratio1->GetMaximum()); 
  if(tmpratio2->GetMinimum() <= tmpratio1->GetMinimum() ) tmpsys2->SetMinimum(0.8*tmpratio2->GetMinimum()); 
  if(tmpratio1->GetMinimum() <= tmpratio2->GetMinimum() ) tmpsys2->SetMinimum(0.8*tmpratio1->GetMinimum()); 
  //tmpsys2->SetMaximum(3.0);
  if(title == "#nu_{#mu} Selection") tmpsys2->SetMaximum(1.5);
  tmpsys2->SetMinimum(0.6);
  tmpsys2->SetMaximum(1.4);
  if(title == "#nu_{#mu} Selection") tmpsys2->SetMinimum(0.5);
  std::cout << "@@@@@@ tmpSys max, min = " << tmpsys2->GetMaximum() << ", " << tmpsys2->GetMinimum() << std::endl;
  tmpsys2->Draw("E2");
  tmpsys->Draw("E2 same");
  //tmpsys2->Draw("E2 same");

  const TAxis *axis = tmpratio1->GetXaxis();
  double lowX  = axis->GetBinLowEdge( axis->GetFirst() );
  double highX = axis->GetBinUpEdge(  axis->GetLast() );

  TLine line;
  line.SetLineStyle(2);
  line.SetLineWidth(2);
  line.SetLineColor(kRed);
  line.DrawLine(lowX, 1., highX, 1.); //creates a new line which is owned by gPad

  tmpratio1->SetMarkerStyle(20);
  //tmpratio1->SetLineColor(kBlue);
  tmpratio1->SetLineColor(kBlack);
  tmpratio1->SetLineWidth(2);
  tmpratio1->Draw("PE same");

  for(int bin=1; bin < tmpratio1->GetNbinsX(); bin++ ) std::cout << "tmpratio1, tmpratio2  = " << tmpratio1->GetBinContent(bin) << tmpratio2->GetBinContent(bin) << std::endl;
  tmpratio2->SetMarkerStyle(20);
  tmpratio2->SetLineColor(kRed);
  tmpratio2->SetLineWidth(2);
  //tmpratio2->Draw("PE same");

  pad2->Update();
  pad2->Modified(); // so it updates the pad with the new changes
  pad2->Draw("");
*/
  //gPad->RedrawAxis();
  //gPad->Update();
  //can.Update();
 
  can.Print(cName+".pdf","pdf");

}
void plot_one(TMatrixD matrix, std::string tag){
    gStyle->SetPalette(kLightTemperature);
    std::vector<TString> channel_names = {"1eNp0pi","1e0p0pi", "nu_mu"};
    int num_channels = 3;
    int num_bins_total = matrix.GetNrows();
    std::vector<int> num_bins = {14,14,14}; 

    TH2D h2_full(matrix);
    h2_full.SetName((tag+"_th2d").c_str());
    TCanvas *c_full = new TCanvas((tag+"_canvas").c_str());
    TPad *p_full = (TPad*)c_full->cd();
    c_full->SetFixedAspectRatio();
    h2_full.Draw("colz");
    h2_full.SetTitle(tag.c_str());
    h2_full.GetXaxis()->SetTitle("Global Bin Number");
    h2_full.GetYaxis()->SetTitle(" ");
    h2_full.GetYaxis()->SetLabelSize(0);

    c_full->SetFrameFillColor(kWhite);
    c_full->SetFillColor(kWhite);
    p_full->SetFillColor(kWhite);

    c_full->SetRightMargin(0.150);
    c_full->SetLeftMargin(0.250);
    c_full->SetTopMargin(0.10);
    int use_full =0;

    double percent_left = 0.15;
    double nice_shift = num_bins_total*0.02;

    for(int ic = 0; ic < num_channels; ic++){

        TText * tmd = new TText(-num_bins_total*percent_left*0.15, use_full+nice_shift*0.5, (channel_names[ic]) );
        tmd->SetTextColor(kBlack);
        tmd->SetTextSize(0.03);
        tmd->SetTextAlign(31);
        tmd->Draw();
       
        TLine *lv = new TLine(-num_bins_total*percent_left, num_bins.at(ic)+use_full, num_bins_total, num_bins.at(ic)+use_full);
        TLine *lh = new TLine(num_bins.at(ic)+use_full,0, num_bins.at(ic)+use_full, num_bins_total*1.045);
        lv->SetLineWidth(3);
        lh->SetLineWidth(3);
        lv->SetLineColor(kRed);
        lh->SetLineColor(kRed);
        use_full+=num_bins.at(ic);
        lv->Draw();
        lh->Draw();

    }
    gStyle->SetOptStat(0);
    std::cout << "print matrix...." << std::endl;
    c_full->Print((tag+"_collapsed.pdf").c_str(),"pdf");
}
