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
#include "TPaveStats.h"
#include "TH1.h"
#include "TH1C.h"
#include "TString.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TMath.h"
#include "TSystem.h"
#include "TMatrixT.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TError.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TObjArray.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TH2F.h"
#include "TImage.h"
#include "TGraph.h"
#include "TLatex.h"

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


void DrawDataMCAndSyst(TString cName, TH1D* MC, TH1D* data, std::string var, std::string title, double chi2);
float CalcChi(float **invert_matrix, float* core, float *sig, int bins);

int main(int argc, char* argv[])
{
  std::string xml = "example.xml";
  std::string tag = "example";
  std::string fakedata = "example";
  int constr_mode = 0;
  double width = 1.0;
  std::string var = "var";
  int iarg = 0;
  opterr=1;
  int index;
  bool sys = false;
  bool detsys = false;
  bool combined = false;
  bool np = false;
  bool zp = false;
  bool mcerr = false;
  bool addext = false;
  bool fc = false;
  
  const struct option longopts[] =
    {
      {"xml", 		required_argument, 	0, 'x'},
      {"tag", 		required_argument, 	0, 't'},
      {"fakedata", 	required_argument, 	0, 'f'},
      {"var", 		required_argument, 	0, 'v'},
      {"constr_mode",	required_argument, 	0, 'l'},
      {"sys",	        no_argument, 		0, 's'},
      {"detsys",	no_argument, 		0, 'd'},
      {"combined",	no_argument, 		0, 'c'},
      {"np",	        no_argument, 		0, 'n'},
      {"zp",	        no_argument, 		0, 'z'},
      {"mcerr",	        no_argument, 		0, 'm'},
      {"addext",	no_argument, 		0, 'a'},
      {"fc",	        no_argument, 		0, 'y'},
      {0,		no_argument, 		0,  0},
    };
  
  while(iarg != -1)
    {
      iarg = getopt_long(argc,argv, "x:t:f:l:v:sdcnzmay", longopts, &index);
      
      switch(iarg)
	{
	case 'x':
	  xml = optarg;
	  break;
	case 't':
	  tag = optarg;
	  break;
	case 'f':
	  fakedata = optarg;
	  break;
	case '?':
	case 'h':
	  std::cout<<"Allowed arguments:"<<std::endl;
	  std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
	  return 0;
	}
    }
  
  //Load up the central value spectra we computed in example 1, to act as a signal spectra
  SBNspec sig_spectra("../bin/"+tag+".SBNspec.root",xml);
  SBNchi chi_h1(sig_spectra);	
  
  //define bins number
  //the xml assumes that we always have numu channel
  //and at least a 1e0p channel or 1eNp channel 
  int num_bins = sig_spectra.num_bins_total;
  int num_bins_coll = sig_spectra.num_bins_total_compressed;
  
  //Load up our covariance matricies we calculated in example1 (we could also load up single variation ones)
  TFile * fsys = new TFile(Form("../bin/%s.SBNcovar.root",tag.c_str()),"read");
  TMatrixD *cov = (TMatrixD*)fsys->Get("collapsed_covariance");
  TMatrixD *full_cov = (TMatrixD*)fsys->Get("full_covariance");
  TMatrixD *fraccov = (TMatrixD*)fsys->Get("collapsed_frac_covariance");
  
  //data 
  std::vector<double> datavec;
  TFile *f_data = new TFile(fakedata.c_str(),"read");
  TH1D* histo_leptonE = (TH1D*)f_data->Get("nu_uBooNE_1eNp_leptonenergy_data");
  TH1D* histo_leptonAngle = (TH1D*)f_data->Get("nu_uBooNE_1eNp_leptonangle_data");
  TH1D* histo_protonE = (TH1D*)f_data->Get("nu_uBooNE_1eNp_protonenergy_data");
  TH1D* histo_protonAngle = (TH1D*)f_data->Get("nu_uBooNE_1eNp_protonangle_data");

  for(int bin = 1; bin <= histo_leptonE->GetNbinsX(); bin++ ) datavec.push_back(histo_leptonE->GetBinContent(bin));
  for(int bin = 1; bin <= histo_leptonAngle->GetNbinsX(); bin++ ) datavec.push_back(histo_leptonAngle->GetBinContent(bin));
  for(int bin = 1; bin <= histo_protonE->GetNbinsX(); bin++ ) datavec.push_back(histo_protonE->GetBinContent(bin));
  for(int bin = 1; bin <= histo_protonAngle->GetNbinsX(); bin++ ) datavec.push_back(histo_protonAngle->GetBinContent(bin));
  
  //collapse vector
  sig_spectra.CollapseVector();
  //collapse matrix

  TH1D *collapsed_mc = new TH1D("1eNp_othogonals_mc","1eNp_othogonals_mc",num_bins_coll, 0., num_bins_coll);
  for( int bin = 1; bin <= collapsed_mc->GetNbinsX(); bin++ ){ collapsed_mc->SetBinContent(bin, sig_spectra.collapsed_vector[bin-1]); collapsed_mc->SetBinError(bin, sqrt(sig_spectra.collapsed_vector[bin-1]));} 
  TH1D *collapsed_data = new TH1D("1eNp_othogonals_data","1eNp_othogonals_data",num_bins_coll, 0., num_bins_coll);
  for( int bin = 1; bin <= collapsed_data->GetNbinsX(); bin++ ){ collapsed_data->SetBinContent(bin, datavec[bin-1]); collapsed_data->SetBinError(bin, sqrt(datavec[bin-1])); } 

  TMatrixD Mcov(num_bins_coll,num_bins_coll);
  for(int i=0; i < num_bins_coll; i++){
    for(int j=0; j < num_bins_coll; j++){
      if(i==j) Mcov(i,j) = (*cov)(i,j) + datavec[i];
    }
  }
  //chi_h1.CollapseModes(Mcov, Mcoll);

  TMatrixT<double> invcovmatrix = chi_h1.InvertMatrix(Mcov);
  std::cout << "invert matrix: " << std::endl;
  invcovmatrix.Print();
  float** invmatrix = new float*[num_bins_coll];
  for(int i=0; i < num_bins_coll; i++){
    invmatrix[i] = new float[num_bins_coll];
  }
  for(int i=0; i < num_bins_coll; i++){
    for(int j=0; j < num_bins_coll; j++){
      invmatrix[i][j] = invcovmatrix[i][j];
    }
  }
  float* mc = new float[num_bins_coll];
  float* data = new float[num_bins_coll];
  std::vector<double> data_vec;
  data_vec.resize(num_bins_coll);
  for(int i=0; i < num_bins_coll; i++){
    mc[i] = sig_spectra.collapsed_vector[i];
    data[i] = datavec[i];
  }
  float chi2 = CalcChi(invmatrix,mc,data,num_bins_coll);
  
  DrawDataMCAndSyst(tag, collapsed_mc, collapsed_data, "shower energy [GeV] | shower angle [rad] | proton energy [GeV] | proton angle [rad]", "1eNp0#pi #nu_{e} Selection", chi2);

  //fill histograms
  TFile * fspec1 = new TFile(Form("collapsed_%s_MC.SBNspec.root",tag.c_str()),"recreate");
  fspec1->cd();
  (TH1D*)collapsed_mc->Write("nu_uBooNE_1eNp_mc");
  fspec1->Close();
  
  TFile * fspec2 = new TFile(Form("collapsed_%s_Data.SBNspec.root",tag.c_str()),"recreate");
  fspec2->cd();
  (TH1D*)collapsed_data->Write("nu_uBooNE_1eNp_data");
  fspec2->Close();

  TFile * fcoll = new TFile(Form("collapsed_%s.SBNcovar.root",tag.c_str()),"recreate");
  fcoll->cd();
  (TMatrixD*)cov->Write("full_covariance");
  (TMatrixD*)fraccov->Write("frac_covariance");
  fcoll->Close();
   
  return 0;
  
}


void DrawDataMCAndSyst(TString cName, TH1D* MC, TH1D* data, std::string var, std::string title, double chi2){
 
  std::cout << "cname: " << cName << std::endl; 
  TCanvas can(cName,cName,1200,800);
  
  TH1D* tmpMC = (TH1D*)MC->Clone("tempCV");
  TH1D* sysErr = (TH1D*)MC->Clone("tempSysCV");
  TH1D* tmpData = (TH1D*)data->Clone("tempData");

  TH1D* tmpMC_staterr = (TH1D*)MC->Clone("tempCV_staterr2");
  for( int bin = 1; bin < tmpMC_staterr->GetNbinsX()+1; bin++){ tmpMC_staterr->SetBinError(bin,sqrt(tmpMC_staterr->GetBinContent(bin))); }

  // Upper plot will be in pad1
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.06, 1.0, 1.0);
  pad1->SetGridx(2);        // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad

  double max = 1.1;
  double max1 = 1.1;
  double maxMC = tmpMC->GetMaximum();
  double maxData = tmpData->GetMaximum();
  if(maxMC > maxData) max *= maxMC;
  else if(maxData > maxMC ) max *= maxData;

  std::cout << "cName = " << std::endl;
  sysErr->SetTitleSize(80);
  sysErr->SetTitleFont(63);
  sysErr->SetTitle(title.c_str());
  sysErr->SetStats(0);
  sysErr->GetYaxis()->SetTitleSize(30);
  sysErr->GetYaxis()->SetTitleFont(43); //apparently changing 43->43 causes the title to not be printed on canvas
  sysErr->GetYaxis()->SetTitle("Events/0.1 GeV");
  sysErr->GetXaxis()->SetTitle(var.c_str());
  sysErr->GetXaxis()->SetTitleFont(43);
  sysErr->GetXaxis()->SetTitleSize(30);
  sysErr->GetXaxis()->SetTitleOffset(1.25);
  sysErr->GetXaxis()->SetLabelSize(0.04);

  tmpMC->SetLineColor(kRed);
  tmpMC->SetLineWidth(3);
  tmpMC->SetLineStyle(1);
  tmpData->SetLineWidth(2);
  tmpData->SetMarkerStyle(20);
  tmpData->SetLineColor(kBlack);

  sysErr->SetLineColor(kRed-10);
  sysErr->SetLineWidth(1);
  sysErr->SetFillColor(kRed-10);
  sysErr->SetFillStyle(1001);

  sysErr->SetMaximum(1.2*max);
  sysErr->SetMinimum(0.);

  sysErr->Draw("E2");
  tmpMC->Draw("hist same"); 
  tmpData->Draw("PE same"); 

  TLatex latex;
  TString chisq = Form("#chi^{2} = %.2f",chi2);
  TString chisq_dof = Form("#chi^{2}/%d = %.2f",tmpMC->GetNbinsX(),chi2/float(tmpMC->GetNbinsX()));
  latex.SetTextSize(0.05);
  latex.SetTextAlign(22);  //align at top
  latex.SetTextFont(22);  //align at top
  double pos1 = 1.4*tmpData->GetMaximum();
  double pos2 = 1.6*tmpData->GetMaximum();
  if(chi2>0) latex.DrawLatex(10.,pos1,chisq);
  latex.SetTextSize(0.05);
  latex.SetTextAlign(22);  //align at top
  latex.SetTextFont(22);  //align at top
  if(chi2>0) latex.DrawLatex(10.,pos2,chisq_dof);

  double x1=0.0; 
  double y1=0.0; 
  double x2=0.0; 
  double y2=0.0; 
  x1=0.6;
  x2=0.9;
  y1=0.6;
  y2=0.9;

  TLegend  *legend = new TLegend(x1,y1,x2,y2); // we need different positions for the legend to not 
  // get the plot titles for the legend
  legend->AddEntry(tmpMC,Form("#nu_{e}: %.2f",tmpMC->Integral()),"le"); 
  legend->AddEntry(tmpData,Form("data: %.2f",tmpData->Integral()),"pe");
  legend->Draw("same");

  gPad->RedrawAxis();
  gPad->Update();
 
  can.Print(cName+".pdf","pdf");

  return;

}

float CalcChi(float **invert_matrix, float* core, float *sig, int bins){
    float tchi = 0;

    for(int i =0; i<bins; i++){
        for(int j =0; j<bins; j++){
            tchi += (core[i]-sig[i])*invert_matrix[i][j]*(core[j]-sig[j] );
            if(i==j)std::cout << "1. invert_matrix: " << invert_matrix[i][j] << ", " << core[i] << ", " << sig[i] << ", " << core[j] << ", " << sig[j] << ", " << tchi << std::endl;
        }
    }
    std::cout << "tchi = " << tchi << std::endl;

   return tchi;
}

