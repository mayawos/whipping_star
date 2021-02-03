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


TH1D* GetTotalError(std::string name, TMatrixD errMatrix, TH1D* histo );
void DrawMCAndSyst(TString cName, TH1D* MC, std::string var, std::string title);
void DrawDataMCAndSyst(TString cName, TH1D* MC, TH1D* data, std::string var, std::string title);
void DrawMCAndSystOverlay(TString cName, TH1D* MC1, TH1D* MC2, std::string var, std::string title);
void CollapseSubchannels(TMatrixD & M, TMatrixD & Mc, std::vector<int> num_bins, int num_channels, std::vector<int> num_subchannels);
void Draw_Stacked(TObjArray histos, TH1D *data, TString samplename, TPad *pad = 0, bool normalised = false, std::string stacktitle = "", std::string var = "" );
void plot_one(TMatrixD matrix, TH1D *h_nue, std::string tag);
double CalcChi(TH1D *h_data, TH1D *h_mc, TMatrixD inv_matrix, TH1D *h_nue);
void make_separate_stackedplots(bool detsys, std::vector<TH1D*> histos, std::string channelname, std::string tag, std::string var, TFile *file, TMatrixD *cov, TMatrixD *fraccov, TMatrixD *corr );
//void make_separate_stackedplots(bool detsys, std::vector<TH1D*> histos, std::string channelname, std::string tag, std::string var, TFile *file, TMatrixD *cov, TMatrixD *fraccov, TMatrixD *corr, TString filefake );
void AddPlotLabel( const char * 	label,
                   const double 	x,
                   const double 	y,
                   const double 	size = 0.05,
                   const int 	color = 1,
                   const int 	font = 43,
                   const int 	align = 22,
                   const double 	angle = 0 );




int main(int argc, char* argv[])
{
  std::string xml = "example.xml";
  std::string tag = "example";
  double width = 1.0;
  std::string var = "var";
  int iarg = 0;
  opterr=1;
  int index;
  bool sys = false;
  bool detsys = false;
  int mass_start = -1;
  std::cout << "A" << std::endl; 
  const struct option longopts[] =
    {
      {"xml", 		required_argument, 	0, 'x'},
      {"tag", 		required_argument, 	0, 't'},
      {"var", 		required_argument, 	0, 'v'},
      {"sys",	no_argument, 0, 's'},
      {"detsys",	no_argument, 0, 'd'},
      {"part", required_argument,0,'p'},
      {0,			no_argument, 		0,  0},
    };
  
  std::cout << "B" << std::endl; 
  while(iarg != -1)
    {
      iarg = getopt_long(argc,argv, "x:t:v:scp:d", longopts, &index);
      
      switch(iarg)
	{
	case 'x':
  std::cout << "B1" << std::endl; 
	  xml = optarg;
	  break;
	case 't':
  std::cout << "B2" << std::endl; 
	  tag = optarg;
	  break;
	case 'v':
  std::cout << "B4" << std::endl; 
	  var = optarg;
	  break;
	case 's':
  std::cout << "B5" << std::endl; 
	  sys = true;
	  break;
	case 'd':
  std::cout << "B6" << std::endl; 
	  detsys = true;
	  break;
	case 'p':
  std::cout << "B7" << std::endl; 
	  mass_start = atoi(optarg);
	  break;
	case '?':
	case 'h':
	  std::cout<<"Allowed arguments:"<<std::endl;
	  std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
	  return 0;
	}
    }
  std::cout << "C" << std::endl; 
  
  //fetch the TFile for the MC spectrum
  TString fname = "../bin/"+tag+".SBNspec.root"; 
  std::cout << "fname = " << fname << std::endl;
  TFile *file = new TFile(fname,"read");
  std::vector<TH1D*> histos_1eNp, histos_1e0p, histos_numu;
  
  //loop counter(s)
  int i=0;
  int j=0;
  TKey *key;
  TObject *obj;
  TIter nextkey(gDirectory->GetListOfKeys());
  while (key = (TKey*)nextkey()) {
    std::cout << "i: " << i << std::endl;
    ++i;
    obj = key->ReadObj();
    std::string name = obj->GetName();
    std::cout << "name: " << name << std::endl;
    TH1D* h = (TH1D*)obj;
    std::cout << h->GetName() << std::endl; 
    std::cout << "j: " << j << std::endl;
    j++;
    if(name.find("nu_uBooNE_nue") != std::string::npos ) histos_1eNp.push_back(h);
    if(name.find("nu_uBooNE_1e0p") != std::string::npos ) histos_1e0p.push_back(h);
    if(name.find("nu_uBooNE_numu") != std::string::npos ) histos_numu.push_back(h);
  }
  
  //Load up our covariance matricies we calculated in example1 (we could also load up single variation ones)
  TString filename = "../bin/"+tag+".SBNcovar.root";
  std::cout << "filename = " << filename << std::endl;
  TFile * fsys = new TFile(filename,"read");
  TMatrixD *cov = (TMatrixD*)fsys->Get("full_covariance");
  TMatrixD *fraccov = (TMatrixD*)fsys->Get("frac_covariance");
  TMatrixD *corr = (TMatrixD*)fsys->Get("full_correlation");

  if(histos_1eNp.size() > 0) make_separate_stackedplots(detsys, histos_1eNp, "1eNp", tag, var, file, cov, fraccov, corr);
  if(histos_1e0p.size() > 0) make_separate_stackedplots(detsys, histos_1e0p, "1e0p", tag, var, file, cov, fraccov, corr);
  if(histos_numu.size() > 0) make_separate_stackedplots(detsys, histos_numu, "numu", tag, var, file, cov, fraccov, corr);
  
  return 0;
  
}

void make_separate_stackedplots(bool detsys, std::vector<TH1D*> histos, std::string channelname, std::string tag, std::string var, TFile *file, TMatrixD *cov, TMatrixD *fraccov, TMatrixD *corr ){
  
  //histo per sample
  TH1D *h_nue_data = (TH1D*)histos[0]->Clone(Form("nuedatalee_%s",channelname.c_str())); //assume that the nue N_data == nue N_mc     
  h_nue_data->Reset();
  TH1D *h_nue_lee = (TH1D*)histos[0]->Clone(Form("nuedataleeonly_%s",channelname.c_str())); //assume that the nue N_data == nue N_mc     
  h_nue_lee->Reset();
  for( int h = 0; h < histos.size(); h++ ){
    std::string name = histos[h]->GetName();
    if(name.find("_lee") != std::string::npos ) h_nue_lee->Add(histos[h]);
    h_nue_data->Add(histos[h]);
  }
  h_nue_data->Sumw2(kFALSE);
  
  std::vector<double> input_vec_mc_total, input_vec;
  for(int b = 1; b < h_nue_data->GetNbinsX()+1; b++) input_vec_mc_total.push_back(h_nue_data->GetBinContent(b));
 
  int leestart = -1; 
  int leeend = -1; 
  int sigstart = -1; 
  int sigend = -1; 
  int extstart = -1; 
  int extend = -1; 
  for( int h = 0; h < histos.size(); h++ ){
    std::string histname = histos[h]->GetName();
    if(histname.find("_lee") != std::string::npos ){
      leestart = h*h_nue_data->GetNbinsX();
      leeend = (h+1)*h_nue_data->GetNbinsX();
    }
    if(histname.find("_intrinsic") != std::string::npos ){
      sigstart = h*h_nue_data->GetNbinsX();
      sigend = (h+1)*h_nue_data->GetNbinsX();
    }
    if(histname.find("_extbnb") != std::string::npos ){
      extstart = h*h_nue_data->GetNbinsX();
      extend = (h+1)*h_nue_data->GetNbinsX();
    }
    for( int b=1; b < h_nue_data->GetNbinsX()+1; b++ ){
      std::cout << channelname << " data bin, bin content, bin error = " << b << ", " << h_nue_data->GetBinContent(b) << ", " << h_nue_data->GetBinError(b) << std::endl;
      input_vec.push_back(histos[h]->GetBinContent(b));
    }
  }

  std::cout << "extbnb startbin, endbin = " << extstart << ", " << extend << std::endl;
  
  //Make TObjArray
  TObjArray hist_nue_array;
  for(int h = 0; h < histos.size(); h++ ){
    histos[h]->SetTitle(histos[h]->GetName());
    hist_nue_array.Add(histos[h]);
  }

  std::cout << "1..." << std::endl;
  TString nuesamplename = tag+"_"+channelname+"_withpoiserr_Stacked"; 
  std::string sampletitle;
  if(channelname=="1eNp") sampletitle = "#nu_e 1eNp0#pi Selection";
  if(channelname=="1e0p") sampletitle = "#nu_{e} 1e0p0#pi Selection";
  if(channelname=="numu") sampletitle = "#nu_#mu Contained Selection";
  /*if(channelname=="numu"){
        //fetch the TFile for the data and MC spectrum
        //TString fname = filefake.c_str(;
        TFile *f_data = new TFile(filefake,"read");
        TH1D *h_numu_data = (TH1D*)f_data->Get("nu_uBooNE_numu_data");
        h_numu_data->SetName("nu_uBooNE_numu_data");
        h_numu_data->SetTitle("nu_uBooNE_numu_data");
        TString cName = "numu_stacked_with_data_overlay";
        Draw_Stacked(hist_nue_array,h_numu_data,nuesamplename, 0, false, sampletitle, var);	
  }
  else if(channelname=="1eNp"){
        //fetch the TFile for the data and MC spectrum
        //TString fname = filefake.c_str();
        TFile *f_data = new TFile(filefake,"read");
        TH1D *h_1eNp_data = (TH1D*)f_data->Get("nu_uBooNE_nue_data");
        h_1eNp_data->SetName("nu_uBooNE_nue_intrinsic");
        h_1eNp_data->SetTitle("nu_uBooNE_nue_data");
        TString cName = "1eNp_stacked_with_data_overlay";
        Draw_Stacked(hist_nue_array,h_1eNp_data,nuesamplename, 0, false, sampletitle, var);	
  }
  else if(channelname=="1e0p"){
        //fetch the TFile for the data and MC spectrum
        //TString fname = "../bin/1e0p_numu_reco_e_H1_fakedata2.SBNspec.root";
        TFile *f_data = new TFile(filefake,"read");
        TH1D *h_1e0p_data = (TH1D*)f_data->Get("nu_uBooNE_1e0p_data");
        h_1e0p_data->SetName("nu_uBooNE_1e0p_intrinsic");
        h_1e0p_data->SetTitle("nu_uBooNE_1e0p_data");
        TString cName = "1e0p_stacked_with_data_overlay";
        Draw_Stacked(hist_nue_array,h_1e0p_data,nuesamplename, 0, false, sampletitle, var);	
  }
  else*/ Draw_Stacked(hist_nue_array, h_nue_data, nuesamplename, 0, false, sampletitle, var );
  //create SBNspec and SBNchi object
  //========================================
  /*TString bgname = "../bin/"+tag+".SBNspec.root";
  std::cout << bgname << std::endl;
  SBNspec bg(bgname.Data(),Form("../bin/%s.xml",tag.c_str()));
  bg.CalcFullVector();

  std::cout << "segfault 1" << std::endl;
  if(channelname=="1eNp"){
    for(auto& h: bg.hist){
      std::string hname = (std::string)h.GetName();
      std::cout << "hname 1eNp = " << hname << std::endl;
      if(hname.find("_lee") != std::string::npos ) bg.Scale(hname.c_str(),0.0);
      if(hname.find("nu_uBooNE_1e0p") != std::string::npos ) bg.Scale(hname.c_str(),0.0);
      if(hname.find("nu_uBooNE_numu") != std::string::npos ) bg.Scale(hname.c_str(),0.0);
    }
  }
  if(channelname=="1e0p"){
    for(auto& h: bg.hist){
      std::string hname = (std::string)h.GetName();
      if(hname.find("_lee") != std::string::npos ) bg.Scale(hname.c_str(),0.0);
      if(hname.find("nu_uBooNE_1eNp") != std::string::npos ) bg.Scale(hname.c_str(),0.0);
      if(hname.find("nu_uBooNE_numu") != std::string::npos ) bg.Scale(hname.c_str(),0.0);
    }
  }
  if(channelname=="numu"){
    for(auto& h: bg.hist){
      std::string hname = (std::string)h.GetName();
      if(hname.find("nu_uBooNE_1eNp") != std::string::npos ) bg.Scale(hname.c_str(),0.0);
      if(hname.find("nu_uBooNE_1e0p") != std::string::npos ) bg.Scale(hname.c_str(),0.0);
    }
  }

  SBNchi chi_h0(bg,fraccov);
 
  TMatrixD Mtotal(cov->GetNcols(), cov->GetNrows());
  TMatrixD Mtotalfraccov(fraccov->GetNcols(), fraccov->GetNrows());
  TMatrixD Mtotalcorr(corr->GetNcols(), corr->GetNrows());
  TMatrixD Mdetsys(cov->GetNcols(), cov->GetNrows());
  Mdetsys.Zero();
 
  Mtotal.Zero();
  Mtotalfraccov.Zero();
  Mtotalcorr.Zero();
 
  bool useBDT=false;
  if(tag.find("BDT") != std::string::npos) useBDT=true;
  if(detsys) chi_h0.FillDetSysMatrix(Mdetsys, bg,useBDT);
  Mtotal = *cov;
  for(int i = 0; i < Mtotal.GetNrows(); i++ ){ if(Mdetsys(i,i) !=0 ) std::cout << "Mdetsys before = " << Mtotal(i,i) << std::endl; }
  Mtotal += Mdetsys;
  for(int i = 0; i < Mtotal.GetNrows(); i++ ){ if(Mdetsys(i,i) !=0 )std::cout << "Mdetsys after = " << Mtotal(i,i) << std::endl; }

  std::vector<double> bg_vector = bg.full_vector;
  std::cout << bg_vector.size() << std::endl;
  for(int i=0; i<cov->GetNcols(); i++){
    for(int j=0; j<cov->GetNrows(); j++){
      if(Mtotal(i,j) != Mtotal(i,j) || std::isnan(Mtotal(i,j)) ) Mtotal = 0.0;
      if(bg_vector[i] != 0 || bg_vector[j] !=0 ) Mtotalfraccov(i,j) = Mtotal(i,j)/(bg_vector[i]*bg_vector[j]);
      if(bg_vector[i] != 0 || bg_vector[j] !=0 ) Mtotalcorr(i,j) = Mtotal(i,j)/(sqrt(Mtotal(i,i))*sqrt(Mtotal(j,j)));
      if(i==j) std::cout << "bin, uncertainties = " << i+1 << ", " << 100*sqrt(Mtotalfraccov(i,j)) << std::endl;
    }
  }

  //add detsys
  TString mname = Form("../bin/%s_nodetsys.SBNcovar.root",tag.c_str());
  if(detsys) mname = Form("../bin/%s_detsys.SBNcovar.root",tag.c_str());
  TFile *covarfile = new TFile(mname.Data(),"recreate");
  covarfile->cd();

  std::cout << "segfault 7" << std::endl;
  (TMatrixD*)Mtotal.Write("full_covariance"); 
  (TMatrixD*)Mtotalfraccov.Write("frac_covariance"); 
  (TMatrixD*)Mtotalcorr.Write("full_correlation");
 
  plot_one(Mtotal, h_nue_data, Form("SBNfit_covariance_matrix_%s_full",tag.c_str()));
  plot_one(Mtotalfraccov, h_nue_data, Form("SBNfit_fractional_covariance_matrix_%s_full",tag.c_str()));
  plot_one(Mtotalcorr, h_nue_data, Form("SBNfit_correlation_matrix_%s_full",tag.c_str()));
  
  covarfile->Write();
  covarfile->Close();
   
  //int collsize = input_vec_mc_total.size();	
  int collsize = histos[0]->GetNbinsX();	
  std::cout << "collsize=" << collsize << std::endl;
  //collapse the covariance matrix to only its channels
  TMatrixD Mcol(collsize,collsize);
  TMatrixD Mcolfrac(collsize,collsize);
  TMatrixD Mcolcorr(collsize,collsize);
  Mcol.Zero();
  Mcolfrac.Zero();
  Mcolcorr.Zero();
  
  std::cout << "segfault 8" << std::endl;
  int num_channels = 1;
  std::vector<int> num_subchannels = {int(histos.size())};
  std::vector<int> num_bins = {h_nue_data->GetNbinsX()};
  CollapseSubchannels(Mtotal, Mcol, num_bins, num_channels, num_subchannels);
  
  double toterror=0;
  double totbins=0;
  std::cout << "h" << std::endl;
  for(int i=0; i<collsize; i++){
    for(int j=0; j<collsize; j++){
      / *if(h_nue_data->GetBinContent(i+1) != 0)* / Mcolfrac(i,j) = Mcol(i,j)/(h_nue_data->GetBinContent(i+1)*h_nue_data->GetBinContent(j+1));
      //else Mcolfrac(i,j) = 0.0;
      / *if(h_nue_data->GetBinContent(i+1) != 0)* / Mcolcorr(i,j) = Mcol(i,j)/(sqrt(Mcol(i,i))*sqrt(Mcol(j,j)));
      //else Mcolfrac(i,j) = 0.0;
      if(i==j) std::cout << "bin, uncertainties = " << i+1 << ", " << 100*sqrt(Mcolfrac(i,j)) << std::endl;
      toterror += Mcolfrac(i,j);
      totbins += h_nue_data->GetBinContent(i+1);
    }
  }

  std::cout << "total error = " << sqrt(toterror) << std::endl;
  std::cout << "total bins = " << totbins << std::endl;
  std::cout << "total error/total bins = " << sqrt(toterror)/totbins << std::endl;
  
  //add detsys
  TString fname = Form("../bin/%s_withpoiserr_bgonly.SBNcovar.root",tag.c_str());
  if(detsys) fname = Form("../bin/%s_withpoiserr_detsys.SBNcovar.root",tag.c_str());
  TFile *newmatrix = new TFile(fname.Data(),"recreate");
  newmatrix->cd();

  (TMatrixD*)Mcol.Write("full_covariance"); 
  (TMatrixD*)Mcolfrac.Write("frac_covariance"); 
  (TMatrixD*)Mcolcorr.Write("full_correlation");
 
  plot_one(Mcol, h_nue_data, Form("SBNfit_covariance_matrix_%s_withpoiserr_mc",tag.c_str()));
  plot_one(Mcolfrac, h_nue_data, Form("SBNfit_fractional_covariance_matrix_%s_withpoiserr_mc",tag.c_str()));
  plot_one(Mcolcorr, h_nue_data, Form("SBNfit_correlation_matrix_%s_withpoiserr_mc",tag.c_str()));
  
  newmatrix->Write();
  newmatrix->Close();
 */
}

//borrow this function from SBNchi
void CollapseSubchannels(TMatrixD & M, TMatrixD & Mc, std::vector<int> num_bins, int num_channels, std::vector<int> num_subchannels){
  bool debug = true;
  if(debug)	std::cout<<"Starting:M "<<M.GetNcols()<<" "<<M.GetNrows()<<" "<<std::endl;
  if(debug)	std::cout<<"Starting:Mc "<<Mc.GetNcols()<<" "<<Mc.GetNrows()<<" "<<std::endl;
  
  std::vector<std::vector<TMatrixD>> Summed(num_channels, std::vector<TMatrixD>(num_channels) );	//Initialise a matrix of matricies, to ZERO.
  for(int ic = 0; ic < num_channels; ic++){
    for(int jc =0; jc < num_channels; jc++){
      Summed[ic][jc].ResizeTo(num_bins[jc],num_bins[ic]) ;// This is CORRECT, do not switch (ie Summed[0][1] = size (num_bins[1], num_bins[0])
      Summed[ic][jc] = 0.0;
    }
  }
  
  int mrow = 0.0;
  int mcol = 0.0;
  
  for(int ic = 0; ic < num_channels; ic++){ 	 //Loop over all rows
    for(int jc =0; jc < num_channels; jc++){ //Loop over all columns
      
      if(debug)std::cout<<"Diagonal! : "<<ic<<" "<<jc<<" mcol is: "<<mcol<<" mrow is: "<<mrow<<std::endl;
      
      for(int m=0; m < num_subchannels[ic]; m++){
	for(int n=0; n< num_subchannels[jc]; n++){ //For each big block, loop over all subchannels summing toGether
	  Summed[ic][jc] +=  M.GetSub(mrow+n*num_bins[jc] ,mrow + n*num_bins[jc]+num_bins[jc]-1, mcol + m*num_bins[ic], mcol+ m*num_bins[ic]+num_bins[ic]-1 );
	  //std::cout << "Summed[" << ic << "][" << jc << "] = " << std::endl;
	  //Summed[ic][jc].Print();
	}
      }
      
      
      mrow += num_subchannels[jc]*num_bins[jc];//As we work our way left in columns, add on that many bins
    }//end of column loop
    
    mrow = 0; // as we end this row, reSet row count, but jump down 1 column
    mcol += num_subchannels[ic]*num_bins[ic];
  }//end of row loop
  
  ///********************************* And put them back toGether! ************************//
  Mc.Zero();
  mrow = 0;
  mcol = 0;
  
  //Repeat again for Contracted matrix
  for(int ic = 0; ic < num_channels; ic++){
    for(int jc =0; jc < num_channels; jc++){
      Mc.SetSub(mrow,mcol,Summed[ic][jc]);
      mrow += num_bins[jc];
    }
    
    mrow = 0;
    mcol +=num_bins[ic];
  }
  
  return;
}

//Pretty up the Plots!
void Draw_Stacked(TObjArray histos, TH1D *data, 
		  TString samplename,
		  TPad *pad,
		  bool normalised,
		  std::string stacktitle,
		  std::string var )
{// this function draws histograms stacked and correctly takes into account the
  // stats boxes for each


  for( int b=1; b < data->GetNbinsX()+1; b++ ) std::cout << "bin, bin error = " << data->GetBinContent(b) << ", " << data->GetBinError(b) << std::endl;
  gStyle->SetOptStat(0); 
  if( histos.GetEntries() == 0 ) return; // nothing to do
  
  // Initial set up
  TObjArray statsboxes;    
  std::vector<std::string> legends_str;
  std::vector<int> colours;
  TPaveStats *st1;
  
  const int n = histos.GetEntries();
  
  // choose the colours
  colours.push_back(kRed-3);
  colours.push_back(kPink+1);
  colours.push_back(kMagenta-6);
  colours.push_back(kAzure+4);
  colours.push_back(kSpring-1);
  colours.push_back(kAzure+1);
  
  // lets open and draw the canvas 
  TCanvas *canvas;
  TString cName = Form("%s",samplename.Data());
  if(pad == 0){
    canvas = new TCanvas(cName,"Stacked Histograms", 1200., 800.);
    pad = (TPad*)canvas->cd();
  }
  pad->cd();
  pad->SetGridx(2);        // Vertical grid
  pad->SetTicks(0,0);
  pad->SetRightMargin(0.1);
  pad->SetBottomMargin(0.13);
  
  // lets take the first histoname and see if we can match a title to its which will be HS stack title
  if( stacktitle == "" ) stacktitle = ((TH1F*)histos[0])->GetTitle();
  std::string plottitle = std::string(((TH1F*)histos[0])->GetTitle());
  if( plottitle.find("uBooNE_nue_") != std::string::npos ) stacktitle = "#nu_{e} 1eNp0#pi selection"; 
  else if( plottitle.find("uBooNE_nue0pNp_") != std::string::npos ) stacktitle = "#nu_{e} 1e0p0#pi selection";
  else if( plottitle.find("uBooNE_numu_") != std::string::npos ) stacktitle = "#nu_{#mu} inclusive selection"; 
  else if( plottitle.find("uBooNE_nueall_") != std::string::npos ) stacktitle = "#nu_{e} inclusive selection"; 
  THStack *Hs = new THStack("hs2",stacktitle.c_str());
  
  //Set Axis Units
  //Hs->GetXaxis()->SetTitle("Number of Hits in Slice");
  
  // Set up the LEGEND
  TLegend  *legend = new TLegend(0.49,0.6,0.9,0.9); // we need different positions for the legend to not 
  legend->SetHeader("MicroBooNE Simulation, Preliminary","L");
  legend->SetNColumns(2);
  // get the plot titles for the legend
  for(int i = 0; i < n ; i++){
    
    TH1D *h =  (TH1D*)histos[i];
    std::cout << "histo = " << h->GetTitle() << std::endl;
    //h->Scale(1,"width");
    h->SetLineWidth(0);
    h->SetLineColor(colours[i]);
    h->SetFillColor(colours[i]);
    std::string histosample = std::string(h->GetTitle());
    if( histosample.find("extbnb") != std::string::npos ){ 
      h->SetFillStyle(3354);
      h->SetFillColor(kBlack);
      h->SetLineColor(kBlack);
      std::string htitle = Form("BNB Off: %2.2f",h->Integral());
      std::cout << "htitle = " << htitle << std::endl;
      h->SetTitle(htitle.c_str());
    }
    else if( histosample.find("lee") != std::string::npos ){
      h->SetFillColor(kOrange+1);
      h->SetLineColor(kOrange+1);
      std::string htitle = Form("#nu_{e} LEE: %2.2f",h->Integral());
      //std::cout << "htitle = " << htitle << std::endl;
      h->SetTitle(htitle.c_str());
    }
    else if( histosample.find("ncpi0") != std::string::npos ){
      h->SetFillColor(kAzure+1);
      h->SetLineColor(kAzure+1);
      std::string htitle = Form("#nu NC #pi^{0}: %2.2f",h->Integral());
      //std::cout << "htitle = " << htitle << std::endl;
      h->SetTitle(htitle.c_str());
    }
    else if( histosample.find("ccpi0") != std::string::npos ){
      h->SetFillColor(kAzure+4);
      h->SetLineColor(kAzure+4);
      std::string htitle = Form("#nu CC #pi^{0}: %2.2f",h->Integral());
      //std::cout << "htitle = " << htitle << std::endl;
      h->SetTitle(htitle.c_str());
    }
    else if( histosample.find("intrinsic") != std::string::npos ){
      h->SetFillColor(kSpring);
      h->SetLineColor(kSpring);
      std::string htitle = Form("#nu_{e} CC: %2.2f",h->Integral());
      //std::cout << "htitle = " << htitle << std::endl;
      h->SetTitle(htitle.c_str());
    }
    else if( histosample.find("ncnopi") != std::string::npos ){
      h->SetFillColor(kAzure+7);
      h->SetLineColor(kAzure+7);
      std::string htitle = Form("NC 0#pi: %2.2f",h->Integral());
      //std::cout << "htitle = " << htitle << std::endl;
      h->SetTitle(htitle.c_str());
    }
    else if( histosample.find("CCmuNoPi") != std::string::npos ){
      h->SetFillColor(kAzure-2);
      h->SetLineColor(kAzure-2);
      std::string htitle = Form("#nu_{#mu} CC 0#pi: %2.2f",h->Integral());
      //std::cout << "htitle = " << htitle << std::endl;
      h->SetTitle(htitle.c_str());
    }
    else if( histosample.find("NCcPiNoPi0") != std::string::npos ){
      h->SetFillColor(kAzure-5);
      h->SetLineColor(kAzure-5);
      std::string htitle = Form("NC #pi 0#pi^{0}: %2.2f",h->Integral());
      //std::cout << "htitle = " << htitle << std::endl;
      h->SetTitle(htitle.c_str());
    }
    else if( histosample.find("CCmuCPiNoPi0") != std::string::npos ){
      h->SetFillColor(kAzure-8);
      h->SetLineColor(kAzure-8);
      std::string htitle = Form("#nu_{#mu} CC #pi 0#pi^{0}: %2.2f",h->Integral());
      //std::cout << "htitle = " << htitle << std::endl;
      h->SetTitle(htitle.c_str());
    }
    else if( histosample.find("dirt") != std::string::npos ){
      h->SetFillColor(kOrange+3);
      h->SetLineColor(kOrange+3);
      std::string htitle = Form("Dirt: %2.2f",h->Integral());
      //std::cout << "htitle = " << htitle << std::endl;
      h->SetTitle(htitle.c_str());
    }
    else if( histosample.find("nue_nu") != std::string::npos || histosample.find("nue0p_nu") != std::string::npos || histosample.find("nueall_nu") != std::string::npos || histosample.find("numu_bnb") != std::string::npos ){
      h->SetFillColor(kRed+2);
      h->SetLineColor(kRed+2);
      std::string htitle = Form("#nu_{#mu} BNB: %2.2f",h->Integral());
      //std::cout << "htitle = " << htitle << std::endl;
      h->SetTitle(htitle.c_str());
    }
    else if( histosample.find("dirt") != std::string::npos ){
      h->SetFillColor(kOrange+3);
      h->SetLineColor(kOrange+3);
      std::string htitle = Form("Dirt: %2.2f",h->Integral());
      std::cout << "htitle = " << htitle << std::endl;
      h->SetTitle(htitle.c_str());
    }
    else if( histosample.find("nue_nu") != std::string::npos || histosample.find("nue0p_nu") != std::string::npos || histosample.find("nueall_nu") != std::string::npos || histosample.find("numu_nu") != std::string::npos || histosample.find("1e0p_nu") != std::string::npos ){
      h->SetFillColor(kCyan);
      h->SetLineColor(kCyan);
      std::string htitle = Form("#nu_{#mu} CC: %2.2f",h->Integral());
      std::cout << "htitle = " << htitle << std::endl;
      h->SetTitle(htitle.c_str());
    }
    gStyle->SetLegendTextSize(0.035);
    legends_str.push_back( h->GetTitle() );
    legend->AddEntry(h,legends_str[i].c_str(),"f"); 
    std::cout << "histo name =  " << h->GetTitle() << std::endl;
    std::cout << "integral =  " << h->Integral() << std::endl;
    Hs->Add(h,"sames");
  }
    legend->AddEntry(data,"MC Stats. Error","lep"); 
  
  float heightboxes;
  // the position of the top corner
  float top_corner,deltay;
  // don't draw errors if it is 1
  int no_error = 1;
  top_corner = 0.9;

  /* 
  // if the array has more than 0 histograms lets specify the stat boxes 
  if (histos.GetEntries() > 0) { 
  heightboxes = (float)0.5/(float)histos.GetEntries();
  }
  else
  heightboxes = (float)0.5;
  */
  
  //HACK -need to draw a smaller histogram first to set the x-axis range correctly 
  //for stacks with variable-width bins
  //TH1D* tmp_mnv = (TH1D*)Hs->GetHists()->At(0);
  Hs->SetMaximum(1.4*data->GetMaximum());
  
  data->SetMarkerStyle(1);
  data->SetMarkerSize(2);
  data->SetMarkerColor(kBlack);
  data->SetLineWidth(2);
  data->SetLineColor(kBlack);
  
  data->Draw("PE1"); 
  // DRAW not stacked to get correct stats boxes
  if (no_error == 1) // do not draw errors
    Hs->Draw("HIST");
  else
    Hs->Draw("HISTE");

  Hs->GetXaxis()->SetTitle( var.c_str() );
  Hs->GetYaxis()->SetTitle( "Events" );
  //Hs->GetYaxis()->SetNdivisions(220, kTRUE);
  //Hs->GetXaxis()->SetNdivisions(220, kTRUE);
  Hs->GetXaxis()->SetTitleFont(43);
  Hs->GetYaxis()->SetTitleFont(43);
  Hs->GetXaxis()->SetTitleSize(40);
  Hs->GetYaxis()->SetTitleSize(40);
  Hs->GetXaxis()->SetLabelFont(43);
  Hs->GetYaxis()->SetLabelFont(43);
  Hs->GetXaxis()->SetLabelSize(25);
  Hs->GetYaxis()->SetLabelSize(25);
  Hs->GetXaxis()->CenterTitle(kTRUE);

  data->DrawCopy("PE1 same");
  legend->Draw("");
  AddPlotLabel("MicroBooNE Preliminary, Run 1 scaled to 1.01e+21 POT",0.4,0.8);
  
  gPad->Modified(); // so it updates the pad with the new changes
 
  canvas->Modified(); 
  canvas->Print(Form("%s.pdf",cName.Data()),"pdf"); 
  
  delete canvas;
  
  return;
}

TH1D* GetTotalError(std::string name, TMatrixD errMatrix, TH1D *histo ){
  // Make a copy of this histogram as a TH1D and rename it
  TString hname = name+"_total_error";
  TH1D *err = (TH1D*)histo->Clone(hname);
  err->Reset();
  
  const int highBin = err->GetNbinsX() + 1;
  const int lowBin = 0;
  
  for( int iBin = lowBin; iBin <= highBin; ++iBin )
    {
      double derr = errMatrix[iBin][iBin];
      err->SetBinContent( iBin+1, ( derr > 0 ) ? sqrt(derr): 0. );
    }
  
  return err;
}

void DrawMCAndSyst(TString cName, TH1D* MC, std::string var, std::string title){

  TCanvas can(cName,cName,1200,800);

  TH1D* tmpMC = (TH1D*)MC->Clone("tempCV");
  TH1D* sysErr = (TH1D*)MC->Clone("tempSysCV");
  TH1D* tmpsys = (TH1D*)MC->Clone("tempSysRatio");

  for( int bin = 1; bin < tmpsys->GetNbinsX()+1; bin++){
    tmpsys->SetBinContent(bin,1.0);
    tmpsys->SetBinError(bin,tmpMC->GetBinError(bin)/tmpMC->GetBinContent(bin));
    sysErr->SetBinContent(bin,tmpMC->GetBinContent(bin));
    sysErr->SetBinError(bin,tmpMC->GetBinError(bin));
  }
  
  // Upper plot will be in pad1
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 1.0);
  pad1->SetBottomMargin(0.05); // Upper and lower plot are joined
  pad1->SetGridx(2);        // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad

  tmpMC->SetLineColor(kRed);
  tmpMC->SetLineWidth(3);
  tmpMC->SetLineStyle(1);
  sysErr->SetTitleSize(80);
  sysErr->SetTitleFont(43);
  sysErr->SetTitle(title.c_str());
  sysErr->SetLineColor(kRed-10);
  sysErr->SetLineWidth(1);
  sysErr->SetFillColor(kRed-10);
  sysErr->SetFillStyle(1001);
  sysErr->SetStats(0);
  sysErr->GetYaxis()->SetTitleSize(25);
  sysErr->GetYaxis()->SetTitleFont(43);
  sysErr->GetYaxis()->SetTitle("Events");
  sysErr->GetXaxis()->SetLabelSize(0);
  sysErr->GetXaxis()->SetTitleFont(43);
  sysErr->GetXaxis()->SetTitle(var.c_str());
  sysErr->SetMaximum(1.2*sysErr->GetMaximum());
  sysErr->Draw("E2");
  tmpMC->Draw("HIST same");

  // lower plot will be in pad
  can.cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.35);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.25);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad

  tmpsys->SetTitle("");
  tmpsys->GetYaxis()->SetTitle("uncertainties");
  tmpsys->GetXaxis()->SetTitle(var.c_str());
  tmpsys->SetStats(0);
  tmpsys->SetLineColor(kRed);
  tmpsys->SetLineWidth(1);
  tmpsys->SetFillColor(kRed-10);
  tmpsys->SetFillStyle(1001);
  tmpsys->GetYaxis()->SetTitleSize(25);
  tmpsys->GetYaxis()->SetTitleFont(43);
  tmpsys->GetXaxis()->SetTitleSize(25);
  tmpsys->GetXaxis()->SetTitleOffset(3.5);
  tmpsys->GetYaxis()->SetTitleOffset(1.0);
  tmpsys->GetXaxis()->SetLabelSize(0.1);
  tmpsys->GetYaxis()->SetLabelSize(0.08);
  tmpsys->GetYaxis()->SetNdivisions(210, kTRUE);
  tmpsys->GetXaxis()->SetTitleFont(43);
  tmpsys->SetMaximum(1.4);
  tmpsys->SetMinimum(0.6);
  tmpsys->Draw("E2");

  const TAxis *axis = tmpsys->GetXaxis();
  double lowX  = axis->GetBinLowEdge( axis->GetFirst() );
  double highX = axis->GetBinUpEdge(  axis->GetLast() );

  TLine line;
  line.SetLineStyle(2);
  line.SetLineWidth(2);
  line.SetLineColor(kRed);
  line.DrawLine(lowX, 1., highX, 1.); //creates a new line which is owned by gPad
 
  
  pad2->Modified(); // so it updates the pad with the new changes
  pad2->Draw("");

  gPad->RedrawAxis();
  gPad->Update();
 
  can.Print(cName+".pdf","pdf");

}

void DrawMCAndSystOverlay(TString cName, TH1D* MC1, TH1D* MC2, std::string var, std::string title){

  TCanvas can(cName,cName,1200,800);

  TH1D* tmpMC = (TH1D*)MC1->Clone("tempCV");
  TH1D* tmpMC2 = (TH1D*)MC2->Clone("tempCV2");
  TH1D* sysErr = (TH1D*)MC1->Clone("tempSysCV");
  TH1D* sysErr2 = (TH1D*)MC2->Clone("tempSysCV2");
  TH1D* tmpsys = (TH1D*)MC1->Clone("tempSysRatio");
  TH1D* tmpsys2 = (TH1D*)MC2->Clone("tempSysRatio2");

  for( int bin = 1; bin < tmpsys->GetNbinsX()+1; bin++){
    tmpsys->SetBinContent(bin,1.0);
    tmpsys->SetBinError(bin,tmpMC->GetBinError(bin)/tmpMC->GetBinContent(bin));
    tmpsys2->SetBinContent(bin,1.0);
    tmpsys2->SetBinError(bin,tmpMC2->GetBinError(bin)/tmpMC2->GetBinContent(bin));
    sysErr->SetBinContent(bin,tmpMC->GetBinContent(bin));
    sysErr->SetBinError(bin,tmpMC->GetBinError(bin));
    sysErr2->SetBinContent(bin,tmpMC2->GetBinContent(bin));
    sysErr2->SetBinError(bin,tmpMC2->GetBinError(bin));
  }
  
  // Upper plot will be in pad1
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 1.0);
  pad1->SetBottomMargin(0.05); // Upper and lower plot are joined
  pad1->SetGridx(2);        // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad

  tmpMC->SetLineColor(kRed);
  tmpMC->SetLineWidth(3);
  tmpMC->SetLineStyle(1);
  tmpMC2->SetLineColor(kMagenta);
  tmpMC2->SetLineWidth(3);
  tmpMC2->SetLineStyle(2);
  sysErr->SetTitleSize(80);
  sysErr->SetTitleFont(43);
  sysErr->SetTitle(title.c_str());
  sysErr->SetLineColor(kRed-10);
  sysErr->SetLineWidth(1);
  sysErr2->SetLineColor(kMagenta-10);
  sysErr2->SetLineWidth(1);
  sysErr2->SetFillColor(kMagenta-10);
  sysErr2->SetFillStyle(3144);
  sysErr->SetStats(0);
  sysErr->GetYaxis()->SetTitleSize(25);
  sysErr->GetYaxis()->SetTitleFont(43);
  sysErr->GetYaxis()->SetTitle("Events");
  sysErr->GetXaxis()->SetLabelSize(0);
  sysErr->GetXaxis()->SetTitleFont(43);
  sysErr->GetXaxis()->SetTitle(var.c_str());
  sysErr->SetMaximum(1.2*sysErr->GetMaximum());
  sysErr->Draw("E2");
  sysErr2->Draw("E2 same");
  tmpMC->Draw("HIST same");
  tmpMC2->Draw("HIST same");

  
  // lower plot will be in pad
  can.cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.35);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.25);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad

  tmpsys->SetTitle("");
  tmpsys->GetYaxis()->SetTitle("uncertainties");
  tmpsys->GetXaxis()->SetTitle(var.c_str());
  tmpsys->SetStats(0);
  tmpsys->SetLineColor(kRed);
  tmpsys->SetLineWidth(3);
  //tmpsys->SetFillColor(kRed-10);
  //tmpsys->SetFillStyle(1001);
  tmpsys2->SetLineColor(kBlue);
  tmpsys2->SetLineWidth(2);
  tmpsys2->SetMarkerStyle(20);
  tmpsys->SetMarkerStyle(21);
  tmpsys2->SetMarkerColor(kBlue);
  tmpsys->SetMarkerColor(kRed);
  //tmpsys2->SetFillColor(kMagenta-10);
  //tmpsys2->SetFillStyle(3144);
  tmpsys->GetYaxis()->SetTitleSize(25);
  tmpsys->GetYaxis()->SetTitleFont(43);
  tmpsys->GetXaxis()->SetTitleSize(25);
  tmpsys->GetXaxis()->SetTitleOffset(3.5);
  tmpsys->GetYaxis()->SetTitleOffset(1.0);
  tmpsys->GetXaxis()->SetLabelSize(0.1);
  tmpsys->GetYaxis()->SetLabelSize(0.08);
  tmpsys->GetYaxis()->SetNdivisions(210, kTRUE);
  tmpsys->GetXaxis()->SetTitleFont(43);
  tmpsys->SetMaximum(1.4);
  tmpsys->SetMinimum(0.6);
  tmpsys->Draw("PE1");
  tmpsys2->Draw("PE1 same");

  const TAxis *axis = tmpsys->GetXaxis();
  double lowX  = axis->GetBinLowEdge( axis->GetFirst() );
  double highX = axis->GetBinUpEdge(  axis->GetLast() );

  TLine line;
  line.SetLineStyle(2);
  line.SetLineWidth(2);
  line.SetLineColor(kRed);
  line.DrawLine(lowX, 1., highX, 1.); //creates a new line which is owned by gPad 
  
  pad2->Modified(); // so it updates the pad with the new changes
  //pad2->Draw("");

  gPad->RedrawAxis();
  gPad->Update();
 
  can.Print(cName+".pdf","pdf");

}
void DrawDataMCAndSyst(TString cName, TH1D* MC, TH1D* data, std::string var, std::string title){

  TCanvas can(cName,cName,1200,800);

  TH1D* tmpMC = (TH1D*)MC->Clone("tempCV");
  TH1D* sysErr = (TH1D*)MC->Clone("tempSysCV");
  TH1D* tmpData = (TH1D*)data->Clone("tempData");

  TH1D* tmpratio = (TH1D*)data->Clone("tempRatio");
  tmpratio->Reset();
  tmpratio->Divide(tmpData,tmpMC,1.0,1.0,"B");
  std::cout << "nbins data, mc = " << tmpData->GetNbinsX() << ", " << tmpMC->GetNbinsX() << std::endl;
  TH1D* tmpsys = (TH1D*)MC->Clone("tempSysRatio");

  for( int bin = 1; bin < tmpsys->GetNbinsX()+1; bin++){
    tmpsys->SetBinContent(bin,1.0);
    tmpsys->SetBinError(bin,tmpMC->GetBinError(bin)/tmpMC->GetBinContent(bin));
    sysErr->SetBinContent(bin,tmpMC->GetBinContent(bin));
    sysErr->SetBinError(bin,tmpMC->GetBinError(bin));
    //std::cout << "temp ratio = " << tmpratio->GetBinContent(bin) << std::endl; 
  }
  
  // Upper plot will be in pad1
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 1.0);
  pad1->SetBottomMargin(0.05); // Upper and lower plot are joined
  pad1->SetGridx(2);        // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad

  tmpMC->SetLineColor(kRed);
  tmpMC->SetLineWidth(3);
  tmpMC->SetLineStyle(1);
  sysErr->SetTitleSize(80);
  sysErr->SetTitleFont(43);
  sysErr->SetTitle(title.c_str());
  tmpData->SetLineWidth(2);
  tmpData->SetLineColor(kBlack);
  sysErr->SetLineColor(kRed-10);
  sysErr->SetLineWidth(1);
  sysErr->SetFillColor(kRed-10);
  sysErr->SetFillStyle(1001);
  sysErr->SetStats(0);
  sysErr->GetYaxis()->SetTitleSize(25);
  sysErr->GetYaxis()->SetTitleFont(43);
  sysErr->GetYaxis()->SetTitle("Events");
  sysErr->GetXaxis()->SetLabelSize(0);
  sysErr->SetMaximum(1.2*tmpData->GetMaximum());
  sysErr->Draw("E2");
  tmpMC->Draw("hist same"); 
  tmpData->Draw("PE same"); 
  
  // lower plot will be in pad
  can.cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.35);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.25);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad

  tmpsys->SetTitle("");
  tmpsys->GetYaxis()->SetTitle("uncertainties");
  tmpsys->GetXaxis()->SetTitle(var.c_str());
  tmpsys->SetStats(0);
  tmpsys->SetLineColor(kRed-10);
  tmpsys->SetLineWidth(1);
  tmpsys->SetFillColor(kRed-10);
  tmpsys->SetFillStyle(1001);
  tmpsys->GetYaxis()->SetTitleSize(25);
  tmpsys->GetYaxis()->SetTitleFont(43);
  tmpsys->GetXaxis()->SetTitleSize(25);
  tmpsys->GetXaxis()->SetTitleOffset(3.5);
  tmpsys->GetYaxis()->SetTitleOffset(1.0);
  tmpsys->GetXaxis()->SetLabelSize(0.1);
  tmpsys->GetYaxis()->SetLabelSize(0.08);
  tmpsys->GetYaxis()->SetNdivisions(210, kTRUE);
  tmpsys->GetXaxis()->SetTitleFont(43);
  tmpsys->SetMaximum(2.0);
  tmpsys->SetMinimum(0.0);
  tmpsys->Draw("E2");

  const TAxis *axis = tmpratio->GetXaxis();
  double lowX  = axis->GetBinLowEdge( axis->GetFirst() );
  double highX = axis->GetBinUpEdge(  axis->GetLast() );

  TLine line;
  line.SetLineStyle(2);
  line.SetLineWidth(2);
  line.SetLineColor(kRed);
  line.DrawLine(lowX, 1., highX, 1.); //creates a new line which is owned by gPad
 
  tmpratio->SetLineColor(kBlack);
  tmpratio->SetLineWidth(2);
  tmpratio->Draw("PE same");

  pad2->Update();
  pad2->Modified(); // so it updates the pad with the new changes
  pad2->Draw("");

  gPad->RedrawAxis();
  gPad->Update();
 
  can.Print(cName+".pdf","pdf");

}

void plot_one(TMatrixD matrix, TH1D *h_nue, std::string tag){

    std::string dir="/uboone/data/users/wospakrk/SBNFitPlots/";
    std::vector<std::string> channel_names = {"nu uBooNE nue intrinsic","nu uBooNE numu inclusive"};
    int num_channels = 1;
    int num_bins_total = matrix.GetNrows();
    std::vector<int> num_bins = {h_nue->GetNbinsX()}; 
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

        TText * tmd = new TText(-num_bins_total*percent_left*0.15, use_full+nice_shift*0.5, (channel_names[ic]).c_str() );
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

    c_full->Print((dir+tag+"_withpoiserr_collapsed.pdf").c_str(),"pdf");
}

double CalcChi(TH1D *h_data, TH1D *h_mc, TMatrixD inv_matrix, TH1D *h_nue){

  double tchi=0.0;

  for( int i = 0; i < inv_matrix.GetNcols(); i++ ){
    for( int j = 0; j < inv_matrix.GetNrows(); j++ ){
      if( i > (h_nue->GetNbinsX()-1) && j > (h_nue->GetNbinsX()-1) ){
        int numubin_i = i-h_nue->GetNbinsX()+1;
        int numubin_j = j-h_nue->GetNbinsX()+1;
        std::cout << "i, j, numubin_i, numubin_j, h_data_i, h_mc_i, cov_matrix, h_data_j,  = " << i << ", " << j << ", " << numubin_i << ", " << numubin_j << ", " << h_data->GetBinContent(numubin_i+1) << ", " << h_mc->GetBinContent(numubin_i+1) << std::endl;
        double tchiperbin = (h_data->GetBinContent(numubin_i+1)-h_mc->GetBinContent(numubin_i+1))*inv_matrix(i,j)*(h_data->GetBinContent(numubin_j+1)-h_mc->GetBinContent(numubin_j+1));
        tchi += tchiperbin;
      }
    }
  }
  std::cout << "tchi = " << tchi << std::endl;
  return tchi;

}

void AddPlotLabel( const char * 	label,
                   const double 	x,
                   const double 	y,
                   const double 	size,
                   const int 	        color,
                   const int 	        font,
                   const int 	        align,
                   const double 	angle ) 
{
  TLatex *latex = new TLatex( x, y, label );

  latex->SetNDC();
  latex->SetTextSize(size);
  latex->SetTextColor(color);
  latex->SetTextFont(font);
  latex->SetTextAlign(align);
  latex->SetTextAngle(angle);
  latex->Draw("same");
}
