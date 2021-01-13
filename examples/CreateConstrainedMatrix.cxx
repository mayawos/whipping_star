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


void DrawMCAndSystOverlay(TString cName, TH1D* MC1, TH1D* MC2, std::string var, std::string title);
void DrawDataMCAndSyst(TString cName, TH1D* MC1, TH1D* MC2, TH1D* data, std::string var, std::string title);
void plot_one(TMatrixD matrix, std::vector<double> spectrum, std::vector<std::string> channel_names, std::vector<int> num_bins, std::string tag, bool collapsed);
bool verbose=true;

int main(int argc, char* argv[])
{
  std::string xml = "example.xml";
  std::string tag = "example";
  std::string fakedata = "example";
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
  
  const struct option longopts[] =
    {
      {"xml", 		required_argument, 	0, 'x'},
      {"tag", 		required_argument, 	0, 't'},
      {"fakedata", 	required_argument, 	0, 'f'},
      {"var", 		required_argument, 	0, 'v'},
      {"sys",	        no_argument, 		0, 's'},
      {"detsys",	no_argument, 		0, 'd'},
      {"combined",	no_argument, 		0, 'c'},
      {"np",	        no_argument, 		0, 'n'},
      {"zp",	        no_argument, 		0, 'z'},
      {"mcerr",	        no_argument, 		0, 'm'},
      {"addext",	no_argument, 		0, 'a'},
      {0,		no_argument, 		0,  0},
    };
  
  while(iarg != -1)
    {
      iarg = getopt_long(argc,argv, "x:t:f:v:sdcnzma", longopts, &index);
      
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
	case 'v':
	  var = optarg;
	  break;
	case 's':
	  sys = true;
	  break;
	case 'd':
	  detsys = true;
	  break;
	case 'c':
	  combined = true;
	  break;
	case 'n':
	  np = true;
	  break;
	case 'z':
	  zp = true;
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
  
  std::string detsyslabel = "";
  if(detsys) detsyslabel = "_with_detsys"; 
  gStyle->SetPalette(kLightTemperature);
  
  //set up text files to write the table
  std::ofstream outfile1, outfile2, outfile3, outfile4, outfile5, outfile6, outfile7, outfile8, outfile9;
  outfile1.open(Form("full_1eNp_lee_systematicstable_beforeconstraint%s.txt",detsyslabel.c_str())); 
  outfile2.open(Form("full_1eNp_lee_systematicstable_afterconstraint%s.txt",detsyslabel.c_str())); 
  outfile3.open(Form("full_1eNp_systematicstable_beforeconstraint%s.txt",detsyslabel.c_str())); 
  outfile4.open(Form("full_1eNp_systematicstable_afterconstraint%s.txt",detsyslabel.c_str())); 
  outfile5.open(Form("full_1e0p_lee_systematicstable_beforeconstraint%s.txt",detsyslabel.c_str())); 
  outfile6.open(Form("full_1e0p_lee_systematicstable_afterconstraint%s.txt",detsyslabel.c_str())); 
  outfile7.open(Form("full_1e0p_systematicstable_beforeconstraint%s.txt",detsyslabel.c_str())); 
  outfile8.open(Form("full_1e0p_systematicstable_afterconstraint%s.txt",detsyslabel.c_str())); 
  outfile9.open(Form("full_numu_systematicstable%s.txt",detsyslabel.c_str())); 

  //Load up the central value spectra we computed in example 1, to act as a signal spectra
  SBNspec sig_spectra("../bin/"+tag+".SBNspec.root",xml);
  
  //Repeat but this time scale the LEE signal subchannel to 0, to act as our background spectra
  SBNspec bkg_spectra("../bin/"+tag+".SBNspec.root",xml);
  bkg_spectra.Scale("lee", 0.0); // Scale() will scale all that match so bkg_spectra.Scale("leesignal",0.0); would work too 
  
  TMatrixD Mdetsys(sig_spectra.num_bins_total,sig_spectra.num_bins_total);
  TMatrixD Mfracdetsys(sig_spectra.num_bins_total,sig_spectra.num_bins_total);
  TMatrixD Mcorrdetsys(sig_spectra.num_bins_total,sig_spectra.num_bins_total);
  TMatrixD Mdetsys0(sig_spectra.num_bins_total,sig_spectra.num_bins_total);
  Mdetsys.Zero();
  Mdetsys0.Zero();
  SBNchi chi_h1(sig_spectra);	
  SBNchi chi_h0(bkg_spectra);	
  if(detsys){ 
    chi_h1.FillDetSysMatrix(Mdetsys,sig_spectra);
    chi_h0.FillDetSysMatrix(Mdetsys0,bkg_spectra);
  }
 
  std::vector<double> sig_fullvec, bkg_fullvec, sig_collvec, bkg_collvec, sig_errfullvec, bkg_errfullvec, sig_exterrfullvec, bkg_exterrfullvec, numu_data_vec, numu_vec, delta_numu; 
  sig_spectra.CalcFullVector();
  sig_fullvec = sig_spectra.full_vector;
  bkg_fullvec = bkg_spectra.full_vector;
  //frac detsys
  for(int i=0; i < Mdetsys.GetNrows(); i++){
  for(int j=0; j < Mdetsys.GetNcols(); j++){
    Mfracdetsys(i,j) = 0.0;
    if(sig_spectra.full_vector.at(i) !=0 && sig_spectra.full_vector.at(j) !=0 ) Mfracdetsys(i,j) = Mdetsys(i,j)/(sig_spectra.full_vector.at(i)*sig_spectra.full_vector.at(j));
    Mcorrdetsys(i,j) = 0.0;
    if(Mdetsys(i,j) != 0 ) Mcorrdetsys(i,j) = Mdetsys(i,j)/sqrt(Mdetsys(i,i)*Mdetsys(j,j));
  }
  }
 
  //collapse vector
  sig_spectra.CollapseVector();
  bkg_spectra.CollapseVector();
  sig_collvec = sig_spectra.collapsed_vector;
  bkg_collvec = bkg_spectra.collapsed_vector;

  //mc intrinsic
  sig_spectra.CalcErrorVector();
  sig_errfullvec = sig_spectra.full_error;
  sig_exterrfullvec = sig_spectra.full_ext_err_vector;
  bkg_errfullvec = bkg_spectra.full_error;
  bkg_exterrfullvec = bkg_spectra.full_ext_err_vector;

  //collapse matrix detsys
  TMatrixD detsyscoll;
  TMatrixD fracdetsyscoll;
  TMatrixD corrdetsyscoll;
  TMatrixD fraccoll;
  TMatrixD fraccoll_noext;
  TMatrixD corrcoll;
  TMatrixD detsyscoll0;
  TMatrixD fraccov;
  TMatrixD corr;
  detsyscoll.ResizeTo(sig_spectra.num_bins_total_compressed, sig_spectra.num_bins_total_compressed);
  fracdetsyscoll.ResizeTo(sig_spectra.num_bins_total_compressed, sig_spectra.num_bins_total_compressed);
  corrdetsyscoll.ResizeTo(sig_spectra.num_bins_total_compressed, sig_spectra.num_bins_total_compressed);
  fraccoll.ResizeTo(sig_spectra.num_bins_total_compressed, sig_spectra.num_bins_total_compressed);
  fraccoll_noext.ResizeTo(sig_spectra.num_bins_total_compressed, sig_spectra.num_bins_total_compressed);
  corrcoll.ResizeTo(sig_spectra.num_bins_total_compressed, sig_spectra.num_bins_total_compressed);
  detsyscoll0.ResizeTo(bkg_spectra.num_bins_total_compressed, bkg_spectra.num_bins_total_compressed);
  fraccov.ResizeTo(sig_spectra.num_bins_total_compressed, sig_spectra.num_bins_total_compressed);
  corr.ResizeTo(sig_spectra.num_bins_total_compressed, sig_spectra.num_bins_total_compressed);
  chi_h1.CollapseModes(Mdetsys, detsyscoll);
  chi_h0.CollapseModes(Mdetsys0, detsyscoll0);

  //channel names
  std::vector<std::string> channel_names;
  std::vector<std::string> channel_names_coll;
  std::vector<int> num_bins;
  std::vector<int> num_bins_coll;

  //get extbnb
  std::vector<double> sig_collvec_noext;
  std::vector<double> collapsed_ext;

  TH1D *h_1eNp_bg_ext, *h_1e0p_bg_ext, *h_numu_ext;
  for(auto& h: sig_spectra.hist){
      channel_names.push_back(h.GetName());
      num_bins.push_back(h.GetNbinsX());
      TString hname = h.GetName();
      if(np || combined ){if(hname=="nu_uBooNE_1eNp_bg_ext") h_1eNp_bg_ext = (TH1D*)h.Clone("nu_uBooNE_1eNp_bg_ext");}
      if(zp || combined ){if(hname=="nu_uBooNE_1e0p_bg_ext") h_1e0p_bg_ext = (TH1D*)h.Clone("nu_uBooNE_1eNp_bg_ext");}
      if(hname=="nu_uBooNE_numu_ext") h_numu_ext = (TH1D*)h.Clone("nu_uBooNE_1eNp_bg_ext");
  }
  if(np || combined ){ for(int bin=1; bin < h_1eNp_bg_ext->GetNbinsX()+1; bin++) collapsed_ext.push_back(0.0);}
  if(np || combined ){ for(int bin=1; bin < h_1eNp_bg_ext->GetNbinsX()+1; bin++) collapsed_ext.push_back(h_1eNp_bg_ext->GetBinContent(bin));}
  if(zp || combined ){ for(int bin=1; bin < h_1e0p_bg_ext->GetNbinsX()+1; bin++) collapsed_ext.push_back(0.0);}
  if(zp || combined ){ for(int bin=1; bin < h_1e0p_bg_ext->GetNbinsX()+1; bin++) collapsed_ext.push_back(h_1e0p_bg_ext->GetBinContent(bin));}

  for(int bin=1; bin < h_numu_ext->GetNbinsX()+1; bin++) collapsed_ext.push_back(h_numu_ext->GetBinContent(bin));

  //define bins number 
  int numu_bins = h_numu_ext->GetNbinsX();
  int nue_bins = sig_spectra.num_bins_total_compressed-numu_bins;
  int np_bins = 0, zp_bins = 0;
  if(np || combined ) np_bins = h_1eNp_bg_ext->GetNbinsX();
  if(zp || combined ) zp_bins = h_1e0p_bg_ext->GetNbinsX();
  
  for(int i=0; i < detsyscoll.GetNrows(); i++){
     sig_collvec_noext.push_back( sig_collvec[i] - collapsed_ext[i] );
  }
  
  //frac detsys
  for(int i=0; i < detsyscoll.GetNrows(); i++){
  for(int j=0; j < detsyscoll.GetNcols(); j++){
    fracdetsyscoll(i,j) = 0.0;
    if(sig_collvec[i] !=0 && sig_collvec[j] !=0 ) fracdetsyscoll(i,j) = detsyscoll(i,j)/(sig_collvec[i]*sig_collvec[j]);
    corrdetsyscoll(i,j) = 0.0;
    if(detsyscoll(i,j) != 0 ) corrdetsyscoll(i,j) = detsyscoll(i,j)/sqrt(detsyscoll(i,i)*detsyscoll(j,j));
  }
  }

  for(auto& h: sig_spectra.hist){
      channel_names.push_back(h.GetName());
      num_bins.push_back(h.GetNbinsX());
  }
  if(combined){
    channel_names_coll.push_back("1eNp LEE");
    channel_names_coll.push_back("1eNp nue + bkg");
    channel_names_coll.push_back("1e0p LEE");
    channel_names_coll.push_back("1e0p nue + bkg");
    channel_names_coll.push_back("numu");
  }
  else if(np){
    channel_names_coll.push_back("1eNp LEE");
    channel_names_coll.push_back("1eNp nue + bkg");
    channel_names_coll.push_back("numu");
  }
  else if(zp){
    channel_names_coll.push_back("1e0p LEE");
    channel_names_coll.push_back("1e0p nue + bkg");
    channel_names_coll.push_back("numu");
  }
  for(int c=0; c < channel_names_coll.size(); c++){ 
    if( channel_names_coll[c].find("1eNp") != std::string::npos ) num_bins_coll.push_back(np_bins);
    if( channel_names_coll[c].find("1e0p") != std::string::npos ) num_bins_coll.push_back(zp_bins);
    if( channel_names_coll[c].find("numu") != std::string::npos ) num_bins_coll.push_back(numu_bins);
  }

  plot_one(Mcorrdetsys, sig_spectra.full_vector, channel_names, num_bins, tag+"_correlation_matrix", false);
  plot_one(Mfracdetsys, sig_spectra.full_vector, channel_names, num_bins, tag+"_fractional_covariance_matrix", false);
  plot_one(corrdetsyscoll, sig_collvec, channel_names_coll, num_bins_coll, tag+"_correlation_matrix", true);
  plot_one(fracdetsyscoll, sig_collvec, channel_names_coll, num_bins_coll, tag+"_fractional_covariance_matrix", true);
  
  //correlation detsys
  TFile * fdetsys = new TFile(Form("DetSys_%s.SBNcovar.root",tag.c_str()),"recreate");
  fdetsys->cd();
  (TMatrixD*)Mdetsys.Write("full_covariance");
  (TMatrixD*)Mfracdetsys.Write("frac_covariance");
  (TMatrixD*)Mcorrdetsys.Write("full_correlation");
  (TMatrixD*)detsyscoll.Write("collapsed_covariance");
  (TMatrixD*)fracdetsyscoll.Write("collapsed_frac_covariance");
  (TMatrixD*)corrdetsyscoll.Write("collapsed_correlation");
  fdetsys->Close(); 
  
 
  //fetch the TFile for the data and MC spectrum
  TFile *f_data = new TFile(fakedata.c_str(),"read");
  TH1D *h_nue_np_fake_data;
  if(tag.find("fakedata") != std::string::npos ) h_nue_np_fake_data = (TH1D*)f_data->Get("nu_uBooNE_1eNp_data");
  TH1D *h_nue_0p_fake_data;
  if(tag.find("fakedata") != std::string::npos ) h_nue_0p_fake_data = (TH1D*)f_data->Get("nu_uBooNE_1e0p_data");
  TH1D *h_numu_data = (TH1D*)f_data->Get("nu_uBooNE_numu_data");
  
  for(int bin=0; bin < numu_bins; bin++) numu_data_vec.push_back(h_numu_data->GetBinContent(bin+1)); 
  for(int bin=0; bin < numu_bins; bin++ ) numu_vec.push_back(sig_collvec[nue_bins+bin]);
  for(int bin=0; bin < numu_bins; bin++ ) delta_numu.push_back(numu_data_vec[bin]-numu_vec[bin]);
  for(int bin=0; bin < numu_bins; bin++ ) std::cout << "numu data, numu mc, delta_numu = " << numu_data_vec[bin] << ", " << numu_vec[bin] << ", " << delta_numu[bin] << std::endl; 
  
  //Load up our covariance matricies we calculated in example1 (we could also load up single variation ones)
  TFile * fsys = new TFile(Form("../bin/%s.SBNcovar.root",tag.c_str()),"read");
  TMatrixD *cov = (TMatrixD*)fsys->Get("collapsed_covariance");
  TMatrixD *fullcov = (TMatrixD*)fsys->Get("full_covariance");
  TMatrixD *cov0 = (TMatrixD*)cov->Clone("coll_covar");
  //Recalculate cov matrix by adding mc intrinsic error
  TMatrixD collapsed;
  collapsed.ResizeTo(sig_spectra.num_bins_total_compressed, sig_spectra.num_bins_total_compressed);
  fullcov->Zero();
  TMatrixT<double> fullsig = chi_h1.CalcCovarianceMatrix(fullcov, sig_fullvec, sig_errfullvec, false);  
  chi_h1.CollapseModes(fullsig, collapsed);
  collapsed.Print();
  if(mcerr){ 
	std::cout << "ADD MC ERR" << std::endl;
	tag=tag+"_with_mc_err";
        for(int i=0; i < detsyscoll.GetNrows(); i++) std::cout << tag << "  Matrix diag mc err, sample: " << sqrt((*cov)(i,i)) << ", " << collapsed(i,i) << ", " << sig_collvec[i] << std::endl; 
  	*cov += collapsed;
  }
  //create the zero bin ext bnb err matrix
  TMatrixD fullext;
  fullext.ResizeTo(sig_spectra.num_bins_total, sig_spectra.num_bins_total);
  for(int i=0; i < sig_spectra.num_bins_total; i++ ) fullext(i,i) += sig_exterrfullvec[i]*sig_exterrfullvec[i];

  TMatrixD collext;
  collext.ResizeTo(sig_spectra.num_bins_total_compressed, sig_spectra.num_bins_total_compressed);
  //collapse...
  chi_h1.CollapseModes(fullext, collext);
  if(addext){
	std::cout << "ADD ZERO BIN ERROR" << std::endl; 
	tag=tag+"_with_zerobin_err";
        //add to full covariance matrix (collapsed)
        *cov += collext;
  }
  if(detsys)  *cov += detsyscoll;
  //frac cov
  for(int i=0; i < detsyscoll.GetNrows(); i++){
  for(int j=0; j < detsyscoll.GetNcols(); j++){
    fraccov(i,j) = 0.0;
    if(sig_collvec[i] !=0 && sig_collvec[j] !=0 ) fraccov(i,j) = (*cov)(i,j)/(sig_collvec[i]*sig_collvec[j]);
    corr(i,j) = 0.0;
    if((*cov)(i,j) != 0 ) corr(i,j) = (*cov)(i,j)/sqrt((*cov)(i,i)*(*cov)(j,j));
    if(i==j) std::cout << tag << "  Matrix diag: " << sqrt((*cov)(i,i)) << ", " << sqrt(fraccov(i,j)) << std::endl; 
  }
  }

  if(!addext) collext.Zero();
  if(!mcerr) collapsed.Zero();

  for(int i=0; i < sig_spectra.num_bins_total_compressed; i++ ){ collapsed(i,i) += collext(i,i); std::cout << "coll ext err: " << collext(i,i) << std::endl; }

  //calculate cov matrix before constraint
  //frac cov
  for(int i=0; i < (*cov).GetNrows(); i++){
  for(int j=0; j < (*cov).GetNcols(); j++){
    fraccoll(i,j) = 0.0;
    if(sig_collvec[i] != 0.0 && sig_collvec[j] != 0.0) fraccoll(i,j) = (*cov)(i,j)/(sig_collvec[i]*sig_collvec[j]);
    fraccoll_noext(i,j) = 0.0;
    if(sig_collvec_noext[i] != 0.0 && sig_collvec_noext[j] != 0.0) fraccoll_noext(i,j) = (*cov)(i,j)/(sig_collvec_noext[i]*sig_collvec_noext[j]);
    corrcoll(i,j) = 0.0;
    if((*cov)(i,j) != 0 ) corrcoll(i,j) = (*cov)(i,j)/sqrt((*cov)(i,i)*(*cov)(j,j));
  }
  }

  //write out the cov matrix before constraint only do this for the combined channels
  if(combined){
    for(int i=0; i < (*cov).GetNrows(); i++){
      for(int j=0; j < (*cov).GetNcols(); j++){
        if( i==j && i < np_bins ) outfile1  << sqrt(fraccoll(i,j)) << " " << sqrt(fraccoll_noext(i,j)) << std::endl;
        else if( i==j && i >= np_bins && i < np_bins+np_bins ) outfile3  << sqrt(fraccoll(i,j)) << " " << sqrt(fraccoll_noext(i,j)) << std::endl;
        else if( i==j && i >= np_bins+np_bins && i < np_bins+np_bins+zp_bins ) outfile5  << sqrt(fraccoll(i,j)) << " " << sqrt(fraccoll_noext(i,j)) << std::endl;
        else if( i==j && i >= np_bins+np_bins+zp_bins && i < np_bins+np_bins+zp_bins+zp_bins ) outfile7  << sqrt(fraccoll(i,j)) << " " << sqrt(fraccoll_noext(i,j)) << std::endl;
        else if( i==j && i >= np_bins+np_bins+zp_bins+zp_bins && i < np_bins+np_bins+zp_bins+zp_bins+numu_bins ) outfile9  << sqrt(fraccoll(i,j)) << " " << sqrt(fraccoll_noext(i,j)) << std::endl;
      }
    }
  }

  plot_one(corrcoll, sig_collvec, channel_names_coll, num_bins_coll, tag+"_correlation_matrix_xsecfluxdetsys", true);
  plot_one(fraccoll, sig_collvec, channel_names_coll, num_bins_coll, tag+"_fractional_covariance_matrix_xsecfluxdetsys", true);


  //Draw before constraint
  TH1D* h_1eNp_lee_before = (TH1D*) h_1eNp_bg_ext->Clone("1eNp_lee_before");
  TH1D* h_1eNp_bg_before  = (TH1D*) h_1eNp_bg_ext->Clone("1eNp_bg_before");
  TH1D* h_1e0p_lee_before = (TH1D*) h_1e0p_bg_ext->Clone("1e0p_lee_before");
  TH1D* h_1e0p_bg_before  = (TH1D*) h_1e0p_bg_ext->Clone("1e0p_bg_before");
  TH1D* h_numu_before     = (TH1D*) h_numu_data->Clone("numu_before"); 
  
  //Draw before constraint
  h_1eNp_lee_before->Reset();
  h_1eNp_bg_before->Reset();
  h_1e0p_lee_before->Reset();
  h_1e0p_bg_before->Reset();
  h_numu_before->Reset();
    
  TString cname;
  
  TH1D* histo_dummy = (TH1D*)h_numu_data->Clone("dummy");
  histo_dummy->Reset();
  
  if(combined){     
    std::cout << "sig_collvec size = " << sig_collvec.size() << std::endl;
    //set bin content
    for(int bin=0; bin < np_bins; bin++ ) h_1eNp_lee_before->SetBinContent(bin+1,sig_collvec[bin]);
    for(int bin=0; bin < np_bins; bin++ ) h_1eNp_bg_before->SetBinContent(bin+1,sig_collvec[bin+np_bins]);
    for(int bin=0; bin < zp_bins; bin++ ) h_1e0p_lee_before->SetBinContent(bin+1,sig_collvec[bin+2*np_bins]);
    for(int bin=0; bin < zp_bins; bin++ ) h_1e0p_bg_before->SetBinContent(bin+1,sig_collvec[bin+2*np_bins+zp_bins]);
    for(int bin=0; bin < numu_bins; bin++ ) h_numu_before->SetBinContent(bin+1,sig_collvec[bin+2*np_bins+2*zp_bins]);
    //set bin error
    for(int bin=0; bin < np_bins; bin++ ) h_1eNp_lee_before->SetBinError(bin+1,sqrt((*cov)[bin][bin]));
    for(int bin=0; bin < np_bins; bin++ ) h_1eNp_bg_before->SetBinError(bin+1,sqrt((*cov)[bin+np_bins][bin+np_bins]));
    for(int bin=0; bin < zp_bins; bin++ ) h_1e0p_lee_before->SetBinError(bin+1,sqrt((*cov)[bin+2*np_bins][bin+2*np_bins]));
    for(int bin=0; bin < zp_bins; bin++ ) h_1e0p_bg_before->SetBinError(bin+1,sqrt((*cov)[bin+2*np_bins+zp_bins][bin+2*np_bins+zp_bins]));
    for(int bin=0; bin < numu_bins; bin++ ) h_numu_before->SetBinError(bin+1,sqrt((*cov)[bin+2*np_bins+2*zp_bins][bin+2*np_bins+2*zp_bins]));
  }else if(np){
    //set bin content
    for(int bin=0; bin < np_bins; bin++ ) h_1eNp_lee_before->SetBinContent(bin+1,sig_collvec[bin]);
    for(int bin=0; bin < np_bins; bin++ ) h_1eNp_bg_before->SetBinContent(bin+1,sig_collvec[bin+np_bins]);
    for(int bin=0; bin < numu_bins; bin++ ) h_numu_before->SetBinContent(bin+1,sig_collvec[bin+2*np_bins]);
    //set bin error
    for(int bin=0; bin < np_bins; bin++ ) h_1eNp_lee_before->SetBinError(bin+1,sqrt((*cov)[bin][bin]));
    for(int bin=0; bin < np_bins; bin++ ) h_1eNp_bg_before->SetBinError(bin+1,sqrt((*cov)[bin+np_bins][bin+np_bins]));
    for(int bin=0; bin < numu_bins; bin++ ) h_numu_before->SetBinError(bin+1,sqrt((*cov)[bin+2*np_bins][bin+2*np_bins]));
  }else if(zp){
    //set bin content
    for(int bin=0; bin < zp_bins; bin++ ) h_1e0p_lee_before->SetBinContent(bin+1,sig_collvec[bin]);
    for(int bin=0; bin < zp_bins; bin++ ) h_1e0p_bg_before->SetBinContent(bin+1,sig_collvec[bin+zp_bins]);
    for(int bin=0; bin < numu_bins; bin++ ) h_numu_before->SetBinContent(bin+1,sig_collvec[bin+2*zp_bins]);
    //set bin error
    for(int bin=0; bin < zp_bins; bin++ ) h_1e0p_lee_before->SetBinError(bin+1,sqrt((*cov)[bin][bin]));
    for(int bin=0; bin < zp_bins; bin++ ) h_1e0p_bg_before->SetBinError(bin+1,sqrt((*cov)[bin+zp_bins][bin+zp_bins]));
    for(int bin=0; bin < numu_bins; bin++ ) h_numu_before->SetBinError(bin+1,sqrt((*cov)[bin+2*zp_bins][bin+2*zp_bins]));
  }

  if(detsys) tag += "_detsys";
  cname = tag + "_before_constraint_numu";
  DrawDataMCAndSyst(cname, h_numu_before, h_numu_before, h_numu_data, "Reconstructed Visible Energy [GeV]", "#nu_{#mu} Selection");
	
  if(combined){
    //numu before constraint
    TH1D* h_1eNp_sum_before = (TH1D*)h_1eNp_lee_before->Clone("Np_sum_lee_before");
    h_1eNp_sum_before->Add(h_1eNp_bg_before);
    TH1D* h_1e0p_sum_before = (TH1D*)h_1e0p_lee_before->Clone("Zp_sum_lee_before");
    h_1e0p_sum_before->Add(h_1e0p_bg_before);
    
    cname = tag + "_before_constraint_1eNp";
    DrawDataMCAndSyst(cname, h_1eNp_sum_before, h_1eNp_bg_before, histo_dummy, "Reconstructed Visible Energy [GeV]", "#nu_{e} 1eNp0#pi Selection");	
    cname = tag + "_before_constraint_1e0p";
    DrawDataMCAndSyst(cname, h_1e0p_sum_before, h_1e0p_bg_before, histo_dummy, "Reconstructed Visible Energy [GeV]", "#nu_{e} 1e0p0#pi Selection");	
    
  }
  
  //Block Matrices
  TMatrixD numumatrix(numu_bins,numu_bins);
  TMatrixD InvertNumu(numu_bins,numu_bins);
  TMatrixD nuenumumatrix(nue_bins,numu_bins);
  TMatrixD numunuematrix(numu_bins,nue_bins);
  TMatrixD nuematrix(nue_bins,nue_bins);
  
  for( int i = 0; i < cov->GetNcols(); i++ ){
    for( int j = 0; j < cov->GetNrows(); j++ ){
      if( i >= nue_bins && j >= nue_bins) numumatrix(i-nue_bins,j-nue_bins) = (*cov)(i,j);
      if( i >= nue_bins && j < nue_bins) numunuematrix(i-nue_bins,j) = (*cov)(i,j);
      if( i < nue_bins && j >= nue_bins) nuenumumatrix(i,j-nue_bins) = (*cov)(i,j);
      if( i < nue_bins && j < nue_bins) nuematrix(i,j) = (*cov)(i,j);
    }
  }
  
  //add the stats error to the cov matrix
  //we are adding the mc stats error which include the mc stats error for the MC numu
  //only adding this for the option where we do not add the mc stats error
  for( int i = 0; i < numumatrix.GetNcols(); i++ ){ 
       numumatrix(i,i) += numu_vec[i];
  }

  //constrained procedure
  InvertNumu = numumatrix;
  InvertNumu.Zero();
  TDecompSVD svdnumu( numumatrix );
  
  if (!svdnumu.Decompose() ) {
    std::cout << "Decomposition failed, matrix singular ?" << std::endl;
  }else{
    InvertNumu = numumatrix.Invert();
    InvertNumu.Print();
  }
  
  TMatrixD A(nuenumumatrix.GetNrows(),nuenumumatrix.GetNcols());
  TMatrixD B(nuematrix.GetNrows(),nuematrix.GetNcols());
  TMatrixD C(nuenumumatrix.GetNrows(),nuenumumatrix.GetNcols());
  
  A.Mult(nuenumumatrix,InvertNumu);
  B.Mult(A,numunuematrix);
  C.Mult(nuenumumatrix,InvertNumu);
  
  //constrained spectrum 
  std::vector<double> input_nue_constrained;
  std::vector<double> input_nue_constrained_noext;
  
  for( int binx = 0; binx < nue_bins; binx++ ){
    double scaled = 0.;
    for( int biny = 0; biny < delta_numu.size(); biny++ ){
      scaled += C(binx,biny) * delta_numu[biny];
      //std::cout << "binx, biny, delta_numu = " << binx << ", " << biny << ", " << delta_numu[biny] << ", " << C(binx,biny) << std::endl;
    }
    double constnue =  sig_collvec[binx] + scaled;
    double constnue_noext =  sig_collvec[binx] + scaled - collapsed_ext[binx];
    std::cout << "check that constrained results make sense: binx, scaled = " << binx << ", " << scaled << std::endl;
    input_nue_constrained.push_back(constnue);
    input_nue_constrained_noext.push_back(constnue_noext);
  }
  TMatrixD constnuematrix(nuematrix.GetNrows(),nuematrix.GetNcols());
  TMatrixD fracnuematrix(nuematrix.GetNrows(),nuematrix.GetNcols());
  TMatrixD fracnuematrix_noext(nuematrix.GetNrows(),nuematrix.GetNcols());
  TMatrixD corrnuematrix(nuematrix.GetNrows(),nuematrix.GetNcols());
  for( int i = 0; i < nuematrix.GetNrows(); i++ ){
    for( int j = 0; j < nuematrix.GetNcols(); j++ ){
      constnuematrix(i,j) = nuematrix(i,j) - B(i,j);
      if( input_nue_constrained[i] != 0 && input_nue_constrained[j] != 0 ) fracnuematrix(i,j) = constnuematrix(i,j)/(input_nue_constrained[i]*input_nue_constrained[j]);
      else fracnuematrix(i,j) = 0.0;
      if( input_nue_constrained_noext[i] != 0 && input_nue_constrained_noext[j] != 0 ) fracnuematrix_noext(i,j) = constnuematrix(i,j)/(input_nue_constrained_noext[i]*input_nue_constrained_noext[j]);
      else fracnuematrix_noext(i,j) = 0.0;
      if( combined && i==j && i < np_bins ) outfile2  << sqrt(fracnuematrix(i,j)) << ", " << fracnuematrix_noext(i,j) << std::endl;
      else if( combined && i==j && i >= np_bins && i < 2*np_bins ) outfile4  << sqrt(fracnuematrix(i,j)) << " " << sqrt(fracnuematrix_noext(i,j)) << std::endl;
      else if( combined && i==j && i >= 2*np_bins && i < 2*np_bins+zp_bins ) outfile6  << sqrt(fracnuematrix(i,j)) << " " << sqrt(fracnuematrix_noext(i,j)) << std::endl;
      else if( combined && i==j && i >= 2*np_bins+zp_bins && i < 2*np_bins+2*zp_bins ) outfile8  << sqrt(fracnuematrix(i,j)) << " " << sqrt(fracnuematrix_noext(i,j)) << std::endl;
      corrnuematrix(i,j) = constnuematrix(i,j)/sqrt(constnuematrix(i,i)*constnuematrix(j,j));
      if(i==j) std::cout << i << " fracnuematrix = " << fracnuematrix(i,j) << std::endl;
    }
  }
   
  //Draw after constraint
  TH1D* h_1eNp_lee;
  TH1D* h_1eNp_bg;
  TH1D* h_1eNp_data;
  TH1D* h_1e0p_lee;
  TH1D* h_1e0p_bg;
  TH1D* h_1e0p_data; 
  
  h_1eNp_lee  = (TH1D*)h_1eNp_bg_ext->Clone("nu_uBooNE_1eNp_lee"); 
  h_1eNp_bg   = (TH1D*)h_1eNp_bg_ext->Clone("nu_uBooNE_1eNp_bg");
  h_1eNp_data = (TH1D*)h_1eNp_bg_ext->Clone("nu_uBooNE_1eNp_data");
  h_1e0p_lee  = (TH1D*)h_1e0p_bg_ext->Clone("nu_uBooNE_1e0p_lee");
  h_1e0p_bg   = (TH1D*)h_1e0p_bg_ext->Clone("nu_uBooNE_1e0p_bg");
  h_1e0p_data = (TH1D*)h_1e0p_bg_ext->Clone("nu_uBooNE_1e0p_data");
  
  //reset the histograms
  h_1eNp_lee->Reset(); 
  h_1eNp_bg->Reset();   
  h_1eNp_data->Reset();
  h_1e0p_lee->Reset();
  h_1e0p_bg->Reset();
  h_1e0p_data->Reset();
 
  for(int b=0; b < collapsed.GetNcols(); b++ ) std::cout << "collapsed, sqrt = " << collapsed[b][b] << ", " << sqrt(collapsed[b][b]) << std::endl; 
  if(combined){
    //set bin content
    for(int bin=0; bin < np_bins; bin++ ) h_1eNp_lee->SetBinContent(bin+1,input_nue_constrained[bin]);
    for(int bin=0; bin < np_bins; bin++ ) h_1eNp_bg->SetBinContent(bin+1,input_nue_constrained[bin+np_bins]);
    for(int bin=0; bin < np_bins; bin++ ) h_1eNp_data->SetBinContent(bin+1,0.);
    for(int bin=0; bin < zp_bins; bin++ ) h_1e0p_lee->SetBinContent(bin+1,input_nue_constrained[bin+2*np_bins]);
    for(int bin=0; bin < zp_bins; bin++ ) h_1e0p_bg->SetBinContent(bin+1,input_nue_constrained[bin+2*np_bins+zp_bins]);
    for(int bin=0; bin < zp_bins; bin++ ) h_1e0p_data->SetBinContent(bin+1,0);
    //set bin content error
    for(int bin=0; bin < np_bins; bin++ ) h_1eNp_lee->SetBinError(bin+1,sqrt(constnuematrix[bin][bin]));
    for(int bin=0; bin < np_bins; bin++ ) h_1eNp_bg->SetBinError(bin+1,sqrt(constnuematrix[bin+np_bins][bin+np_bins]));
    for(int bin=0; bin < np_bins; bin++ ) h_1eNp_data->SetBinError(bin+1,0.);
    for(int bin=0; bin < zp_bins; bin++ ) h_1e0p_lee->SetBinError(bin+1,sqrt(constnuematrix[bin+2*np_bins][bin+2*np_bins]));
    for(int bin=0; bin < zp_bins; bin++ ) h_1e0p_bg->SetBinError(bin+1,sqrt(constnuematrix[bin+2*np_bins+zp_bins][2*np_bins+zp_bins]));
    for(int bin=0; bin < zp_bins; bin++ ) h_1e0p_data->SetBinError(bin+1,0.);
  }else if(np){
    for(int bin=0; bin < np_bins; bin++ ) h_1eNp_lee->SetBinContent(bin+1,input_nue_constrained[bin]);
    for(int bin=0; bin < np_bins; bin++ ) h_1eNp_bg->SetBinContent(bin+1,input_nue_constrained[bin+np_bins]);
    for(int bin=0; bin < np_bins; bin++ ) h_1eNp_data->SetBinContent(bin+1,0.);
    for(int bin=0; bin < np_bins; bin++ ) h_1eNp_lee->SetBinError(bin+1,sqrt(constnuematrix[bin][bin]));
    for(int bin=0; bin < np_bins; bin++ ) h_1eNp_bg->SetBinError(bin+1,sqrt(constnuematrix[bin+np_bins][bin+1*np_bins]));
    for(int bin=0; bin < np_bins; bin++ ) h_1eNp_data->SetBinError(bin+1,0.);
  }else if(zp){
    for(int bin=0; bin < zp_bins; bin++ ) h_1e0p_lee->SetBinContent(bin+1,input_nue_constrained[bin]);
    for(int bin=0; bin < zp_bins; bin++ ) h_1e0p_bg->SetBinContent(bin+1,input_nue_constrained[bin+zp_bins]);
    for(int bin=0; bin < zp_bins; bin++ ) h_1e0p_data->SetBinContent(bin+1,0.);
    for(int bin=0; bin < zp_bins; bin++ ) h_1e0p_lee->SetBinError(bin+1,sqrt(constnuematrix[bin][bin]));
    for(int bin=0; bin < zp_bins; bin++ ) h_1e0p_bg->SetBinError(bin+1,sqrt(constnuematrix[bin+zp_bins][bin+zp_bins]));
    for(int bin=0; bin < zp_bins; bin++ ) h_1e0p_data->SetBinError(bin+1,0.);
  }
 
  /*if(combined){ 
    for(int bin=0; bin < 14; bin++ ){ 
      std::cout << "@@ before constraint " << tag << std::endl;
      std::cout << "@@ 1eNp lee, nue, ratio: " <<  h_1eNp_lee_before->GetBinContent(bin+1) << ", " << h_1eNp_bg_before->GetBinContent(bin+1) << ", " << -h_1eNp_lee_before->GetBinContent(bin+1)+h_1eNp_bg_before->GetBinContent(bin+1) << std::endl;     
      std::cout << "@@ 1e0p lee, nue, ratio: " <<  h_1e0p_lee_before->GetBinContent(bin+1) << ", " << h_1e0p_bg_before->GetBinContent(bin+1) << ", " << h_1e0p_lee_before->GetBinContent(bin+1)/h_1e0p_bg_before->GetBinContent(bin+1) << std::endl; 
    }
    for(int bin=0; bin < 14; bin++ ){ 
      std::cout << "@@ after constraint " << tag << std::endl;
      std::cout << "@@ 1eNp lee, nue, ratio: " <<  h_1eNp_lee->GetBinContent(bin+1) << ", " << h_1eNp_bg->GetBinContent(bin+1) << ", " << -h_1eNp_lee->GetBinContent(bin+1)+h_1eNp_bg->GetBinContent(bin+1) << std::endl; 
      std::cout << "@@ 1e0p lee, nue, ratio: " <<  h_1e0p_lee->GetBinContent(bin+1) << ", " << h_1e0p_bg->GetBinContent(bin+1) << ", " << h_1e0p_lee->GetBinContent(bin+1)/h_1e0p_bg->GetBinContent(bin+1) << std::endl;
    }
  }*/

  //set bin content error
  if(combined){
    //nue before constraint
    TH1D* h_1eNp_sum_before = (TH1D*)h_1eNp_lee_before->Clone("Np_sum_lee_before");
    h_1eNp_sum_before->Add(h_1eNp_bg_before);
    TH1D* h_1e0p_sum_before = (TH1D*)h_1e0p_lee_before->Clone("Zp_sum_lee_before");
    h_1e0p_sum_before->Add(h_1e0p_bg_before);
    
    TH1D* h_1eNp_sum = (TH1D*)h_1eNp_bg->Clone("1eNp_sum_bg");
    h_1eNp_sum->Add(h_1eNp_lee);
    TH1D* h_1e0p_sum = (TH1D*)h_1e0p_bg->Clone("1e0p_sum_bg");
    h_1e0p_sum->Add(h_1e0p_lee);
    
    //Draw before and after constraint
    cname = tag + "_after_constraint_1e0p";
    DrawDataMCAndSyst(cname, h_1e0p_sum, h_1e0p_bg, histo_dummy, "Reconstructed Visible Energy [GeV]", "#nu_{e} Selection");	
    cname = tag + "_after_constraint_1eNp";
    DrawDataMCAndSyst(cname, h_1eNp_sum, h_1eNp_bg, histo_dummy, "Reconstructed Visible Energy [GeV]", "#nu_{e} Selection");	
    
    //Draw before and after constraint
    cname = tag + "_before_after_constraint_1e0p_SUM";
    DrawDataMCAndSyst(cname, h_1e0p_sum_before, h_1e0p_sum, histo_dummy, "Reconstructed Visible Energy [GeV]", "#nu_{e} 1eNp0#pi Selection");	
    cname = tag + "_before_after_constraint_1eNp_SUM";
    DrawDataMCAndSyst(cname, h_1eNp_sum_before, h_1eNp_sum, histo_dummy, "Reconstructed Visible Energy [GeV]", "#nu_{e} 1eNp0#pi Selection");	
    cname = tag + "_before_after_constraint_1e0p_BKG";
    DrawDataMCAndSyst(cname, h_1e0p_sum_before, h_1e0p_sum, histo_dummy, "Reconstructed Visible Energy [GeV]", "#nu_{e} 1e0p0#pi Selection");	
    cname = tag + "_before_after_constraint_1eNp_BKG";
    DrawDataMCAndSyst(cname, h_1eNp_sum_before, h_1eNp_sum, histo_dummy, "Reconstructed Visible Energy [GeV]", "#nu_{e} 1eNp0#pi Selection");	
  }
  
  TFile * fspec = new TFile(Form("constrained_%s.SBNspec.root",tag.c_str()),"recreate");
  fspec->cd();
  if(combined){
    (TH1D*)h_1eNp_lee->Write("nu_uBooNE_1eNp_lee");
    (TH1D*)h_1eNp_bg->Write("nu_uBooNE_1eNp_bg");
    (TH1D*)h_1e0p_lee->Write("nu_uBooNE_1e0p_lee");
    (TH1D*)h_1e0p_bg->Write("nu_uBooNE_1e0p_bg");
  }else if(np){
    std::cout << "filling np 3" << std::endl;
    (TH1D*)h_1eNp_lee->Write("nu_uBooNE_1eNp_lee");
    (TH1D*)h_1eNp_bg->Write("nu_uBooNE_1eNp_bg");
  }else if(zp){
    (TH1D*)h_1e0p_lee->Write("nu_uBooNE_1e0p_lee");
    (TH1D*)h_1e0p_bg->Write("nu_uBooNE_1e0p_bg");
  }
  fspec->Close();
  
  TFile * fspec2 = new TFile(Form("unconstrained_%s.SBNspec.root",tag.c_str()),"recreate");
  fspec2->cd();
  if(combined){
    h_1eNp_lee_before->SetName("nu_uBooNE_1eNp_lee");
    h_1eNp_lee_before->SetTitle("nu_uBooNE_1eNp_lee");
    (TH1D*)h_1eNp_lee_before->Write("nu_uBooNE_1eNp_lee");
    h_1eNp_bg_before->SetName("nu_uBooNE_1eNp_bg");
    h_1eNp_bg_before->SetTitle("nu_uBooNE_1eNp_bg");
    (TH1D*)h_1eNp_bg_before->Write("nu_uBooNE_1eNp_bg");
    h_1e0p_lee_before->SetName("nu_uBooNE_1e0p_lee");
    h_1e0p_lee_before->SetTitle("nu_uBooNE_1e0p_lee");
    (TH1D*)h_1e0p_lee_before->Write("nu_uBooNE_1e0p_lee");
    h_1e0p_bg_before->SetName("nu_uBooNE_1e0p_bg");
    h_1e0p_bg_before->SetTitle("nu_uBooNE_1e0p_bg");
    (TH1D*)h_1e0p_bg_before->Write("nu_uBooNE_1e0p_bg");
  }else if(np){
    h_1eNp_lee_before->SetName("nu_uBooNE_1eNp_lee");
    h_1eNp_lee_before->SetTitle("nu_uBooNE_1eNp_lee");
    (TH1D*)h_1eNp_lee_before->Write("nu_uBooNE_1eNp_lee");
    h_1eNp_bg_before->SetName("nu_uBooNE_1eNp_bg");
    h_1eNp_bg_before->SetTitle("nu_uBooNE_1eNp_bg");
    (TH1D*)h_1eNp_bg_before->Write("nu_uBooNE_1eNp_bg");
  }else if(zp){
    h_1e0p_lee_before->SetName("nu_uBooNE_1e0p_lee");
    h_1e0p_lee_before->SetTitle("nu_uBooNE_1e0p_lee");
    (TH1D*)h_1e0p_lee_before->Write("nu_uBooNE_1e0p_lee");
    h_1e0p_bg_before->SetName("nu_uBooNE_1e0p_bg");
    h_1e0p_bg_before->SetTitle("nu_uBooNE_1e0p_bg");
    (TH1D*)h_1e0p_bg_before->Write("nu_uBooNE_1e0p_bg");
  }
  fspec2->Close();

  std::cout << "constnuematrix size = " << constnuematrix.GetNcols() << ", " << constnuematrix.GetNrows() << std::endl;
  TFile * fcovar = new TFile(Form("constrained_%s.SBNcovar.root",tag.c_str()),"recreate");
  fcovar->cd();
  (TMatrixD*)constnuematrix.Write("full_covariance");
  (TMatrixD*)fracnuematrix.Write("frac_covariance");
  (TMatrixD*)corrnuematrix.Write("full_correlation");
  fcovar->Close();
  
  TFile * fcovar2 = new TFile(Form("unconstrained_%s.SBNcovar.root",tag.c_str()),"recreate");
  fcovar2->cd();
  cov->Write("full_covariance");
  (TMatrixD*)fraccov.Write("frac_covariance");
  (TMatrixD*)corr.Write("full_correlation");
  fcovar2->Close();

  TFile * fcollext = new TFile(Form("collapsed_zero_ext_bin_%s.SBNcovar.root",tag.c_str()),"recreate");
  fcollext->cd();
  (TMatrixD*)collext.Write("full_covariance");
  fcollext->Close();

  return 0;
  
}

void DrawDataMCAndSyst(TString cName, TH1D* MC1, TH1D* MC2, TH1D* data, std::string var, std::string title){
  
  TCanvas can(cName,cName,1200,800);
  
  TH1D* tmpMC1 = (TH1D*)MC1->Clone("tempCV");
  TH1D* tmpMC2 = (TH1D*)MC2->Clone("tempCVnolee");
  TH1D* sysErr1 = (TH1D*)MC1->Clone("tempSysCV");
  TH1D* sysErr2 = (TH1D*)MC2->Clone("tempSysCVnolee");
  TH1D* tmpData = (TH1D*)data->Clone("tempData");

  TH1D* tmpratio1 = (TH1D*)data->Clone("tempRatio1");
  TH1D* tmpMC1_staterr = (TH1D*)MC1->Clone("tempCV_staterr2");
  for( int bin = 1; bin < tmpMC1_staterr->GetNbinsX()+1; bin++){ tmpMC1_staterr->SetBinError(bin,sqrt(tmpMC1_staterr->GetBinContent(bin))); }
  tmpratio1->Reset();
  tmpratio1->Divide(tmpData,tmpMC1_staterr,1.0,1.0,"B");
  if(verbose) std::cout << "nbins data, mc = " << tmpData->Integral() << ", " << tmpMC1->Integral() << std::endl;

  TH1D* tmpratio2 = (TH1D*)data->Clone("tempRatio2");
  TH1D* tmpMC2_staterr = (TH1D*)MC2->Clone("tempCV_staterr2");
  for( int bin = 1; bin < tmpMC2_staterr->GetNbinsX()+1; bin++){ tmpMC2_staterr->SetBinError(bin,sqrt(tmpMC2_staterr->GetBinContent(bin))); }
  tmpratio2->Reset();
  tmpratio2->Divide(tmpData,tmpMC2_staterr,1.0,1.0,"B");
  if(verbose) std::cout << "nbins data, mc nolee = " << tmpData->Integral() << ", " << tmpMC2->Integral() << std::endl;

  TH1D* tmpsys = (TH1D*)MC2->Clone("tempSysRatio");
  TH1D* tmpsys2 = (TH1D*)MC1->Clone("tempSysRatio1");

  for( int bin = 1; bin < tmpsys->GetNbinsX()+1; bin++){
    tmpsys->SetBinContent(bin,1.0);
    tmpsys->SetBinError(bin,tmpMC2->GetBinError(bin)/tmpMC2->GetBinContent(bin));
    if(verbose) std::cout << "tempSys MC2 " << tmpsys->GetBinError(bin) << std::endl;
    tmpsys2->SetBinContent(bin,1.0);
    tmpsys2->SetBinError(bin,tmpMC1->GetBinError(bin)/tmpMC1->GetBinContent(bin));
    if(verbose) std::cout << "tempSys MC1 " << tmpsys2->GetBinError(bin) << std::endl;
    sysErr1->SetBinContent(bin,tmpMC1->GetBinContent(bin));
    sysErr1->SetBinError(bin,tmpMC1->GetBinError(bin));
    sysErr2->SetBinContent(bin,tmpMC2->GetBinContent(bin));
    sysErr2->SetBinError(bin,tmpMC2->GetBinError(bin));
  }
  
  // Upper plot will be in pad1
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.06, 1.0, 1.0);
  //TPad *pad1 = new TPad("pad1", "pad1", 0, 0.45, 1, 1.0); //if we want to split the plot again
  //pad1->SetBottomMargin(0.05); // Upper and lower plot are joined
  pad1->SetGridx(2);        // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad

  double max = 1.1;
  double max1 = 1.1;
  double maxMC1 = tmpMC1->GetMaximum();
  double maxMC2 = tmpMC2->GetMaximum();
  double maxData = tmpData->GetMaximum();
  if(maxMC1 > maxData) max1 *= maxMC1;
  if(maxData > maxMC1 ) max1 *= maxData;
  if(maxMC2 > max1) max *= maxMC2;
  else max = max1;

  std::cout << "cName = " << std::endl;
  if(std::string(cName).find("1eNp") != std::string::npos) max = 12.0;
  if(std::string(cName).find("1e0p") != std::string::npos) max = 8.0;

  sysErr1->SetTitleSize(80);
  sysErr1->SetTitleFont(63);
  sysErr1->SetTitle(title.c_str());
  sysErr1->SetStats(0);
  sysErr1->GetYaxis()->SetTitleSize(30);
  sysErr1->GetYaxis()->SetTitleFont(43); //apparently changing 43->43 causes the title to not be printed on canvas
  sysErr1->GetYaxis()->SetTitle("Events/0.1 GeV");
  sysErr1->GetXaxis()->SetTitle("Reconstructed Visible Energy [GeV]");
  sysErr1->GetXaxis()->SetTitleFont(43);
  sysErr1->GetXaxis()->SetTitleSize(30);
  sysErr1->GetXaxis()->SetTitleOffset(1.25);
  sysErr1->GetXaxis()->SetLabelSize(0.04);

  tmpMC1->SetLineColor(kBlue);
  tmpMC1->SetLineWidth(3);
  tmpMC1->SetLineStyle(1);
  tmpMC2->SetLineColor(kRed);
  tmpMC2->SetLineWidth(3);
  tmpMC2->SetLineStyle(1);
  tmpData->SetLineWidth(2);
  tmpData->SetMarkerStyle(20);
  tmpData->SetLineColor(kBlack);
  if(std::string(cName).find("nowgt") != std::string::npos && std::string(cName).find("nue") != std::string::npos ) tmpData->SetLineColor(kBlack);

  sysErr1->SetLineColor(kBlue-10);
  sysErr1->SetLineWidth(1);
  sysErr1->SetFillColor(kBlue-10);
  sysErr1->SetFillStyle(1001);

  sysErr2->SetLineColor(kRed-10);
  sysErr2->SetLineWidth(1);
  sysErr2->SetFillColor(kRed-10);
  sysErr2->SetFillStyle(1001);

  sysErr1->SetMaximum(1.2*max);
  sysErr1->SetMinimum(0.);

  sysErr1->Draw("E2");
  sysErr2->Draw("E2 same");
  tmpMC1->Draw("hist same"); 
  tmpMC2->Draw("hist same"); 
  if(title == "#nu_{#mu} Selection")tmpData->Draw("PE same"); 
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
    y1=0.6;
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
      legend->AddEntry(tmpMC2,"#nu_{#mu} CC","l"); 
      legend->AddEntry(tmpData,"data","l");
    }
  }else{
    if(std::string(cName).find("nowgt") != std::string::npos ){
      legend->AddEntry(tmpMC1,"#nu_e, Genie Wgt","l"); 
      legend->AddEntry(tmpMC2,"#nu_e, no Genie wgt","l"); 
      legend->AddEntry(tmpData,"#nu_e, post-constraint","l");
    }else if(std::string(cName).find("before_after") != std::string::npos ){
      legend->AddEntry(tmpMC1,"before constraint","l"); 
      legend->AddEntry(tmpMC2,"after constraint","l"); 
    }else{
      legend->AddEntry(tmpMC1,Form("LEE: %.2f",tmpMC1->Integral()-tmpMC2->Integral()),"l"); 
      legend->AddEntry(tmpMC2,Form("intrinsic: %.2f",tmpMC2->Integral()),"le");
      //legend->AddEntry(tmpData,"fake data","le");
    }
  } 
  legend->Draw("same");

  //remove the ratio
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
  tmpsys2->GetYaxis()->SetTitle("(BNB/MC+EXT)");
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
  tmpsys2->GetYaxis()->SetTitleOffset(1.0);
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
  tmpsys2->SetMinimum(0.);
  tmpsys2->SetMaximum(2.0);
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
  tmpratio1->SetLineColor(kBlue);
  tmpratio1->SetLineWidth(2);
  tmpratio1->Draw("PE same");

  tmpratio2->SetMarkerStyle(20);
  tmpratio2->SetLineColor(kRed);
  tmpratio2->SetLineWidth(2);
  tmpratio2->Draw("PE same");

  pad2->Update();
  pad2->Modified(); // so it updates the pad with the new changes
  pad2->Draw("");
*/
  gPad->RedrawAxis();
  gPad->Update();
 
  can.Print(cName+".pdf","pdf");

}

void plot_one(TMatrixD matrix, std::vector<double> spectrum, std::vector<std::string> channel_names, std::vector<int> num_bins, std::string tag, bool collapsed){
    gStyle->SetPalette(kRedBlue);
    std::string title;
    if(collapsed) title = "Collapsed";
    else title = "Full";
    if(tag.find("correlation") != std::string::npos ) title = title + " Correlation Matrix";
    else if(tag.find("frac") != std::string::npos ) title = title + " Fractional Covariance Matrix";
    else title = title + " Covariance Matrix";

    int num_bins_total = matrix.GetNrows();
    TH2D h2_full(matrix);
    //h2_full.SetName((tag+"_th2d").c_str());
    h2_full.SetTitle(title.c_str());
    h2_full.SetName(title.c_str());
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

    for(int ic = 0; ic < channel_names.size(); ic++){

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
    gStyle->SetOptStat(0);
    if(collapsed) tag = tag+"_collapsed.pdf";
    else tag = tag+"_full.pdf"; 
    c_full->Print(tag.c_str(),"pdf");
}
