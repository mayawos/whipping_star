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
  
  const struct option longopts[] =
    {
      {"xml", 		required_argument, 	0, 'x'},
      {"tag", 		required_argument, 	0, 't'},
      {"fakedata", 		required_argument, 	0, 'f'},
      {"var", 		required_argument, 	0, 'v'},
      {"sys",	no_argument, 0, 's'},
      {"detsys",	no_argument, 0, 'd'},
      {"combined",	no_argument, 0, 'c'},
      {"np",	no_argument, 0, 'n'},
      {"zp",	no_argument, 0, 'z'},
      {0,			no_argument, 		0,  0},
    };
  
  while(iarg != -1)
    {
      iarg = getopt_long(argc,argv, "x:t:f:v:sdcnz", longopts, &index);
      
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
	case '?':
	case 'h':
	  std::cout<<"Allowed arguments:"<<std::endl;
	  std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
	  return 0;
	}
    }
  
  //Load up the central value spectra we computed in example 1, to act as a signal spectra
  SBNspec sig_spectra("../bin/"+tag+".SBNspec.root",xml);
  
  //Repeat but this time scale the LEE signal subchannel to 0, to act as our background spectra
  SBNspec bkg_spectra("../bin/"+tag+".SBNspec.root",xml);
  bkg_spectra.Scale("lee", 0.0); // Scale() will scale all that match so bkg_spectra.Scale("leesignal",0.0); would work too 
 
  //find bin numbers for numu, high nue, 
  int numu_bins =  sig_spectra.GetIndiciesFromSubchannel("nu_uBooNE_numu_bnb").size();
  int nue_bins = sig_spectra.GetIndiciesFromSubchannel("nu_uBooNE_1eNp_bg_intrinsic").size(); 
  int lee_bins = sig_spectra.GetIndiciesFromSubchannel("nu_uBooNE_1eNp_sig_lee").size(); 
  int lownue_bins = 7; 
  int highnue_bins = nue_bins-lownue_bins;
  int numuhighnue_bins = numu_bins+highnue_bins; 

  //fetch the TFile for the data and MC spectrum
  TFile *f_data = new TFile(fakedata.c_str(),"read");
  TH1D *h_numu_data = (TH1D*)f_data->Get("nu_uBooNE_numu_data");
  TH1D *h_1eNp_data = (TH1D*)f_data->Get("nu_uBooNE_1eNp_data");
  SBNspec specdata(fakedata,"../bin/np_numu_reco_e_H1_DATA_mc_collabOct.xml");
  specdata.CalcFullVector();
  std::vector<double> numu_highnue_spec = specdata.full_vector;
  std::vector<double> numu_highnue_data;
  
  for(int i=lownue_bins; i<numu_highnue_spec.size(); i++){ std::cout << "i, numu_highnue_spec[i] = " << i << ", " << numu_highnue_spec[i] << std::endl; numu_highnue_data.push_back(numu_highnue_spec[i]);}

  TMatrixD Mdetsys(sig_spectra.num_bins_total,sig_spectra.num_bins_total);
  TMatrixD Mfracdetsys(sig_spectra.num_bins_total,sig_spectra.num_bins_total);
  TMatrixD Mcorrdetsys(sig_spectra.num_bins_total,sig_spectra.num_bins_total);
  TMatrixD Mdetsys0(sig_spectra.num_bins_total,sig_spectra.num_bins_total);
  Mdetsys.Zero();
  Mdetsys0.Zero();
  SBNchi chi_h1(sig_spectra);	
  SBNchi chi_h0(bkg_spectra);	
  if(detsys){ 
    std::cout << "***FILL DETSYS*****"<< std::endl; 
    chi_h1.FillDetSysMatrix(Mdetsys,sig_spectra);
    chi_h0.FillDetSysMatrix(Mdetsys0,bkg_spectra);
  }
 
  sig_spectra.CalcFullVector();
  //frac detsys
  for(int i=0; i < Mdetsys.GetNrows(); i++){
  for(int j=0; j < Mdetsys.GetNcols(); j++){
    Mfracdetsys(i,j) = 0.0;
    if(sig_spectra.full_vector.at(i) !=0 && sig_spectra.full_vector.at(j) !=0 ) Mfracdetsys(i,j) = Mdetsys(i,j)/(sig_spectra.full_vector.at(i)*sig_spectra.full_vector.at(j));
    Mcorrdetsys(i,j) = 0.0;
    if(Mdetsys(i,j) != 0 ) Mcorrdetsys(i,j) = Mdetsys(i,j)/sqrt(Mdetsys(i,i)*Mdetsys(j,j));
  }
  }
 
  std::vector<double> sig_collvec, bkg_collvec, numu_data_vec, numu_vec, delta_numu; 
  //collapse vector
  sig_spectra.CollapseVector();
  bkg_spectra.CollapseVector();
  sig_collvec = sig_spectra.collapsed_vector;
  bkg_collvec = bkg_spectra.collapsed_vector;

  //collapse matrix detsys
  TMatrixD detsyscoll;
  TMatrixD fracdetsyscoll;
  TMatrixD corrdetsyscoll;
  TMatrixD detsyscoll0;
  detsyscoll.ResizeTo(sig_spectra.num_bins_total_compressed, sig_spectra.num_bins_total_compressed);
  fracdetsyscoll.ResizeTo(sig_spectra.num_bins_total_compressed, sig_spectra.num_bins_total_compressed);
  corrdetsyscoll.ResizeTo(sig_spectra.num_bins_total_compressed, sig_spectra.num_bins_total_compressed);
  detsyscoll0.ResizeTo(bkg_spectra.num_bins_total_compressed, bkg_spectra.num_bins_total_compressed);
  chi_h1.CollapseModes(Mdetsys, detsyscoll);
  chi_h0.CollapseModes(Mdetsys0, detsyscoll0);

  //frac detsys
  for(int i=0; i < detsyscoll.GetNrows(); i++){
  for(int j=0; j < detsyscoll.GetNcols(); j++){
    fracdetsyscoll(i,j) = 0.0;
    if(sig_collvec[i] !=0 && sig_collvec[j] !=0 ) fracdetsyscoll(i,j) = detsyscoll(i,j)/(sig_collvec[i]*sig_collvec[j]);
    corrdetsyscoll(i,j) = 0.0;
    if(detsyscoll(i,j) != 0 ) corrdetsyscoll(i,j) = detsyscoll(i,j)/sqrt(detsyscoll(i,i)*detsyscoll(j,j));
  }
  }

  std::vector<std::string> channel_names;
  std::vector<std::string> channel_names_coll;
  std::vector<int> num_bins;
  std::vector<int> num_bins_coll;

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
  for(int c=0; c < 3; c++){ 
    num_bins_coll.push_back(lee_bins);
    num_bins_coll.push_back(nue_bins);
    num_bins_coll.push_back(numu_bins);
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
  
   
  //Load up our covariance matricies we calculated in example1 (we could also load up single variation ones)
  TFile * fsys = new TFile(Form("../bin/%s.SBNcovar.root",tag.c_str()),"read");
  TMatrixD * cov = (TMatrixD*)fsys->Get("collapsed_covariance");
  TMatrixD *cov0 = (TMatrixD*)cov->Clone("coll_covar");
  *cov = *cov + detsyscoll;

  //get histos to represent the 3 main channels: lee, nue, numu
  TH1D *h_1eNp_lee_before, *h_1eNp_bg_before, *h_1e0p_lee_before, *h_1e0p_bg_before, *h_numu_before;
  for(auto& h: sig_spectra.hist){
    std::string hname = h.GetName();
    if(hname=="nu_uBooNE_1eNp_sig_lee") h_1eNp_lee_before = (TH1D*)h.Clone("nu_uBooNE_1eNp_lee_before");
    if(hname=="nu_uBooNE_1eNp_sig_lee") h_1eNp_bg_before = (TH1D*)h.Clone("nu_uBooNE_1eNp_bg_before");
    if(hname=="nu_uBooNE_1eNp_sig_lee") h_1e0p_lee_before = (TH1D*)h.Clone("nu_uBooNE_1e0p_bg_before");
    if(hname=="nu_uBooNE_1eNp_sig_lee") h_1e0p_bg_before = (TH1D*)h.Clone("nu_uBooNE_1e0p_bg_before");
    if(hname=="nu_uBooNE_numu_bnb") h_numu_before = (TH1D*)h.Clone("nu_uBooNE_numu_before");
  }
  h_1eNp_lee_before->Reset();
  h_1eNp_bg_before->Reset();
  h_1e0p_lee_before->Reset();
  h_1e0p_bg_before->Reset();
  h_numu_before->Reset();

  TString cname;
  
  TH1D* histo_dummy = (TH1D*)h_1eNp_data->Clone("dummy");
  histo_dummy->Reset();
 
  std::vector<double> numuhighnuevec;
  if(combined){     
    //set bin content
    for(int bin=0; bin < 14; bin++ ) h_1eNp_lee_before->SetBinContent(bin+1,sig_collvec[bin]);
    for(int bin=0; bin < 14; bin++ ) h_1eNp_bg_before->SetBinContent(bin+1,sig_collvec[bin+14]);
    for(int bin=0; bin < 14; bin++ ) h_1e0p_lee_before->SetBinContent(bin+1,sig_collvec[bin+28]);
    for(int bin=0; bin < 14; bin++ ) h_1e0p_bg_before->SetBinContent(bin+1,sig_collvec[bin+42]);
    for(int bin=0; bin < 14; bin++ ) h_numu_before->SetBinContent(bin+1,sig_collvec[bin+56]);
    //set bin error
    for(int bin=0; bin < 14; bin++ ) h_1eNp_lee_before->SetBinError(bin+1,sqrt((*cov)[bin][bin]));
    for(int bin=0; bin < 14; bin++ ) h_1eNp_bg_before->SetBinError(bin+1,sqrt((*cov)[bin+14][bin+14]));
    for(int bin=0; bin < 14; bin++ ) h_1e0p_lee_before->SetBinError(bin+1,sqrt((*cov)[bin+28][bin+28]));
    for(int bin=0; bin < 14; bin++ ) h_1e0p_bg_before->SetBinError(bin+1,sqrt((*cov)[bin+42][bin+42]));
    for(int bin=0; bin < 14; bin++ ) h_numu_before->SetBinError(bin+1,sqrt((*cov)[bin+56][bin+56]));
  }else if(np){
    //set bin content
    for(int bin=0; bin < lee_bins; bin++ ) h_1eNp_lee_before->SetBinContent(bin+1,sig_collvec[bin]);
    for(int bin=lownue_bins; bin < nue_bins; bin++ ) numuhighnuevec.push_back(sig_collvec[bin+lee_bins]);
    for(int bin=lownue_bins; bin < nue_bins; bin++ ) h_1eNp_data->SetBinContent(bin+1,numu_highnue_data[bin-lownue_bins]);
    for(int bin=0; bin < nue_bins; bin++ ) h_1eNp_bg_before->SetBinContent(bin+1,sig_collvec[bin+lee_bins]);
    for(int bin=0; bin < numu_bins; bin++ ) h_numu_before->SetBinContent(bin+1,sig_collvec[bin+lee_bins+nue_bins]);
    for(int bin=0; bin < numu_bins; bin++ ) numuhighnuevec.push_back(sig_collvec[bin+lee_bins+nue_bins]);
    //set bin error
    for(int bin=0; bin < lee_bins; bin++ ) h_1eNp_lee_before->SetBinError(bin+1,sqrt((*cov)[bin][bin]));
    for(int bin=0; bin < nue_bins; bin++ ) h_1eNp_bg_before->SetBinError(bin+1,sqrt((*cov)[bin+lee_bins][bin+lee_bins]));
    for(int bin=0; bin < numu_bins; bin++ ) h_numu_before->SetBinError(bin+1,sqrt((*cov)[bin+numu_bins][bin+numu_bins]));
    for(int bin=lownue_bins; bin < nue_bins; bin++ ) h_1eNp_data->SetBinError(bin+1,sqrt(numu_highnue_data[bin-lownue_bins]));
  }else if(zp){
    //set bin content
    for(int bin=0; bin < 14; bin++ ) h_1e0p_lee_before->SetBinContent(bin+1,sig_collvec[bin]);
    for(int bin=0; bin < 14; bin++ ) h_1e0p_bg_before->SetBinContent(bin+1,sig_collvec[bin+14]);
    for(int bin=0; bin < 14; bin++ ) h_numu_before->SetBinContent(bin+1,sig_collvec[bin+28]);
    //set bin error
    for(int bin=0; bin < 14; bin++ ) h_1e0p_lee_before->SetBinError(bin+1,sqrt((*cov)[bin][bin]));
    for(int bin=0; bin < 14; bin++ ) h_1e0p_bg_before->SetBinError(bin+1,sqrt((*cov)[bin+14][bin+14]));
    for(int bin=0; bin < 14; bin++ ) h_numu_before->SetBinError(bin+1,sqrt((*cov)[bin+28][bin+28]));
  }

  if(detsys) tag += "_detsys";
  cname = tag + "before_constraint_numu";
  DrawDataMCAndSyst(cname, h_numu_before, h_numu_before, h_numu_data, "Reconstructed Visible Energy [GeV]", "#nu_{#mu} Selection");

  TH1D* h_1e0p_sum_before;
  TH1D* h_1eNp_sum_before;
  	
  if(combined){

    //numu before constraint
    h_1eNp_sum_before = (TH1D*)h_1eNp_lee_before->Clone("Np_sum_lee_before");
    h_1eNp_sum_before->Add(h_1eNp_bg_before);
    std::cout << "Nbins h_1eNp_sum_before: " << h_1eNp_sum_before->GetNbinsX() << std::endl;
    h_1e0p_sum_before = (TH1D*)h_1e0p_lee_before->Clone("Zp_sum_lee_before");
    h_1e0p_sum_before->Add(h_1e0p_bg_before);
    
    cname = tag + "_before_constraint_1eNp";
    DrawDataMCAndSyst(cname, h_1eNp_sum_before, h_1eNp_bg_before, h_1eNp_data, "Reconstructed Visible Energy [GeV]", "#nu_{e} 1eNp0#pi Selection");	
    cname = tag + "_before_constraint_1e0p";
    DrawDataMCAndSyst(cname, h_1e0p_sum_before, h_1e0p_bg_before, histo_dummy, "Reconstructed Visible Energy [GeV]", "#nu_{e} 1e0p0#pi Selection");	
    
  }

  //numu constraint
  for(int bin=0; bin < numu_bins; bin++) numu_data_vec.push_back(h_numu_data->GetBinContent(bin+1)); 
  for(int bin=0; bin < numu_bins; bin++ ) numu_vec.push_back(sig_collvec[nue_bins+bin]);
  for(int bin=0; bin < numu_bins; bin++ ) delta_numu.push_back(numu_data_vec[bin]-numu_vec[bin]);
  std::vector<double> delta_numu_highnue;
  for(int bin=0; bin < numu_highnue_data.size(); bin++ ){
    //std::cout << "bin, numu_highnue_data[bin], numuhighnuevec[bin] = " << bin << ", " << numu_highnue_data[bin] << ", " << numuhighnuevec[bin] << std::endl; 
    delta_numu_highnue.push_back(numu_highnue_data[bin]-numuhighnuevec[bin]);
  }

  //Block Matrices -- numu+HE sidebands
  lownue_bins += lee_bins;
  //std::cout << "lownue_bins = " << lownue_bins << std::endl;
  //std::cout << "numuhighnue_bins = " << numuhighnue_bins << std::endl;
  TMatrixD numuhighnuematrix(numuhighnue_bins,numuhighnue_bins);
  TMatrixD InvertNumuHighNue(numuhighnue_bins,numuhighnue_bins);
  TMatrixD lownuenumuhighnuematrix(lownue_bins,numuhighnue_bins);
  TMatrixD numuhighnuelownuematrix(numuhighnue_bins,lownue_bins);
  TMatrixD lownuematrix(lownue_bins,lownue_bins);
  
  //std::cout << "(*cov).GetNrows(),(*cov).GetNcols() = " << (*cov).GetNrows() << ", " << (*cov).GetNcols() << std::endl;  
  //std::cout << "(numuhighnuelownuematrix.GetNrows(),numuhighnuelownuematrix.GetNcols() = " << numuhighnuelownuematrix.GetNrows() << ", " << numuhighnuelownuematrix.GetNcols() << std::endl;  
  //std::cout << "(lownuenumuhighnuematrix.GetNrows(),lownuenumuhighnuematrix.GetNcols() = " << lownuenumuhighnuematrix.GetNrows() << ", " << lownuenumuhighnuematrix.GetNcols() << std::endl;  
  //std::cout << "(numuhighnuematrix.GetNrows(),numuhighnuematrix.GetNcols() = " << numuhighnuematrix.GetNrows() << ", " << numuhighnuematrix.GetNcols() << std::endl;  
  for( int i = 0; i < cov->GetNcols(); i++ ){
    for( int j = 0; j < cov->GetNrows(); j++ ){
      if( i >= numuhighnue_bins && j >= numuhighnue_bins) numuhighnuematrix(i-lownue_bins,j-lownue_bins) = (*cov)(i,j);
      if( i >= numuhighnue_bins && j < lownue_bins) numuhighnuelownuematrix(i-lownue_bins,j) = (*cov)(i,j);
      if( i < lownue_bins && j >= numuhighnue_bins) lownuenumuhighnuematrix(i,j-lownue_bins) = (*cov)(i,j);
      if( i < lownue_bins && j < lownue_bins) lownuematrix(i,j) = (*cov)(i,j);
    }
  }
  
  for( int i = 0; i < numuhighnuematrix.GetNcols(); i++ ) numuhighnuematrix(i,i) += numuhighnuevec[i]; 
  
  //constrained procedure
  InvertNumuHighNue = numuhighnuematrix;
  InvertNumuHighNue.Zero();
  TDecompSVD svdnumu( numuhighnuematrix );
  
  if (!svdnumu.Decompose() ) {
    std::cout << "Decomposition failed, matrix singular ?" << std::endl;
  }else{
    InvertNumuHighNue = numuhighnuematrix.Invert();
  }
  
  TMatrixD A(lownuenumuhighnuematrix.GetNrows(),lownuenumuhighnuematrix.GetNcols());
  TMatrixD B(lownuematrix.GetNrows(),lownuematrix.GetNcols());
  TMatrixD C(lownuenumuhighnuematrix.GetNrows(),lownuenumuhighnuematrix.GetNcols());
  
  A.Mult(lownuenumuhighnuematrix,InvertNumuHighNue);
  B.Mult(A,numuhighnuelownuematrix);
  C.Mult(lownuenumuhighnuematrix,InvertNumuHighNue);
  
  //constrained spectrum 
  std::vector<double> input_lownue_constrained;
 
  for( int binx = 0; binx < lownue_bins; binx++ ){
    double scaled = 0.;
    for( int biny = 0; biny < delta_numu_highnue.size(); biny++ ){
      scaled += C(binx,biny) * delta_numu_highnue[biny];
      //std::cout << "binx, biny, delta_numu = " << binx << ", " << biny << ", " << delta_numu[biny] << std::endl;
    }
    double constnue =  sig_collvec[binx] + scaled;
    std::cout << "binx, scaled = " << binx << ", " << scaled << std::endl;
    input_lownue_constrained.push_back(constnue);
  }
  
  TMatrixD constnuematrix(lownuematrix.GetNrows(),lownuematrix.GetNcols());
  TMatrixD fracnuematrix(lownuematrix.GetNrows(),lownuematrix.GetNcols());
  TMatrixD corrnuematrix(lownuematrix.GetNrows(),lownuematrix.GetNcols());
  TMatrixD covnue(lownuematrix.GetNrows(),lownuematrix.GetNcols());
  TMatrixD fraccovnue(lownuematrix.GetNrows(),lownuematrix.GetNcols());
  TMatrixD corrnue(lownuematrix.GetNrows(),lownuematrix.GetNcols());
  for( int i = 0; i < lownuematrix.GetNrows(); i++ ){
    for( int j = 0; j < lownuematrix.GetNcols(); j++ ){
      constnuematrix(i,j) = lownuematrix(i,j) - B(i,j);
      if( input_lownue_constrained[i] != 0 && input_lownue_constrained[j] != 0 ) fracnuematrix(i,j) = constnuematrix(i,j)/(input_lownue_constrained[i]*input_lownue_constrained[j]);
      else fracnuematrix(i,j) = 0.0;
      corrnuematrix(i,j) = constnuematrix(i,j)/sqrt(constnuematrix(i,i)*constnuematrix(j,j));
      if(i==j) std::cout << i << " lownuematrix(i,j), D(i,j), constnuematrix = " << lownuematrix(i,j) << ", " << B(i,j) << ", " << constnuematrix(i,j) << std::endl;
      //if(i==j) std::cout << i << " fracnuematrix = " << fracnuematrix(i,j) << std::endl;
      covnue(i,j) = lownuematrix(i,j) - B(i,j);
      if( sig_collvec[i] != 0 && sig_collvec[j] != 0 ) fraccovnue(i,j) = covnue(i,j)/(sig_collvec[i]*sig_collvec[j]);
      else fraccovnue(i,j) = 0.0;
      corrnue(i,j) = covnue(i,j)/sqrt(covnue(i,i)*covnue(j,j));
    }
  }
  
  
  //Draw after constraint
  TH1D* h_1eNp_lee;
  TH1D* h_1eNp_bg;
  //TH1D* h_1eNp_data;
  TH1D* h_1e0p_lee;
  TH1D* h_1e0p_bg;
  TH1D* h_1e0p_data; 
  
  h_1eNp_lee  = (TH1D*)h_1eNp_bg_before->Clone("nu_uBooNE_1eNp_lee"); 
  h_1eNp_bg   = (TH1D*)h_1eNp_bg_before->Clone("nu_uBooNE_1eNp_bg");
  h_1eNp_data = (TH1D*)h_1eNp_bg_before->Clone("nu_uBooNE_1eNp_data");
  h_1e0p_lee  = (TH1D*)h_1e0p_lee_before->Clone("nu_uBooNE_1e0p_lee");
  h_1e0p_bg   = (TH1D*)h_1e0p_bg_before->Clone("nu_uBooNE_1e0p_bg");
  h_1e0p_data = (TH1D*)h_1e0p_bg_before->Clone("nu_uBooNE_1e0p_data");
  
  //reset the histograms
  h_1eNp_lee->Reset(); 
  h_1eNp_bg->Reset();   
  h_1eNp_data->Reset();
  h_1e0p_lee->Reset();
  h_1e0p_bg->Reset();
  h_1e0p_data->Reset();
 
  lownue_bins = 7; 
  if(combined){
    //set bin content
    for(int bin=0; bin < 14; bin++ ) h_1eNp_lee->SetBinContent(bin+1,input_lownue_constrained[bin+0*numu_bins]);
    for(int bin=0; bin < 14; bin++ ) h_1eNp_bg->SetBinContent(bin+1,input_lownue_constrained[bin+1*numu_bins]);
    for(int bin=0; bin < 14; bin++ ) h_1eNp_data->SetBinContent(bin+1,0.);
    for(int bin=0; bin < 14; bin++ ) h_1e0p_lee->SetBinContent(bin+1,input_lownue_constrained[bin+2*numu_bins]);
    for(int bin=0; bin < 14; bin++ ) h_1e0p_bg->SetBinContent(bin+1,input_lownue_constrained[bin+3*numu_bins]);
    for(int bin=0; bin < 14; bin++ ) h_1e0p_data->SetBinContent(bin+1,0);
    //set bin content error
    for(int bin=0; bin < 14; bin++ ) h_1eNp_lee->SetBinError(bin+1,sqrt(constnuematrix[bin+0*numu_bins][bin+0*numu_bins]));
    for(int bin=0; bin < 14; bin++ ) h_1eNp_bg->SetBinError(bin+1,sqrt(constnuematrix[bin+1*numu_bins][bin+1*numu_bins]));
    for(int bin=0; bin < 14; bin++ ) h_1eNp_data->SetBinError(bin+1,0.);
    for(int bin=0; bin < 14; bin++ ) h_1e0p_lee->SetBinError(bin+1,sqrt(constnuematrix[bin+2*numu_bins][bin+2*numu_bins]));
    for(int bin=0; bin < 14; bin++ ) h_1e0p_bg->SetBinError(bin+1,sqrt(constnuematrix[bin+3*numu_bins][bin+3*numu_bins]));
    for(int bin=0; bin < 14; bin++ ) h_1e0p_data->SetBinError(bin+1,0.);
  }else if(np){
    for(int bin=0; bin < lee_bins; bin++ ) h_1eNp_lee->SetBinContent(bin+1,input_lownue_constrained[bin]);
    for(int bin=0; bin < lownue_bins; bin++ ) h_1eNp_bg->SetBinContent(bin+1,input_lownue_constrained[bin+lee_bins]);
    for(int bin=lownue_bins; bin < nue_bins; bin++ ) h_1eNp_bg->SetBinContent(bin+1,sig_collvec[bin+lee_bins]);
    for(int bin=lownue_bins; bin < nue_bins; bin++ ) h_1eNp_data->SetBinContent(bin+1,numu_highnue_data[bin-lownue_bins]);
    for(int bin=0; bin < lee_bins; bin++ ) h_1eNp_lee->SetBinError(bin+1,sqrt(constnuematrix[bin][bin]));
    for(int bin=0; bin < lownue_bins; bin++ ){ std::cout << "bin+lee_bins = " << bin+lee_bins << std::endl; h_1eNp_bg->SetBinError(bin+1,sqrt(constnuematrix[bin+lee_bins][bin+lee_bins]));}
    for(int bin=lownue_bins; bin < nue_bins; bin++ ){ std::cout << "bin+lee_bins = " << bin+lee_bins << std::endl; h_1eNp_bg->SetBinError(bin+1,sqrt((*cov)[bin+lee_bins][bin+lee_bins]));}
    for(int bin=lownue_bins; bin < nue_bins; bin++ ) h_1eNp_data->SetBinError(bin+1,sqrt(numu_highnue_data[bin-lownue_bins]));
  }else if(zp){
    for(int bin=0; bin < 14; bin++ ) h_1e0p_lee->SetBinContent(bin+1,input_lownue_constrained[bin+0*numu_bins]);
    for(int bin=0; bin < 14; bin++ ) h_1e0p_bg->SetBinContent(bin+1,input_lownue_constrained[bin+1*numu_bins]);
    for(int bin=0; bin < 14; bin++ ) h_1e0p_data->SetBinContent(bin+1,0.);
    for(int bin=0; bin < 14; bin++ ) h_1e0p_lee->SetBinError(bin+1,sqrt(constnuematrix[bin+0*numu_bins][bin+0*numu_bins]));
    for(int bin=0; bin < 14; bin++ ) h_1e0p_bg->SetBinError(bin+1,sqrt(constnuematrix[bin+1*numu_bins][bin+1*numu_bins]));
    for(int bin=0; bin < 14; bin++ ) h_1e0p_data->SetBinError(bin+1,0.);
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
    //nue after constraint
    TH1D* h_1eNp_sum = (TH1D*)h_1eNp_bg->Clone("1eNp_sum_bg");
    h_1eNp_sum->Add(h_1eNp_lee);
    TH1D* h_1e0p_sum = (TH1D*)h_1e0p_bg->Clone("1e0p_sum_bg");
    h_1e0p_sum->Add(h_1e0p_lee);
    
    //Draw after constraint
    cname = tag + "_after_constraint_1e0p";
    DrawDataMCAndSyst(cname, h_1e0p_sum, h_1e0p_bg, histo_dummy, "Reconstructed Visible Energy [GeV]", "#nu_{e} Selection");	
    cname = tag + "_after_constraint_1eNp";
    DrawDataMCAndSyst(cname, h_1eNp_sum, h_1eNp_bg, h_1eNp_data, "Reconstructed Visible Energy [GeV]", "#nu_{e} Selection");	
    
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
  (TMatrixD*)covnue.Write("full_covariance");
  (TMatrixD*)fraccovnue.Write("frac_covariance");
  (TMatrixD*)corrnue.Write("full_correlation");
  fcovar2->Close();
 
  std::cout <<" numuhighnuematrix " << numuhighnuematrix.GetNrows() << ", " << numuhighnuematrix.GetNcols() << std::endl;
  std::cout <<" InvertNumuHighNue " << InvertNumuHighNue.GetNrows() << ", " << InvertNumuHighNue.GetNcols() << std::endl;
  std::cout <<" lownuenumuhighnuematrix " << lownuenumuhighnuematrix.GetNrows() << ", " << lownuenumuhighnuematrix.GetNcols() << std::endl;
  std::cout <<" numuhighnuelownuematrix " << numuhighnuelownuematrix.GetNrows() << ", " << numuhighnuelownuematrix.GetNcols() << std::endl;
  std::cout <<" lownuematrix " << lownuematrix.GetNrows() << ", " << lownuematrix.GetNcols() << std::endl;
 
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
  std::cout << "nbins data, mc lee = " << tmpData->GetNbinsX() << ", " << tmpMC1->GetNbinsX() << std::endl;
  std::cout << "nbins data, mc = " << tmpData->Integral() << ", " << tmpMC1->Integral() << std::endl;

  TH1D* tmpratio2 = (TH1D*)data->Clone("tempRatio2");
  TH1D* tmpMC2_staterr = (TH1D*)MC2->Clone("tempCV_staterr2");
  for( int bin = 1; bin < tmpMC2_staterr->GetNbinsX()+1; bin++){ tmpMC2_staterr->SetBinError(bin,sqrt(tmpMC2_staterr->GetBinContent(bin))); }
  tmpratio2->Reset();
  tmpratio2->Divide(tmpData,tmpMC2_staterr,1.0,1.0,"B");
  std::cout << "nbins data, mc nolee = " << tmpData->GetNbinsX() << ", " << tmpMC2->GetNbinsX() << std::endl;

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
  //TPad *pad1 = new TPad("pad1", "pad1", 0, 0.45, 1, 1.0);
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
  sysErr1->GetYaxis()->SetTitleFont(43);
  sysErr1->GetYaxis()->SetTitle("Events");
  sysErr1->GetXaxis()->SetTitleSize(30);
  sysErr1->GetXaxis()->SetTitleFont(43);
  sysErr1->GetXaxis()->SetTitle("Reconstructed Visible Energy [GeV]");
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

  for( int bin = 1; bin < tmpMC1->GetNbinsX()+1; bin++){ tmpMC1->SetBinError(bin,sqrt(tmpMC1->GetBinContent(bin))); }
  for( int bin = 1; bin < tmpMC2->GetNbinsX()+1; bin++){ tmpMC2->SetBinError(bin,sqrt(tmpMC2->GetBinContent(bin))); }
  sysErr1->Draw("E2");
  sysErr2->Draw("E2 same");
  tmpMC1->Draw("histE same"); 
  tmpMC2->Draw("histE same"); 
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
      legend->AddEntry(tmpData,Form("data: %.2f",data->Integral()),"le");
    }
  } 
  legend->Draw("same");

  //remove ratio
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

  const TAxis *axis = tmpsys->GetXaxis();
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
  std::cout << "tmpratio2, tmpratio1 = " << tmpratio2->GetBinContent(8) << ", " << tmpratio1->GetBinContent(8) << std::endl;

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
