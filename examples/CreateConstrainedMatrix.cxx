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
bool verbose=false;

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
  bool fc = false;
  
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
      {"fc",	        no_argument, 		0, 'y'},
      {0,		no_argument, 		0,  0},
    };
  
  while(iarg != -1)
    {
      iarg = getopt_long(argc,argv, "x:t:f:v:sdcnzmay", longopts, &index);
      
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
	case 'y':
	  fc = true;
	  break;
	case '?':
	case 'h':
	  std::cout<<"Allowed arguments:"<<std::endl;
	  std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
	  return 0;
	}
    }
  
  bool doFakedata = false;
  std::string detsyslabel = "";
  std::string mcerrlabel = "";
  std::string zerobinlabel = "";
  if( mcerr ) mcerrlabel = "_with_mc_err"; 
  if( addext ) zerobinlabel = "_with_zerobin_err"; 
  if( detsys ) detsyslabel = "_with_detsys"; 
  if( tag.find("fakedata") != std::string::npos ) doFakedata = true; 

  gStyle->SetPalette(kLightTemperature);
  
  //set up text files to write the systematics table
  //ONLY WRITE TABLES FOR THE COMBINED CHANNEL
  std::vector<std::ofstream> outfiles(9);//, outfile2, outfile3, outfile4, outfile5, outfile6, outfile7, outfile8, outfile9;
  if(combined){
    outfiles[0].open(Form("full_1eNp_lee_systematicstable_beforeconstraint_%s%s%s%s.txt",tag.c_str(),mcerrlabel.c_str(),zerobinlabel.c_str(),detsyslabel.c_str())); 
    outfiles[1].open(Form("full_1eNp_lee_systematicstable_afterconstraint_%s%s%s%s.txt",tag.c_str(),mcerrlabel.c_str(),zerobinlabel.c_str(),detsyslabel.c_str())); 
    outfiles[2].open(Form("full_1eNp_systematicstable_beforeconstraint_%s%s%s%s.txt",tag.c_str(),mcerrlabel.c_str(),zerobinlabel.c_str(),detsyslabel.c_str())); 
    outfiles[3].open(Form("full_1eNp_systematicstable_afterconstraint_%s%s%s%s.txt",tag.c_str(),mcerrlabel.c_str(),zerobinlabel.c_str(),detsyslabel.c_str())); 
    outfiles[4].open(Form("full_1e0p_lee_systematicstable_beforeconstraint_%s%s%s%s.txt",tag.c_str(),mcerrlabel.c_str(),zerobinlabel.c_str(),detsyslabel.c_str())); 
    outfiles[5].open(Form("full_1e0p_lee_systematicstable_afterconstraint_%s%s%s%s.txt",tag.c_str(),mcerrlabel.c_str(),zerobinlabel.c_str(),detsyslabel.c_str())); 
    outfiles[6].open(Form("full_1e0p_systematicstable_beforeconstraint_%s%s%s%s.txt",tag.c_str(),mcerrlabel.c_str(),zerobinlabel.c_str(),detsyslabel.c_str())); 
    outfiles[7].open(Form("full_1e0p_systematicstable_afterconstraint_%s%s%s%s.txt",tag.c_str(),mcerrlabel.c_str(),zerobinlabel.c_str(),detsyslabel.c_str())); 
    outfiles[8].open(Form("full_numu_systematicstable_beforeconstraint_%s%s%s%s.txt",tag.c_str(),mcerrlabel.c_str(),zerobinlabel.c_str(),detsyslabel.c_str())); 
  }
  
  std::vector<double> sig_fullvec, bkg_fullvec, sig_collvec, bkg_collvec, sig_errfullvec, bkg_errfullvec, sig_exterrfullvec, bkg_exterrfullvec, numu_data_vec, numu_vec, delta_numu; 
  
  //Load up the central value spectra we computed in example 1, to act as a signal spectra
  SBNspec sig_spectra("../bin/"+tag+".SBNspec.root",xml);
  //sig_spectra.CalcFullVector();
  sig_fullvec = sig_spectra.full_vector;
  //collapse vector
  sig_spectra.CollapseVector();
  sig_collvec = sig_spectra.collapsed_vector;
  if(verbose){ std::cout << "print the signal collapsed vector: " << sig_spectra.PrintCollapsedVector() << std::endl; }

  //Repeat but this time scale the LEE signal subchannel to 0, to act as our background spectra
  SBNspec bkg_spectra("../bin/"+tag+".SBNspec.root",xml);
  bkg_spectra.Scale("lee", 0.0); // Scale() will scale all that match so bkg_spectra.Scale("leesignal",0.0); would work too 
  bkg_fullvec = bkg_spectra.full_vector;

  //collapse vector
  bkg_spectra.CollapseVector();
  bkg_collvec = bkg_spectra.collapsed_vector;

  //define bins number
  //the xml assumes that we always have numu channel
  //and at least a 1e0p channel or 1eNp channel 
  int num_bins = sig_spectra.num_bins_total;
  int num_bins_coll = sig_spectra.num_bins_total_compressed;
  int np_bins = 0, zp_bins = 0, numu_bins = 0;
  //create vector of channels names
  //assume there is a nue and numu channels
  //check xml if run into issues
  std::vector<std::string> channel_names_coll;

  //get hist based on order of the channels
  std::vector<TH1D*> channel_hists;
  std::vector<TH1D*> channel_hists_after;
  std::vector<TH1D*> channel_hists_after_signalbins;
  std::vector<TH1D*> channel_hists_data_signalbins;
  std::vector<int> channel_bins;
  std::vector<int> channel_signalbins;
  std::vector<int> channel_edges;
  std::vector<int> num_bins_coll_vec;
  int edge = 0;
  channel_edges.push_back(edge);

  for(int ic=0; ic<sig_spectra.num_channels; ic++){
    std::string subchannel_name = Form("%s_%s_%s_%s",sig_spectra.mode_names[0].c_str(),sig_spectra.detector_names[0].c_str(),sig_spectra.channel_names[ic].c_str(),sig_spectra.subchannel_names[ic][0].c_str());
    if(verbose) std::cout << "channel name: " << subchannel_name << std::endl;
    int wb = sig_spectra.GetGlobalBinNumber(1,subchannel_name);
    if(verbose) std::cout << "which global bin: " << wb << std::endl;
    int wh = sig_spectra.GetHistNumber(wb);
    if(verbose) std::cout << "which histo: " << wh << std::endl;
    int bins = sig_spectra.GetIndiciesFromSubchannel(subchannel_name).size();
    if(verbose) std::cout << "bins: " << bins << std::endl;
    channel_bins.push_back(bins);
    edge += bins;
    channel_edges.push_back(edge);
    TH1D* histo1 = (TH1D*)sig_spectra.hist[wh].Clone(Form("%s_before",sig_spectra.channel_names[ic].c_str()));
    TH1D* histo2 = (TH1D*)sig_spectra.hist[wh].Clone(Form("%s_after",sig_spectra.channel_names[ic].c_str()));
    TH1D* histo3 = (TH1D*)sig_spectra.hist[wh].Clone(Form("%s_after_signalbins",sig_spectra.channel_names[ic].c_str()));
    TH1D* histo4 = (TH1D*)sig_spectra.hist[wh].Clone(Form("%s_data_signalbins",sig_spectra.channel_names[ic].c_str()));
    channel_hists.push_back(histo1);
    channel_hists_after.push_back(histo2);
    histo3->Reset();
    channel_hists_after_signalbins.push_back(histo3);
    channel_signalbins.push_back(histo4->FindBin(0.85) - 1);
    if(verbose) std::cout << "#### ic, channel bins = " << channel_signalbins.back() << std::endl;
    histo4->Reset();
    channel_hists_data_signalbins.push_back(histo4);
  }

  //use the vector of channel histograms to define 
  //binnings and name of channels of collapsed spectra
  for(int ic=0; ic<sig_spectra.num_channels; ic++){
    std::string name = channel_hists[ic]->GetName();
    if(verbose) std::cout << "channel name: " << name << std::endl;
    if(name.find("1eNp_sig") != std::string::npos ) np_bins = channel_hists[ic]->GetNbinsX();
    if(name.find("1e0p_sig") != std::string::npos ) zp_bins = channel_hists[ic]->GetNbinsX();
    if(name.find("numu")     != std::string::npos ) numu_bins = channel_hists[ic]->GetNbinsX();
    if(name.find("1eNp_sig") != std::string::npos ) channel_names_coll.push_back("1eNp LEE");
    if(name.find("1eNp_sig") != std::string::npos ) num_bins_coll_vec.push_back(channel_hists[ic]->GetNbinsX());
    if(name.find("1eNp_bg")  != std::string::npos ) channel_names_coll.push_back("1eNp nue + bkg");
    if(name.find("1eNp_bg")  != std::string::npos ) num_bins_coll_vec.push_back(channel_hists[ic]->GetNbinsX());
    if(name.find("1e0p_sig") != std::string::npos ) channel_names_coll.push_back("1e0p LEE");
    if(name.find("1e0p_sig") != std::string::npos ) num_bins_coll_vec.push_back(channel_hists[ic]->GetNbinsX());
    if(name.find("1e0p_bg")  != std::string::npos ) channel_names_coll.push_back("1e0p nue + bkg");
    if(name.find("1e0p_bg")  != std::string::npos ) num_bins_coll_vec.push_back(channel_hists[ic]->GetNbinsX());
    if(name.find("numu")     != std::string::npos ) channel_names_coll.push_back("numu");
    if(name.find("numu")     != std::string::npos ) num_bins_coll_vec.push_back(channel_hists[ic]->GetNbinsX());;
  }
  int nue_bins = num_bins_coll-numu_bins;
  
  //create vector of subchanels/histos names
  std::vector<std::string> subchannel_names;
  std::vector<int> num_bins_full;
  for(auto& h: sig_spectra.hist){
    subchannel_names.push_back(h.GetName());
    num_bins_full.push_back(h.GetNbinsX());
  }
   
  //store the mc stats error and the zero bin error
  sig_spectra.CalcErrorVector();
  sig_errfullvec = sig_spectra.full_error;
  sig_exterrfullvec = sig_spectra.full_ext_err_vector;
  bkg_errfullvec = bkg_spectra.full_error;
  bkg_exterrfullvec = bkg_spectra.full_ext_err_vector;

  //get extbnb to create vector of spectra with no bnb
  //this will be used to calculate the systematics for the table
  std::vector<double> sig_collvec_noext;
  std::vector<double> collapsed_ext;
  TH1D *h_1eNp_bg_ext, *h_1e0p_bg_ext, *h_numu_ext;
  TH1D *h_1eNp_bg_nue, *h_1e0p_bg_nue, *h_numu_bnb;
  for(auto& h: sig_spectra.hist){
    TString hname = h.GetName();
    if(np || combined ){ if(hname=="nu_uBooNE_1eNp_bg_ext") h_1eNp_bg_ext = (TH1D*)h.Clone("nu_uBooNE_1eNp_bg_ext");}
    if(zp || combined ){ if(hname=="nu_uBooNE_1e0p_bg_ext") h_1e0p_bg_ext = (TH1D*)h.Clone("nu_uBooNE_1e0p_bg_ext");}
    if(hname=="nu_uBooNE_numu_ext") h_numu_ext = (TH1D*)h.Clone("nu_uBooNE_1eNp_bg_ext");
  }
  if( !doFakedata && ( np || combined ) ){ for(int bin=1; bin < np_bins+1; bin++) collapsed_ext.push_back(0.0);}
  if( !doFakedata && ( np || combined ) ){ for(int bin=1; bin < np_bins+1; bin++) collapsed_ext.push_back(h_1eNp_bg_ext->GetBinContent(bin));}
  if( !doFakedata && ( zp || combined ) ){ for(int bin=1; bin < zp_bins+1; bin++) collapsed_ext.push_back(0.0);}
  if( !doFakedata && ( zp || combined ) ){ for(int bin=1; bin < zp_bins+1; bin++) collapsed_ext.push_back(h_1e0p_bg_ext->GetBinContent(bin));}
  
  if(!doFakedata){
    for(int bin=1; bin < h_numu_ext->GetNbinsX()+1; bin++) collapsed_ext.push_back(h_numu_ext->GetBinContent(bin));
  }

  //spectra without the ext bins   
  for(int i=0; i < num_bins_coll; i++){
    if(!doFakedata) sig_collvec_noext.push_back( sig_collvec[i] - collapsed_ext[i] );
    else sig_collvec_noext.push_back( sig_collvec[i] );
  }

  //prepare matrices for the detector systematics 
  TMatrixD Mdetsys(num_bins,num_bins);
  TMatrixD Mfracdetsys(num_bins,num_bins);
  TMatrixD Mcorrdetsys(num_bins,num_bins);
  TMatrixD Mdetsys0(num_bins,num_bins);
  Mdetsys.Zero();
  Mdetsys0.Zero();
  SBNchi chi_h1(sig_spectra);	
  SBNchi chi_h0(bkg_spectra);	
  
  //fill in the full covariance detsys based on the values stored in SBNchi::FillDetSysMatrix
  if(detsys){ 
    chi_h1.FillDetSysMatrix(Mdetsys,sig_spectra,false); // !false! means this function will return a full covariance matrix,
    chi_h0.FillDetSysMatrix(Mdetsys0,bkg_spectra,false);// NOT frac covariance matrix!
  }
  
  //calculate frac detsys and correlation matrices
  for(int i=0; i < Mdetsys.GetNrows(); i++){
    for(int j=0; j < Mdetsys.GetNcols(); j++){
      Mfracdetsys(i,j) = 0.0;
      if(sig_fullvec.at(i) !=0 && sig_fullvec.at(j) !=0 ) Mfracdetsys(i,j) = Mdetsys(i,j)/(sig_fullvec.at(i)*sig_fullvec.at(j));
      Mcorrdetsys(i,j) = 0.0;
      if(Mdetsys(i,j) != 0 ) Mcorrdetsys(i,j) = Mdetsys(i,j)/sqrt(Mdetsys(i,i)*Mdetsys(j,j));
    }
  }
  
  //prepare collapsed matrix detsys
  TMatrixD detsyscoll;
  TMatrixD fracdetsyscoll;
  TMatrixD fracdetsyscoll_noext;
  TMatrixD corrdetsyscoll;
  TMatrixD fraccoll;
  TMatrixD fraccoll_noext;
  TMatrixD corrcoll;
  TMatrixD detsyscoll0;
  TMatrixD fraccov;
  TMatrixD corr;
  //resize matrix to num_bins_coll
  detsyscoll.ResizeTo(num_bins_coll, num_bins_coll);
  fracdetsyscoll.ResizeTo(num_bins_coll, num_bins_coll);
  fracdetsyscoll_noext.ResizeTo(num_bins_coll, num_bins_coll);
  corrdetsyscoll.ResizeTo(num_bins_coll, num_bins_coll);
  fraccoll.ResizeTo(num_bins_coll, num_bins_coll);
  fraccoll_noext.ResizeTo(num_bins_coll, num_bins_coll);
  corrcoll.ResizeTo(num_bins_coll, num_bins_coll);
  detsyscoll0.ResizeTo(bkg_spectra.num_bins_total_compressed, bkg_spectra.num_bins_total_compressed);
  fraccov.ResizeTo(num_bins_coll, num_bins_coll);
  corr.ResizeTo(num_bins_coll, num_bins_coll);
  chi_h1.CollapseModes(Mdetsys, detsyscoll);
  chi_h0.CollapseModes(Mdetsys0, detsyscoll0);

  
  //frac detsys
  for(int i=0; i < detsyscoll.GetNrows(); i++){
    for(int j=0; j < detsyscoll.GetNcols(); j++){
      fracdetsyscoll(i,j) = 0.0;
      if(sig_collvec[i] !=0 && sig_collvec[j] !=0 ) fracdetsyscoll(i,j) = detsyscoll(i,j)/(sig_collvec[i]*sig_collvec[j]);
      fracdetsyscoll_noext(i,j) = 0.0;
      if(sig_collvec_noext[i] !=0 && sig_collvec_noext[j] !=0 ) fracdetsyscoll(i,j) = detsyscoll(i,j)/(sig_collvec_noext[i]*sig_collvec_noext[j]);
      corrdetsyscoll(i,j) = 0.0;
      if(detsyscoll(i,j) != 0 ) corrdetsyscoll(i,j) = detsyscoll(i,j)/sqrt(detsyscoll(i,i)*detsyscoll(j,j));
    }
  }
  
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
  
  //print the collapsed matrix plots
  plot_one(Mcorrdetsys, sig_fullvec, subchannel_names, num_bins_full, tag+"_correlation_matrix", false);
  plot_one(Mfracdetsys, sig_fullvec, subchannel_names, num_bins_full, tag+"_fractional_covariance_matrix", false);
  plot_one(corrdetsyscoll, sig_collvec, channel_names_coll, num_bins_coll_vec, tag+"_correlation_matrix", true);
  plot_one(fracdetsyscoll, sig_collvec, channel_names_coll, num_bins_coll_vec, tag+"_fractional_covariance_matrix", true);
  
  //Load up our covariance matricies we calculated in example1 (we could also load up single variation ones)
  TFile * fsys = new TFile(Form("../bin/%s.SBNcovar.root",tag.c_str()),"read");
  TMatrixD *cov = (TMatrixD*)fsys->Get("collapsed_covariance");
  TMatrixD *fraccov0 = (TMatrixD*)fsys->Get("collapsed_frac_covariance");
  TMatrixD *fullfrac = (TMatrixD*)fsys->Get("frac_covariance");
   
  //Recalculate cov matrix by adding mc intrinsic error
  TMatrixD dummy;
  dummy.ResizeTo(num_bins, num_bins);
  TMatrixD collapsed;
  collapsed.ResizeTo(num_bins_coll, num_bins_coll);
  TMatrixD mc_stats_collapsed;
  mc_stats_collapsed.ResizeTo(num_bins_coll, num_bins_coll);
  TMatrixT<double> fullsig = chi_h1.CalcCovarianceMatrix(fullfrac, sig_fullvec, sig_errfullvec, false);  
  TMatrixT<double> mc_stats = chi_h1.CalcCovarianceMatrix(&dummy, sig_fullvec, sig_errfullvec, false);  
  chi_h1.CollapseModes(fullsig, collapsed);
  chi_h1.CollapseModes(mc_stats, mc_stats_collapsed);
  if(mcerr){ 
    tag=tag+"_with_mc_err";
    for(int i=0; i < detsyscoll.GetNrows(); i++){ 
      if(verbose) std::cout << tag << "  Matrix diag mc err: " << sqrt((*cov)(i,i)) << " ==> " << sqrt(collapsed(i,i)) << std::endl; 
      (*cov)(i,i) = collapsed(i,i);
      if(verbose) std::cout << " Matrix diag = " << sqrt((*cov)(i,i)) << std::endl;
    }
  }

  //create the zero bin ext bnb err matrix
  TMatrixD fullext;
  fullext.ResizeTo(num_bins, num_bins);
  if(doFakedata) std::fill(sig_exterrfullvec.begin(), sig_exterrfullvec.end(), 0.);
  for(int i=0; i < num_bins; i++ ) fullext(i,i) += sig_exterrfullvec[i]*sig_exterrfullvec[i];
  
  TMatrixD collext;
  collext.ResizeTo(num_bins_coll, num_bins_coll);
  //collapse...
  chi_h1.CollapseModes(fullext, collext);

  //add zero bins error
  if(addext){
    tag=tag+"_with_zerobin_err";
    //add to full covariance matrix (collapsed)
    for(int i=0; i < detsyscoll.GetNrows(); i++){
      if(verbose) std::cout << tag << " adding zero bin error Matrix diag: " << sqrt(((*cov))(i,i)) << " + " << sqrt(collext(i,i)) << " = "; 
      (*cov)(i,i) += collext(i,i);
      if(verbose) std::cout << sqrt((*cov)(i,i)) << std::endl;
    }
  }

  //add the detsys error
  if(detsys){
    for(int i=0; i < detsyscoll.GetNrows(); i++){
      if(verbose) std::cout << tag << " adding detsys Matrix diag: " << sqrt(((*cov))(i,i)) << " + " << sqrt(detsyscoll(i,i)) << " = ";   
      (*cov)(i,i) += detsyscoll(i,i);
      if(verbose) std::cout << sqrt((*cov)(i,i)) << std::endl;
    }
  }
  
  //maybe store a different errors that will work with current FC code?
  /*for(int i=0; i < detsyscoll.GetNrows(); i++){
    if(verbose) std::cout << "---mc err = " << mc_stats_collapsed(i,i);
    if(fc) mc_stats_collapsed += collext;
    if(verbose) std::cout << "+++mc err = " << mc_stats_collapsed(i,i);
  }*/

  //calculate the full systematics frac cov and correlation
  for(int i=0; i < detsyscoll.GetNrows(); i++){
    for(int j=0; j < detsyscoll.GetNcols(); j++){
      fraccov(i,j) = 0.0;
      if(sig_collvec[i] !=0 && sig_collvec[j] !=0 ) fraccov(i,j) = (*cov)(i,j)/(sig_collvec[i]*sig_collvec[j]);
      corr(i,j) = 0.0;
      if((*cov)(i,j) != 0 ) corr(i,j) = (*cov)(i,j)/sqrt((*cov)(i,i)*(*cov)(j,j));
      if(verbose){ if(i==j) std::cout << tag << "  Matrix diag: " << sqrt((*cov)(i,i)) << ", " << sqrt(fraccov(i,j)) << std::endl;} 
    }
  }
  
  //zero out these if the flags are not turned on
  if(!addext) collext.Zero();
  if(!mcerr) collapsed.Zero();
   
  //calculate frac cov matrix and correlation matrix before constraint
  for(int i=0; i < (*cov).GetNrows(); i++){
    for(int j=0; j < (*cov).GetNcols(); j++){
      fraccoll(i,j) = 0.0;
      if(sig_collvec[i] != 0 && sig_collvec[j] != 0) fraccoll(i,j) = (*cov)(i,j)/(sig_collvec[i]*sig_collvec[j]);
      fraccoll_noext(i,j) = 0.0;
      if(sig_collvec_noext[i] != 0 && sig_collvec_noext[j] != 0) fraccoll_noext(i,j) = (*cov)(i,j)/(sig_collvec_noext[i]*sig_collvec_noext[j]);
      corrcoll(i,j) = 0.0;
      if((*cov)(i,j) != 0 ) corrcoll(i,j) = (*cov)(i,j)/sqrt((*cov)(i,i)*(*cov)(j,j));
    }
  }
  
  //write out the diagonals of frac cov matrix before constraint only do this for the combined channels
  if(combined && !doFakedata){ 
    for(int ic=0; ic<sig_spectra.num_channels; ic++){
      for( int i = 0; i < (*cov).GetNrows(); i++ ){
	for( int j = 0; j < (*cov).GetNcols(); j++ ){
	  if( i==j){
            int is = 2*ic; //0-->0, 1-->2, 2-->4, 3-->6
            if( i >= channel_edges[ic] && i < channel_edges[ic+1] ) outfiles[is]  << sqrt(fraccoll(i,j)) << " " << sqrt(fraccoll_noext(i,j)) << " " << sig_collvec[i] << " " << sig_collvec_noext[i] << std::endl;
	  }// end if i==j 
	}// end loop of matrix rows
      }// end loop of matrix cols
    }// end loop of channels 
  }// end if we are computing combined channels and not fakedata 

  //print the collapsed matrix plot
  plot_one(corrcoll, sig_collvec, channel_names_coll, num_bins_coll_vec, tag+"_correlation_matrix_xsecfluxdetsys", true);
  plot_one(fraccoll, sig_collvec, channel_names_coll, num_bins_coll_vec, tag+"_fractional_covariance_matrix_xsecfluxdetsys", true);
  
  //loop over all channels
  for(int ic=0; ic<sig_spectra.num_channels; ic++){
    channel_hists[ic]->Reset();
    //set bin content
    for(int bin=0; bin < channel_bins[ic]; bin++) channel_hists[ic]->SetBinContent(bin+1,sig_collvec[bin+channel_edges[ic]]);
    //set bin error
    for(int bin=0; bin < channel_bins[ic]; bin++) channel_hists[ic]->SetBinError(bin+1,sqrt((*cov)[bin+channel_edges[ic]][bin+channel_edges[ic]]));
  }
  

  //prepare the tags to identify plots 
  if(fc) tag += "_fc"; 
  if(detsys) tag += "_detsys";
  
  //prepare the canvas and dummy histogram
  //to be passed to the DrawDataMCAndSyst
  TString cname;
  
  //fetch the TFile for the data spectrum
  //to be passed into the DrawDataMCAndSyst
  //and to use the numu data spectrum in the constraint
  TFile *f_data = new TFile(fakedata.c_str(),"read");
  //1eNp data
  TH1D *h_nue_np_fake_data;
  if( doFakedata && (np || combined) ) h_nue_np_fake_data = (TH1D*)f_data->Get("nu_uBooNE_1eNp_data");
  //1e0p data
  TH1D *h_nue_0p_fake_data;
  if( doFakedata && (zp || combined) ) h_nue_0p_fake_data = (TH1D*)f_data->Get("nu_uBooNE_1e0p_data");
  //numu data 
  TH1D *h_numu_data = (TH1D*)f_data->Get("nu_uBooNE_numu_data");

  //prepare the sum histograms (LEE+BG)
  TH1D* h_1eNp_sum_before, *h_1eNp_bg_before; 
  TH1D* h_1e0p_sum_before, *h_1e0p_bg_before; 
  TH1D* h_1eNp_dummy, *h_1e0p_dummy; 
  TH1D* h_numu_before; 
  if(combined){
    for(int ic=0; ic<sig_spectra.num_channels; ic++){
      std::cout << channel_hists[ic]->GetName() << std::endl;
      std::string name = channel_hists[ic]->GetName();
      if(name.find("1eNp_sig") != std::string::npos){ std::cout << "found: " << channel_hists[ic]->GetName() << std::endl; h_1eNp_sum_before = (TH1D*)channel_hists[ic]->Clone("Np_sum_lee");}
      if(name.find("1eNp_bg") != std::string::npos ){ std::cout << "found: " << channel_hists[ic]->GetName() << std::endl; h_1eNp_bg_before  = (TH1D*)channel_hists[ic]->Clone("Np_bg");}
      if(name.find("1eNp_bg") != std::string::npos ){ std::cout << "found: " << channel_hists[ic]->GetName() << std::endl; h_1eNp_dummy = (TH1D*)channel_hists[ic]->Clone("Np_dummy"); h_1eNp_dummy->Reset();}
      if(name.find("1e0p_sig") != std::string::npos){ std::cout << "found: " << channel_hists[ic]->GetName() << std::endl; h_1e0p_sum_before = (TH1D*)channel_hists[ic]->Clone("Zp_sum_lee");}
      if(name.find("1e0p_bg") != std::string::npos ){ std::cout << "found: " << channel_hists[ic]->GetName() << std::endl; h_1e0p_bg_before  = (TH1D*)channel_hists[ic]->Clone("Zp_bg");}
      if(name.find("1e0p_bg") != std::string::npos ){ std::cout << "found: " << channel_hists[ic]->GetName() << std::endl; h_1e0p_dummy = (TH1D*)channel_hists[ic]->Clone("Zp_dummy"); h_1e0p_dummy->Reset();}
      if(name.find("numu") != std::string::npos ){ std::cout << "found: " << channel_hists[ic]->GetName() << std::endl; h_numu_before = (TH1D*)channel_hists[ic]->Clone("numu_before"); }
    }
    
    //1eNp before constraint
    h_1eNp_sum_before->Add(h_1eNp_bg_before);
    //1e0p before constraint
    h_1e0p_sum_before->Add(h_1e0p_bg_before);
    
    //Draw the before constraint plots
    
    cname = tag + "_before_constraint_1eNp";
    if(doFakedata) DrawDataMCAndSyst(cname, h_1eNp_sum_before, h_1eNp_bg_before, h_nue_np_fake_data, "Reconstructed Visible Energy [GeV]", "#nu_{e} 1eNp0#pi Selection");	
    else DrawDataMCAndSyst(cname, h_1eNp_sum_before, h_1eNp_bg_before, h_1eNp_dummy, "Reconstructed Visible Energy [GeV]", "#nu_{e} 1eNp0#pi Selection");	
    cname = tag + "_before_constraint_1e0p";
    if(doFakedata) DrawDataMCAndSyst(cname, h_1e0p_sum_before, h_1e0p_bg_before, h_nue_0p_fake_data, "Reconstructed Visible Energy [GeV]", "#nu_{e} 1e0p0#pi Selection");	
    else DrawDataMCAndSyst(cname, h_1e0p_sum_before, h_1e0p_bg_before, h_1e0p_dummy, "Reconstructed Visible Energy [GeV]", "#nu_{e} 1e0p0#pi Selection");	
    cname = tag + "_before_constraint_numu";
    DrawDataMCAndSyst(cname, h_numu_before, h_numu_before, h_numu_data, "Reconstructed Visible Energy [GeV]", "#nu_{#mu} Selection");	
  }
 
  //======================================================
  // THIS IS THE START OF THE CONSTRAINT PROCEDURE
  //======================================================
 
  //calculate the difference betweeen numu data and mc 
  //to be used in the constraint procedure 
  for(int bin=0; bin < numu_bins; bin++) numu_data_vec.push_back(h_numu_data->GetBinContent(bin+1)); 
  for(int bin=0; bin < numu_bins; bin++ ) numu_vec.push_back(sig_collvec[nue_bins+bin]);
  for(int bin=0; bin < numu_bins; bin++ ) delta_numu.push_back(numu_data_vec[bin]-numu_vec[bin]);

  if(verbose){
    for(int bin=0; bin < numu_bins; bin++ ) std::cout << "numu data, numu mc, delta_numu = " << numu_data_vec[bin] << ", " << numu_vec[bin] << ", " << delta_numu[bin] << std::endl; 
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
  for( int i = 0; i < numumatrix.GetNcols(); i++ ){ 
    numumatrix(i,i) += numu_vec[i];
  }
  
  //Take the numu block of the matrix and invert it
  InvertNumu = numumatrix;
  InvertNumu.Zero();
  TDecompSVD svdnumu( numumatrix );
  if (!svdnumu.Decompose() ) {
    std::cout << "Decomposition failed, matrix singular ?" << std::endl;
  }else{
    InvertNumu = numumatrix.Invert();
    if(verbose) InvertNumu.Print();
  }
  
  //multiply the block matrices to get the matrices used to
  //constraint the systematics: Matrix B
  //constraint the spectrum: Matrix C
  std::cout << "nuenumumatrix.GetNrows(), nuematrix.GetNrows(),nuenumumatrix.GetNrows() = " << nuenumumatrix.GetNrows() << ", " << nuematrix.GetNrows() << ", " << nuenumumatrix.GetNrows() << std::endl; 
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
      if(verbose) std::cout << "binx, biny, delta_numu = " << binx << ", " << biny << ", " << delta_numu[biny] << ", " << C(binx,biny) << std::endl;
    }
    double constnue =  sig_collvec[binx] + scaled;
    double constnue_noext;
    if(!doFakedata) constnue_noext =  sig_collvec[binx] + scaled - collapsed_ext[binx];
    else constnue_noext =  sig_collvec[binx] + scaled;
    if(verbose) std::cout << "check that constrained results make sense: binx, scaled = " << binx << ", " << scaled << std::endl;
    input_nue_constrained.push_back(constnue);
    if(!doFakedata) input_nue_constrained_noext.push_back(constnue_noext);
    else input_nue_constrained_noext.push_back(constnue);
  }
 
  //======================================================
  // THIS IS THE END OF THE CONSTRAINT PROCEDURE
  //======================================================

  //constrained systematics 
  TMatrixD constnuematrix(nuematrix.GetNrows(),nuematrix.GetNcols());
  TMatrixD fracnuematrix(nuematrix.GetNrows(),nuematrix.GetNcols());
  TMatrixD fracnuematrix_noext(nuematrix.GetNrows(),nuematrix.GetNcols());
  TMatrixD corrnuematrix(nuematrix.GetNrows(),nuematrix.GetNcols());
  TMatrixD covnue_signal(nuematrix.GetNrows(),nuematrix.GetNcols());
  TMatrixD fracnue_signal(nuematrix.GetNrows(),nuematrix.GetNcols());
  TMatrixD fracnue_noext_signal(nuematrix.GetNrows(),nuematrix.GetNcols());
  TMatrixD corrnue_signal(nuematrix.GetNrows(),nuematrix.GetNcols());
  for( int i = 0; i < nuematrix.GetNrows(); i++ ){
    for( int j = 0; j < nuematrix.GetNcols(); j++ ){
      constnuematrix(i,j) = nuematrix(i,j) - B(i,j);
      if( input_nue_constrained[i] != 0 && input_nue_constrained[j] != 0 ) fracnuematrix(i,j) = constnuematrix(i,j)/(input_nue_constrained[i]*input_nue_constrained[j]);
      else fracnuematrix(i,j) = 0.0;
      if( input_nue_constrained_noext[i] != 0 && input_nue_constrained_noext[j] != 0) fracnuematrix_noext(i,j) = constnuematrix(i,j)/(input_nue_constrained_noext[i]*input_nue_constrained_noext[j]);
      else fracnuematrix_noext(i,j) = 0.0;
      corrnuematrix(i,j) = constnuematrix(i,j)/sqrt(constnuematrix(i,i)*constnuematrix(j,j));
      if(verbose){ if(i==j) std::cout << i << " fracnuematrix = " << fracnuematrix(i,j) << std::endl; }
      covnue_signal(i,j) = nuematrix(i,j) - B(i,j);
      if( input_nue_constrained[i] != 0 && input_nue_constrained[j] != 0 ) fracnue_signal(i,j) = constnuematrix(i,j)/(input_nue_constrained[i]*input_nue_constrained[j]);
      else fracnue_signal(i,j) = 0.0;
      if( input_nue_constrained_noext[i] != 0 && input_nue_constrained_noext[j] != 0) fracnue_noext_signal(i,j) = constnuematrix(i,j)/(input_nue_constrained_noext[i]*input_nue_constrained_noext[j]);
      else fracnue_noext_signal(i,j) = 0.0;
      corrnue_signal(i,j) = constnuematrix(i,j)/sqrt(constnuematrix(i,i)*constnuematrix(j,j));
      if(channel_bins.size() == 4 && (i >= channel_signalbins[0] && i < channel_bins[0]) || (i >= channel_signalbins[1] && i < channel_bins[1]) || (j >= channel_signalbins[0] && j < channel_bins[0]) || (j >= channel_signalbins[1] && j < channel_bins[1]) ){ fracnue_signal(i,j) = 0.0; covnue_signal(i,j) = 0.0; corrnue_signal = 0.0; fracnue_noext_signal = 0.0; }
      if(channel_bins.size() == 7 && (i >= channel_signalbins[0] && i < channel_bins[0]) || (i >= channel_signalbins[1] && i < channel_bins[1]) || (j >= channel_signalbins[0] && j < channel_bins[0]) || (j >= channel_signalbins[1] && j < channel_bins[1]) ){ fracnue_signal(i,j) = 0.0; covnue_signal(i,j) = 0.0; corrnue_signal = 0.0; fracnue_noext_signal = 0.0;}
      if(channel_bins.size() == 7 && (i >= channel_signalbins[3] && i < channel_bins[3]) || (i >= channel_signalbins[4] && i < channel_bins[4]) || (j >= channel_signalbins[3] && j < channel_bins[3]) || (j >= channel_signalbins[4] && j < channel_bins[4]) ){ fracnue_signal(i,j) = 0.0; covnue_signal(i,j) = 0.0; corrnue_signal = 0.0; fracnue_noext_signal = 0.0;}
    }
  }
  //write out the diagonals of the frac covariance matrix to the tex format
  if(combined && !doFakedata){ 
    for(int ic=0; ic<sig_spectra.num_channels; ic++){
      if(verbose) std::cout << "ic, channel_edges[ic] = " << ic << ", " << channel_edges[ic] << std::endl;
      for( int i = 0; i < nuematrix.GetNrows(); i++ ){
	for( int j = 0; j < nuematrix.GetNcols(); j++ ){
	  if( i==j){
	    int is = 2*ic+1;//0-->1,1-->3,2-->5,3-->7
	    if( i >= channel_edges[ic] && i < channel_edges[ic+1] ) outfiles[is]  << sqrt(fracnuematrix(i,j)) << " " << sqrt(fracnuematrix_noext(i,j)) << " " << input_nue_constrained[i] << " " << input_nue_constrained_noext[i] << std::endl;
	  }// end if i==j 
	}// end loop of matrix rows
      }// end loop of matrix cols
    }// end loop of channels 
  }// end if we are computing combined channels and not fakedata
  
  //loop over all channels
  for(int ic=0; ic<sig_spectra.num_channels-1; ic++){
    //set bin content
    for(int bin=0; bin < channel_bins[ic]; bin++){ 
      channel_hists_after[ic]->SetBinContent(bin+1,input_nue_constrained[bin+channel_edges[ic]]);
      channel_hists_after[ic]->SetBinError(bin+1,sqrt((constnuematrix)[bin+channel_edges[ic]][bin+channel_edges[ic]]));
      if(channel_hists_after[ic]->GetBinContent(bin+1) < 0 ) channel_hists_after[ic]->SetBinContent(bin+1,0.);
    }
    //set bin content
    for(int bin=0; bin < channel_signalbins[0]; bin++){ 
      channel_hists_after_signalbins[ic]->SetBinContent(bin+1,input_nue_constrained[bin+channel_edges[ic]]);
      channel_hists_after_signalbins[ic]->SetBinError(bin+1,sqrt((constnuematrix)[bin+channel_edges[ic]][bin+channel_edges[ic]]));
      if(channel_hists_after_signalbins[ic]->GetBinContent(bin+1) < 0 ) channel_hists_after_signalbins[ic]->SetBinContent(bin+1,0.);
    }
    //set bin error -- note that the new FC code adds the mc stats separately
    //in our case the MC stats is already added to the resulting covariance matrix
    //the fc flag is still in testing mode and need further validation
    if(fc){  for(int bin=0; bin < channel_bins[ic]; bin++) channel_hists_after[ic]->SetBinError(bin+1,sqrt(mc_stats_collapsed[bin+channel_edges[ic]][bin+channel_edges[ic]])); }
    else{ for(int bin=0; bin < channel_bins[ic]; bin++) channel_hists_after[ic]->SetBinError(bin+1,sqrt(constnuematrix[bin+channel_edges[ic]][bin+channel_edges[ic]])); }
  }

  TH1D *h_1eNp_data_signal, *h_1e0p_data_signal;
  //1eNp data
  if( doFakedata && (np || combined) ){ h_1eNp_data_signal = (TH1D*)h_nue_np_fake_data->Clone("data_1eNp");
  h_1eNp_data_signal->Reset();
  for(int bin=0; bin < h_1eNp_data_signal->FindBin(0.85)-1; bin++) h_1eNp_data_signal->SetBinContent(bin+1,h_nue_np_fake_data->GetBinContent(bin+1));
  for(int bin=0; bin < h_1eNp_data_signal->FindBin(0.85)-1; bin++) h_1eNp_data_signal->SetBinError(bin+1,h_nue_np_fake_data->GetBinError(bin+1));}
  //1e0p data
  if( doFakedata && (zp || combined) ){ h_1e0p_data_signal = (TH1D*)h_nue_0p_fake_data->Clone("data_1e0p");
  h_1e0p_data_signal->Reset();
  for(int bin=0; bin < h_1e0p_data_signal->FindBin(0.85)-1; bin++) h_1e0p_data_signal->SetBinContent(bin+1,h_nue_0p_fake_data->GetBinContent(bin+1));
  for(int bin=0; bin < h_1e0p_data_signal->FindBin(0.85)-1; bin++) h_1e0p_data_signal->SetBinError(bin+1,h_nue_0p_fake_data->GetBinError(bin+1)); }
  

  //Draw the after constraint
  TH1D* h_1eNp_sum, *h_1eNp_bg; 
  TH1D* h_1e0p_sum, *h_1e0p_bg; 
  if(combined){
    for(int ic=0; ic<sig_spectra.num_channels; ic++){
      std::cout << channel_hists_after[ic]->GetName() << std::endl;
      std::string name = channel_hists_after[ic]->GetName();
      if(name.find("1eNp_sig_after") != std::string::npos ){ std::cout << "found: " << channel_hists_after[ic]->GetName() << std::endl; h_1eNp_sum = (TH1D*)channel_hists_after[ic]->Clone("Np_sum_lee_after");}
      if(name.find("1eNp_bg_after") != std::string::npos  ){ std::cout << "found: " << channel_hists_after[ic]->GetName() << std::endl; h_1eNp_bg  = (TH1D*)channel_hists_after[ic]->Clone("Np_bg_after");}
      if(name.find("1e0p_sig_after") != std::string::npos ){ std::cout << "found: " << channel_hists_after[ic]->GetName() << std::endl; h_1e0p_sum = (TH1D*)channel_hists_after[ic]->Clone("Zp_sum_lee_after");}
      if(name.find("1e0p_bg_after") != std::string::npos  ){ std::cout << "found: " << channel_hists_after[ic]->GetName() << std::endl; h_1e0p_bg  = (TH1D*)channel_hists_after[ic]->Clone("Zp_bg_after");}
    }
    
    //1eNp after constraint
    h_1eNp_sum->Add(h_1eNp_bg);
    //1e0p after constraint
    h_1e0p_sum->Add(h_1e0p_bg);
    
    //Draw the before constraint plots
    cname = tag + "after_constraint_1eNp";
    if(doFakedata) DrawDataMCAndSyst(cname, h_1eNp_sum, h_1eNp_bg, h_nue_np_fake_data, "Reconstructed Visible Energy [GeV]", "#nu_{e} 1eNp0#pi Selection");	
    else DrawDataMCAndSyst(cname, h_1eNp_sum, h_1eNp_bg, h_1eNp_dummy, "Reconstructed Visible Energy [GeV]", "#nu_{e} 1eNp0#pi Selection");	
    cname = tag + "after_constraint_1e0p";
    if(doFakedata) DrawDataMCAndSyst(cname, h_1e0p_sum, h_1e0p_bg, h_nue_0p_fake_data, "Reconstructed Visible Energy [GeV]", "#nu_{e} 1e0p0#pi Selection");	
    else DrawDataMCAndSyst(cname, h_1e0p_sum, h_1e0p_bg, h_1e0p_dummy, "Reconstructed Visible Energy [GeV]", "#nu_{e} 1e0p0#pi Selection");	
  }
  
  //For the unconstrained matrix, we need to make the number of columns to be consistent with the number of the nue bins
  //Basically we are making a copy of the covariance matrix without the numu part
  TMatrixD covnue(nuematrix.GetNrows(),nuematrix.GetNcols());
  TMatrixD corrnue(nuematrix.GetNrows(),nuematrix.GetNcols());
  TMatrixD fracnue(nuematrix.GetNrows(),nuematrix.GetNcols());
  TMatrixD fracnue_noext(nuematrix.GetNrows(),nuematrix.GetNcols());
  for( int i = 0; i < nuematrix.GetNrows(); i++ ){
    for( int j = 0; j < nuematrix.GetNcols(); j++ ){
      covnue(i,j) = (*cov)(i,j);
      fracnue(i,j) = 0.0;
      fracnue_signal(i,j) = 0.0;
      if(sig_collvec[i] != 0 && sig_collvec[j] != 0) fracnue(i,j) = (*cov)(i,j)/(sig_collvec[i]*sig_collvec[j]);
      fracnue_noext(i,j) = 0.0;
      if(sig_collvec_noext[i] != 0 && sig_collvec_noext[j] != 0) fracnue_noext(i,j) = (*cov)(i,j)/(sig_collvec_noext[i]*sig_collvec_noext[j]);
      corrnue(i,j) = 0.0;
      if((*cov)(i,j) != 0 ) corrnue(i,j) = (*cov)(i,j)/sqrt((*cov)(i,i)*(*cov)(j,j));
      if(verbose){ if(i==j) std::cout << "before const, after const: " << fracnue(i,j) << ", " << fracnuematrix(i,j) << std::endl; }
    }
  }
 
  //write out the output 
  TFile * fspec = new TFile(Form("constrained_%s.SBNspec.root",tag.c_str()),"recreate");
  fspec->cd();
  //loop over all channels
  for(int ic=0; ic<sig_spectra.num_channels-1; ic++){
    std::string name = channel_hists_after[ic]->GetName();
    std::string title = "";
    if(name.find("1eNp_sig_after") != std::string::npos ) title = "nu_uBooNE_1eNp_lee";
    if(name.find("1eNp_bg_after") != std::string::npos ) title = "nu_uBooNE_1eNp_bg";
    if(name.find("1eNp_data") != std::string::npos ) title = "nu_uBooNE_1eNp_data";
    if(name.find("1e0p_sig_after") != std::string::npos ) title = "nu_uBooNE_1e0p_lee";
    if(name.find("1e0p_bg_after") != std::string::npos ) title = "nu_uBooNE_1e0p_bg";
    if(name.find("1e0p_data") != std::string::npos ) title = "nu_uBooNE_1e0p_data";
    std::cout << "write out " <<title<< std::endl;
    channel_hists_after[ic]->SetName(title.c_str());
    channel_hists_after[ic]->SetTitle(title.c_str());
    (TH1D*)channel_hists_after[ic]->Write(title.c_str());
  }
  fspec->Close();
  
  TFile * fspec2 = new TFile(Form("unconstrained_%s.SBNspec.root",tag.c_str()),"recreate");
  fspec2->cd();
  //loop over all channels
  for(int ic=0; ic<sig_spectra.num_channels-1; ic++){
    std::string name = channel_hists[ic]->GetName();
    std::string title = "";
    if(name.find("1eNp_sig") != std::string::npos ) title = "nu_uBooNE_1eNp_lee";
    if(name.find("1eNp_bg") != std::string::npos ) title = "nu_uBooNE_1eNp_bg";
    if(name.find("1eNp_data") != std::string::npos ) title = "nu_uBooNE_1eNp_data";
    if(name.find("1e0p_sig") != std::string::npos ) title = "nu_uBooNE_1e0p_lee";
    if(name.find("1e0p_bg") != std::string::npos ) title = "nu_uBooNE_1e0p_bg";
    if(name.find("1e0p_data") != std::string::npos ) title = "nu_uBooNE_1e0p_data";
    std::cout << "write out " <<title<< std::endl;
    channel_hists[ic]->SetName(title.c_str());
    channel_hists[ic]->SetTitle(title.c_str());
    (TH1D*)channel_hists[ic]->Write(title.c_str());
  }
  fspec2->Close();
  
  //write out the output 
  TFile * fspec3 = new TFile(Form("constrained_%s_signal.SBNspec.root",tag.c_str()),"recreate");
  fspec3->cd();
  //loop over all channels
  for(int ic=0; ic<sig_spectra.num_channels-1; ic++){
    std::string name = channel_hists_after_signalbins[ic]->GetName();
    std::string title = "";
    if(name.find("1eNp_sig_after") != std::string::npos ) title = "nu_uBooNE_1eNp_lee";
    if(name.find("1eNp_bg_after") != std::string::npos ) title = "nu_uBooNE_1eNp_bg";
    if(name.find("1eNp_data") != std::string::npos ) title = "nu_uBooNE_1eNp_data";
    if(name.find("1e0p_sig_after") != std::string::npos ) title = "nu_uBooNE_1e0p_lee";
    if(name.find("1e0p_bg_after") != std::string::npos ) title = "nu_uBooNE_1e0p_bg";
    if(name.find("1e0p_data") != std::string::npos ) title = "nu_uBooNE_1e0p_data";
    std::cout << "write out " <<title<< std::endl;
    channel_hists_after_signalbins[ic]->SetName(title.c_str());
    channel_hists_after_signalbins[ic]->SetTitle(title.c_str());
    (TH1D*)channel_hists_after_signalbins[ic]->Write(title.c_str());
  }
  fspec3->Close();
  
  //write out the output 
  TFile * fspec4 = new TFile(Form("constrained_%s_signal_DATA.SBNspec.root",tag.c_str()),"recreate");
  fspec4->cd();
  if( combined || np ) (TH1D*)h_1eNp_data_signal->Write("nu_uBooNE_1eNp_data");
  if( combined || zp ) (TH1D*)h_1e0p_data_signal->Write("nu_uBooNE_1e0p_data");
  fspec4->Close();
  
  TFile * fcovar = new TFile(Form("constrained_%s.SBNcovar.root",tag.c_str()),"recreate");
  fcovar->cd();
  (TMatrixD*)constnuematrix.Write("full_covariance");
  (TMatrixD*)fracnuematrix.Write("frac_covariance");
  (TMatrixD*)corrnuematrix.Write("full_correlation");
  fcovar->Close();
  
  TFile * fcovar2 = new TFile(Form("unconstrained_%s.SBNcovar.root",tag.c_str()),"recreate");
  fcovar2->cd();
  (TMatrixD*)covnue.Write("full_covariance");
  (TMatrixD*)fracnue.Write("frac_covariance");
  (TMatrixD*)corrnue.Write("full_correlation");
  fcovar2->Close();
  
  TFile * fcovar3 = new TFile(Form("constrained_%s_signal.SBNcovar.root",tag.c_str()),"recreate");
  fcovar3->cd();
  (TMatrixD*)covnue_signal.Write("full_covariance");
  (TMatrixD*)fracnue_signal.Write("frac_covariance");
  (TMatrixD*)corrnue_signal.Write("full_correlation");
  fcovar3->Close();
  
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
  //if(verbose) std::cout << "nbins data, mc = " << tmpData->Integral() << ", " << tmpMC1->Integral() << std::endl;

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
  if(title == "#nu_{#mu} Selection" || std::string(cName).find("fakedata") != std::string::npos ) tmpData->Draw("PE same"); 

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
      if( std::string(cName).find("fakedata") != std::string::npos ) legend->AddEntry(tmpData,"fake data","le");
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
    if( tag.find("frac") != std::string::npos || tag.find("corr") != std::string::npos ) h2_full.SetMaximum(1.0);
    if( tag.find("corr") != std::string::npos ) h2_full.SetMinimum(-1.0);
    if( tag.find("frac") != std::string::npos ) h2_full.SetMinimum(0.0);

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
