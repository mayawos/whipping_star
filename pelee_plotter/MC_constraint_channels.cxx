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
#include "TLatex.h"
#include "TLegendEntry.h"

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


TH1D* GetTotalError(std::string name, TMatrixD errMatrix, TH1D* histo,TH1D* &ratio);
void DrawMCAndSyst(TString cName,TH1D* tmpMC, TH1D* sys, TH1D* ratio, std::string var, std::string title );
void CollapseSubchannels(TMatrixD & M, TMatrixD & Mc, std::vector<int> num_bins, int num_channels, std::vector<int> num_subchannels);
void plot_one(TMatrixD matrix, TH1D *h_nue, TH1D *h_constrain, std::string tag);
void plot_one_noconstraint(TMatrixD matrix, TH1D *h_nue, std::string tag);
void Draw_Stacked(TObjArray histos, TH1D *data, TPad *pad = 0, bool normalised = false, std::string stacktitle = "", std::string var = "" );
void Draw_Stacked2(TObjArray histos, TH1D *data, TString samplename, TPad *pad = 0, bool normalised = false, std::string stacktitle = "", std::string var = "" );
void Draw_Stacked3(TObjArray histos, TH1D *data, TString samplename, TPad *pad = 0, bool normalised = false, std::string stacktitle = "", std::string var = "" );
void DrawStackedMCAndSyst(TString cName, TH1D* MC, TObjArray hist_nue_array, TH1D* sys, TH1D* ratio, std::string var, std::string title);
void DrawStackedMCAndSyst2(TString cName, TH1D* MC, TObjArray hist_nue_array, TH1D* sys, TH1D* ratio1, TH1D* ratio2, std::string var, std::string title);
void DrawSyst(TString cName, TH1D* MC, TObjArray hist_nue_array, TH1D* sys, TH1D* ratio1, TH1D* ratio2, std::string var, std::string title);
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
	bool sys = true;
	bool detsys = false;
	int mass_start = -1;
	
	const struct option longopts[] =
	  {
	    {"xml", 		required_argument, 	0, 'x'},
	    {"tag", 		required_argument, 	0, 't'},
	    {"width", 		required_argument, 	0, 'w'},
	    {"var", 		required_argument, 	0, 'v'},
	    {"sys",	no_argument, 0, 's'},
	    {"detsys",	no_argument, 0, 'd'},
	    {"part", required_argument,0,'p'},
	    {0,			no_argument, 		0,  0},
	  };
	
	while(iarg != -1)
	  {
	    iarg = getopt_long(argc,argv, "x:t:w:v:d:scp:g", longopts, &index);
	    
	    switch(iarg)
	      {
	      case 'x':
		xml = optarg;
		break;
	      case 't':
		tag = optarg;
		break;
	      case 'w':
		width = atoi(optarg);
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
	      case 'p':
		mass_start = atoi(optarg);
		break;
	      case '?':
	      case 'h':
		//std::cout<<"Allowed arguments:"<<std::endl;
		//std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
		return 0;
	      }
	  }
	
	
        std::string dir="../bin/";

	
	//get the H1
	std::vector<TH1D*> histos_h1;
	TH1D* histo;
	TString fnameh1;
        fnameh1 = "../bin/"+tag+".SBNspec.root";
	//if(tag.find("nue_numu") != std::string::npos ) fnameh1 = "../bin/"+tag+".SBNspec.root"; 
	//if(tag.find("1e0p_numu") != std::string::npos ) fnameh1 = "../bin/"+tag+".SBNspec.root"; 
	TFile *fileh1 = new TFile(fnameh1,"read");
	
	//loop counter(s)
	
	int l=0;
	int m=0;
	TKey *key1;
	TObject *obj1;
	TIter nextkey1(gDirectory->GetListOfKeys());
	while (key1 = (TKey*)nextkey1()) {
	  obj1 = key1->ReadObj();
          std::string name = obj1->GetName();
	  histo = (TH1D*)obj1;
	  if(name.find("nu_uBooNE_nue") != std::string::npos || name.find("nu_uBooNE_1e0p") != std::string::npos ) histos_h1.push_back(histo);
	  //if(name.find("nu_uBooNE_nue") != std::string::npos ) histos_h1.push_back(histo);
          std::cout << histos_h1.back()->GetName() << std::endl;
	}
       	
	//Make TObjArray
	TObjArray hist_nueH1_array;
	for(int h = 0; h < histos_h1.size(); h++ ){
	  histos_h1[h]->SetTitle(histos_h1[h]->GetName());
	  hist_nueH1_array.Add(histos_h1[h]);
	}
	//fileh1->Close();
	
        TString fname = dir+tag+"_BKG_noext.SBNspec.root";
        TFile file(fname,"read");
	
        
        std::vector<TH1D*> histos_nue;
        std::vector<TH1D*> histos_constrain;
	
        //loop counter(s)
        int i=0;
        int j=0;
        TKey *key;
	TObject *obj;
        TIter nextkey(gDirectory->GetListOfKeys());
        while (key = (TKey*)nextkey()) {
          obj = key->ReadObj();
          std::string name = obj->GetName();
          TH1D* h = (TH1D*)obj;
	  if( name.find("nu_uBooNE_nue_") != std::string::npos || name.find("nu_uBooNE_1e0p_") != std::string::npos ) histos_nue.push_back(h);
	  //if( name.find("nu_uBooNE_nue_") != std::string::npos ) histos_nue.push_back(h);
	  else histos_constrain.push_back(h);
	}

        SBNspec bg(fname.Data(),xml);
        //make the text files
        std::ofstream nue_before;
        std::ofstream nue_after; 
        std::ofstream constrain_before;
        nue_before << std::fixed << std::setprecision( 2 );
        nue_after << std::fixed << std::setprecision( 2 );
        constrain_before << std::fixed << std::setprecision( 2 );
        std::string labeldetsys = "";
        if(detsys) labeldetsys = "detsys_";
        nue_before.open(Form("nue_channel_%s_%sbeforeconstraint.txt",tag.c_str(),labeldetsys.c_str()));
        constrain_before.open(Form("constrain_channel_%s_%sbeforeconstraint.txt",tag.c_str(),labeldetsys.c_str()));
        nue_after.open(Form("nue_channel_%s_%safterconstraint.txt",tag.c_str(),labeldetsys.c_str()));
        
        if(tag.find("H1") != std::string::npos ){	
        //fill the spectra to fool SBNfit
        TString fsig = dir+"constrained_"+tag+"_signal.SBNspec.root";
        TFile filesig(fsig,"recreate");
	
        TH1D *h_signal = (TH1D*)histos_nue[0]->Clone("nu_uBooNE_nue_signal");
        std::cout << "adding hist nue " << histos_nue[0]->Integral() << std::endl;
        for(int h = 1; h < histos_nue.size(); h++){ 
	  std::cout << "adding " << histos_nue[h]->Integral() << std::endl;
	  h_signal->Add(histos_nue[h]);
	  std::cout << "summed " << h_signal->Integral() << std::endl;
        }
        h_signal->Write();
        filesig.Close();
	
        
        //fill the spectra to fool SBNfit
        TString fbkg = dir+"constrained_"+tag+"_BKG.SBNspec.root";
        TFile filebkg(fbkg,"recreate");
	
        TH1D *h_bkg = (TH1D*)histos_nue[0]->Clone("nu_uBooNE_nue_signal");
        std::cout << "adding hist nue " << histos_nue[0]->Integral() << std::endl;
        for(int h = 1; h < histos_nue.size(); h++){ 
	  std::cout << "bkg adding " << histos_nue[h]->Integral() << std::endl;
          std::string name = (std::string)histos_nue[h]->GetName();
          if(name.find("nu_uBooNE_nue_lee") == std::string::npos) h_bkg->Add(histos_nue[h]);
	  std::cout << "summed " << h_bkg->Integral() << std::endl;
        }
	
        h_bkg->Write();
        filebkg.Close();
        }

        //populate vector with the MC spectrum of the nue and constrain channels          
        std::vector<double> input_vec;
	for( int h = 0; h < histos_nue.size(); h++ ){
          std::cout << "h, h_nue = " << h << ", " << histos_nue[h]->Integral() << std::endl;
          for(int b = 1; b < histos_nue[h]->GetNbinsX()+1; b++){ 
	    std::cout << "b, h_nue = " << b << ", " << histos_nue[h]->GetBinContent(b) << std::endl;
	    input_vec.push_back(histos_nue[h]->GetBinContent(b));
            
          }
	}
        std::cout << "constraint size = " << histos_constrain.size() << std::endl;
	for( int h = 0; h < histos_constrain.size(); h++ ){
          std::cout << "h_constrain = " << histos_constrain[h]->Integral() << std::endl;
          for(int b = 1; b < histos_constrain[h]->GetNbinsX()+1; b++) input_vec.push_back(histos_constrain[h]->GetBinContent(b));
	}
        std::cout << "a" << std::endl;
        //also create histogram summing up the events in each channels
        TH1D *h_nue_total = (TH1D*)histos_nue[0]->Clone("nu_uBooNE_nue"); //assume that the nue N_data == nue N_mc     
        std::cout << "b" << std::endl;
	h_nue_total->Reset();
        std::cout << "c" << std::endl;
	for( int h = 0; h < histos_nue.size(); h++ ){
	  h_nue_total->Add(histos_nue[h]);
	}
        std::cout << "d" << std::endl;
 	h_nue_total->Sumw2(kFALSE);
        std::cout << "e" << std::endl;
        //also create histogram summing up the events in each channels
        std::cout << "histos h1 size = " << histos_h1.size() << std::endl;
        TH1D *h_nueH1_total = (TH1D*)histos_h1[0]->Clone("nu_uBooNE_nue_H1"); //assume that the nue N_data == nue N_mc     
        std::cout << "f" << std::endl;
	h_nueH1_total->Reset();
        std::cout << "try here! " <<  std::endl;
	for( int h = 0; h < histos_h1.size(); h++ ){
	  h_nueH1_total->Add(histos_h1[h]);
	}
 	h_nueH1_total->Sumw2(kFALSE);
	//constrain
        TH1D *h_constrain_total = (TH1D*)histos_constrain[0]->Clone("nu_uBooNE_constrain"); //assume that the nue N_data == nue N_mc     
        h_constrain_total->Reset();
        for( int h = 0; h < histos_constrain.size(); h++ ){
          h_constrain_total->Add(histos_constrain[h]);
        }
        h_constrain_total->Sumw2(kFALSE);
	
        std::vector<double> input_vec_tot;
        for(int b = 1; b < h_nue_total->GetNbinsX()+1; b++) input_vec_tot.push_back(h_nue_total->GetBinContent(b));
        for(int b = 1; b < h_constrain_total->GetNbinsX()+1; b++) input_vec_tot.push_back(h_constrain_total->GetBinContent(b));
	
        std::cout << "size of input_vec, input_vec_total = " << input_vec.size() << ", " << input_vec_tot.size() << std::endl;
	
	//Load up our covariance matricies we calculated in example1 (we could also load up single variation ones)
	TString filename = "../bin/"+tag+"_nolee.SBNcovar.root";
        std::cout << "covar_filename = " << filename << std::endl;
	TFile * fsys = new TFile(filename,"read");
	TString mname = "full_covariance";
	TMatrixD *cov = (TMatrixD*)fsys->Get(mname);
	mname = "frac_covariance";
	TMatrixD *fraccov = (TMatrixD*)fsys->Get(mname);
	mname = "full_correlation";
	TMatrixD *corr = (TMatrixD*)fsys->Get(mname);
	
        SBNchi chi_h0(bg,fraccov);	
        TFile * fconstr = new TFile(Form("constrained_%s.SBNcovar.root",tag.c_str()),"recreate");
	//Now we have all the necessary files, we are going to start the constraint
	//First, add the stat errors from the MC spectrum!
	
	//Create covariance matrix to store the stats error
	TMatrixD Mstat(cov->GetNcols(), cov->GetNrows());
	Mstat.Zero();
	for( int i=0; i < Mstat.GetNcols(); i++ ){
          Mstat(i,i) = input_vec[i];
        } 
        std::cout << "C***" << std::endl;
	//Check that matrix is symmetrix and then add them up
	//And then define the total covariance matrix in all its glory
	TFile *fsystdetsys = new TFile(Form("../bin/%s_detsys.SBNcov.root", tag.c_str()),"recreate");

	TMatrixD Mtotal(cov->GetNcols(), cov->GetNrows());
	TMatrixD Mtotalfraccov(fraccov->GetNcols(), fraccov->GetNrows());
	TMatrixD Mtotalcorr(corr->GetNcols(), corr->GetNrows());

	//Mstat.Print();
	//Mtotal.Print();
        int collsize = input_vec_tot.size();	
	//collapse the covariance matrix to only its channels
	TMatrixD Mcol(collsize,collsize);
	TMatrixD Mcolfrac(collsize,collsize);
	TMatrixD Mcolcorr(collsize,collsize);
        TMatrixD Mdetsys(cov->GetNcols(), cov->GetNrows());
        TMatrixD Mpoisstat(cov->GetNcols(), cov->GetNrows());
        Mpoisstat.Zero();
        Mdetsys.Zero();

        bool useBDT = true;
        if(tag.find("BDT") != std::string::npos ) useBDT = true;
        if(detsys){ std::cout << "***FILL DETSYS*****"<< std::endl; 
        chi_h0.FillDetSysMatrix(Mdetsys,bg);}
	Mtotal.Zero();
	Mtotalfraccov.Zero();
	Mtotalcorr.Zero();
	
	//Mtotalfraccov = *fraccov;
        sys=true;	
        //chi_h0.ZeroOutLEEMatrix(*cov,bg);	
	//if(sys){
	  std::cout<<"Using sys only in covariance matrix"<<std::endl;
	  Mtotalfraccov = *fraccov + Mdetsys;
	/*}else{
	  std::cout<<" Using stats+sys in covariance matrix"<<std::endl;
	  Mtotal = Mstat + Mpoisstat + *cov + Mdetsys;
	}*/
        for(int i=0; i<Mtotal.GetNcols(); i++){
          for(int j=0; j<Mtotal.GetNcols(); j++){
          //if(i==j && input_vec[i] != 0 ) std::cout << "Mdetsys, diagonal = " << Mdetsys(i,j)/(input_vec[i]*input_vec[i]) << std::endl;
          //if(i==j && input_vec[i] != 0 ) std::cout << "cov, diagonal = " << Mcov(i,j)/(input_vec[i]*input_vec[i]) << std::endl;
          //if(i==j && input_vec[i] != 0 ){
          Mtotal(i,j) = Mtotalfraccov(i,j)*(input_vec[i]*input_vec[j]);
          if( i==j && i%14 == 0 )std::cout << "i, Mtotal, Mtotalfraccov, input_vec, error  = " << i << ", " << Mtotal(i,j) << ", " << Mtotalfraccov(i,j) << ", " << input_vec[i]*input_vec[j] << std::endl;
          //}
        }
        }
	//Mstat.Print();
	//Mdetsys.Print();
        	
	//collapse the covariance matrix to only its channels
        //int collsize = input_vec_tot.size();	
	//TMatrixD Mcol(collsize,collsize);
	//TMatrixD Mcolfrac(collsize,collsize);
	//TMatrixD Mcolcorr(collsize,collsize);
	
	Mcol.Zero();
	Mcolfrac.Zero();
	Mcolcorr.Zero();
	
	int num_channels = 2;
        std::cout << "size: " << histos_nue.size() << ", " << histos_constrain.size() << std::endl;
	std::vector<int> num_subchannels = {int(histos_nue.size()),int(histos_constrain.size())};
	//std::vector<int> num_subchannels = {int(histos_nue.size())};
	std::vector<int> num_bins = {h_nue_total->GetNbinsX(),h_constrain_total->GetNbinsX()};
	//std::vector<int> num_bins = {h_nue_total->GetNbinsX()};
        std::cout << "numbins = " << num_bins[0] << ", " << num_bins[1] << std::endl;	
	CollapseSubchannels(Mtotal, Mcol, num_bins, num_channels, num_subchannels);
	
        double nue_err_tot = 0.; 
        double nue_bin_tot = 0.; 
        double constrain_err_tot = 0.; 
        double constrain_bin_tot = 0.;
 
        for(int i=0; i<collsize; i++){
	  for(int j=0; j<collsize; j++){
            if(i==j) std::cout << "cov diagonal = " << Mcol(i,j) << std::endl;
	    Mcolfrac(i,j) = Mcol(i,j)/(input_vec_tot[i]*input_vec_tot[j]);
            if(i==j) std::cout << "error diagonal = " << Mcolfrac(i,j) << std::endl;
	    if( (i==j) && i < (h_nue_total->GetNbinsX()) ){ 
              nue_before << i+1 << " \t " << 100.*sqrt(Mcolfrac(i,j)) << " \t " << input_vec_tot[i] << "\n"; 
	      nue_err_tot += Mcolfrac(i,j); 
	      nue_bin_tot += input_vec_tot[i]*input_vec_tot[j]; 
            }
	    if( (i==j) && i >= (h_nue_total->GetNbinsX()) ){ 
              constrain_before << i-h_nue_total->GetNbinsX()+1 << " \t " << 100.*sqrt(Mcolfrac(i,j)) <<  " \t " << input_vec_tot[i] <<  "\n"; 
	      constrain_err_tot += Mcolfrac(i,j); 
	      constrain_bin_tot += input_vec_tot[i]*input_vec_tot[j];
            } 
	    Mcolcorr(i,j) = Mcol(i,j)/(sqrt(Mcol(i,i))*sqrt(Mcol(j,j)));
	  }
        }
        //nue_before << sqrt(nue_bin_tot) << " \t" << sqrt(nue_err_tot) << "\n";
        //constrain_before << sqrt(constrain_bin_tot) << " \t" << sqrt(constrain_err_tot) << "\n";
        //nue_before << "total \t" << sqrt(nue_err_tot/nue_bin_tot) << "\n";
        //constrain_before << "total \t" << sqrt(constrain_err_tot/constrain_bin_tot) << "\n";

        if(detsys){	
        plot_one(Mcol, h_nue_total, h_constrain_total, "SBNfit_covariance_matrix_"+tag+"_detsys");
        plot_one(Mcolfrac, h_nue_total, h_constrain_total, "SBNfit_fractional_covariance_matrix_"+tag+"_detsys");
        plot_one(Mcolcorr, h_nue_total, h_constrain_total, "SBNfit_correlation_matrix_"+tag+"_detsys");
	}else{
        plot_one(Mcol, h_nue_total, h_constrain_total, "SBNfit_covariance_matrix_"+tag);
        plot_one(Mcolfrac, h_nue_total, h_constrain_total, "SBNfit_fractional_covariance_matrix_"+tag);
        plot_one(Mcolcorr, h_nue_total, h_constrain_total, "SBNfit_correlation_matrix_"+tag);
        }
        TFile * f_nue = new TFile(Form("%s_signal.SBNcovar.root",tag.c_str()),"recreate");
        //now do nue only
        TMatrixD fullcov_nue(h_nue_total->GetNbinsX(),h_nue_total->GetNbinsX()); 
        TMatrixD fraccov_nue(h_nue_total->GetNbinsX(),h_nue_total->GetNbinsX()); 
        TMatrixD fullcorr_nue(h_nue_total->GetNbinsX(),h_nue_total->GetNbinsX()); 
        std::cout << "bin \t frac uncertainties (after)" << std::endl;
        for( int i = 0; i < fullcov_nue.GetNcols(); i++ ){
          for( int j = 0; j < fullcov_nue.GetNrows(); j++ ){
	    fullcov_nue(i,j) = Mcol(i,j);
	    fraccov_nue(i,j) = Mcolfrac(i,j);
	    fullcorr_nue(i,j) = Mcolcorr(i,j);
          }
        }
        //nue_after << "total error \t" << nue_after_err_tot*nue_after_bin_tot << "\n";
        
        (TMatrixD*)fullcov_nue.Write("full_covariance");
        (TMatrixD*)fraccov_nue.Write("frac_covariance");
        (TMatrixD*)fullcorr_nue.Write("full_correlation");
        f_nue->Close();
        //Mstat.Print();
        //Mtotal.Print();
	
	//Make TObjArray
	TObjArray hist_nue_array;
	for(int h = 0; h < histos_nue.size(); h++ ){
	  histos_nue[h]->SetTitle(histos_nue[h]->GetName());
	  hist_nue_array.Add(histos_nue[h]);
	}
	
        //Get error matrix as TH1D and overlay the central value
        TH1D *ratio = (TH1D*)h_nue_total->Clone("rat_before_constraint");
        TH1D *err_before = GetTotalError("before_constraint", Mcol, h_nue_total, ratio);
        TString cname = tag+"_before_constraint";
        //DrawMCAndSyst(cname,h_nue_total,err_before,ratio,var,"#nu_{e} Selection");
        DrawStackedMCAndSyst(cname,h_nueH1_total,hist_nueH1_array,err_before,ratio,var,"#nu_{e} Selection");
	
        //invert the covariance matrix 
        TMatrixD InvertedMatrix = Mcol;
        InvertedMatrix.Zero();
        Mcol.SetTol(1.e-25);
        TDecompSVD svd( Mcol );
        if (!svd.Decompose()) {
          std::cout << "Decomposition failed, matrix singular ?" << std::endl;
          //InvertedMatrix = Mcol.Invert();
        }else{
          std::cout << "Decomposition not failed" << std::endl; 
          InvertedMatrix = Mcol.Invert();
        }
        //InvertedMatrix.Print();
        //Adding the pull term contributions
        for( int i = 0; i < InvertedMatrix.GetNcols(); i++ ){
          for( int j = 0; j < InvertedMatrix.GetNrows(); j++ ){
            if(i==j)std::cout << "BEFORE -- Inverted Matrix i, j = " << i << ", " << j << " = " << (InvertedMatrix)(i,j) << std::endl;
            if( (i==j) && i > (h_nue_total->GetNbinsX()-1) ){ std::cout << "i = " << i << ", " << (1/input_vec_tot[i]) << std::endl; (InvertedMatrix)(i,j) += (1/input_vec_tot[i]);}
            if(i==j)std::cout << "Mcol i, j = " << i << ", " << j << " = " << (Mcol)(i,j) << std::endl;
            if(i==j)std::cout << "AFTER -- Inverted Matrix i, j = " << i << ", " << j << " = " << (InvertedMatrix)(i,j) << std::endl;
          }
        }
        //InvertedMatrix.Print();
        //Invert back the matrix
        TMatrixD ExpandedMatrix = InvertedMatrix;
        ExpandedMatrix.Zero();
        InvertedMatrix.SetTol(1.e-25);
        TDecompSVD svd2( InvertedMatrix );
        if (!svd2.Decompose()) {
          //std::cout << "Decomposition failed, matrix singular ?" << std::endl;
        }else{
          ExpandedMatrix = InvertedMatrix.Invert();
        }
        //ExpandedMatrix.Print();
	
        //create root file containing the constrained matrix
        (TMatrixD*)ExpandedMatrix.Write("full_covariance");
	
        TMatrixD FracExpandedMatrix(ExpandedMatrix.GetNcols(),ExpandedMatrix.GetNrows()); 
        TMatrixD CorrExpandedMatrix(ExpandedMatrix.GetNcols(),ExpandedMatrix.GetNrows()); 
        for( int i = 0; i < ExpandedMatrix.GetNcols(); i++ ){
          for( int j = 0; j < ExpandedMatrix.GetNrows(); j++ ){
            if(i==j)std::cout << "ExpandedMatrix " << i << ", " << j << " = " << ExpandedMatrix(i,j) << std::endl; 
            FracExpandedMatrix(i,j) = ExpandedMatrix(i,j)/(input_vec_tot[i]*input_vec_tot[j]);
            if(i==j)std::cout << "FracExpandedMatrix " << i << ", " << j << " = " << FracExpandedMatrix(i,j) << std::endl;
            CorrExpandedMatrix(i,j) = ExpandedMatrix(i,j)/(sqrt(ExpandedMatrix(i,i))*sqrt(ExpandedMatrix(j,j)));
          }
        }
        (TMatrixD*)FracExpandedMatrix.Write("frac_covariance");
        (TMatrixD*)CorrExpandedMatrix.Write("full_correlation");
	
        plot_one(ExpandedMatrix, h_nue_total, h_constrain_total, "SBNfit_constrained_covariance_matrix_"+tag+"_mc");
        plot_one(FracExpandedMatrix, h_nue_total, h_constrain_total, "SBNfit_constrained_fractional_covariance_matrix_"+tag+"_mc");
        plot_one(CorrExpandedMatrix, h_nue_total, h_constrain_total, "SBNfit_constrained_correlation_matrix_"+tag+"_mc");

        fconstr->Close();
	
        double nue_after_err_tot = 0.; 
        double nue_after_bin_tot = 0.; 
        TFile * fconstr_nue = new TFile(Form("constrained_%s_signal.SBNcovar.root",tag.c_str()),"recreate");
        //now do nue only
        TMatrixD const_fullcov_nue(h_nue_total->GetNbinsX(),h_nue_total->GetNbinsX()); 
        TMatrixD const_fraccov_nue(h_nue_total->GetNbinsX(),h_nue_total->GetNbinsX()); 
        TMatrixD const_fullcorr_nue(h_nue_total->GetNbinsX(),h_nue_total->GetNbinsX()); 
        std::cout << "bin \t frac uncertainties (after)" << std::endl;
        for( int i = 0; i < const_fullcov_nue.GetNcols(); i++ ){
          for( int j = 0; j < const_fullcov_nue.GetNrows(); j++ ){
	    const_fullcov_nue(i,j) = ExpandedMatrix(i,j);
	    const_fraccov_nue(i,j) = FracExpandedMatrix(i,j);
	    if(i==j){ 
              nue_after << i+1 << " \t " << 100.*sqrt(FracExpandedMatrix(i,j)) << " \t " << input_vec_tot[i] <<  "\n";
	      nue_after_err_tot += sqrt(FracExpandedMatrix(i,j))*input_vec_tot[i]; 
	      nue_after_bin_tot += input_vec_tot[i];
            } 
	    const_fullcorr_nue(i,j) = CorrExpandedMatrix(i,j);
          }
        }
        //nue_after << "total error \t" << nue_after_err_tot*nue_after_bin_tot << "\n";
        
        (TMatrixD*)const_fullcov_nue.Write("frac_covariance");
        (TMatrixD*)const_fraccov_nue.Write("frac_covariance");
        (TMatrixD*)const_fullcorr_nue.Write("full_correlation");
        fconstr_nue->Close();
	
        double sum = 0.0;
        TH1D *ratio2 = (TH1D*)h_nue_total->Clone("rat_after_constraint");
        TH1D *err_after = GetTotalError("after_constraint", ExpandedMatrix, h_nue_total, ratio2);
		for( int bin = 1; bin < ratio2->GetNbinsX()+1; bin++ ){
		  std::cout << "bin " << bin << "   ratio = " << ratio->GetBinError(bin) << "     ratio2 = " << ratio2->GetBinError(bin) << "  ratio2/ratio = " << (ratio->GetBinError(bin)-ratio2->GetBinError(bin))/ratio->GetBinError(bin) << std::endl; sum += (ratio->GetBinError(bin)-ratio2->GetBinError(bin))/ratio->GetBinError(bin);
		}
		//std::cout << "average reduction = " << sum/ratio2->GetNbinsX() << std::endl;
		cname = tag+"_after_constraint";
		//DrawMCAndSyst(cname,h_nue_total,err_after,ratio2,var,"#nu_{e} Selection");
		DrawStackedMCAndSyst(cname,h_nueH1_total,hist_nueH1_array,err_after,ratio2,var,"#nu_{e} Selection");
		cname = tag+"_beforeafter_constraint";
		DrawStackedMCAndSyst2(cname,h_nueH1_total,hist_nueH1_array,err_after,ratio,ratio2,var,"#nu_{e} Selection");
		DrawSyst(cname,h_nueH1_total,hist_nueH1_array,err_after,ratio,ratio2,var,"#nu_{e} Selection");
		TString cName2 = tag+"_stacked";
		TCanvas can2(cName2,cName2,2400,1800);
		Draw_Stacked3(hist_nueH1_array, h_nueH1_total, cName2, 0, false, "#nu_{e} Inclusive Selection", var );
		//gPad->Modified(); // so it updates the pad with the new changes
		can2.Print(cName2+".pdf","pdf");

		nue_before.close();
		constrain_before.close();
		nue_after.close();

		return 0;
		
		
	}

	TH1D* GetTotalError(std::string name, TMatrixD errMatrix, TH1D *histo, TH1D* &ratio ){
	  // Make a copy of this histogram as a TH1D and rename it
	  TString hname = name+"_total_error";
	  TH1D *err = (TH1D*)histo->Clone(hname);
	  err->Reset();
	  ratio->Reset();

	  const int highBin = err->GetNbinsX() + 1;
	  const int lowBin = 0;

	  for( int iBin = lowBin; iBin <= highBin; ++iBin ){
	    double derr = errMatrix[iBin][iBin];
	    err->SetBinContent( iBin+1, ( derr > 0 ) ? sqrt(derr): 0. );
	    //std::cout << "cv, err = " << histo->GetBinContent(iBin+1) << ", " << err->GetBinContent(iBin+1) << ", " << err->GetBinContent(iBin+1)/histo->GetBinContent(iBin+1) << std::endl;
	    ratio->SetBinContent(iBin+1, 0);
	    ratio->SetBinError(iBin+1, err->GetBinContent(iBin+1)/histo->GetBinContent(iBin+1));
	  }

	  return err;
	}

	void DrawMCAndSyst(TString cName, TH1D* MC, TH1D* sys, TH1D* ratio, std::string var, std::string title){

	  TCanvas can(cName,cName,2400,1600);
	  TH1D* tmpMC = (TH1D*)MC->Clone("tempCV");
	  TH1D* tmpratio = (TH1D*)ratio->Clone("tempRatio");
	  TH1D* sysErr = (TH1D*)sys->Clone("totalError");
	  sysErr->Reset();
	  //sysErr->Scale(1,"width");
	  for( int bin = 1; bin < tmpMC->GetNbinsX()+1; bin++){
	    sysErr->SetBinContent(bin,tmpMC->GetBinContent(bin));
	    sysErr->SetBinError(bin,sys->GetBinContent(bin));
	    tmpratio->SetBinContent(bin,1);
	    tmpratio->SetBinError(bin,ratio->GetBinError(bin));
	    //std::cout << "bin, cv, error = " << bin << ", " << tmpratio->GetBinContent(bin) << ", " << tmpratio->GetBinError(bin) << std::endl;
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
	  sysErr->GetXaxis()->SetTitleSize(15);
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

	  tmpratio->SetTitle("");
	  tmpratio->GetYaxis()->SetTitle("uncertainties");
	  tmpratio->GetXaxis()->SetTitle(var.c_str());
	  tmpratio->SetStats(0);
	  tmpratio->SetLineColor(kRed);
	  tmpratio->SetLineWidth(1);
	  tmpratio->SetFillColor(kRed-10);
	  tmpratio->SetFillStyle(1001);
	  tmpratio->GetYaxis()->SetTitleSize(25);
	  tmpratio->GetYaxis()->SetTitleFont(43);
	  tmpratio->GetXaxis()->SetTitleSize(25);
	  tmpratio->GetXaxis()->SetTitleOffset(1.0);
	  tmpratio->GetYaxis()->SetTitleOffset(1.0);
	  tmpratio->GetXaxis()->SetLabelSize(0.1);
	  tmpratio->GetYaxis()->SetLabelSize(0.08);
	  tmpratio->GetYaxis()->SetNdivisions(210, kTRUE);
	  tmpratio->GetXaxis()->SetTitleFont(43);
	  tmpratio->SetMaximum(1.6);
	  tmpratio->SetMinimum(0.4);
	  tmpratio->Draw("E2");

	  const TAxis *axis = tmpratio->GetXaxis();
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

	void DrawStackedMCAndSyst(TString cName, TH1D* MC, TObjArray hist_nue_array, TH1D* sys, TH1D* ratio, std::string var, std::string title){

	  TCanvas can(cName,cName,2400,1600);
	  TH1D* tmpMC = (TH1D*)MC->Clone("tempCV");
	  TH1D* tmpratio = (TH1D*)ratio->Clone("tempRatio");
	  TH1D* sysErr = (TH1D*)sys->Clone("totalError");
	  sysErr->Reset();
	  //sysErr->Scale(1,"width");
	  for( int bin = 1; bin < tmpMC->GetNbinsX()+1; bin++){
	    sysErr->SetBinContent(bin,tmpMC->GetBinContent(bin));
	    sysErr->SetBinError(bin,sys->GetBinContent(bin));
	    tmpratio->SetBinContent(bin,1);
	    tmpratio->SetBinError(bin,ratio->GetBinError(bin));
	    //std::cout << "bin, cv, error = " << bin << ", " << tmpratio->GetBinContent(bin) << ", " << tmpratio->GetBinError(bin) << std::endl;
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
	  sysErr->SetTitleSize(120);
	  sysErr->SetTitleFont(43);
	  sysErr->SetTitle(title.c_str());
	  sysErr->SetLineColor(kBlack);
	  sysErr->SetLineWidth(1);
	  //sysErr->SetFillColor(kGray-10);
	  //sysErr->SetFillStyle(3244);
	  sysErr->SetStats(0);
	  sysErr->GetYaxis()->SetTitleSize(25);
	  sysErr->GetYaxis()->SetTitleFont(43);
	  sysErr->GetYaxis()->SetTitle("Events");
	  sysErr->GetXaxis()->SetLabelSize(0);
	  sysErr->GetXaxis()->SetTitleFont(43);
	  sysErr->GetXaxis()->SetTitle(var.c_str());
	  sysErr->GetXaxis()->SetTitleSize(30);
	  sysErr->SetMaximum(1.2*sysErr->GetMaximum());
	  std::cout << "call draw stacked" <<std::endl; 
	  Draw_Stacked2(hist_nue_array, tmpMC, cName, 0, false, "#nu_{e} Inclusive Selection", var );
	  gPad->Modified(); // so it updates the pad with the new changes
	 
	  // lower plot will be in pad
	  can.cd();          // Go back to the main canvas before defining pad2
	  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.35);
	  pad2->SetTopMargin(0.05);
	  pad2->SetBottomMargin(0.25);
	  pad2->SetGridx(); // vertical grid
	  pad2->Draw();
	  pad2->cd();       // pad2 becomes the current pad

	  tmpratio->SetTitle("");
	  tmpratio->GetYaxis()->SetTitle("uncertainties");
	  tmpratio->GetXaxis()->SetTitle(var.c_str());
	  tmpratio->SetStats(0);
	  tmpratio->SetLineColor(kBlack);
	  tmpratio->SetLineWidth(1);
	  tmpratio->SetFillColor(kGray);
	  tmpratio->SetFillStyle(1001);
	  tmpratio->GetYaxis()->SetTitleSize(60);
	  tmpratio->GetYaxis()->SetTitleFont(43);
	  tmpratio->GetXaxis()->SetTitleSize(60);
	  tmpratio->GetXaxis()->SetTitleOffset(3.);
	  tmpratio->GetYaxis()->SetTitleOffset(1.0);
	  tmpratio->GetXaxis()->SetLabelSize(0.1);
	  tmpratio->GetYaxis()->SetLabelSize(0.1);
	  tmpratio->GetYaxis()->SetNdivisions(210, kTRUE);
	  tmpratio->GetXaxis()->SetTitleFont(43);
	  tmpratio->SetMaximum(1.6);
	  tmpratio->SetMinimum(0.4);
	  tmpratio->Draw("E2");

	  const TAxis *axis = tmpratio->GetXaxis();
	  double lowX  = axis->GetBinLowEdge( axis->GetFirst() );
	  double highX = axis->GetBinUpEdge(  axis->GetLast() );

	  TLine line;
	  line.SetLineStyle(2);
	  line.SetLineWidth(2);
	  line.SetLineColor(kBlack);
	  line.DrawLine(lowX, 1., highX, 1.); //creates a new line which is owned by gPad
	  
	  gPad->Modified(); // so it updates the pad with the new changes

	  gPad->RedrawAxis();
	  gPad->Update();
	 
	  can.Print(cName+".pdf","pdf");

	}

	void DrawSyst(TString cName, TH1D* MC, TObjArray hist_nue_array, TH1D* sys, TH1D* ratio1, TH1D* ratio2, std::string var, std::string title){

	  TCanvas can(cName,cName,2400,1800);
	  TH1D* tmpMC = (TH1D*)MC->Clone("tempCV");
	  TH1D* tmpratio1 = (TH1D*)ratio1->Clone("tempRatio1");
	  TH1D* tmpratio2 = (TH1D*)ratio2->Clone("tempRatio2");
	  TH1D* sysErr = (TH1D*)sys->Clone("totalError");
	  sysErr->Reset();
	  //sysErr->Scale(1,"width");
	  for( int bin = 1; bin < tmpMC->GetNbinsX()+1; bin++){
	    sysErr->SetBinContent(bin,tmpMC->GetBinContent(bin));
	    sysErr->SetBinError(bin,sys->GetBinContent(bin));
	    tmpratio1->SetBinContent(bin,1);
	    tmpratio1->SetBinError(bin,ratio1->GetBinError(bin));
	    tmpratio2->SetBinContent(bin,1);
	    tmpratio2->SetBinError(bin,ratio2->GetBinError(bin));
	    //std::cout << "bin, cv, error = " << bin << ", " << tmpratio->GetBinContent(bin) << ", " << tmpratio->GetBinError(bin) << std::endl;
	  }
	  
	  std::string plottitle = std::string(((TH1F*)hist_nue_array[0])->GetName());
	  std::cout << "plottitle " << plottitle << std::endl;
	  std::string stacktitle;
	  if( plottitle.find("uBooNE_nue_") != std::string::npos ) stacktitle = "Systematics Uncertainties, #nu_{e} 1eNp0#pi selection"; 
	  else if( plottitle.find("uBooNE_1e0p_") != std::string::npos ) stacktitle = "Systematics Uncertainties, #nu_{e} 1e0p0#pi selection";
	  else if( plottitle.find("uBooNE_numu_") != std::string::npos ) stacktitle = "#nu_{#mu} inclusive selection"; 
	  //else if( plottitle.find("uBooNE_nueall_") != std::string::npos ) stacktitle = "#nu_{e} inclusive selection"; 

	  tmpratio1->SetTitle(stacktitle.c_str());
	  tmpratio1->GetYaxis()->SetTitle("uncertainties");
	  tmpratio1->GetXaxis()->SetTitle( "Reconstructed Energy [GeV]" );
	  tmpratio1->SetStats(0);
	  tmpratio1->SetLineColor(kBlack);
	  tmpratio1->SetLineWidth(1);
	  tmpratio1->SetFillColor(kGray);
	  tmpratio2->SetFillColorAlpha(kRed, 0.65);
	  tmpratio1->SetFillStyle(1001);
	  tmpratio2->SetFillStyle(1001);
	  tmpratio1->GetYaxis()->SetTitleSize(80);
	  tmpratio1->GetYaxis()->SetTitleFont(43);
	  tmpratio1->GetXaxis()->SetTitleSize(80);
	  tmpratio1->GetXaxis()->SetTitleOffset(0.95);
	  tmpratio1->GetYaxis()->SetTitleOffset(0.95);
	  tmpratio1->GetXaxis()->SetLabelSize(0.045);
	  tmpratio1->GetYaxis()->SetLabelSize(0.045);
	  tmpratio1->GetYaxis()->SetNdivisions(210, kTRUE);
	  tmpratio1->GetXaxis()->SetTitleFont(43);
	  tmpratio1->GetXaxis()->CenterTitle(kTRUE);
	  tmpratio1->SetMaximum(1.6);
	  tmpratio1->SetMinimum(0.4);
	  tmpratio1->Draw("E2");
	  tmpratio2->Draw("E2 same");

	  const TAxis *axis = tmpratio1->GetXaxis();
	  double lowX  = axis->GetBinLowEdge( axis->GetFirst() );
	  double highX = axis->GetBinUpEdge(  axis->GetLast() );

	  TLine line;
	  line.SetLineStyle(2);
	  line.SetLineWidth(2);
	  line.SetLineColor(kBlack);
	  TLegend  *leg2 = new TLegend(0.25,0.78,0.75,0.88);
	  leg2->SetNColumns(2);
	  leg2->SetTextSize(0.045);
	  leg2->AddEntry(tmpratio1,"before const.","f"); 
	  leg2->AddEntry(tmpratio2,"after const.","f"); 
	  leg2->Draw("same");

	  gPad->Modified(); // so it updates the pad with the new changes
	  //pad2->Draw("");

	  //gPad->RedrawAxis();
	  gPad->Update();
	 
	  can.Print(cName+"_systonly.pdf","pdf");

	}
	void DrawStackedMCAndSyst2(TString cName, TH1D* MC, TObjArray hist_nue_array, TH1D* sys, TH1D* ratio1, TH1D* ratio2, std::string var, std::string title){

	  TCanvas can(cName,cName,2400,1800);
	  TH1D* tmpMC = (TH1D*)MC->Clone("tempCV");
	  TH1D* tmpratio1 = (TH1D*)ratio1->Clone("tempRatio1");
	  TH1D* tmpratio2 = (TH1D*)ratio2->Clone("tempRatio2");
	  TH1D* sysErr = (TH1D*)sys->Clone("totalError");
	  sysErr->Reset();
	  //sysErr->Scale(1,"width");
	  for( int bin = 1; bin < tmpMC->GetNbinsX()+1; bin++){
	    sysErr->SetBinContent(bin,tmpMC->GetBinContent(bin));
	    sysErr->SetBinError(bin,sys->GetBinContent(bin));
	    tmpratio1->SetBinContent(bin,1);
	    tmpratio1->SetBinError(bin,ratio1->GetBinError(bin));
	    tmpratio2->SetBinContent(bin,1);
	    tmpratio2->SetBinError(bin,ratio2->GetBinError(bin));
	    //std::cout << "bin, cv, error = " << bin << ", " << tmpratio->GetBinContent(bin) << ", " << tmpratio->GetBinError(bin) << std::endl;
	  }
	  
	  // Upper plot will be in pad1
	  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.5, 1, 1.0);
	  pad1->SetBottomMargin(0.025); // Upper and lower plot are joined
	  pad1->SetGridx(2);        // Vertical grid
	  pad1->Draw();             // Draw the upper pad: pad1
	  pad1->cd();               // pad1 becomes the current pad

	  //tmpMC->SetLineColor(kRed);
	  tmpMC->SetLineWidth(3);
	  tmpMC->SetLineStyle(1);
	  /*sysErr->SetTitleSize(120);
	  sysErr->SetTitleFont(43);
	  sysErr->SetTitle(title.c_str());
	  sysErr->SetLineColor(kBlack);
	  sysErr->SetLineWidth(1);
	  //sysErr->SetFillColor(kGray-10);
	  //sysErr->SetFillStyle(3244);
	  sysErr->SetStats(0);
	  sysErr->GetYaxis()->SetTitleSize(25);
	  sysErr->GetYaxis()->SetTitleFont(43);
	  sysErr->GetYaxis()->SetTitle("Events");
	  sysErr->GetXaxis()->SetLabelSize(0);
	  sysErr->GetXaxis()->SetTitleFont(43);
	  sysErr->GetXaxis()->SetTitle(var.c_str());
	  sysErr->GetXaxis()->SetTitleSize(30);
	  sysErr->SetMaximum(1.2*sysErr->GetMaximum());*/
	  std::cout << "call draw stacked" <<std::endl; 
	  Draw_Stacked2(hist_nue_array, tmpMC, cName, 0, false, "#nu_{e} Inclusive Selection", var );
	  //sysErr->Draw("E1 same");
	  gPad->Modified(); // so it updates the pad with the new changes
	 
	  // lower plot will be in pad
	  can.cd();          // Go back to the main canvas before defining pad2
	  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.5);
	  pad2->SetTopMargin(0.01);
	  pad2->SetBottomMargin(0.2);
	  pad2->SetGridx(); // vertical grid
	  pad2->Draw();
	  pad2->cd();       // pad2 becomes the current pad

	  tmpratio1->SetTitle("");
	  tmpratio1->GetYaxis()->SetTitle("uncertainties");
	  tmpratio1->GetXaxis()->SetTitle(var.c_str());
	  tmpratio1->SetStats(0);
	  tmpratio1->SetLineColor(kBlack);
	  tmpratio1->SetLineWidth(1);
	  tmpratio1->SetFillColor(kGray);
	  tmpratio2->SetFillColorAlpha(kRed, 0.65);
	  tmpratio1->SetFillStyle(1001);
	  tmpratio2->SetFillStyle(1001);
	  tmpratio1->GetYaxis()->SetTitleSize(60);
	  tmpratio1->GetYaxis()->SetTitleFont(43);
	  tmpratio1->GetXaxis()->SetTitleSize(65);
	  tmpratio1->GetXaxis()->SetTitleOffset(2.0);
	  tmpratio1->GetYaxis()->SetTitleOffset(0.85);
	  tmpratio1->GetXaxis()->SetLabelSize(0.065);
	  tmpratio1->GetYaxis()->SetLabelSize(0.065);
	  tmpratio1->GetYaxis()->SetNdivisions(210, kTRUE);
	  tmpratio1->GetXaxis()->SetTitleFont(43);
	  tmpratio1->SetMaximum(1.6);
	  tmpratio1->SetMinimum(0.4);
	  tmpratio1->Draw("E2");
	  tmpratio2->Draw("E2 same");

	  const TAxis *axis = tmpratio1->GetXaxis();
	  double lowX  = axis->GetBinLowEdge( axis->GetFirst() );
	  double highX = axis->GetBinUpEdge(  axis->GetLast() );

	  TLine line;
	  line.SetLineStyle(2);
	  line.SetLineWidth(2);
	  line.SetLineColor(kBlack);
	  TLegend  *leg2 = new TLegend(0.28,0.78,0.6,0.93);
	  leg2->SetNColumns(2);
	  leg2->SetTextSize(0.07);
	  leg2->AddEntry(tmpratio1,"before const.","f"); 
	  leg2->AddEntry(tmpratio2,"after const.","f"); 
	  leg2->Draw("same");

	  gPad->Modified(); // so it updates the pad with the new changes
	  //pad2->Draw("");

	  //gPad->RedrawAxis();
	  gPad->Update();
	 
	  can.Print(cName+".pdf","pdf");

	}


	//borrow this function from SBNchi
	void CollapseSubchannels(TMatrixD & M, TMatrixD & Mc, std::vector<int> num_bins, int num_channels, std::vector<int> num_subchannels){
	  bool debug = true;
	  if(debug) std::cout<<"Starting:M "<<M.GetNcols()<<" "<<M.GetNrows()<<" "<<std::endl;
	  if(debug) std::cout<<"Starting:Mc "<<Mc.GetNcols()<<" "<<Mc.GetNrows()<<" "<<std::endl;
	 
	  //std::cout << "num_channels, num_subchannels = " << num_channels << " , " << num_subchannels[0] << ", " << num_subchannels[1] << std::endl;
	  std::vector<std::vector<TMatrixD>> Summed(num_channels, std::vector<TMatrixD>(num_channels) );	//Initialise a matrix of matricies, to ZERO.
	  for(int ic = 0; ic < num_channels; ic++){
	    for(int jc =0; jc < num_channels; jc++){
	      //std::cout << "num_bins[jc],num_bins[ic] = " << num_bins[jc] << ", " << num_bins[ic] << std::endl;
	      Summed[ic][jc].ResizeTo(num_bins[jc],num_bins[ic]) ;// This is CORRECT, do not switch (ie Summed[0][1] = size (num_bins[1], num_bins[0])
	      Summed[ic][jc] = 0.0;
	    }
	  }
	  
	  int mrow = 0.0;
	  int mcol = 0.0;

	  for(int i=0; i < M.GetNcols(); i++ ){
	  for(int j=0; j < M.GetNrows(); j++ ){
	    if(isnan( M(i,j) )) M(i,j) = 0.0;
	  }
	  }

	  //std::cout << "****M.Print()****" << std::endl;
	  //M.Print();

	  for(int ic = 0; ic < num_channels; ic++){ 	 //Loop over all rows
	    for(int jc =0; jc < num_channels; jc++){ //Loop over all columns
	      
	      if(debug) std::cout<<"Diagonal! : "<<ic<<" "<<jc<<" mcol is: "<<mcol<<" mrow is: "<<mrow<<std::endl;
	      
	      for(int m=0; m < num_subchannels[ic]; m++){
		for(int n=0; n< num_subchannels[jc]; n++){ //For each big block, loop over all subchannels summing toGether
		  ////std::cout << mrow << ", " << n << ", " << num_bins[jc] << ", " << ", " << mcol << ", " << m << ", " << m << ", " << num_bins[ic] << std::endl;
		  Summed[ic][jc] +=  M.GetSub(mrow+n*num_bins[jc] ,mrow + n*num_bins[jc]+num_bins[jc]-1, mcol + m*num_bins[ic], mcol+ m*num_bins[ic]+num_bins[ic]-1 );
		  ////std::cout << "Summed[" << ic << "][" << jc << "] = " << std::endl;
		  //Summed[ic][jc].Print();
		}
	      }
	      mrow += num_subchannels[jc]*num_bins[jc];//As we work our way left in columns, add on that many bins
	    }//end of column loop
	    
	    //std::cout << "mrow  " << mrow << std::endl;  
	    mrow = 0; // as we end this row, reSet row count, but jump down 1 column
	    mcol += num_subchannels[ic]*num_bins[ic];
	    //std::cout << "collapsed mrow,mcol  = " << mrow << ", " << mcol << std::endl;
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
	    //std::cout << "mrow  " << mrow << std::endl;  
	    mrow = 0;
	    mcol +=num_bins[ic];
	  }

	  //std::cout << "collapsed mrow,mcol  = " << mrow << ", " << mcol << std::endl;
	  return;
	}

	void plot_one_noconstraint(TMatrixD matrix, TH1D *h_nue, std::string tag){
            gStyle->SetPalette(kLightTemperature);
	    std::string dir="/uboone/data/users/wospakrk/SBNFitPlots/";
	    //std::vector<std::string> channel_names = {"nu uBooNE nue intrinsic","nu uBooNE constraint"};
	    std::vector<std::string> channel_names = {"1eNp","1e0p"};
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
	    gStyle->SetOptStat(0);
	    c_full->Print((tag+"_collapsed.pdf").c_str(),"pdf");

	}
	void plot_one(TMatrixD matrix, TH1D *h_nue, TH1D *h_constrain, std::string tag){
            gStyle->SetPalette(kLightTemperature);
	    std::string dir="/uboone/data/users/wospakrk/SBNFitPlots/";
	    std::vector<std::string> channel_names;
            if(tag.find("nue_numu") != std::string::npos ) channel_names.push_back("1eNp");
            else channel_names.push_back("1e0p");
            channel_names.push_back("numu");
	    int num_channels = 2;
	    int num_bins_total = matrix.GetNrows();
	    std::vector<int> num_bins = {h_nue->GetNbinsX(), h_constrain->GetNbinsX()}; 
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
	    gStyle->SetOptStat(0);
	    c_full->Print((tag+"_collapsed.pdf").c_str(),"pdf");
	}

	void Draw_Stacked(TObjArray histos, TH1D *data,
			  TPad *pad,
			  bool normalised,
			  std::string stacktitle,
			  std::string var )
	{// this function draws histograms stacked and correctly takes into account the
	  // stats boxes for each

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
	  
	  // lets take the first histoname and see if we can match a title to its which will be HS stack title
	  if( stacktitle == "" ) stacktitle = ((TH1F*)histos[0])->GetTitle();
	  std::string plottitle = std::string(((TH1F*)histos[0])->GetTitle());
	  //if( plottitle.find("uBooNE_1e0p_") != std::string::npos ) stacktitle = "#nu_{e} 1e0p0#pi selection";
	  //else if( plottitle.find("uBooNE_constrain_inclusive") != std::string::npos ) stacktitle = "#nu_{#mu} inclusive selection"; 
	  if( plottitle.find("uBooNE_nue_") != std::string::npos ) stacktitle = "#nu_{e} 1eNp0#pi selection"; 
	  if( plottitle.find("uBooNE_1e0p_") != std::string::npos ) stacktitle = "#nu_{e} 1e0p0#pi selection"; 
	  if( plottitle.find("uBooNE_numu_") != std::string::npos ) stacktitle = "#nu_{#mu} Contained Selection"; 
	  stacktitle = "#nu_{e} 1eNp selection";  
	  THStack *Hs = new THStack("hs2",stacktitle.c_str());


	  // Set up the LEGEND
	  TLegend  *legend = new TLegend(0.60,0.50,0.98,0.9); // we need different positions for the legend to not 
	  // get the plot titles for the legend
	  for(int i = 0; i < n ; i++){
	    
	    TH1D *h =  (TH1D*)histos[i];
	    //std::cout << "histo = " << h->GetTitle() << std::endl;
	    //h->Scale(1,"width");
	    h->SetLineWidth(0);
	    //h->SetLineColor(colours[i]);
	    //h->SetFillColor(colours[i]);
	    std::string histosample = std::string(h->GetTitle());
	    if( histosample.find("ext") != std::string::npos ){ 
	      h->SetFillStyle(3354);
	      h->SetFillColor(kBlack);
	      h->SetLineColor(kBlack);
	      std::string htitle = Form("BNB Off: %2.2f",h->Integral());
	      //std::cout << "htitle = " << htitle << std::endl;
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
	      std::string htitle = Form("NC 0 #pi: %2.2f",h->Integral());
	      //std::cout << "htitle = " << htitle << std::endl;
	      h->SetTitle(htitle.c_str());
	    }
	    else if( histosample.find("CCmuNoPi") != std::string::npos ){
	      h->SetFillColor(kAzure-2);
	      h->SetLineColor(kAzure-2);
	      std::string htitle = Form("#nu_#mu CC 0 #pi: %2.2f",h->Integral());
	      //std::cout << "htitle = " << htitle << std::endl;
	      h->SetTitle(htitle.c_str());
	    }
	    else if( histosample.find("NCcPiNoPi0") != std::string::npos ){
	      h->SetFillColor(kAzure-5);
	      h->SetLineColor(kAzure-5);
	      std::string htitle = Form("NC c#pi 0#pi^{0}: %2.2f",h->Integral());
	      //std::cout << "htitle = " << htitle << std::endl;
	      h->SetTitle(htitle.c_str());
	    }
	    else if( histosample.find("CCmuCPiNoPi0") != std::string::npos ){
	      h->SetFillColor(kAzure-8);
	      h->SetLineColor(kAzure-8);
	      std::string htitle = Form("#nu_#mu CC c#pi 0#pi^{0}: %2.2f",h->Integral());
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
	    else if( histosample.find("nue_nu") != std::string::npos || histosample.find("1e0p_nu") != std::string::npos || histosample.find("nueall_nu") != std::string::npos || histosample.find("numu_nu") != std::string::npos ){
	      h->SetFillColor(kCyan);
	      h->SetLineColor(kCyan);
	      std::string htitle = Form("#nu_{#mu} BNB: %2.2f",h->Integral());
	      //std::cout << "htitle = " << htitle << std::endl;
	      h->SetTitle(htitle.c_str());
	    }
	    legends_str.push_back( h->GetTitle() );
	    legend->AddEntry(h,legends_str[i].c_str(),"f"); 

	    //std::cout << "histo name =  " << h->GetTitle() << std::endl;
	    //std::cout << "integral =  " << h->Integral() << std::endl;
	    Hs->Add(h,"sames");
	  }
	  
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
	  Hs->GetXaxis()->SetTitleSize(80);
	  Hs->GetYaxis()->SetTitleSize(80);
	  Hs->GetXaxis()->SetLabelFont(43);
	  Hs->GetYaxis()->SetLabelFont(43);
	  Hs->GetXaxis()->SetLabelSize(0.1);
	  Hs->GetYaxis()->SetLabelSize(0.1);
	  Hs->GetXaxis()->CenterTitle(kTRUE);

	  //legend->Draw("");
	  
	  gPad->Modified(); // so it updates the pad with the new changes
	 
	}


	void Draw_Stacked2(TObjArray histos, TH1D *data, 
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
	  TPaveStats *st1;
	  
	  const int n = histos.GetEntries();
	  
	  // lets take the first histoname and see if we can match a title to its which will be HS stack title
	  if( stacktitle == "" ) stacktitle = ((TH1F*)histos[0])->GetTitle();
	  std::string plottitle = std::string(((TH1F*)histos[0])->GetName());
	  std::cout << "plottitle " << plottitle << std::endl;
	  if( plottitle.find("uBooNE_nue_") != std::string::npos ) stacktitle = "#nu_{e} 1eNp0#pi selection"; 
	  else if( plottitle.find("uBooNE_1e0") != std::string::npos ) stacktitle = "#nu_{e} 1e0p0#pi selection";
	  //else if( plottitle.find("uBooNE_numu_") != std::string::npos ) stacktitle = "#nu_{#mu} inclusive selection"; 
	  //else if( plottitle.find("uBooNE_nueall_") != std::string::npos ) stacktitle = "#nu_{e} inclusive selection"; 
	  THStack *Hs = new THStack("hs2",stacktitle.c_str());
	  
	  //Set Axis Units
	  //Hs->GetXaxis()->SetTitle("Number of Hits in Slice");
	  
	  // Set up the LEGEND
	  TLegend  *legend = new TLegend(0.49,0.45,0.9,0.9); // we need different positions for the legend to not 
	  legend->SetHeader("MicroBooNE Simulation, Preliminary 6.95x10^{20}POT","L");
	  legend->SetNColumns(2);
	  //legend->SetTextSize(0.052);
	  // get the plot titles for the legend
	  for(int i = 0; i < n ; i++){
	    
	    TH1D *h =  (TH1D*)histos[i];
	    std::cout << "histo = " << h->GetTitle() << std::endl;
	    //h->Scale(1,"width");
	    h->SetLineWidth(0);
	    std::string histosample = std::string(h->GetTitle());
	    if( histosample.find("ext") != std::string::npos ){ 
	      h->SetFillStyle(3354);
	      h->SetFillColor(kBlack);
	      h->SetLineColor(kBlack);
	      std::string htitle = Form("BNB Off: %2.2f",h->Integral());
	      std::cout << "cek htitle = " << htitle << std::endl;
	      h->SetTitle(htitle.c_str());
	    }
	    else if( histosample.find("lee") != std::string::npos ){
	      h->SetFillColor(kOrange+1);
	      h->SetLineColor(kOrange+1);
	      std::string htitle = Form("#nu_{e} LEE: %2.2f",h->Integral());
	      std::cout << "cek htitle = " << htitle << std::endl;
	      h->SetTitle(htitle.c_str());
	    }
	    else if( histosample.find("ncpi0") != std::string::npos ){
	      h->SetFillColor(kAzure+1);
	      h->SetLineColor(kAzure+1);
	      std::string htitle = Form("#nu NC #pi^{0}: %2.2f",h->Integral());
	      std::cout << "cek htitle = " << htitle << std::endl;
	      h->SetTitle(htitle.c_str());
	    }
	    else if( histosample.find("ccpi0") != std::string::npos ){
	      h->SetFillColor(kAzure+4);
	      h->SetLineColor(kAzure+4);
	      std::string htitle = Form("#nu CC #pi^{0}: %2.2f",h->Integral());
	      std::cout << "cek htitle = " << htitle << std::endl;
	      h->SetTitle(htitle.c_str());
	    }
	    else if( histosample.find("intrinsic") != std::string::npos ){
	      h->SetFillColor(kSpring);
	      h->SetLineColor(kSpring);
	      std::string htitle = Form("#nu_{e} CC: %2.2f",h->Integral());
	      std::cout << "cek htitle = " << htitle << std::endl;
	      h->SetTitle(htitle.c_str());
	    }
	    else if( histosample.find("ncnopi") != std::string::npos ){
	      h->SetFillColor(kAzure+7);
	      h->SetLineColor(kAzure+7);
	      std::string htitle = Form("NC 0#pi: %2.2f",h->Integral());
	      std::cout << "htitle = " << htitle << std::endl;
	      h->SetTitle(htitle.c_str());
	    }
	    else if( histosample.find("CCmuNoPi") != std::string::npos ){
	      h->SetFillColor(kAzure-2);
	      h->SetLineColor(kAzure-2);
	      std::string htitle = Form("#nu_{#mu} CC 0#pi: %2.2f",h->Integral());
	      std::cout << "htitle = " << htitle << std::endl;
	      h->SetTitle(htitle.c_str());
	    }
	    else if( histosample.find("NCcPiNoPi0") != std::string::npos ){
	      h->SetFillColor(kAzure-5);
	      h->SetLineColor(kAzure-5);
	      std::string htitle = Form("NC #pi 0#pi^{0}: %2.2f",h->Integral());
	      std::cout << "htitle = " << htitle << std::endl;
	      h->SetTitle(htitle.c_str());
	    }
	    else if( histosample.find("CCmuCPiNoPi0") != std::string::npos ){
	      h->SetFillColor(kAzure-8);
	      h->SetLineColor(kAzure-8);
	      std::string htitle = Form("#nu_{#mu} CC #pi 0#pi^{0}: %2.2f",h->Integral());
	      std::cout << "htitle = " << htitle << std::endl;
	      h->SetTitle(htitle.c_str());
	    }
	    else if( histosample.find("dirt") != std::string::npos ){
	      h->SetFillColor(kOrange+3);
	      h->SetLineColor(kOrange+3);
	      std::string htitle = Form("Dirt: %2.2f",h->Integral());
	      std::cout << "cek htitle = " << htitle << std::endl;
	      h->SetTitle(htitle.c_str());
	    }
	    else if( histosample.find("nue_nu") != std::string::npos || histosample.find("nue0p_nu") != std::string::npos || histosample.find("nueall_nu") != std::string::npos || histosample.find("numu_nu") != std::string::npos || histosample.find("1e0p_nu") != std::string::npos ){
	      h->SetFillColor(kCyan);
	      h->SetLineColor(kCyan);
	      std::string htitle = Form("#nu_{#mu} CC: %2.2f",h->Integral());
	      std::cout << "cek htitle = " << htitle << std::endl;
	      h->SetTitle(htitle.c_str());
	    }
	    gStyle->SetLegendTextSize(0.05);
	    legends_str.push_back( h->GetTitle() );
	    legend->AddEntry(h,legends_str[i].c_str(),"f"); 
	    TLegendEntry *header = (TLegendEntry*)legend->GetListOfPrimitives()->First();
	    header->SetTextSize(.048);
	    std::cout << "histo name =  " << h->GetTitle() << std::endl;
	    std::cout << "histosample =  " << histosample << std::endl;
	    Hs->Add(h,"same");
	  }
	  legend->AddEntry(data,"Stats. Error","lep"); 
	  
	  float heightboxes;
	  // the position of the top corner
	  float top_corner,deltay;
	  // don't draw errors if it is 1
	  int no_error = 1;
	  top_corner = 0.9;

	  
	  //HACK -need to draw a smaller histogram first to set the x-axis range correctly 
	  //for stacks with variable-width bins
	  //TH1D* tmp_mnv = (TH1D*)Hs->GetHists()->At(0);
	  Hs->SetMaximum(1.7*data->GetMaximum());
	  
	  // DRAW not stacked to get correct stats boxes
	  if (no_error == 1) // do not draw errors
	    Hs->Draw("HIST");
	  else
	    Hs->Draw("HISTE");

	  //Hs->SetTitle(50);
	  Hs->GetXaxis()->SetTitle( "Reconstructed Energy [GeV]" );
	  Hs->GetYaxis()->SetTitle( "Events" );
	  Hs->GetXaxis()->SetTitleFont(43);
	  Hs->GetYaxis()->SetTitleFont(43);
	  Hs->GetXaxis()->SetTitleSize(60);
	  Hs->GetYaxis()->SetTitleSize(60);
	  Hs->GetXaxis()->SetLabelFont(43);
	  Hs->GetYaxis()->SetLabelFont(43);
	  Hs->GetXaxis()->SetLabelSize(60);
	  Hs->GetXaxis()->SetLabelOffset(3);
	  Hs->GetYaxis()->SetLabelSize(60);
	  Hs->GetXaxis()->CenterTitle(kTRUE);

	  data->SetMarkerStyle(1);
	  data->SetMarkerSize(2);
	  data->SetMarkerColor(kBlack);
	  data->SetLineWidth(2);
	  data->SetLineColor(kBlack);
	  
	  data->Draw("PE1 same");
	  // DRAW not stacked to get correct stats boxes
	  if (no_error == 1) // do not draw errors
	    Hs->Draw("HIST");
	  else
	    Hs->Draw("HISTE");

	  data->DrawCopy("PE1 same");
	  legend->Draw("same");
	  AddPlotLabel("MicroBooNE Preliminary, Run 1 scaled to 1.01e+21 POT",0.4,0.8);
	  
	  gPad->Modified(); // so it updates the pad with the new changes
	}

	void Draw_Stacked3(TObjArray histos, TH1D *data, 
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
	  TPaveStats *st1;
	  
	  const int n = histos.GetEntries();
	  
	  // lets take the first histoname and see if we can match a title to its which will be HS stack title
	  if( stacktitle == "" ) stacktitle = ((TH1F*)histos[0])->GetTitle();
	  std::string plottitle = std::string(((TH1F*)histos[0])->GetName());
	  std::cout << "plottitle " << plottitle << std::endl;
	  if( plottitle.find("uBooNE_nue") != std::string::npos ) stacktitle = "#nu_{e} 1eNp0#pi selection"; 
	  else if( plottitle.find("uBooNE_1e0p") != std::string::npos ) stacktitle = "#nu_{e} 1e0p0#pi selection";
	  else if( plottitle.find("uBooNE_numu_") != std::string::npos ) stacktitle = "#nu_{#mu} inclusive selection"; 
	  //else if( plottitle.find("uBooNE_nueall_") != std::string::npos ) stacktitle = "#nu_{e} inclusive selection"; 
	  THStack *Hs = new THStack("hs2",stacktitle.c_str());
	  
	  //Set Axis Units
	  //Hs->GetXaxis()->SetTitle("Number of Hits in Slice");
	  
	  // Set up the LEGEND
	  TLegend  *legend = new TLegend(0.43,0.57,0.9,0.9); // we need different positions for the legend to not 
	  legend->SetHeader("MicroBooNE Simulation, Preliminary 6.95x10^{20}POT","L");
	  legend->SetNColumns(2);
	  //legend->SetTextSize(0.052);
	  // get the plot titles for the legend
	  for(int i = 0; i < n ; i++){
	    
	    TH1D *h =  (TH1D*)histos[i];
	    std::cout << "histo = " << h->GetTitle() << std::endl;
	    //h->Scale(1,"width");
	    h->SetLineWidth(0);
	    std::string histosample = std::string(h->GetTitle());
	    if( histosample.find("ext") != std::string::npos ){ 
	      h->SetFillStyle(3354);
	      h->SetFillColor(kBlack);
	      h->SetLineColor(kBlack);
	      std::string htitle = Form("BNB Off: %2.2f",h->Integral());
	      std::cout << "cek htitle = " << htitle << std::endl;
	      h->SetTitle(htitle.c_str());
	    }
	    else if( histosample.find("lee") != std::string::npos ){
	      h->SetFillColor(kOrange+1);
	      h->SetLineColor(kOrange+1);
	      std::string htitle = Form("#nu_{e} LEE: %2.2f",h->Integral());
	      std::cout << "cek htitle = " << htitle << std::endl;
	      h->SetTitle(htitle.c_str());
	    }
	    else if( histosample.find("ncpi0") != std::string::npos ){
	      h->SetFillColor(kAzure+1);
	      h->SetLineColor(kAzure+1);
	      std::string htitle = Form("#nu NC #pi^{0}: %2.2f",h->Integral());
	      std::cout << "cek htitle = " << htitle << std::endl;
	      h->SetTitle(htitle.c_str());
	    }
	    else if( histosample.find("ccpi0") != std::string::npos ){
	      h->SetFillColor(kAzure+4);
	      h->SetLineColor(kAzure+4);
	      std::string htitle = Form("#nu CC #pi^{0}: %2.2f",h->Integral());
	      std::cout << "cek htitle = " << htitle << std::endl;
	      h->SetTitle(htitle.c_str());
	    }
	    else if( histosample.find("intrinsic") != std::string::npos ){
	      h->SetFillColor(kSpring);
	      h->SetLineColor(kSpring);
	      std::string htitle = Form("#nu_{e} CC: %2.2f",h->Integral());
	      std::cout << "cek htitle = " << htitle << std::endl;
	      h->SetTitle(htitle.c_str());
	    }
	    else if( histosample.find("ncnopi") != std::string::npos ){
	      h->SetFillColor(kAzure+7);
	      h->SetLineColor(kAzure+7);
	      std::string htitle = Form("NC 0#pi: %2.2f",h->Integral());
	      std::cout << "htitle = " << htitle << std::endl;
	      h->SetTitle(htitle.c_str());
	    }
	    else if( histosample.find("CCmuNoPi") != std::string::npos ){
	      h->SetFillColor(kAzure-2);
	      h->SetLineColor(kAzure-2);
	      std::string htitle = Form("#nu_{#mu} CC 0#pi: %2.2f",h->Integral());
	      std::cout << "htitle = " << htitle << std::endl;
	      h->SetTitle(htitle.c_str());
	    }
	    else if( histosample.find("NCcPiNoPi0") != std::string::npos ){
	      h->SetFillColor(kAzure-5);
	      h->SetLineColor(kAzure-5);
	      std::string htitle = Form("NC #pi 0#pi^{0}: %2.2f",h->Integral());
	      std::cout << "htitle = " << htitle << std::endl;
	      h->SetTitle(htitle.c_str());
	    }
	    else if( histosample.find("CCmuCPiNoPi0") != std::string::npos ){
	      h->SetFillColor(kAzure-8);
	      h->SetLineColor(kAzure-8);
	      std::string htitle = Form("#nu_{#mu} CC #pi 0#pi^{0}: %2.2f",h->Integral());
	      std::cout << "htitle = " << htitle << std::endl;
	      h->SetTitle(htitle.c_str());
	    }
	    else if( histosample.find("dirt") != std::string::npos ){
	      h->SetFillColor(kOrange+3);
	      h->SetLineColor(kOrange+3);
	      std::string htitle = Form("Dirt: %2.2f",h->Integral());
	      std::cout << "cek htitle = " << htitle << std::endl;
	      h->SetTitle(htitle.c_str());
	    }
	    else if( histosample.find("nue_nu") != std::string::npos || histosample.find("nue0p_nu") != std::string::npos || histosample.find("nueall_nu") != std::string::npos || histosample.find("numu_nu") != std::string::npos || histosample.find("1e0p_nu") != std::string::npos ){
	      h->SetFillColor(kCyan);
	      h->SetLineColor(kCyan);
	      std::string htitle = Form("#nu_{#mu} CC: %2.2f",h->Integral());
	      std::cout << "cek htitle = " << htitle << std::endl;
	      h->SetTitle(htitle.c_str());
	    }
	    gStyle->SetLegendTextSize(0.032);
	    legends_str.push_back( h->GetTitle() );
	    legend->AddEntry(h,legends_str[i].c_str(),"f"); 
	    TLegendEntry *header = (TLegendEntry*)legend->GetListOfPrimitives()->First();
	    header->SetTextSize(.03);
	    //header->SetTextFont(62);
	    std::cout << "histo name =  " << h->GetTitle() << std::endl;
	    std::cout << "histosample =  " << histosample << std::endl;
	    Hs->Add(h,"same");
	  }
	  legend->AddEntry(data,"Stats. Error","lep"); 
	  
	  float heightboxes;
	  // the position of the top corner
	  float top_corner,deltay;
	  // don't draw errors if it is 1
	  int no_error = 1;
	  top_corner = 0.9;

	  
	  //HACK -need to draw a smaller histogram first to set the x-axis range correctly 
	  //for stacks with variable-width bins
	  //TH1D* tmp_mnv = (TH1D*)Hs->GetHists()->At(0);
	  Hs->SetMaximum(1.7*data->GetMaximum());
	  
	  // DRAW not stacked to get correct stats boxes
	  if (no_error == 1) // do not draw errors
	    Hs->Draw("HIST");
	  else
	    Hs->Draw("HISTE");

	  //Hs->SetTitle(50);
	  Hs->GetXaxis()->SetTitle( "Reconstructed Energy [GeV]" );
	  Hs->GetYaxis()->SetTitle( "Events" );
	  Hs->GetXaxis()->SetTitleFont(43);
	  Hs->GetYaxis()->SetTitleFont(43);
	  Hs->GetXaxis()->SetTitleSize(80);
	  Hs->GetYaxis()->SetTitleSize(80);
	  Hs->GetXaxis()->SetLabelFont(43);
	  Hs->GetYaxis()->SetLabelFont(43);
	  Hs->GetXaxis()->SetLabelSize(80);
	  //Hs->GetXaxis()->SetLabelOffset(3);
	  Hs->GetYaxis()->SetLabelSize(80);
	  Hs->GetXaxis()->CenterTitle(kTRUE);

	  data->SetMarkerStyle(1);
	  data->SetMarkerSize(2);
	  data->SetMarkerColor(kBlack);
	  data->SetLineWidth(2);
	  data->SetLineColor(kBlack);
	  
	  data->Draw("PE1 same");
	  // DRAW not stacked to get correct stats boxes
	  if (no_error == 1) // do not draw errors
	    Hs->Draw("HIST");
	  else
	    Hs->Draw("HISTE");

	  data->DrawCopy("PE1 same");
	  legend->Draw("same");
	  AddPlotLabel("MicroBooNE Preliminary, Run 1 scaled to 1.01e+21 POT",0.4,0.8);
	  
	  gPad->Modified(); // so it updates the pad with the new changes
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
