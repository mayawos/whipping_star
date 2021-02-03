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

int main(int argc, char* argv[])
{
	std::string xml = "example.xml";
	std::string tag1 = "example1";
	std::string tag2 = "example2";
	double width = 1.0;
	std::string var = "Reconstructed Energy [GeV]";
	int iarg = 0;
	opterr=1;
	int index;
	bool sys = false;
        bool detsys = true;
	int mass_start = -1;

	const struct option longopts[] =
	{
		{"xml", 		required_argument, 	0, 'x'},
		{"tag1", 		required_argument, 	0, 't'},
		{"tag2", 		required_argument, 	0, 'u'},
		{"width", 		required_argument, 	0, 'w'},
		{"var", 		required_argument, 	0, 'v'},
		{"sys",	no_argument, 0, 's'},
		{"detsys",	        no_argument, 0, 'd'},
		{"part", required_argument,0,'p'},
		{0,			no_argument, 		0,  0},
	};

	while(iarg != -1)
	{
		iarg = getopt_long(argc,argv, "x:t:u:w:v:s:d:cp:g", longopts, &index);

		switch(iarg)
		{
			case 'x':
				xml = optarg;
				break;
			case 't':
				tag1 = optarg;
				break;
			case 'u':
				tag2 = optarg;
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
				std::cout<<"Allowed arguments:"<<std::endl;
				std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
				return 0;
		}
	}

        std::string varlabel = "";
        std::string chanlabel = "";
        if(tag1.find("reco_e") != std::string::npos) varlabel = "Enu_reco";	
        else if(tag1.find("true_e") != std::string::npos) varlabel = "Enu_true";	
        else if(tag1.find("reco_elep") != std::string::npos) varlabel = "Elepton_reco";	
        else if(tag1.find("true_elep") != std::string::npos) varlabel = "Elepton_true";	
        if(tag1.find("1e0p") != std::string::npos) chanlabel += "_1e0p";
        else chanlabel += "_1eNp"; 	
        std::vector<TH1D*> numus, nues, nues_nolee;
        //fetch the TFile for the MC spectrum
        TString fname = "/uboone/data/users/wospakrk/09-07-2020/"+tag1+".SBNspec.root";
        std::cout << "fname tag1 = " << fname << std::endl; 
        TFile *file = new TFile(fname,"read");
        //TH1D *h_nue_total, *h_numu;

        //loop counter(s)
        int i=0;
        int j=0;
        std::vector<std::string> names;
        TKey *key;
	TObject *obj;
        TIter nextkey(gDirectory->GetListOfKeys());
        while (key = (TKey*)nextkey()) {
          obj = key->ReadObj();
          std::string name = obj->GetName();
          std::cout << "name = " << name << std::endl;
          TH1D* h = (TH1D*)obj;
          names.push_back(name);
	  if( name.find("nu_uBooNE_nue_") != std::string::npos || name.find("nu_uBooNE_1e0p_") != std::string::npos ){ 
            if( name.find("lee") == std::string::npos ) nues_nolee.push_back(h);
            nues.push_back(h);
          }else{ 
            numus.push_back(h);
          }
	}
        std::cout << "add the MC histos" << std::endl;
        //now add the MC histos
        std::cout << "add the MC histos 1" << std::endl;

        TH1D *h_nue_total = (TH1D*)nues[0]->Clone("nue_total");
        h_nue_total->Reset();
        for(int h = 0; h < nues.size(); h++){ 
           h_nue_total->Add(nues[h]);
           std::cout << "h_nue_total : " << nues[h]->GetName() << ", " << nues[h]->Integral() << ", " << h_nue_total->Integral() << std::endl;
        }
        TH1D *h_nue_total_nolee = (TH1D*)nues_nolee[0]->Clone("nue_total_nolee");
        h_nue_total_nolee->Reset();
        for(int h = 0; h < nues_nolee.size(); h++){   
           h_nue_total_nolee->Add(nues_nolee[h]);
           std::cout << "h_nue_total_nolee : " << nues_nolee[h]->GetName() << ", " << nues_nolee[h]->Integral() << ", " << h_nue_total_nolee->Integral() << std::endl;
        }

        std::cout << "nues tot size, nues no lee tot size = " << nues.size() << ", " << nues_nolee.size() << std::endl;
        std::cout << "h_nue_total, h_nue_total_nolee  = " << h_nue_total->Integral() << ", " << h_nue_total_nolee->Integral() << std::endl;
 
        TH1D *h_numu = (TH1D*)numus[0]->Clone("numu_total");
        std::cout << "add the MC histos 2" << std::endl;
        for(int h = 1; h < numus.size(); h++) h_numu->Add(numus[h]);


        std::cout << "A" << std::endl;
        std::vector<double> input_vec_mc_total, input_vec_mc_total_nolee;
        for(int b = 1; b < h_nue_total->GetNbinsX()+1; b++) input_vec_mc_total.push_back(h_nue_total->GetBinContent(b));
        for(int b = 1; b < h_nue_total_nolee->GetNbinsX()+1; b++) input_vec_mc_total_nolee.push_back(h_nue_total_nolee->GetBinContent(b));
        for(int b = 1; b < h_numu->GetNbinsX()+1; b++){
	  input_vec_mc_total.push_back(h_numu->GetBinContent(b));
	  input_vec_mc_total_nolee.push_back(h_numu->GetBinContent(b));
        }
        
        //populate vector with the MC spectrum of the nue and numu channels          
        std::vector<double> input_vec_mc, input_vec_mc_nolee;
        for(int h=0; h < nues.size(); h++ ){
          for(int b = 1; b < h_nue_total->GetNbinsX()+1; b++) input_vec_mc.push_back(nues[h]->GetBinContent(b));
        }
        for(int h=0; h < nues_nolee.size(); h++ ){
          for(int b = 1; b < h_nue_total_nolee->GetNbinsX()+1; b++) input_vec_mc_nolee.push_back(nues_nolee[h]->GetBinContent(b));
        }

        for(int h=0; h < numus.size(); h++ ){
        for(int b = 1; b < h_numu->GetNbinsX()+1; b++){ 
          input_vec_mc.push_back(numus[h]->GetBinContent(b));
          input_vec_mc_nolee.push_back(numus[h]->GetBinContent(b));
	}
        }
	
        std::cout << "B" << std::endl;
        //fetch the TFile for the data and MC spectrum
        fname = "../bin/"+tag2+".SBNspec.root";
        //fname = "../bin/SBNfit_variation_plots_"+tag1+".root";
        std::cout << "B1" << std::endl;

        std::cout << "fname = " <<fname << std::endl;

        TFile *f_data = new TFile(fname,"read");
        TH1D *h_nue_fake_data;
        if(tag1.find("nue_numu") != std::string::npos ){ 
            h_nue_fake_data = (TH1D*)f_data->Get("nu_uBooNE_nue_data");
            //h_nue_fake_data = (TH1D*)f_data->Get("All_UBGenie/nu_uBooNE_nue_intrinsic_75"); 
            //std::cout << "h_nue_fake_data intrinsic = " << h_nue_fake_data->GetBinContent(2) << std::endl;
            //h_nue_fake_data->Add((TH1D*)f_data->Get("All_UBGenie/nu_uBooNE_nue_lee_75"));  
            //std::cout << "h_nue_fake_data intrinsic+lee = " << h_nue_fake_data->GetBinContent(2) << std::endl;
        }else{
            h_nue_fake_data = (TH1D*)f_data->Get("nu_uBooNE_1e0p_data");
            //h_nue_fake_data = (TH1D*)f_data->Get("All_UBGenie/nu_uBooNE_1e0p_intrinsic_75");  
            //std::cout << "h_nue_fake_data intrinsic = " << h_nue_fake_data->GetBinContent(2) << std::endl;
            //h_nue_fake_data->Add((TH1D*)f_data->Get("All_UBGenie/nu_uBooNE_1e0p_intrinsic_75"));  
            //std::cout << "h_nue_fake_data intrinsic+lee = " << h_nue_fake_data->GetBinContent(2) << std::endl;
        }
        //h_nue_fake_data->SetName("nu_uBooNE_nue_fakedata");
        //h_nue_fake_data->SetTitle("nu_uBooNE_nue_fakedata");
        TH1D *h_numu_data = (TH1D*)f_data->Get("nu_uBooNE_numu_data");
        //TH1D *h_numu_data = (TH1D*)f_data->Get("All_UBGenie/nu_uBooNE_numu_bnb_75");
        //h_numu_data->SetName("nu_uBooNE_numu_bnb_data");
        //h_numu_data->SetTitle("nu_uBooNE_numu_bnb_data");

        //fetch the TFile for the data and MC spectrum
        std::cout << "C" << std::endl;
        //populate vector with the data spectrum of the nue and numu channels          
        std::vector<double> input_vec_data;
        for(int b = 1; b < h_nue_fake_data->GetNbinsX()+1; b++) input_vec_data.push_back(h_nue_fake_data->GetBinContent(b));
        for(int b = 1; b < (h_numu_data->GetNbinsX()+1); b++) input_vec_data.push_back(h_numu_data->GetBinContent(b));
        std::cout << "D" << std::endl;
	
	//Load up our covariance matricies we calculated in example1 (we could also load up single variation ones)
	TString filename = "/uboone/data/users/wospakrk/09-07-2020/"+tag1+".SBNcovar.root";
	TString filename0 = "/uboone/data/users/wospakrk/09-07-2020/"+tag1+"_nolee.SBNcovar.root";
        std::cout << filename << std::endl;
	TFile *fsys = new TFile(filename,"read");
	TFile *fsys0 = new TFile(filename0,"read");
	TMatrixD *cov = (TMatrixD*)fsys->Get("full_covariance");
	TMatrixD *cov0 = (TMatrixD*)fsys0->Get("full_covariance");
	TMatrixD *fraccov = (TMatrixD*)fsys->Get("frac_covariance");
	TMatrixD *fraccov0 = (TMatrixD*)fsys0->Get("frac_covariance");
	TMatrixD *corr = (TMatrixD*)fsys->Get("full_correlation");
		
        TString bgname = "/uboone/data/users/wospakrk/09-07-2020/"+tag1+"_nolee.SBNspec.root";
        SBNspec bg(bgname.Data(),xml);
        TString signame = "/uboone/data/users/wospakrk/09-07-2020/"+tag1+".SBNspec.root";
        SBNspec sig(signame.Data(),xml);

        SBNchi chi_h0(bg,fraccov);	
        SBNchi chi_h1(sig,fraccov);	
	//Now we have all the necessary files, we are going to start the constraint
	//First, add the stat errors from the MC spectrum!
	
	//Create covariance matrix to store the stats error
	TMatrixD Mstat(cov->GetNcols(), cov->GetNrows());
	Mstat.Zero();
	for( int i=0; i < Mstat.GetNcols(); i++ ) Mstat(i,i) = input_vec_mc_total[i]; 
	
	//Check that matrix is symmetrix and then add them up
	//And then define the total covariance matrix in all its glory
	TMatrixD Mtotal(cov->GetNcols(), cov->GetNrows());
	TMatrixD Mtotal_nolee(cov->GetNcols(), cov->GetNrows());
	TMatrixD Mtotalfraccov(fraccov->GetNcols(), fraccov->GetNrows());
	TMatrixD Mtotalcorr(corr->GetNcols(), corr->GetNrows());
        TMatrixD Mdetsys(cov->GetNcols(), cov->GetNrows());
        TMatrixD Mdetsys0(cov->GetNcols(), cov->GetNrows());
        Mdetsys.Zero();
        Mdetsys0.Zero();
        if(detsys){ std::cout << "***FILL DETSYS*****"<< std::endl; chi_h0.FillDetSysMatrix(Mdetsys0,bg,true); chi_h1.FillDetSysMatrix(Mdetsys,sig,true);}
	
	Mtotal.Zero();
	Mtotal_nolee.Zero();
	Mtotalfraccov.Zero();
	Mtotalcorr.Zero();
	
	Mtotalfraccov = *fraccov;
	Mtotalcorr = *corr;

       	
	std::cout<<"Using syst only in covariance matrix"<<std::endl;
	//Mtotalfraccov = *fraccov + Mdetsys;
	Mtotal = *cov + Mdetsys;
	Mtotal_nolee = *cov0 + Mdetsys0;
	//Mtotal = Mdetsys;

        //get the core spectrum
        /*bg.CalcFullVector();
        std::vector<double> spectrum = bg.full_vector;
      
        for(int i=0; i < Mtotal.GetNrows(); i++){
        for(int j=0; j < Mtotal.GetNcols(); j++){
          Mtotal(i,j) = Mtotalfraccov(i,j)*spectrum[i]*spectrum[j]; 
        }
        }*/
        std::cout << "###### Mtotal ######" << std::endl;
	//Mtotal.Print();
	
	//collapse the covariance matrix to only its channels
	TMatrixD Mcol(input_vec_mc_total.size(),input_vec_mc_total.size());
	TMatrixD Mcol0(input_vec_mc_total.size(),input_vec_mc_total.size());
	TMatrixD Mcoldetsys(input_vec_mc_total.size(),input_vec_mc_total.size());
	TMatrixD Mcoldetsys0(input_vec_mc_total.size(),input_vec_mc_total.size());
	TMatrixD Mcolfrac(input_vec_mc_total.size(),input_vec_mc_total.size());
	TMatrixD Mcolcorr(input_vec_mc_total.size(),input_vec_mc_total.size());

	Mcol.Zero();
	Mcol0.Zero();
	Mcoldetsys.Zero();
	Mcoldetsys0.Zero();
	Mcolfrac.Zero();
	Mcolcorr.Zero();

	int num_channels = 2;
	std::vector<int> num_subchannels = {nues.size(),numus.size()};
	std::vector<int> num_bins = {h_nue_total_nolee->GetNbinsX(),h_numu->GetNbinsX()};
	
	CollapseSubchannels(Mtotal, Mcol, num_bins, num_channels, num_subchannels);
	CollapseSubchannels(Mdetsys, Mcoldetsys, num_bins, num_channels, num_subchannels);
	CollapseSubchannels(Mtotal_nolee, Mcol0, num_bins, num_channels, num_subchannels);
	CollapseSubchannels(Mdetsys0, Mcoldetsys0, num_bins, num_channels, num_subchannels);
	CollapseSubchannels(Mtotalfraccov, Mcolfrac, num_bins, num_channels, num_subchannels);
	CollapseSubchannels(Mtotalcorr, Mcolcorr, num_bins, num_channels, num_subchannels);
        
        for(int i=0; i < Mcol0.GetNrows(); i++){
          for(int j=0; j < Mcol0.GetNcols(); j++){
            if(i != j || i > 0) continue;
            std::cout << "i, j, Mcol0 before = " << i << ", " << j << " = " << Mcol0(i,j) << std::endl; 
            std::cout << "i, j, Mcoldetsys0 before = " << i << ", " << j << " = " << Mcoldetsys0(i,j) << std::endl; 
            std::cout << "i, j, Mcol before = " << i << ", " << j << " = " << Mcol(i,j) << std::endl; 
            std::cout << "i, j, Mcoldetsys before = " << i << ", " << j << " = " << Mcoldetsys(i,j) << std::endl; 
          }
        }

        std::cout << "====================================================== " << std::endl;
        std::cout << "====================================================== " << std::endl;

        //Mcol += Mcoldetsys;
        //Mcol0 += Mcoldetsys0;
        for(int i=0; i < Mcol0.GetNrows(); i++){
          for(int j=0; j < Mcol0.GetNcols(); j++){
            if(i != j || i > 0) continue;
            std::cout << "i, j, Mcol0 after = " << i << ", " << j << " = " << Mcol0(i,j) << std::endl; 
            std::cout << "i, j, Mcol after = " << i << ", " << j << " = " << Mcol(i,j) << std::endl; 
          }
        }

        plot_one(Mcol, h_nue_total, h_numu, "SBNfit_covariance_matrix_"+tag1+"_mc");
        plot_one(Mcol0, h_nue_total_nolee, h_numu, "SBNfit_covariance_matrix_"+tag1+"_mc");
        plot_one(Mcolfrac, h_nue_total, h_numu, "SBNfit_fractional_covariance_matrix_"+tag1+"_mc");
        plot_one(Mcolcorr, h_nue_total, h_numu, "SBNfit_correlation_matrix_"+tag1+"_mc");
        
        TMatrixD* Mcol2 = (TMatrixD*)Mcol.Clone("Mcol2");
        TMatrixD* Mcol3 = (TMatrixD*)Mcol0.Clone("Mcol3");
       	
        //Draw before constraint	
        for( int bin = 0; bin < Mcol.GetNcols(); bin++ ){
           if( bin > (h_nue_total_nolee->GetNbinsX()-1) ){ 
             int numubin = bin-h_nue_total_nolee->GetNbinsX()+1; 
             double derr = Mcol(bin,bin) > 0 ? sqrt(Mcol(bin,bin)): 0.;
             h_numu->SetBinError(numubin, derr);
            }else{
             double derr = Mcol(bin,bin) > 0 ? sqrt(Mcol(bin,bin)): 0.; 
             double derr0 = Mcol0(bin,bin) > 0 ? sqrt(Mcol0(bin,bin)): 0.; 
             h_nue_total_nolee->SetBinError(bin+1, derr0);
             h_nue_total->SetBinError(bin+1, derr);
             std::cout << "before const nue, err  = " << derr/h_nue_total->GetBinContent(bin+1) << std::endl;
             std::cout << "before const nue nolee, err  = " << derr0/h_nue_total_nolee->GetBinContent(bin+1) << std::endl;
           }	
        }

        TString cname = "variancemethod_"+tag1+"_nue_before_data_constraint";
        if(!sys)cname = tag1+"_nue_before_data_constraint_noleesyst";
        std::cout << "DrawDataAndMCSyst nue before const" << std::endl;
        DrawDataMCAndSyst(cname, h_nue_total, h_nue_total_nolee, h_nue_fake_data, var, "#nu_{e} Selection");
	
        cname = "variancemethod_"+tag1+"_numu_before_data_constraint";
        if(!sys)cname = tag1+"_numu_before_data_constraint_noleesyst";
        DrawDataMCAndSyst(cname, h_numu, h_numu, h_numu_data, var, "#nu_{#mu} Selection");

        //add numu stats error	
        for( int bin = 0; bin < Mcol.GetNcols(); bin++ ){
           if( bin > (h_nue_total_nolee->GetNbinsX()-1) ){ 
             int numubin = bin-h_nue_total_nolee->GetNbinsX()+1; 
             Mcol(bin,bin) += input_vec_data[bin];
             Mcol0(bin,bin) += input_vec_data[bin];
             (*Mcol2)(bin,bin) += input_vec_data[bin];
             (*Mcol3)(bin,bin) += input_vec_data[bin];
            }
        }
        //prep the covariance matrix textfiles
        std::ofstream covmatrixnue(Form("cov_matrix_%s%s.txt",varlabel.c_str(),chanlabel.c_str()));
	std::ofstream fraccovmatrixnue(Form("fraccov_matrix_%s%s.txt",varlabel.c_str(),chanlabel.c_str()));
        std::ofstream covmatrixnumu(Form("cov_matrix_%s_numu.txt",varlabel.c_str()));
	std::ofstream fraccovmatrixnumu(Form("fraccov_matrix_%s_numu.txt",varlabel.c_str()));

        //variance method:
        TMatrixD nuematrix(h_nue_total_nolee->GetNbinsX(), h_nue_total_nolee->GetNbinsX());
        TMatrixD nuenumumatrix(h_nue_total_nolee->GetNbinsX(), h_numu->GetNbinsX());
        TMatrixD numunuematrix(h_numu->GetNbinsX(), h_nue_total_nolee->GetNbinsX());
        TMatrixD numumatrix(h_numu->GetNbinsX(), h_numu->GetNbinsX());
        TMatrixD numumatrix2(h_numu->GetNbinsX(), h_numu->GetNbinsX());
        TMatrixD InvertedNumumatrix(h_numu->GetNbinsX(), h_numu->GetNbinsX());

        //zero out the elements
        nuematrix.Zero();
        nuenumumatrix.Zero();
        numunuematrix.Zero();
        numumatrix.Zero();
        numumatrix2.Zero();
        TMatrixD nuematrix_nolee(h_nue_total_nolee->GetNbinsX(), h_nue_total_nolee->GetNbinsX());
        TMatrixD nuenumumatrix_nolee(h_nue_total_nolee->GetNbinsX(), h_numu->GetNbinsX());
        TMatrixD numunuematrix_nolee(h_numu->GetNbinsX(), h_nue_total_nolee->GetNbinsX());
        nuematrix_nolee.Zero();
        nuenumumatrix_nolee.Zero();
        numunuematrix_nolee.Zero();
        std::vector<double> numu_vec;

	for( int i = 0; i < Mcol2->GetNrows(); i++ ){
	  for( int j = 0; j < Mcol2->GetNcols(); j++ ){
	    if( i >= (h_nue_total_nolee->GetNbinsX()) && j >= (h_nue_total_nolee->GetNbinsX()) ){ 
              //if(i==j) std::cout << "numu i, j = " << i << ", " << j << ", " << (i - (h_nue_total_nolee->GetNbinsX())) << ", " << (j - (h_nue_total_nolee->GetNbinsX())) << std::endl;
              numumatrix( i - (h_nue_total_nolee->GetNbinsX()), j - (h_nue_total_nolee->GetNbinsX())) = (*Mcol2)(i,j);
              numumatrix2( i - (h_nue_total_nolee->GetNbinsX()), j - (h_nue_total_nolee->GetNbinsX())) = (*Mcol2)(i,j);
              //if(i==j) numumatrix( i - (h_nue_total_nolee->GetNbinsX()), j - (h_nue_total_nolee->GetNbinsX())) += input_vec_data[i];
              //if(i==j) numumatrix2( i - (h_nue_total_nolee->GetNbinsX()), j - (h_nue_total_nolee->GetNbinsX())) += input_vec_data[i];
              numumatrix2( i - (h_nue_total_nolee->GetNbinsX()), j - (h_nue_total_nolee->GetNbinsX())) = (*Mcol2)(i,j);
              if(j < (Mcol2->GetNcols()-1)) covmatrixnumu << (*Mcol2)(i,j) << ", ";
              else covmatrixnumu << (*Mcol2)(i,j) << "\n";
              if(j < (Mcol2->GetNcols()-1)) fraccovmatrixnumu << (*Mcol2)(i,j)/(input_vec_mc_total[i]*input_vec_mc_total[j]) << ", ";
              else fraccovmatrixnumu << (*Mcol2)(i,j)/(input_vec_mc_total[i]*input_vec_mc_total[j]) << "\n";
              if(i==j) numu_vec.push_back(h_numu_data->GetBinContent((i - (h_nue_total_nolee->GetNbinsX()) + 1)) - h_numu->GetBinContent((i - h_nue_total_nolee->GetNbinsX()) + 1));
            } 
	    else if( i >= (h_nue_total_nolee->GetNbinsX()) && j < (h_nue_total_nolee->GetNbinsX()) ){ 
              numunuematrix( i - (h_nue_total_nolee->GetNbinsX()),j) = (*Mcol2)(i,j);
              //if(i>=j+14) std::cout << "numunue no lee i, j = " << i << ", " << j << " = " << Mcol0(i,j) << std::endl;
              numunuematrix_nolee( i - (h_nue_total_nolee->GetNbinsX()),j) = (*Mcol3)(i,j);
            }
	    else if( i < (h_nue_total_nolee->GetNbinsX()) && j >= (h_nue_total_nolee->GetNbinsX()) ){ 
              nuenumumatrix(i,(j-(h_nue_total_nolee->GetNbinsX()))) = (*Mcol2)(i,j);
              //if(i<j+14) std::cout << "nuenumu i, j = " << i << ", " << j << " = " << Mcol0(i,j) << std::endl;
              nuenumumatrix_nolee(i,(j-(h_nue_total_nolee->GetNbinsX()))) = (*Mcol3)(i,j);
            }
	    else{  
              if(j < (Mcol2->GetNcols()-h_nue_total_nolee->GetNbinsX()-1)) covmatrixnue << (*Mcol2)(i,j) << ", ";
              else covmatrixnue << (*Mcol2)(i,j) << "\n";
              if(j < (Mcol2->GetNcols()-h_nue_total_nolee->GetNbinsX()-1)) fraccovmatrixnue << (*Mcol2)(i,j)/(input_vec_mc_total[i]*input_vec_mc_total[j]) << ", ";
              else fraccovmatrixnue << (*Mcol2)(i,j)/(input_vec_mc_total[i]*input_vec_mc_total[j]) << "\n";
              nuematrix(i,j) = (*Mcol2)(i,j);
              nuematrix_nolee(i,j) = (*Mcol3)(i,j);
              if(i==j) std::cout << "nue i, j = " << i << ", " << j << ", " << nuematrix(i,j) << std::endl;
              if(i==j) std::cout << "nue nolee i, j = " << i << ", " << j << ", " << nuematrix_nolee(i,j) << std::endl;
            }
	  }
	}
        
        covmatrixnue.close(); 
        fraccovmatrixnue.close(); 
        covmatrixnumu.close(); 
        fraccovmatrixnumu.close(); 

	InvertedNumumatrix = numumatrix2;
	InvertedNumumatrix.Zero();
	TDecompSVD svdnumu( numumatrix2 );

	if (!svdnumu.Decompose() ) {
	  std::cout << "Decomposition failed, matrix singular ?" << std::endl;
	}else{
	  InvertedNumumatrix = numumatrix2.Invert();
	}

        TMatrixD C(h_nue_total->GetNbinsX(),h_nue_total->GetNbinsX());
        C.Mult(nuenumumatrix,InvertedNumumatrix);
        TMatrixD C_nolee(h_nue_total->GetNbinsX(),h_nue_total->GetNbinsX());
        C_nolee.Mult(nuenumumatrix_nolee,InvertedNumumatrix);
        TMatrixD D(h_nue_total->GetNbinsX(),h_nue_total->GetNbinsX());
        D.Mult(C,numunuematrix);
        TMatrixD D_nolee(h_nue_total->GetNbinsX(),h_nue_total->GetNbinsX());
        D_nolee.Mult(C_nolee,numunuematrix_nolee);
        TMatrixD E(h_nue_total->GetNbinsX(),h_nue_total->GetNbinsX());
        E.Mult(nuenumumatrix,InvertedNumumatrix);
        TMatrixD E_nolee(h_nue_total->GetNbinsX(),h_nue_total->GetNbinsX());
        E_nolee.Mult(nuenumumatrix_nolee,InvertedNumumatrix);
   
        //constrained spectrum 
        TH1D *h_nue_constrained = (TH1D*)h_nue_total->Clone("constrained_nue"); 
        TH1D *h_nue_constrained_nolee = (TH1D*)h_nue_total_nolee->Clone("constrained_nue_nolee");
        std::vector<double> input_nue_constrained;
        std::vector<double> input_nue_constrained_nolee;
 
        for( int binx = 0; binx < h_nue_total_nolee->GetNbinsX(); binx++ ){
        double scaled = 0.;
        double scaled_nolee = 0.;
        for( int biny = 0; biny < h_nue_total_nolee->GetNbinsX(); biny++ ){
          scaled += E(binx,biny) * numu_vec[biny];
          scaled_nolee += E_nolee(binx,biny) * numu_vec[biny];
          //std::cout << "E_"<<binx<<","<<biny<<"*numu_"<<biny<<" = " << E(binx,biny) << " * " << numu_vec[biny+1] << " = " << scaled << std::endl;
          //std::cout << "E_nolee"<<binx<<","<<biny<<"*numu_"<<biny<<" = " << E_nolee(binx,biny) << " * " << numu_vec[biny+1] << " = " << scaled_nolee << std::endl;
        }
        if(binx==0) std::cout << "scaled = " << scaled << std::endl;
        if(binx==0) std::cout << "scaled no lee= " << scaled_nolee << std::endl;
        double constnue =  h_nue_total->GetBinContent(binx+1) + scaled;
        double constnue_nolee =  h_nue_total_nolee->GetBinContent(binx+1) + scaled_nolee;
        h_nue_constrained->SetBinContent(binx+1,constnue);
        h_nue_constrained_nolee->SetBinContent(binx+1,constnue_nolee);
        input_nue_constrained_nolee.push_back(constnue_nolee);
        input_nue_constrained.push_back(constnue);
        }

        TMatrixD constnuematrix(h_nue_total_nolee->GetNbinsX(), h_nue_total_nolee->GetNbinsX());
        TMatrixD constnuematrix_nolee(h_nue_total_nolee->GetNbinsX(), h_nue_total_nolee->GetNbinsX());
	for( int i = 0; i < h_nue_total->GetNbinsX(); i++ ){
	  for( int j = 0; j < h_nue_total->GetNbinsX(); j++ ){
            constnuematrix(i,j) = nuematrix(i,j) - D(i,j);
            constnuematrix_nolee(i,j) = nuematrix_nolee(i,j) - D_nolee(i,j);
            //if(i==0 && i==j) std::cout << i << " constnuematrix = " << scaled << std::endl;
            //if(i==0 && i==j) std::cout << i << " constnuematrix_nolee  = " << scaled_nolee << std::endl;
          }
        }
        for( int bin = 0; bin < constnuematrix.GetNrows(); bin++ ){
             double derr = constnuematrix(bin,bin) > 0 ? sqrt(constnuematrix(bin,bin)): 0.;
             double derr_nolee = constnuematrix_nolee(bin,bin) > 0 ? sqrt(constnuematrix_nolee(bin,bin)): 0.;
             double onlysyserr = sqrt(fabs(derr*derr - h_nue_constrained->GetBinContent(bin+1))/(h_nue_constrained->GetBinContent(bin+1)*h_nue_constrained->GetBinContent(bin+1)));
             h_nue_constrained->SetBinError(bin+1,derr);
             double onlysyserr_nolee = sqrt(fabs(derr*derr - h_nue_constrained_nolee->GetBinContent(bin+1))/(h_nue_constrained_nolee->GetBinContent(bin+1)*h_nue_constrained_nolee->GetBinContent(bin+1)));
             h_nue_constrained_nolee->SetBinError(bin+1,derr_nolee);
             std::cout << "nue scaled bin, err  = " << bin+1 << ", " << derr/h_nue_constrained->GetBinContent(bin+1) << std::endl;
             std::cout << "bin, h_nue_constrained, h_nue_constrained_err = " << bin+1 << ", " << h_nue_constrained->GetBinContent(bin+1) << ", " << h_nue_constrained->GetBinError(bin+1) << std::endl;
             std::cout << "bin, h_nue_constrained_nolee, h_nue_constrained_err_nolee = " << bin+1 << ", " << h_nue_constrained_nolee->GetBinContent(bin+1) << ", " << h_nue_constrained_nolee->GetBinError(bin+1) << std::endl;
        }
       
        std::vector<double> input_vec_scaled, input_vec_scaled_nolee;
        for(int b = 1; b < h_nue_constrained->GetNbinsX()+1; b++) input_vec_scaled.push_back(h_nue_constrained->GetBinContent(b));
        for(int b = 1; b < h_nue_constrained_nolee->GetNbinsX()+1; b++) input_vec_scaled_nolee.push_back(h_nue_constrained_nolee->GetBinContent(b));

        cname = "variancemethod_"+tag1+"_scaled_nue";
        if(!sys)cname = "variancemethod_"+tag1+"_scaled_nue_noleesyst";
        DrawMCAndSyst(cname, h_nue_constrained, var, "#nu_{e} Selection");	
        cname = "variancemethod_"+tag1+"_univ_overlay_nue";
        if(!sys)cname = "variancemethod_"+tag1+"_univ_overlay_nue_noleesyst";
        std::cout << "Draw Data and MC Syst nue" << std::endl;
        DrawDataMCAndSyst(cname, h_nue_constrained, h_nue_constrained_nolee, h_nue_fake_data, var, "#nu_{e} Selection");	
        cname = "variancemethod_"+tag1+"_overlay_nue";
        if(!sys)cname = "variancemethod_"+tag1+"_overlay_nue_noleesyst";
        DrawMCAndSystOverlay(cname, h_nue_total, h_nue_constrained, var, "#nu_{e} Selection");	
        cname = "variancemethod_"+tag1+"_overlay_nue_nolee";
        if(!sys)cname = "variancemethod_"+tag1+"_overlay_nue_nolee_noleesyst";
        DrawMCAndSystOverlay(cname, h_nue_total_nolee, h_nue_constrained_nolee, var, "#nu_{e} Selection");
	
        //create root file containing the constrained matrix
        //also fill in the txt files with the results
	TFile * fconstr = new TFile(Form("variancemethod_constrained_%s.SBNcovar.root",tag1.c_str()),"recreate");	
        TMatrixD ConstMatrix(h_nue_total->GetNbinsX(),h_nue_total->GetNbinsX()); 
        TMatrixD FracConstMatrix(h_nue_total->GetNbinsX(),h_nue_total->GetNbinsX()); 
        TMatrixD CorrConstMatrix(h_nue_total->GetNbinsX(),h_nue_total->GetNbinsX());
        std::ofstream covmatrix_constrained(Form("constrained_cov_matrix_%s%s.txt",varlabel.c_str(),chanlabel.c_str()));
	std::ofstream fraccovmatrix_constrained(Form("constrained_fraccov_matrix_%s%s.txt",varlabel.c_str(),chanlabel.c_str()));

        for( int i = 0; i < ConstMatrix.GetNrows(); i++ ){
          for( int j = 0; j < ConstMatrix.GetNcols(); j++ ){
            ConstMatrix(i,j) = constnuematrix(i,j);
            if(j<ConstMatrix.GetNcols()-1) covmatrix_constrained << constnuematrix(i,j) << ", ";
            else covmatrix_constrained << constnuematrix(i,j) << "\n";
            FracConstMatrix(i,j) = constnuematrix(i,j)/(input_vec_scaled[i]*input_vec_scaled[j]);
            if(j<ConstMatrix.GetNcols()-1) fraccovmatrix_constrained << constnuematrix(i,j)/(input_vec_scaled[i]*input_vec_scaled[j]) << ", ";
            else fraccovmatrix_constrained << constnuematrix(i,j)/(input_vec_scaled[i]*input_vec_scaled[j]) << "\n";
            CorrConstMatrix(i,j) = constnuematrix(i,j)/(sqrt(constnuematrix(i,i))*sqrt(constnuematrix(j,j)));
            //if(j<ConstMatrix.GetNcols()) corrcovmatrix_constrained << constnuematrix(i,j)/(input_vec_scaled[i]*input_vec_scaled[j]) << ", ";
            //else corrcovmatrix_constrained << constnuematrix(i,j)/(input_vec_scaled[i]*input_vec_scaled[j]) << "\n";
          }
        }
        covmatrix_constrained.close();
        fraccovmatrix_constrained.close();

        (TMatrixD*)ConstMatrix.Write("full_covariance");
        (TMatrixD*)FracConstMatrix.Write("frac_covariance");
        (TMatrixD*)CorrConstMatrix.Write("full_correlation");
    
        plot_one(ConstMatrix, h_nue_total, h_numu, "SBNfit_constrained_covariance_matrix_variancemethod_"+tag1+"_mc");
        plot_one(FracConstMatrix, h_nue_total, h_numu, "SBNfit_constrained_fractional_covariance_matrix_variancemethod_"+tag1+"_mc");
        plot_one(CorrConstMatrix, h_nue_total, h_numu, "SBNfit_constrained_correlation_matrix_variancemethod_"+tag1+"_mc");

        fconstr->Close();	
        //create root file containing the lee
        std::cout << "write TH1D h_nue_scaled" << std::endl;

	TFile * fsignal = new TFile(Form("variancemethod_constrained_%s_signal.SBNspec.root",tag1.c_str()),"recreate");	
        (TH1D*)h_nue_constrained->Write("nu_uBooNE_nue_sig"); 
        fsignal->Close();
	
	TFile * fbkg = new TFile(Form("variancemethod_constrained_%s_bkg.SBNspec.root",tag1.c_str()),"recreate");
        (TH1D*)h_nue_constrained_nolee->Write("nu_uBooNE_nue_sig"); 
        fbkg->Close();
	
	return 0;
}
