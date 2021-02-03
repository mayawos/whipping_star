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


TH1D* GetTotalError(std::string name, TMatrixD errMatrix, TH1D* histo );
void DrawMCAndSyst(TString cName, TH1D* MC, std::string var, std::string title);
void DrawDataMCAndSyst(TString cName, TH1D* MC1, TH1D* MC2, TH1D* data, std::string var, std::string title);
void DrawMCAndSystOverlay(TString cName, TH1D* MC1, TH1D* MC2, std::string var, std::string title);
void CollapseSubchannels(TMatrixD & M, TMatrixD & Mc, std::vector<int> num_bins, int num_channels, std::vector<int> num_subchannels);
void Draw_Stacked(TObjArray histos, TH1D *data, TString samplename, TPad *pad = 0, bool normalised = false, std::string stacktitle = "", std::string var = "" );
void plot_one(TMatrixD matrix, TH1D *h_nue, TH1D *h_numu, std::string tag);
double CalcChiNumu(TH1D *h_data, TH1D *h_mc, TMatrixD inv_matrix, TH1D *h_nue);
double CalcChiNue(TH1D *h_data, TH1D *h_mc, TMatrixD inv_matrix);
double CalcChi(std::vector<double> data, std::vector<double> mc, TMatrixD inv_matrix);
void DrawUniverses(std::vector<TH1D*> histos, std::string name, std::string title, std::string var );

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
        bool detsys = false;
	int mass_start = -1;

	const struct option longopts[] =
	{
		{"xml", 		required_argument, 	0, 'x'},
		{"tag1", 		required_argument, 	0, 't'},
		{"tag2", 		required_argument, 	0, 'u'},
		{"width", 		required_argument, 	0, 'w'},
		{"var", 		required_argument, 	0, 'v'},
		{"sys",	                no_argument,            0, 's'},
		{"detsys",	        no_argument,            0, 'd'},
		{"part", required_argument,0,'p'},
		{0,			no_argument, 		0,  0},
	};

	while(iarg != -1)
	{
		iarg = getopt_long(argc,argv, "x:t:u:w:v:d:s:cp:g", longopts, &index);

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
	
	
        std::vector<TH1D*> numus, nuesnp, nuesnp_nolee, nues0p, nues0p_nolee;
        //fetch the TFile for the MC spectrum
        TString fname = "/uboone/data/users/wospakrk/09-07-2020/"+tag1+".SBNspec.root";
        std::cout << fname << std::endl; 
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
	  if( name.find("nu_uBooNE_nue_") != std::string::npos ){ 
            if( name.find("lee") == std::string::npos ) nuesnp_nolee.push_back(h);
            nuesnp.push_back(h);
	  }else if( name.find("nu_uBooNE_1e0p") != std::string::npos /*|| name.find("nu_uBooNE_1e0p_") != std::string::npos*/ ){ 
            if( name.find("lee") == std::string::npos ) nues0p_nolee.push_back(h);
            nues0p.push_back(h);
	  }else if( name.find("nu_uBooNE_numu") != std::string::npos ){ 
            numus.push_back(h);
	  }/*else if( name.find("nu_uBooNE_numu2") != std::string::npos ){ 
            numus2.push_back(h);
          }*/
	}
        std::cout << "add the MC histos" << std::endl;
        //now add the MC histos
        std::cout << "add the MC histos 1" << std::endl;

        TH1D *h_nue_total_np = (TH1D*)nuesnp[0]->Clone("nue_total_np");
        h_nue_total_np->Reset();
        TH1D *h_nue_total_0p = (TH1D*)nues0p[0]->Clone("nue_total_0p");
        h_nue_total_0p->Reset();
        for(int h = 0; h < nuesnp.size(); h++){ 
           h_nue_total_np->Add(nuesnp[h]);
           h_nue_total_0p->Add(nues0p[h]);
           std::cout << "h_nue_total : " << nuesnp[h]->GetName() << ", " << nuesnp[h]->Integral() << ", " << h_nue_total_np->Integral() << std::endl;
        }
        TH1D *h_nue_total_np_nolee = (TH1D*)nuesnp_nolee[0]->Clone("nue_total_np_nolee");
        h_nue_total_np_nolee->Reset();
        TH1D *h_nue_total_0p_nolee = (TH1D*)nues0p_nolee[0]->Clone("nue_total_0p_nolee");
        h_nue_total_0p_nolee->Reset();
        for(int h = 0; h < nuesnp_nolee.size(); h++){   
	  h_nue_total_np_nolee->Add(nuesnp_nolee[h]);
	  h_nue_total_0p_nolee->Add(nues0p_nolee[h]);
	  std::cout << "h_nue_total_nolee : " << nuesnp_nolee[h]->GetName() << ", " << nues0p_nolee[h]->Integral() << ", " << h_nue_total_np_nolee->Integral() << std::endl;
        }
	
        //std::cout << "nues tot size, nues no lee tot size = " << nuesnp.size() << ", " << nuesnp_nolee.size() << std::endl;
        //std::cout << "numus, numus2 = " << numus.size() << ", " << numus2.size() << std::endl;
	
        TH1D *h_numu_total = (TH1D*)numus[0]->Clone("numu1_total");
        for(int h = 1; h < numus.size(); h++) h_numu_total->Add(numus[h]);
        //TH1D *h_numu2_total = (TH1D*)numus2[0]->Clone("numu2_total");
        //for(int h = 1; h < numus2.size(); h++) h_numu2_total->Add(numus2[h]);
        std::cout << "add the MC histos 2, " << std::endl;
	
        std::cout << "A" << std::endl;
        std::vector<double> input_vec_mc_total, input_vec_mc_total_nolee;
        for(int b = 1; b < h_nue_total_np->GetNbinsX()+1; b++) input_vec_mc_total.push_back(h_nue_total_np->GetBinContent(b));
        for(int b = 1; b < h_nue_total_0p->GetNbinsX()+1; b++) input_vec_mc_total.push_back(h_nue_total_0p->GetBinContent(b));
        for(int b = 1; b < h_numu_total->GetNbinsX()+1; b++) input_vec_mc_total.push_back(h_numu_total->GetBinContent(b));
        //for(int b = 1; b < h_numu2_total->GetNbinsX()+1; b++) input_vec_mc_total.push_back(h_numu2_total->GetBinContent(b));
	
        for(int b = 1; b < h_nue_total_np_nolee->GetNbinsX()+1; b++) input_vec_mc_total_nolee.push_back(h_nue_total_np_nolee->GetBinContent(b));
        for(int b = 1; b < h_nue_total_0p_nolee->GetNbinsX()+1; b++) input_vec_mc_total_nolee.push_back(h_nue_total_0p_nolee->GetBinContent(b));
        for(int b = 1; b < h_numu_total->GetNbinsX()+1; b++) input_vec_mc_total_nolee.push_back(h_numu_total->GetBinContent(b));
        //for(int b = 1; b < h_numu2_total->GetNbinsX()+1; b++) input_vec_mc_total_nolee.push_back(h_numu2_total->GetBinContent(b));
	
        //std::cout << "h_numu_total->GetNbinsX(), h_numu2_total->GetNbinsX() = " << h_numu_total->GetNbinsX() << ", " << h_numu2_total->GetNbinsX() << std::endl;
        for(int b = 1; b < h_numu_total->GetNbinsX()+1; b++) std::cout << "h_numu_total->GetBinContent(b)" << h_numu_total->GetBinContent(b) << std::endl;
        //for(int b = 1; b < h_numu2_total->GetNbinsX()+1; b++) std::cout << "h_numu2_total->GetBinContent(b)" << h_numu2_total->GetBinContent(b) << std::endl;
        
        //populate vector with the MC spectrum of the nue and numu channels          
        std::vector<double> input_vec_mc, input_vec_mc_nolee;
        for(int h=0; h < nuesnp.size(); h++ ){
          for(int b = 1; b < h_nue_total_np->GetNbinsX()+1; b++) input_vec_mc.push_back(nuesnp[h]->GetBinContent(b));
        }
        for(int h=0; h < nues0p.size(); h++ ){
          for(int b = 1; b < h_nue_total_0p->GetNbinsX()+1; b++) input_vec_mc.push_back(nues0p[h]->GetBinContent(b));
        }
        for(int h=0; h < nuesnp_nolee.size(); h++ ){
          for(int b = 1; b < h_nue_total_np_nolee->GetNbinsX()+1; b++) input_vec_mc_nolee.push_back(nuesnp_nolee[h]->GetBinContent(b));
        }
        for(int h=0; h < nues0p_nolee.size(); h++ ){
          for(int b = 1; b < h_nue_total_0p_nolee->GetNbinsX()+1; b++) input_vec_mc_nolee.push_back(nues0p_nolee[h]->GetBinContent(b));
        }
	
        for(int h=0; h < numus.size(); h++ ){
          for(int b = 1; b < h_numu_total->GetNbinsX()+1; b++){ 
            input_vec_mc.push_back(numus[h]->GetBinContent(b));
            input_vec_mc_nolee.push_back(numus[h]->GetBinContent(b));
	  }
        }
	
	//collect all nue
	TH1D *h_nue_total = new TH1D("all_nue","all_nue",h_nue_total_np->GetNbinsX()*2,0,h_nue_total_np->GetNbinsX()*2);
	TH1D *h_nue_total_nolee = new TH1D("all_nue_nolee","all_nue_nolee",h_nue_total_np->GetNbinsX()*2,0,h_nue_total_np->GetNbinsX()*2);
	TH1D *h_numu = new TH1D("nu_uBooNE_numu1_bnb","nu_uBooNE_numu1_bnb",h_numu_total->GetNbinsX(),0.15,1.55);
        std::cout << "B" << std::endl;
        //fetch the TFile for the data and MC spectrum
        fname = "/uboone/data/users/wospakrk/09-07-2020/"+tag2+".SBNspec.root";
        //fname = "../bin/SBNfit_variation_plots_nue_1e0p_numu_reco_e_H1_mc_fakedata5_sig_allgenie.root";
        TFile *f_data = new TFile(fname,"read");
        TH1D *h_nue_fake_data = new TH1D("nu_uBooNE_nue_fakedata","nu_uBooNE_nue_fakedata",2*h_nue_total_np->GetNbinsX(),0,2*h_nue_total_np->GetNbinsX());
        //TH1D *h_numu_data = new TH1D("nu_uBooNE_numu_fakedata","nu_uBooNE_numu_fakedata",h_numu_total->GetNbinsX(),0.15,1.55);
        TH1D *h_nue_np_fake_data = (TH1D*)f_data->Get("nu_uBooNE_nue_data");
        TH1D *h_nue_0p_fake_data = (TH1D*)f_data->Get("nu_uBooNE_1e0p_data");
        TH1D *h_numu_data = (TH1D*)f_data->Get("nu_uBooNE_numu_data");
        //TH1D *h_nue_fake_data = new TH1D("nu_uBooNE_nue_univ","nu_uBooNE_nue_univ",2*h_nue_total_np->GetNbinsX(),0,2*h_nue_total_np->GetNbinsX());
        //TH1D *h_numu_data = new TH1D("nu_uBooNE_numu_univ","nu_uBooNE_numu_univ",h_numu_total->GetNbinsX(),0.15,1.55);
        //TH1D *h_nue_np_fake_data = (TH1D*)f_data->Get("All_Genie/nu_uBooNE_nue_intrinsic_72");
        //TH1D *h_nue_0p_fake_data = (TH1D*)f_data->Get("All_Genie/nu_uBooNE_1e0p_intrinsic_72");
        //TH1D *h_numu_data = (TH1D*)f_data->Get("All_Genie/nu_uBooNE_numu_bnb_72");
        h_numu_data->SetName("nu_uBooNE_numu1_bnb_data");
        h_numu_data->SetTitle("nu_uBooNE_numu1_bnb_data");
        //populate vector with the data spectrum of the nue and numu channels          
        std::vector<double> input_vec_data;
        for(int b = 1; b < h_nue_np_fake_data->GetNbinsX()+1; b++) input_vec_data.push_back(h_nue_np_fake_data->GetBinContent(b));
        for(int b = 1; b < (h_nue_0p_fake_data->GetNbinsX()+1); b++) input_vec_data.push_back(h_nue_0p_fake_data->GetBinContent(b));
        for(int b = 1; b < (h_numu_data->GetNbinsX()+1); b++) input_vec_data.push_back(h_numu_data->GetBinContent(b));
        for(int b=h_nue_total->GetNbinsX(); b < input_vec_data.size(); b++ ){
          std::cout << "h_numu data: " << h_numu_data->GetBinContent((b-h_nue_total->GetNbinsX())+1) << ", " << input_vec_data[b] << std::endl; 
        }
        //for(int b = 1; b < (h_numu2_data->GetNbinsX()+1); b++) input_vec_data.push_back(h_numu2_data->GetBinContent(b));
        //for(int b = 1; b < (h_numu2_data->GetNbinsX()+1); b++) std::cout << "h_numu2_data->GetBinContent(b)) " << h_numu2_data->GetBinContent(b) << std::endl;
	
        //fetch the TFile for the data and MC spectrum
        std::cout << "C" << std::endl;
	
        for(int b=0; b < h_nue_total->GetNbinsX(); b++ ){
          h_nue_total->SetBinContent(b+1,input_vec_mc_total[b]);
          std::cout << "bin, nue data = " << b+1 << ", " << h_nue_total->GetBinContent(b+1) << std::endl;
        }
        for(int b=0; b < h_nue_total->GetNbinsX(); b++ ){
          h_nue_total_nolee->SetBinContent(b+1,input_vec_mc_total_nolee[b]);
          std::cout << "bin, nue no lee data = " << b+1 << ", " << h_nue_total_nolee->GetBinContent(b+1) << std::endl;
        }
        for(int b=h_nue_total->GetNbinsX(); b < input_vec_mc_total.size(); b++ ){
          h_numu->SetBinContent((b-h_nue_total->GetNbinsX())+1,input_vec_mc_total[b]);
          std::cout << "h_numu: " << h_numu->GetBinContent((b-h_nue_total->GetNbinsX())+1) << ", " << input_vec_mc_total[b] << std::endl; 
        }
        for(int b=0; b < h_nue_total->GetNbinsX(); b++ ){
          h_nue_fake_data->SetBinContent(b+1,input_vec_data[b]);
        }
        for(int b=h_nue_total->GetNbinsX(); b < input_vec_data.size(); b++ ){
          h_numu_data->SetBinContent((b-h_nue_total->GetNbinsX())+1,input_vec_data[b]);
          std::cout << "h_numu data: " << h_numu_data->GetBinContent((b-h_nue_total->GetNbinsX())+1) << ", " << input_vec_data[b] << std::endl; 
        }
        std::cout << "input_vec_data.size() = " << input_vec_data.size() << std::endl;
	
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
	TMatrixD *colcov = (TMatrixD*)fsys->Get("collapsed_covariance");
	TMatrixD *colfraccov = (TMatrixD*)fsys->Get("collapsed_frac_covariance");
	TMatrixD *colcorr = (TMatrixD*)fsys->Get("collapsed_correlation");
	
        TString bgname = "/uboone/data/users/wospakrk/09-07-2020/"+tag1+"_nolee.SBNspec.root";
        SBNspec bg(bgname.Data(),xml);
        TString signame = "/uboone/data/users/wospakrk/09-07-2020/"+tag1+"_nolee.SBNspec.root";
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
	
	//collapse the covariance matrix to only its channels
	TMatrixD Mcol(input_vec_mc_total.size(),input_vec_mc_total.size());
	TMatrixD Mcol0(input_vec_mc_total.size(),input_vec_mc_total.size());
	TMatrixD Mcolfrac(input_vec_mc_total.size(),input_vec_mc_total.size());
	TMatrixD Mcolcorr(input_vec_mc_total.size(),input_vec_mc_total.size());
	
	Mcol.Zero();
	Mcol0.Zero();
	Mcolfrac.Zero();
	Mcolcorr.Zero();
	
        Mcol = *colcov;
        Mcolfrac = *colfraccov;
        Mcolcorr = *colcorr;
	
	int num_channels = 3;
	std::vector<int> num_subchannels = {nuesnp.size(),nues0p.size(),numus.size()};
	std::vector<int> num_bins = {h_nue_total_np_nolee->GetNbinsX(),h_nue_total_0p_nolee->GetNbinsX(),h_numu_total->GetNbinsX()};
	
	CollapseSubchannels(Mtotal, Mcol, num_bins, num_channels, num_subchannels);
	CollapseSubchannels(Mtotal_nolee, Mcol0, num_bins, num_channels, num_subchannels);
	CollapseSubchannels(Mtotalfraccov, Mcolfrac, num_bins, num_channels, num_subchannels);
        CollapseSubchannels(Mtotalcorr, Mcolcorr, num_bins, num_channels, num_subchannels);
        
        plot_one(Mcol, h_nue_total, h_numu, "SBNfit_covariance_matrix_variancemethod_"+tag1+"_mc");
        plot_one(Mcol0, h_nue_total_nolee, h_numu, "SBNfit_covariance_matrix_noleevariancemethod_"+tag1+"_mc");
        plot_one(Mcolfrac, h_nue_total_nolee, h_numu, "SBNfit_fractional_covariance_matrix_variancemethod_"+tag1+"_mc");
        plot_one(Mcolcorr, h_nue_total_nolee, h_numu, "SBNfit_correlation_matrix_variancemethod_"+tag1+"_mc");
        
        //for( int i=0; i < Mcol.GetNcols(); i++ ){
	// Mcol(i,i) += input_vec_mc_total[i];
        //}
	
        TMatrixD* Mcol2 = (TMatrixD*)Mcol.Clone("Mcol2");
        TMatrixD* Mcol3 = (TMatrixD*)Mcol0.Clone("Mcol3");
       	
        for(int b=h_nue_total->GetNbinsX(); b < input_vec_data.size(); b++ ){
          std::cout << "###h_numu: " << h_numu->GetBinContent((b-h_nue_total->GetNbinsX())+1) << ", " << input_vec_mc_total[b] << std::endl; 
          std::cout << "###h_numu data: " << h_numu_data->GetBinContent((b-h_nue_total->GetNbinsX())+1) << ", " << input_vec_data[b] << std::endl; 
        }
        //Draw before constraint	
        for( int bin = 0; bin < Mcol2->GetNcols(); bin++ ){
	  if( bin > (h_nue_total_nolee->GetNbinsX()-1) ){ 
	    int numubin = bin-h_nue_total_nolee->GetNbinsX()+1; 
	    double derr = (*Mcol2)[bin][bin] > 0 ? sqrt((*Mcol2)[bin][bin]): 0.;
	    h_numu->SetBinError(numubin, derr);
	    std::cout << "###numu, err  = " << (*Mcol2)[bin][bin] << ", " << h_numu->GetBinContent(numubin+1) << ", " << derr/h_numu->GetBinContent(numubin+1) << std::endl;
	  }else{
	    double derr = (*Mcol2)[bin][bin] > 0 ? sqrt((*Mcol2)[bin][bin]): 0.; 
	    double derr_nolee = (*Mcol3)[bin][bin] > 0 ? sqrt((*Mcol3)[bin][bin]): 0.; 
	    h_nue_total_nolee->SetBinError(bin+1, derr_nolee);
	    h_nue_total->SetBinError(bin+1, derr);
	    std::cout << "###nue, err  = " << (*Mcol2)[bin][bin] << ", " << derr/h_nue_total->GetBinContent(bin+1) << std::endl;
	    std::cout << "###nue nolee, err  = " << (*Mcol3)[bin][bin] << ", " << derr_nolee/h_nue_total_nolee->GetBinContent(bin+1) << std::endl;
	  }	
	}
	
	
	TString cname = "variancemethod_"+tag1+"_nue_before_data_constraint";
	if(!sys)cname = "variancemethod_"+tag1+"_nue_before_data_constraint_noleesyst";
	std::cout << "DrawDataAndMCSyst nue before const" << std::endl;
	DrawDataMCAndSyst(cname, h_nue_total, h_nue_total_nolee, h_nue_fake_data, var, "#nu_{e} Selection");
	
	cname = "variancemethod_"+tag1+"_numu_before_data_constraint";
	if(!sys)cname = "variancemethod_"+tag1+"_numu_before_data_constraint_noleesyst";
	DrawDataMCAndSyst(cname, h_numu, h_numu, h_numu_data, var, "#nu_{#mu} Selection");
	
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
        std::cout << "h_nue_total_nolee->GetNbinsX() = " << h_nue_total_nolee->GetNbinsX() << std::endl;
	for( int i = 0; i < Mcol2->GetNcols(); i++ ){
	  for( int j = 0; j < Mcol2->GetNrows(); j++ ){
	    if( i >= (h_nue_total_nolee->GetNbinsX()) && j >= (h_nue_total_nolee->GetNbinsX()) ){ 
              if(i==j) std::cout << "numu i, j = " << i << ", " << j << ", " << (i - (h_nue_total_nolee->GetNbinsX())) << ", " << (j - (h_nue_total_nolee->GetNbinsX())) << std::endl;
              numumatrix( i - (h_nue_total_nolee->GetNbinsX()), j - (h_nue_total_nolee->GetNbinsX())) = (*Mcol2)(i,j);
              numumatrix2( i - (h_nue_total_nolee->GetNbinsX()), j - (h_nue_total_nolee->GetNbinsX())) = (*Mcol2)(i,j);
              //if(i==j) std::cout << "numu i, j, bin i, bin j = " << i << ", " << j << ", " << (i - h_nue_total_nolee->GetNbinsX()) + 1 << ", " << (j - h_nue_total_nolee->GetNbinsX()) + 1 << " = " << numumatrix( i - (h_nue_total_nolee->GetNbinsX()), j - (h_nue_total_nolee->GetNbinsX())) <<std::endl;
              if(i==j) numu_vec.push_back(h_numu_data->GetBinContent((i - (h_nue_total_nolee->GetNbinsX()) + 1)) - h_numu->GetBinContent((i - h_nue_total_nolee->GetNbinsX()) + 1));
            } 
	    else if( i >= (h_nue_total_nolee->GetNbinsX()) && j < (h_nue_total_nolee->GetNbinsX()) ){ 
              //if(i>=j+14) std::cout << "numunue i, j = " << i << ", " << j << ", " << (i - (h_nue_total_nolee->GetNbinsX())) << ", " << j << std::endl;
              if(i>=j+14) std::cout << "numunue i, j = " << i << ", " << j << " = " << (*Mcol2)(i,j)  << std::endl;
              numunuematrix( i - (h_nue_total_nolee->GetNbinsX()),j) = (*Mcol2)(i,j);
              //if(i>=j+14) std::cout << "numunue no lee i, j = " << i << ", " << j << " = " << (*Mcol3)(i,j) << std::endl;
              numunuematrix_nolee( i - (h_nue_total_nolee->GetNbinsX()),j) = (*Mcol3)(i,j);
            }
	    else if( i < (h_nue_total_nolee->GetNbinsX()) && j >= (h_nue_total_nolee->GetNbinsX()) ){ 
              //if(i<j+14) std::cout << "nuenumu i, j = " << i << ", " << j << ", " << (i) << ", " << (j - (h_nue_total_nolee->GetNbinsX())) << std::endl;
              if(i<j+14) std::cout << "nuenumu i, j = " << i << ", " << j << " = " << (*Mcol2)(i,j) << std::endl;
              nuenumumatrix(i,(j-(h_nue_total_nolee->GetNbinsX()))) = (*Mcol2)(i,j);
              //if(i<j+14) std::cout << "nuenumu nolee i, j = " << i << ", " << j << " = " << (*Mcol3)(i,j) << std::endl;
              nuenumumatrix_nolee(i,(j-(h_nue_total_nolee->GetNbinsX()))) = (*Mcol3)(i,j);
            }
	    else{  
              nuematrix(i,j) = (*Mcol2)(i,j);
              nuematrix_nolee(i,j) = (*Mcol3)(i,j);
              //if(i==j) std::cout << "nue i, j = " << i << ", " << j << ", " << nuematrix(i,j) << std::endl;
              //if(i==j) std::cout << "nue nolee i, j = " << i << ", " << j << ", " << nuematrix_nolee(i,j) << std::endl;
            }
	  }
	}
        std::cout << "============ nuenumumatrix ================" << std::endl;
        nuenumumatrix.Print();
        std::cout << "============ numunuematrix ================" << std::endl;
        numunuematrix.Print();
        std::cout << "===========================================" << std::endl;
        
        //std::cout << "A" << std::endl;
	InvertedNumumatrix = numumatrix2;
	InvertedNumumatrix.Zero();
	TDecompSVD svdnumu( numumatrix2 );

        //std::cout << "B" << std::endl;
	if (!svdnumu.Decompose() ) {
	  std::cout << "Decomposition failed, matrix singular ?" << std::endl;
	}else{
	  InvertedNumumatrix = numumatrix2.Invert();
	}

        TMatrixD C(28,14);
        TMatrixD C_nolee(28,14);
        TMatrixD D(28,28);
        TMatrixD D_nolee(28,28);
        TMatrixD E(28,14);
        TMatrixD E_nolee(28,14);

        C.Mult(nuenumumatrix,InvertedNumumatrix);
        C_nolee.Mult(nuenumumatrix_nolee,InvertedNumumatrix);
        D.Mult(C,numunuematrix);
        D_nolee.Mult(C_nolee,numunuematrix_nolee);
        E.Mult(nuenumumatrix,InvertedNumumatrix);
        E_nolee.Mult(nuenumumatrix,InvertedNumumatrix);

        std::cout << "number of columns nuenumumatrix == number of rows InvertedNumumatrix :  " << nuenumumatrix.GetNcols() << " == " << InvertedNumumatrix.GetNrows() << "?" << std::endl;
        //std::cout << "============ InvertedNumumatrix ================" << std::endl;
        //InvertedNumumatrix.Print();
        std::cout << "============ E ================" << std::endl;
        E.Print();
        //nuenumumatrix_nolee.Print();
        //std::cout << "====================================" << std::endl;
        E_nolee.Print();
        std::cout << "=============== D ===================" << std::endl;
        D.Print();
        //std::cout << "====================================" << std::endl;
        D_nolee.Print();
        //std::cout << "===============%%%%%=================" << std::endl;
      
        //constrained spectrum 
        TH1D *h_nue_constrained = (TH1D*)h_nue_total->Clone("constrained_nue"); 
        TH1D *h_nue_constrained_nolee = (TH1D*)h_nue_total_nolee->Clone("constrained_nue_nolee");
        std::vector<double> input_nue_constrained;
        std::vector<double> input_nue_constrained_nolee;
 
        std::cout << "E" << std::endl;
        for( int binx = 0; binx < h_nue_total_nolee->GetNbinsX(); binx++ ){
        double scaled = 0.;
        double scaled_nolee = 0.;
        for( int biny = 0; biny < h_numu->GetNbinsX(); biny++ ){
          scaled += E(binx,biny) * numu_vec[biny];
          scaled_nolee += E_nolee(binx,biny) * numu_vec[biny];
          //std::cout << "E_"<<binx<<","<<biny<<"*numu_"<<biny<<" = " << E(binx,biny) << " * " << numu_vec[biny+1] << " = " << scaled << std::endl;
          //std::cout << "E_nolee"<<binx<<","<<biny<<"*numu_"<<biny<<" = " << E_nolee(binx,biny) << " * " << numu_vec[biny+1] << " = " << scaled_nolee << std::endl;
        }
        std::cout << "scaled = " << scaled << std::endl;
        std::cout << "scaled no lee= " << scaled_nolee << std::endl;
        double constnue =  h_nue_total->GetBinContent(binx+1) + scaled;
        double constnue_nolee =  h_nue_total_nolee->GetBinContent(binx+1) + scaled_nolee;
        h_nue_constrained->SetBinContent(binx+1,constnue);
        h_nue_constrained_nolee->SetBinContent(binx+1,constnue_nolee);
        input_nue_constrained_nolee.push_back(constnue_nolee);
        input_nue_constrained.push_back(constnue);
        }
        std::cout << "F" << std::endl;

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
             //std::cout << "nue scaled bin, err  = " << bin+1 << ", " << derr/h_nue_constrained->GetBinContent(bin+1) << std::endl;
             //std::cout << "bin, h_nue_constrained, h_nue_constrained_err = " << bin+1 << ", " << h_nue_constrained->GetBinContent(bin+1) << ", " << h_nue_constrained->GetBinError(bin+1) << std::endl;
             //std::cout << "bin, h_nue_constrained_nolee, h_nue_constrained_err_nolee = " << bin+1 << ", " << h_nue_constrained_nolee->GetBinContent(bin+1) << ", " << h_nue_constrained_nolee->GetBinError(bin+1) << std::endl;
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
	TFile * fconstr = new TFile(Form("variancemethod_constrained_%s.SBNcovar.root",tag1.c_str()),"recreate");	
        TMatrixD ConstMatrix(h_nue_total->GetNbinsX(),h_nue_total->GetNbinsX()); 
        TMatrixD FracConstMatrix(h_nue_total->GetNbinsX(),h_nue_total->GetNbinsX()); 
        TMatrixD CorrConstMatrix(h_nue_total->GetNbinsX(),h_nue_total->GetNbinsX()); 
        for( int i = 0; i < ConstMatrix.GetNrows(); i++ ){
          for( int j = 0; j < ConstMatrix.GetNcols(); j++ ){
            ConstMatrix(i,j) = constnuematrix(i,j);
            FracConstMatrix(i,j) = constnuematrix(i,j)/(input_vec_scaled[i]*input_vec_scaled[j]);
            CorrConstMatrix(i,j) = constnuematrix(i,j)/(sqrt(constnuematrix(i,i))*sqrt(constnuematrix(j,j)));
          }
        }
        (TMatrixD*)ConstMatrix.Write("full_covariance");
        (TMatrixD*)FracConstMatrix.Write("frac_covariance");
        (TMatrixD*)CorrConstMatrix.Write("full_correlation");
    
        plot_one(ConstMatrix, h_nue_total, h_numu, "SBNfit_constrained_covariance_matrix_variancemethod_"+tag1+"_mc");
        plot_one(FracConstMatrix, h_nue_total, h_numu, "SBNfit_constrained_fractional_covariance_matrix_variancemethod_"+tag1+"_mc");
        plot_one(CorrConstMatrix, h_nue_total, h_numu, "SBNfit_constrained_correlation_matrix_variancemethod_"+tag1+"_mc");

        fconstr->Close();
	
        //create root file containing the constrained matrix no lee
	TFile * fconstr_nolee = new TFile(Form("variancemethod_constrained_%s_nolee.SBNcovar.root",tag1.c_str()),"recreate");	
        TMatrixD ConstMatrix_nolee(h_nue_total->GetNbinsX(),h_nue_total->GetNbinsX()); 
        TMatrixD FracConstMatrix_nolee(h_nue_total->GetNbinsX(),h_nue_total->GetNbinsX()); 
        TMatrixD CorrConstMatrix_nolee(h_nue_total->GetNbinsX(),h_nue_total->GetNbinsX()); 
        for( int i = 0; i < ConstMatrix_nolee.GetNrows(); i++ ){
          for( int j = 0; j < ConstMatrix_nolee.GetNcols(); j++ ){
            ConstMatrix_nolee(i,j) = constnuematrix_nolee(i,j);
            FracConstMatrix_nolee(i,j) = constnuematrix_nolee(i,j)/(input_vec_scaled_nolee[i]*input_vec_scaled_nolee[j]);
            CorrConstMatrix_nolee(i,j) = constnuematrix_nolee(i,j)/(sqrt(constnuematrix_nolee(i,i))*sqrt(constnuematrix_nolee(j,j)));
          }
        }
        (TMatrixD*)ConstMatrix_nolee.Write("full_covariance");
        (TMatrixD*)FracConstMatrix_nolee.Write("frac_covariance");
        (TMatrixD*)CorrConstMatrix_nolee.Write("full_correlation");
    
        plot_one(ConstMatrix_nolee, h_nue_total, h_numu, "SBNfit_constrained_covariance_matrix_variancemethod_"+tag1+"_nolee_mc");
        plot_one(FracConstMatrix_nolee, h_nue_total, h_numu, "SBNfit_constrained_fractional_covariance_matrix_variancemethod_"+tag1+"_nolee_mc");
        plot_one(CorrConstMatrix_nolee, h_nue_total, h_numu, "SBNfit_constrained_correlation_matrix_variancemethod_"+tag1+"_nolee_mc");

        fconstr_nolee->Close();	
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

void DrawUniverses(std::vector<TH1D*> histos, std::string name, std::string title, std::string var ){
  
  TCanvas can(name.c_str(),name.c_str());
  histos[0]->SetTitle(title.c_str());
  histos[0]->GetXaxis()->SetTitle(var.c_str());
  histos[0]->GetYaxis()->SetTitle("Events / 0.1 GeV");
  histos[0]->SetMaximum(1.2*histos[0]->GetMaximum());
  histos[0]->SetLineColor(1);
  histos[1]->SetLineColor(2);
  histos[2]->SetLineColor(2);
  //histos[1]->SetLineStyle(10);
  histos[2]->SetLineStyle(2);
  for(int i=0; i < 3; i++){
    histos[i]->SetLineWidth(2);
    if(i==0)histos[i]->Draw("hist");
    histos[i]->Draw("histsame");
  }
  TLegend  *legend = new TLegend(0.70,0.70,0.98,0.95); // we need different positions for the legend to not 
  // get the plot titles for the legend
  std::vector<std::string> legends_str = {"CV","Univ 0","Univ 1"};
  for(int i = 0; i < 3 ; i++){
    legend->AddEntry(histos[i],legends_str[i].c_str(),"l"); 
  }
  legend->Draw("same");
  can.Print(Form("%s_universes.pdf",name.c_str()),"pdf");
  
}

//borrow this function from SBNchi
void CollapseSubchannels(TMatrixD & M, TMatrixD & Mc, std::vector<int> num_bins, int num_channels, std::vector<int> num_subchannels){
  bool debug = true;
  //if(debug)	std::cout<<"Starting:M "<<M.GetNcols()<<" "<<M.GetNrows()<<" "<<std::endl;
  //if(debug)	std::cout<<"Starting:Mc "<<Mc.GetNcols()<<" "<<Mc.GetNrows()<<" "<<std::endl;
  
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
      
      //if(debug)std::cout<<"Diagonal! : "<<ic<<" "<<jc<<" mcol is: "<<mcol<<" mrow is: "<<mrow<<std::endl;
      
      for(int m=0; m < num_subchannels[ic]; m++){
	for(int n=0; n< num_subchannels[jc]; n++){ //For each big block, loop over all subchannels summing toGether
	  //std::cout << mrow << ", " << n << ", " << num_bins[jc] << ", " << ", " << mcol << ", " << m << ", " << m << ", " << num_bins[ic] << std::endl;
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
  
  /// ********************************* And put them back toGether! ************************ //
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
  colours.push_back(kOrange-5);
  colours.push_back(kMagenta-6);
  
  // lets open and draw the canvas 
  TCanvas *canvas;
  TString cName = Form("%s",samplename.Data());
  if(pad == 0){
    canvas = new TCanvas(cName,"Stacked Histograms", 1200., 800.);
    pad = (TPad*)canvas->cd();
	  }
  pad->cd();
  pad->SetTicks(0,0);
  pad->SetRightMargin(0.1);
  pad->SetBottomMargin(0.13);
  
  // lets take the first histoname and see if we can match a title to its which will be HS stack title
  if( stacktitle == "" ) stacktitle = ((TH1F*)histos[0])->GetTitle();
  THStack *Hs = new THStack("hs2",stacktitle.c_str());
  
  //Set Axis Units
  //Hs->GetXaxis()->SetTitle("Number of Hits in Slice");
  
  // Set up the LEGEND
  TLegend  *legend = new TLegend(0.70,0.5,0.98,0.7); // we need different positions for the legend to not 
  // get the plot titles for the legend
  for(int i = 0; i < n ; i++){
    
    TH1D *h =  (TH1D*)histos[i];
    std::cout << "histo = " << h->GetTitle() << std::endl;
    legends_str.push_back( h->GetTitle() );
    //h->Scale(1,"width");
    h->SetLineWidth(0);
    h->SetLineColor(colours[i]);
    h->SetFillColor(colours[i]);
    legend->AddEntry(h,legends_str[i].c_str(),"f"); 
    
    std::cout << "histo name =  " << h->GetTitle() << std::endl;
    std::cout << "integral =  " << h->Integral() << std::endl;
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
  
  data->SetMarkerStyle(20);
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
  pad1->SetGridx(2);// Vertical grid
  pad1->Draw(); // Draw the upper pad: pad1
  pad1->cd();   // pad1 becomes the current pad
  
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
  can.cd();  // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.35);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.25);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();   // pad2 becomes the current pad
  
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
  pad1->SetGridx(2);// Vertical grid
  pad1->Draw(); // Draw the upper pad: pad1
  pad1->cd();   // pad1 becomes the current pad
  
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
  sysErr->SetMaximum(1.5*sysErr->GetMaximum());
  sysErr->Draw("E2");
  sysErr2->Draw("E2 same");
  tmpMC->Draw("HIST same");
  tmpMC2->Draw("HIST same");
  
  
  // lower plot will be in pad
  can.cd();  // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.35);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.25);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();   // pad2 becomes the current pad
  
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
  tmpsys->SetMaximum(2.0);
  tmpsys->SetMinimum(0.0);
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
  std::cout << "nbins data, mc = " << tmpData->Integral() << ", " << tmpMC1->Integral() << std::endl;

  TH1D* tmpratio2 = (TH1D*)data->Clone("tempRatio2");
  TH1D* tmpMC2_staterr = (TH1D*)MC2->Clone("tempCV_staterr2");
  for( int bin = 1; bin < tmpMC2_staterr->GetNbinsX()+1; bin++){ tmpMC2_staterr->SetBinError(bin,sqrt(tmpMC2_staterr->GetBinContent(bin))); }
  tmpratio2->Reset();
  tmpratio2->Divide(tmpData,tmpMC2_staterr,1.0,1.0,"B");
  std::cout << "nbins data, mc nolee = " << tmpData->Integral() << ", " << tmpMC2->Integral() << std::endl;

  TH1D* tmpsys = (TH1D*)MC2->Clone("tempSysRatio");
  TH1D* tmpsys2 = (TH1D*)MC1->Clone("tempSysRatio1");

  for( int bin = 1; bin < tmpsys->GetNbinsX()+1; bin++){
    tmpsys->SetBinContent(bin,1.0);
    tmpsys->SetBinError(bin,tmpMC2->GetBinError(bin)/tmpMC2->GetBinContent(bin));
    std::cout << "tempSys MC2 " << tmpsys->GetBinError(bin) << std::endl;
    if(bin==1) std::cout << "tempSys MC2 " << tmpMC2->GetBinError(bin) << "/" << tmpMC2->GetBinContent(bin) << " = " << tmpMC2->GetBinError(bin)/tmpMC2->GetBinContent(bin) << std::endl;
    tmpsys2->SetBinContent(bin,1.0);
    tmpsys2->SetBinError(bin,tmpMC1->GetBinError(bin)/tmpMC1->GetBinContent(bin));
    std::cout << "tempSys MC1 " << tmpsys2->GetBinError(bin) << std::endl;
    //tmpratio->SetBinError(bin,0.01);
    sysErr1->SetBinContent(bin,tmpMC1->GetBinContent(bin));
    sysErr1->SetBinError(bin,tmpMC1->GetBinError(bin));
    sysErr2->SetBinContent(bin,tmpMC2->GetBinContent(bin));
    sysErr2->SetBinError(bin,tmpMC2->GetBinError(bin));
    //std::cout << "temp ratio = " << tmpratio->GetBinContent(bin) << std::endl; 
  }
  
  // Upper plot will be in pad1
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.45, 1, 1.0);
  pad1->SetBottomMargin(0.05); // Upper and lower plot are joined
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

  sysErr1->SetTitleSize(80);
  sysErr1->SetTitleFont(43);
  sysErr1->SetTitle(title.c_str());
  sysErr1->SetStats(0);
  sysErr1->GetYaxis()->SetTitleSize(25);
  sysErr1->GetYaxis()->SetTitleFont(43);
  sysErr1->GetYaxis()->SetTitle("Events");
  sysErr1->GetXaxis()->SetLabelSize(0);

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
      legend->AddEntry(tmpData,"fake data","l");
    }
  }else{
    if(std::string(cName).find("nowgt") != std::string::npos ){
      legend->AddEntry(tmpMC1,"#nu_e, Genie Wgt","l"); 
      legend->AddEntry(tmpMC2,"#nu_e, no Genie wgt","l"); 
      legend->AddEntry(tmpData,"#nu_e, post-constraint","l");
    }else{
      legend->AddEntry(tmpMC1,"#nu_e intrinsic + LEE","l"); 
      legend->AddEntry(tmpMC2,"#nu_e intrinsic","l"); 
      legend->AddEntry(tmpData,"fake data","l");
    }
  } 
  legend->Draw("same");

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
  tmpsys2->Draw("E2 same");

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

  gPad->RedrawAxis();
  gPad->Update();
 
  can.Print(cName+".pdf","pdf");

}

void plot_one(TMatrixD matrix, TH1D *h_nue, TH1D *h_numu, std::string tag){
  
  std::string dir="/uboone/data/users/wospakrk/SBNFitPlots/";
  std::vector<std::string> channel_names = {"nu uBooNE nue intrinsic","nu uBooNE numu BNB"};
  int num_channels = 2;
  int num_bins_total = matrix.GetNrows();
  std::vector<int> num_bins = {h_nue->GetNbinsX(), h_numu->GetNbinsX()}; 
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

double CalcChiNumu(TH1D *h_data, TH1D *h_mc, TMatrixD inv_matrix, TH1D *h_nue){
  
  double tchi=0.0;
  
  for( int i = 0; i < inv_matrix.GetNcols(); i++ ){
    for( int j = 0; j < inv_matrix.GetNrows(); j++ ){
      if( i >= (h_nue->GetNbinsX()) && j >= (h_nue->GetNbinsX()) ){
	int numubin_i = i-h_nue->GetNbinsX();
	int numubin_j = j-h_nue->GetNbinsX();
	//std::cout << "i, j, numubin_i, numubin_j, h_data_i, h_mc_i, cov_matrix, h_data_j,  = " << i << ", " << j << ", " << numubin_i << ", " << numubin_j << ", " << h_data->GetBinContent(numubin_i+1)-h_mc->GetBinContent(numubin_j+1) << ", " << inv_matrix(i,j) << ", " << h_data->GetBinContent(numubin_j+1)-h_mc->GetBinContent(numubin_i+1) <<std::endl;
	double tchiperbin = (h_data->GetBinContent(numubin_i+1)-h_mc->GetBinContent(numubin_i+1))*inv_matrix(i,j)*(h_data->GetBinContent(numubin_j+1)-h_mc->GetBinContent(numubin_j+1));
	tchi += tchiperbin;
      }
    }
  }
  return tchi;
  
}

double CalcChiNue(TH1D *h_data, TH1D *h_mc, TMatrixD inv_matrix){
  
  double tchi=0.0;
  
  for( int i = 0; i < inv_matrix.GetNcols(); i++ ){
    for( int j = 0; j < inv_matrix.GetNrows(); j++ ){
      //if(i==j) std::cout << "i,j = " << i << ", " << j << std::endl;
      if( i < (h_data->GetNbinsX()) && j < (h_data->GetNbinsX()) ){
	//std::cout << "i, j, i, j, h_data_i, h_mc_i, cov_matrix, h_data_j,  = " << i << ", " << j << ", " << i << ", " << j << ", " << h_data->GetBinContent(i+1)-h_mc->GetBinContent(j+1) << ", " << inv_matrix(i,j) << ", " << h_data->GetBinContent(j+1)-h_mc->GetBinContent(i+1) <<std::endl;
	double tchiperbin = (h_data->GetBinContent(i+1)-h_mc->GetBinContent(i+1))*inv_matrix(i,j)*(h_data->GetBinContent(j+1)-h_mc->GetBinContent(j+1));
	tchi += tchiperbin;
      }
    }
  }
  return tchi;
  
}

double CalcChi(std::vector<double> data, std::vector<double> mc, TMatrixD inv_matrix){
  
  double tchi=0.0;
  
  for( int i = 0; i < inv_matrix.GetNcols(); i++ ){
    for( int j = 0; j < inv_matrix.GetNrows(); j++ ){
      double tchiperbin = (data[i]-mc[i])*inv_matrix(i,j)*(data[j]-mc[j]);
      tchi += tchiperbin;
    }
  }
  return tchi;
  
}
