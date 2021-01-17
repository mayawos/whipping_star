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


void plot_one(TMatrixD matrix, std::string tag);

int main(int argc, char* argv[])
{
	std::string xml = "example.xml";
	std::string tag = "example";

	const struct option longopts[] =
	{
		{"xml", 		required_argument, 	0, 'x'},
		{"tag", 		required_argument, 	0, 't'},
		{0,   		        no_argument,     	0,  0},
	};

	int iarg = 0;
	opterr=1;
	int index;

	while(iarg != -1)
	{
		iarg = getopt_long(argc,argv, "x:t:", longopts, &index);

		switch(iarg)
		{
			case 'x':
				xml = optarg;
				break;
			case 't':
				tag = optarg;
                                break;
			case '?':
			case 'h':
				std::cout<<"Allowed arguments:"<<std::endl;
				std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
				return 0;
		}
	}

        gStyle->SetPalette(87);
	//Load up our covariance matricies we calculated in example1 (we could also load up single variation ones)
	TString filename = "../bin/"+tag+".SBNcovar.root";
	TFile * fsys = new TFile(filename,"read");
	TMatrixD *cov = (TMatrixD*)fsys->Get("collapsed_covariance");
std::cout << "4" << std::endl;
	TMatrixD *fraccov = (TMatrixD*)fsys->Get("collapsed_frac_covariance");
std::cout << "5" << std::endl;
	TMatrixD *corr = (TMatrixD*)fsys->Get("collapsed_correlation");

        plot_one(*cov, "SBNfit_covariance_matrix_"+tag);
        plot_one(*fraccov, "SBNfit_fractional_covariance_matrix_"+tag);
        plot_one(*corr, "SBNfit_correlation_matrix_"+tag);

	return 0;
	
	
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
