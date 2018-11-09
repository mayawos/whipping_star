#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <cstring>
#include <array>


#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
/*
#include "TH1F.h"
#include "TString.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TMath.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TH2F.h"
*/


int main(){

	int nBins = 10;
	double unfoldedWgts[] = {0, 5.03093, 4.50515, 3.50515, 2.31959, 1.31959, 0.64948, 0.27835, 0.11340, 0};
	double binedges[] = {0,200,250,300,350,400,450,500,600,800,3000};
		
	double binedgessmall[] = {100,500,900,1300,1800};

	// Read in sample file
	std::cout << "Reading in file" << std::endl;
	TFile *infile = new TFile("/uboone/data/users/yatesla/othersys/inputs_to_sbnfit/arxiv/input_to_sbnfit_August2018.root","READ");
	TTree *TIntrinsic = (TTree*)infile->Get("nue_intrinsic_tree");
	
	// Grab true energy branch
	double etrue,ereco;
	TIntrinsic->SetBranchAddress("true_nu_energy",&etrue);
	TIntrinsic->SetBranchAddress("reco_energy",&ereco);

	// Make outputs
	std::cout << "Creating output file unfolder.root" << std::endl;
	TFile *f = new TFile("unfolder.root","RECREATE");
	TTree *t = new TTree("nue_intrinsic_tree","Tree of the weights to unfold the signal");
	double unfoldedWgt;
	t->Branch("unfoldedWgt",&unfoldedWgt);
	
	TH1D * tru = new TH1D("bkg","",4,binedgessmall);
	TH1D * reco = new TH1D("lee","",4,binedgessmall);

	// Loop through branches and get weights.
	std::cout << "Looping over " << TIntrinsic->GetEntries() << " events" << std::endl;
	for(int i = 0; i < TIntrinsic->GetEntries(); i++){
		TIntrinsic->GetEntry(i);

		for(int iB = 0; iB < nBins; iB++){
			if(etrue < binedges[iB+1]){
				unfoldedWgt = unfoldedWgts[iB];
				break;
			}
		}

		tru->Fill(ereco,.13);
		reco->Fill(ereco,unfoldedWgt*.13);
		
		if(ereco > 800.0){
			std::cout << i << ": " << "RECO: " << ereco << " TRUE: " << etrue << " WEIGHT: " << unfoldedWgt << std::endl;
		}

		t->Fill();
	}

	

	tru->Write();
	reco->Write();


	// Write it all and save it!
	t->Write();
	f->Close();
	infile->Close();

	return 1;
}







