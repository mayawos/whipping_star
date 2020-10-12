#ifndef PELEEPLOTTER_H_
#define PELEEPLOTTER_H_

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
#include "TMatrixD.h"
#include "TPaveStats.h"
#include "TH1.h"
#include "THStack.h"
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

namespace sbn{


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


};
#endif
