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
#include "TRandom.h"
#include "TError.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TLatex.h"

#include "params.h"
#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBNfit.h"
#include "SBNfit3pN.h"
#include "SBNcovariance.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;

/*************************************************************
 *************************************************************
 *		BEGIN example.cxx
 ************************************************************
 ************************************************************/

void CollapseSubchannels(TMatrixD & M, TMatrixD & Mc, std::vector<int> num_bins, int num_channels, std::vector<int> num_subchannels);
void plot_one(TMatrixD matrix, std::vector<int> num_bins, std::string tag);

int main(int argc, char* argv[])
{

  std::string xml = "example.xml";
  bool print_mode = false;
  bool nolee = false;
  bool nocorr = false;
  
  /*************************************************************
   *************************************************************
   *		Command Line Argument Reading
   ************************************************************
   ************************************************************/
  const struct option longopts[] =
    {
      {"xml", 		required_argument, 	0, 'x'},
      {"print", 		no_argument, 	0, 'p'},
      {"tag", 		required_argument,	0, 't'},
      {"nolee", 		no_argument,	0, 'l'},
      {"nocorr", 		no_argument,	0, 'c'},
      {0,			no_argument, 	0,  0},
    };
  
  int iarg = 0;
  opterr=1;
  int index;
  
  //a tag to identify outputs and this specific run. defaults to EXAMPLE1
  std::string tag = "EXAMPLE1";
  
  while(iarg != -1)
    {
      iarg = getopt_long(argc,argv, "x:t:plc", longopts, &index);
      
      switch(iarg)
	{
	case 'x':
	  xml = optarg;
	  break;
	case 'p':
	  print_mode=true;
	  break;
        case 't':
	  tag = optarg;
	  break;
	case 'l':
	  nolee=true;
	  break;
	case 'c':
	  nocorr=true;
	  break;
        case '?':
	case 'h':
	  std::cout<<"Allowed arguments:"<<std::endl;
	  std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
	  std::cout<<"\t-t\t--tag\t\tA unique tag to identify the outputs."<<std::endl;
	  std::cout<<"\t-p\t--print\t\tRuns in print mode, making a lot more plots and Variations. (warning can take a while!) "<<std::endl;
	  return 0;
	}
    }
  
  //std::string dict_location = "../libio/libEventWeight.so";
  //std::cout<<"Trying to load dictionary: "<<dict_location<<std::endl;
  //gSystem->Load(  (dict_location).c_str());
  
  /*************************************************************
   *************************************************************
   *			Main Program Flow
   ************************************************************
   ************************************************************/
  
  std::vector< std::vector<double> > detsysmatrix;
  
  detsysmatrix.push_back({0.0403, 0.0085, -0.0022, 0.0088, -0.0040, 0.0098, 0.0040, -0.0023, -0.0042, 0.0100, -0.0394, 0.0021, -0.0142, 0.0127, 0.0122, -0.0068, -0.0262, -0.0235, -0.0452, -0.0321});
  detsysmatrix.push_back({0.0085, 0.0116, 0.0081, 0.0049, 0.0035, 0.0000, -0.0040, 0.0054, -0.0071, 0.0076, 0.0009, 0.0018, -0.0130, 0.0042, 0.0061, -0.0003, 0.0058, 0.0090, 0.0281, 0.0013});
  detsysmatrix.push_back({-0.0022, 0.0081, 0.0085, 0.0022, 0.0041, -0.0031, -0.0046, 0.0059, -0.0045, 0.0042, 0.0123, 0.0011, -0.0056, -0.0022, 0.0007, 0.0017, 0.0113, 0.0137, 0.0314, 0.0101});
  detsysmatrix.push_back({0.0088, 0.0049, 0.0022, 0.0033, 0.0003, 0.0021, -0.0002, 0.0011, -0.0023, 0.0049, -0.0064, 0.0008, -0.0061, 0.0053, 0.0036, -0.0011, -0.0022, -0.0005, 0.0033, -0.0048});
  detsysmatrix.push_back({-0.0040, 0.0035, 0.0041, 0.0003, 0.0044, -0.0036, -0.0037, 0.0046, -0.0040, -0.0008, 0.0130, -0.0010, -0.0035, -0.0049, 0.0018, 0.0047, 0.0102, 0.0137, 0.0201, 0.0095});
  detsysmatrix.push_back({0.0098, 0.0000, -0.0031, 0.0021, -0.0036, 0.0057, 0.0041, -0.0023, 0.0027, 0.0046, -0.0156, 0.0024, -0.0023, 0.0114, 0.0014, -0.0053, -0.0094, -0.0108, -0.0122, -0.0115});
  detsysmatrix.push_back({0.0040, -0.0040, -0.0046, -0.0002, -0.0037, 0.0041, 0.0052, -0.0020, 0.0057, 0.0015, -0.0088, 0.0026, 0.0040, 0.0071, -0.0024, -0.0042, -0.0076, -0.0092, -0.0175, -0.0068});
  detsysmatrix.push_back({-0.0023, 0.0054, 0.0059, 0.0011, 0.0046, -0.0023, -0.0020, 0.0076, -0.0019, 0.0029, 0.0148, 0.0027, -0.0047, -0.0004, 0.0001, 0.0025, 0.0123, 0.0159, 0.0275, 0.0113});
  detsysmatrix.push_back({-0.0042, -0.0071, -0.0045, -0.0023, -0.0040, 0.0027, 0.0057, -0.0019, 0.0095, -0.0002, -0.0014, 0.0024, 0.0105, 0.0046, -0.0076, -0.0041, -0.0050, -0.0073, -0.0151, -0.0010});
  detsysmatrix.push_back({0.0100, 0.0076, 0.0042, 0.0049, -0.0008, 0.0046, 0.0015, 0.0029, -0.0002, 0.0107, -0.0075, 0.0048, -0.0082, 0.0136, 0.0022, -0.0054, -0.0007, -0.0000, 0.0155, -0.0046});
  detsysmatrix.push_back({-0.0394, 0.0009, 0.0123, -0.0064, 0.0130, -0.0156, -0.0088, 0.0148, -0.0014, -0.0075, 0.0651, -0.0010, 0.0070, -0.0202, -0.0094, 0.0148, 0.0448, 0.0521, 0.0782, 0.0495});
  detsysmatrix.push_back({0.0021, 0.0018, 0.0011, 0.0008, -0.0010, 0.0024, 0.0026, 0.0027, 0.0024, 0.0048, -0.0010, 0.0072, -0.0010, 0.0059, -0.0022, -0.0054, 0.0012, -0.0018, 0.0035, 0.0002});
  detsysmatrix.push_back({-0.0142, -0.0130, -0.0056, -0.0061, -0.0035, -0.0023, 0.0040, -0.0047, 0.0105, -0.0082, 0.0070, -0.0010, 0.0204, -0.0096, -0.0119, 0.0005, -0.0035, -0.0066, -0.0301, 0.0049});
  detsysmatrix.push_back({0.0127, 0.0042, -0.0022, 0.0053, -0.0049, 0.0114, 0.0071, -0.0004, 0.0046, 0.0136, -0.0202, 0.0059, -0.0096, 0.0311, 0.0031, -0.0098, -0.0072, -0.0079, 0.0153, -0.0137});

  detsysmatrix.push_back({0.0122, 0.0061, 0.0007, 0.0036, 0.0018, 0.0014, -0.0024, 0.0001, -0.0076, 0.0022, -0.0094, -0.0022, -0.0119, 0.0031, 0.0097, 0.0022, -0.0026, 0.0007, 0.0024, -0.0077});
  detsysmatrix.push_back({-0.0068, -0.0003, 0.0017, -0.0011, 0.0047, -0.0053, -0.0042, 0.0025, -0.0041, -0.0054, 0.0148, -0.0054, 0.0005, -0.0098, 0.0022, 0.0103, 0.0100, 0.0157, 0.0099, 0.0102});
  //detsysmatrix.push_back({-0.0068, -0.0003, 0.0017, -0.0011, 0.0047, -0.0053, -0.0042, 0.0025, -0.0041, -0.0054, 0.0148, -0.0054, 0.0005, -0.0098, 0.0022, 0.0103, 0.0100, 0.0157, 0.0099, 0.0102});
  detsysmatrix.push_back({-0.0262, 0.0058, 0.0113, -0.0022, 0.0102, -0.0094, -0.0076, 0.0123, -0.0050, -0.0007, 0.0448, 0.0012, -0.0035, -0.0072, -0.0026, 0.0100, 0.0372, 0.0420, 0.0761, 0.0358});
  //detsysmatrix.push_back({-0.0262, 0.0058, 0.0113, -0.0022, 0.0102, -0.0094, -0.0076, 0.0123, -0.0050, -0.0007, 0.0448, 0.0012, -0.0035, -0.0072, -0.0026, 0.0100, 0.0372, 0.0420, 0.0761, 0.0358});
  detsysmatrix.push_back({-0.0235, 0.0090, 0.0137, -0.0005, 0.0137, -0.0108, -0.0092, 0.0159, -0.0073, -0.0000, 0.0521, -0.0018, -0.0066, -0.0079, 0.0007, 0.0157, 0.0420, 0.0541, 0.0870, 0.0397});
  //detsysmatrix.push_back({-0.0235, 0.0090, 0.0137, -0.0005, 0.0137, -0.0108, -0.0092, 0.0159, -0.0073, -0.0000, 0.0521, -0.0018, -0.0066, -0.0079, 0.0007, 0.0157, 0.0420, 0.0541, 0.0870, 0.0397});
  detsysmatrix.push_back({-0.0452, 0.0281, 0.0314, 0.0033, 0.0201, -0.0122, -0.0175, 0.0275, -0.0151, 0.0155, 0.0782, 0.0035, -0.0301, 0.0153, 0.0024, 0.0099, 0.0761, 0.0870, 0.2282, 0.0658});
  //detsysmatrix.push_back({-0.0452, 0.0281, 0.0314, 0.0033, 0.0201, -0.0122, -0.0175, 0.0275, -0.0151, 0.0155, 0.0782, 0.0035, -0.0301, 0.0153, 0.0024, 0.0099, 0.0761, 0.0870, 0.2282, 0.0658});
  detsysmatrix.push_back({-0.0321, 0.0013, 0.0101, -0.0048, 0.0095, -0.0115, -0.0068, 0.0113, -0.0010, -0.0046, 0.0495, 0.0002, 0.0049, -0.0137, -0.0077, 0.0102, 0.0358, 0.0397, 0.0658, 0.0390});
  //detsysmatrix.push_back({-0.0321, 0.0013, 0.0101, -0.0048, 0.0095, -0.0115, -0.0068, 0.0113, -0.0010, -0.0046, 0.0495, 0.0002, 0.0049, -0.0137, -0.0077, 0.0102, 0.0358, 0.0397, 0.0658, 0.0390});
  //detsysmatrix.push_back({-0.0321, 0.0013, 0.0101, -0.0048, 0.0095, -0.0115, -0.0068, 0.0113, -0.0010, -0.0046, 0.0495, 0.0002, 0.0049, -0.0137, -0.0077, 0.0102, 0.0358, 0.0397, 0.0658, 0.0390});
  //detsysmatrix.push_back({-0.0321, 0.0013, 0.0101, -0.0048, 0.0095, -0.0115, -0.0068, 0.0113, -0.0010, -0.0046, 0.0495, 0.0002, 0.0049, -0.0137, -0.0077, 0.0102, 0.0358, 0.0397, 0.0658, 0.0390});
  //detsysmatrix.push_back({-0.0321, 0.0013, 0.0101, -0.0048, 0.0095, -0.0115, -0.0068, 0.0113, -0.0010, -0.0046, 0.0495, 0.0002, 0.0049, -0.0137, -0.0077, 0.0102, 0.0358, 0.0397, 0.0658, 0.0390});
  //detsysmatrix.push_back({-0.0321, 0.0013, 0.0101, -0.0048, 0.0095, -0.0115, -0.0068, 0.0113, -0.0010, -0.0046, 0.0495, 0.0002, 0.0049, -0.0137, -0.0077, 0.0102, 0.0358, 0.0397, 0.0658, 0.0390});
  
  std::vector<double> numu_detsys = {0.096,0.097,0.066,0.051,0.065,0.093,0.081,0.07,0.109,0.122,0.142,0.158,0.18,0.261};
  
  TString bgname = tag+".SBNspec.root "; 
  SBNspec bkg(bgname.Data(),xml);
  std::string uniquestring;
  if(!nolee && !nocorr) uniquestring = "corr_lee";
  else if(nolee && !nocorr) uniquestring = "corr_nolee";
  else if(nocorr && !nolee) uniquestring = "nocorr_lee";
  else if(nocorr && nolee) uniquestring = "nocorr_nolee";
  TString matrixname = "detsys_cov_matrix_"+uniquestring+".root";
  TFile * out = new TFile(matrixname,"recreate");
  TMatrixD fullfracmatrix(bkg.num_bins_total,bkg.num_bins_total);
  TMatrixD fullcorrmatrix(bkg.num_bins_total,bkg.num_bins_total);
  TMatrixD fullcovmatrix(bkg.num_bins_total,bkg.num_bins_total);
  fullfracmatrix.Zero();
  fullcovmatrix.Zero();
  fullcorrmatrix.Zero();
  TMatrixD zpmatrix(6,6);
  zpmatrix.Zero();
  
  int j=0;
  std::vector<int> bins_1eNp_nue, bins_1eNp_lee, bins_1e0p_nue, bins_1e0p_lee, bins_numu_bnb, bins_bg;   
  for(auto& h: bkg.hist){
    std::string hname = h.GetName();
    for( int i=1; i < h.GetNbinsX()+1; i++ ){
      if( hname.find("nue_intrinsic") != std::string::npos ){
	bins_1eNp_nue.push_back(j);
      }
      else if( hname.find("nue_lee") != std::string::npos ){
	bins_1eNp_lee.push_back(j);
      }
      else if( hname.find("1e0p_intrinsic") != std::string::npos ){
        //std::cout << "bins_1e0p_nue = " << j << std::endl;
	bins_1e0p_nue.push_back(j);
      }
      else if( hname.find("1e0p_lee") != std::string::npos ){
        //std::cout << "bins_1e0p_lee = " << j << std::endl;
	bins_1e0p_lee.push_back(j);
      }
      else if( hname.find("numu_bnb") != std::string::npos ){
	bins_numu_bnb.push_back(j);
      }
      else{ 
        if( hname.find("ext") == std::string::npos ){
	  bins_bg.push_back(j);
        }
      }
      j++; 
    }
  }
  
  //now fill in the matrix
  //fill in 1eNp
  for(int i = 0; i < bins_1eNp_nue.size(); i++){
    for(int j = 0; j < bins_1eNp_nue.size(); j++){
      //std::cout << "bins 1eNp " << bins_1eNp_nue[i] << ", " << bins_1eNp_nue[j] << std::endl;
      //std::cout << "bins 1eNp lee " << bins_1eNp_lee[i] << ", " << bins_1eNp_lee[j] << std::endl;
      if(i==j) fullfracmatrix[bins_1eNp_nue[i]][bins_1eNp_nue[j]] = detsysmatrix[i][j];
      else{ 
        if(!nocorr) fullfracmatrix[bins_1eNp_nue[i]][bins_1eNp_nue[j]] = detsysmatrix[i][j];
      }
      if(!nolee){
        if(i==j) fullfracmatrix[bins_1eNp_lee[i]][bins_1eNp_lee[j]] = detsysmatrix[i][j];
        else{ 
          if(!nocorr) fullfracmatrix[bins_1eNp_lee[i]][bins_1eNp_lee[j]] = detsysmatrix[i][j];
        }
      }
    }
  }
  
  //fill in 1e0p
  std::cout << "bins_1e0p_nue size, fullfracmatrix = " << bins_1e0p_nue.size() << ", " << fullfracmatrix.GetNrows() << "x" << fullfracmatrix.GetNcols() <<std::endl;
  for(int i = 0; i < bins_1e0p_nue.size(); i++){
    for(int j = 0; j < bins_1e0p_nue.size(); j++){
      if(i==j) fullfracmatrix[bins_1e0p_nue[i]][bins_1e0p_nue[j]] = detsysmatrix[i+14][j+14];
      else{ 
        if(!nocorr) fullfracmatrix[bins_1e0p_nue[i]][bins_1e0p_nue[j]] = detsysmatrix[i+14][j+14];
      }
      if(!nolee){
        if(i==j) fullfracmatrix[bins_1e0p_lee[i]][bins_1e0p_lee[j]] = detsysmatrix[i+14][j+14];
        else{ 
          if(!nocorr) fullfracmatrix[bins_1e0p_lee[i]][bins_1e0p_lee[j]] = detsysmatrix[i+14][j+14];
        }
      }
      zpmatrix[i][j] = detsysmatrix[i+14][j+14];
      //std::cout << "row,cols = " << i+14 << ", " << j+14 << " = " << detsysmatrix[i+14][j+14] << std::endl;
    }
  }

  //plot_one(zpmatrix, {14,0,0}, "SBNfit_collapsed_frac_zpdetsys_"+tag+"_mc");  
  //fill in 1eNp-1e0p
  for(int i = 0; i < bins_1eNp_nue.size(); i++){
    for(int j = 0; j < bins_1e0p_nue.size(); j++){
      if(!nocorr) fullfracmatrix[bins_1eNp_nue[i]][bins_1e0p_nue[j]] = detsysmatrix[i][j+14];
      if(!nocorr && !nolee) fullfracmatrix[bins_1eNp_lee[i]][bins_1e0p_lee[j]] = detsysmatrix[i][j+14];
    }
  }
  
  //fill in 1e0p-1eNp
  for(int i = 0; i < bins_1e0p_nue.size(); i++){
    for(int j = 0; j < bins_1eNp_nue.size(); j++){
      if(!nocorr) fullfracmatrix[bins_1e0p_nue[i]][bins_1eNp_nue[j]] = detsysmatrix[i+14][j];
      if(!nocorr&&!nolee) fullfracmatrix[bins_1e0p_lee[i]][bins_1eNp_lee[j]] = detsysmatrix[i+14][j];
    }
  }
  
  //fill in numu diagonals
  for(int i = 0; i < bins_numu_bnb.size(); i++){
    std::cout << "numu bins = " << bins_numu_bnb[i] << ", " << bins_numu_bnb[i] << std::endl;
    fullfracmatrix[bins_numu_bnb[i]][bins_numu_bnb[i]] += numu_detsys[i]*numu_detsys[i];
  }
  
  //fill in bg diagonals
  for(int i = 0; i < bins_bg.size(); i++){
    for(int j = 0; j < bins_bg.size(); j++){
      //std::cout << "bg bins = " << bins_bg[i] << ", " << bins_bg[i] << std::endl;
      if(i==j) fullfracmatrix[bins_bg[i]][bins_bg[j]] += 0.2*0.2;
      else{
        if(!nocorr) fullfracmatrix[bins_bg[i]][bins_bg[j]] += 0.2*0.2;
      }
    }
  }

  bkg.CalcFullVector();
  std::cout << "rows, cols = " << fullfracmatrix.GetNrows() << ", " << fullfracmatrix.GetNcols() << std::endl;
  for(int i=0; i < fullfracmatrix.GetNrows(); i++){
    for(int j=0; j < fullfracmatrix.GetNcols(); j++){
      //std::cout << "fullfracmatrix(i,j), bkg.full_vector[i], i, j = " << fullfracmatrix(i,j) << ", " << bkg.full_vector[i] << ", " << i << ", " << j << std::endl;
      fullcovmatrix(i,j) = fullfracmatrix(i,j)*bkg.full_vector[i]*bkg.full_vector[j]; 	 
      fullcorrmatrix(i,j) = fullcovmatrix(i,j)/(sqrt(fullcovmatrix(i,i))*sqrt(fullcovmatrix(j,j))); 	
    }
  }
  //std::cout << "1" << std::endl;
  TMatrixD collmatrix(bkg.num_bins_total_compressed,bkg.num_bins_total_compressed);
  TMatrixD covmatrix(bkg.num_bins_total_compressed,bkg.num_bins_total_compressed);
  TMatrixD corrmatrix(bkg.num_bins_total_compressed,bkg.num_bins_total_compressed);
  collmatrix.Zero();
  covmatrix.Zero();
  corrmatrix.Zero();

  int num_channels = 3;
  std::vector<int> num_subchannels = {11,11,3};
  std::vector<int> num_bins = {14,6,14};
  CollapseSubchannels(fullcovmatrix,covmatrix,num_bins, num_channels, num_subchannels);  

  bkg.CollapseVector();
  std::cout << "size of collapsed vector = " << bkg.collapsed_vector.size() << std::endl;
  //collapsed vector
  //calculate cov matrix and corr matrix
  for(int i=0; i < collmatrix.GetNrows(); i++){
    for(int j=0; j < collmatrix.GetNcols(); j++){
      collmatrix(i,j) = covmatrix(i,j)/(bkg.collapsed_vector[i]*bkg.collapsed_vector[j]); 	
      corrmatrix(i,j) = covmatrix(i,j)/(sqrt(covmatrix(i,i))*sqrt(covmatrix(j,j))); 
    }
  }
  //collmatrix.Print();

  plot_one(collmatrix, num_bins, "SBNfit_collapsed_frac_detsys_"+uniquestring+"_"+tag+"_mc");
  plot_one(covmatrix, num_bins, "SBNfit_collapsed_cov_detsys_"+uniquestring+"_"+tag+"_mc");
  plot_one(corrmatrix, num_bins, "SBNfit_collapsed_corr_detsys_"+uniquestring+"_"+tag+"_mc");

  //corrmatrix.Print(); 
  fullfracmatrix.Write("full_frac_covariance_matrix");
  fullcovmatrix.Write("full_covariance_matrix");
  fullcorrmatrix.Write("full_correlation_matrix");
  collmatrix.Write("collapsed_frac_covariance_matrix");
  covmatrix.Write("collapsed_covariance_matrix");
  corrmatrix.Write("collapsed_correlation_matrix");

  out->Close();
  
  return 0;  
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

void plot_one(TMatrixD matrix, std::vector<int> num_bins, std::string tag){

    std::string dir="/uboone/data/users/wospakrk/SBNFitPlots/";
    std::vector<TString> channel_names = {"#nu_{e} 1eNp","#nu_{e} 1e0p","#nu_{#mu} BNB"};
    int num_channels = num_bins.size();
    int num_bins_total = matrix.GetNrows();
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

        TLatex * tmd = new TLatex(-num_bins_total*percent_left*0.15, use_full+nice_shift*0.8, (channel_names[ic]).Data() );
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

