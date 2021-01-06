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
#include "TH2F.h"
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
#include "TGraphAsymmErrors.h"
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


std::vector<double> SubtractVectors(std::vector<double> vector1, std::vector<double> vector2){
  
  std::vector<double> difference;
  
  for(int i=0; i<vector1.size(); i++)
    difference.push_back( fabs(vector1[i]-vector2[i]));
  
  return difference;
  
}


int main(int argc, char* argv[])
{
  
  std::string xml = "example.xml";
  std::string tag = "";
  std::string var = "chi2";
  std::string channel = "combined";
  float bin = 0.25;
  bool print_mode = false;
  bool debug = false;
  int set=0;
  
  /*************************************************************
   *************************************************************
   *		Command Line Argument Reading
   ************************************************************
   ************************************************************/
  const struct option longopts[] =
    {
      {"xml", 		required_argument, 	0, 'x'},
      {"channel", 	required_argument, 	0, 'c'},
      {"tag", 		required_argument,	0, 't'},
      {"set", 		required_argument,	0, 's'},
      {"bin", 		required_argument,	0, 'b'},
      {"var", 		required_argument,	0, 'v'},
      {"print", 	no_argument, 		0, 'p'},
      {"debug", 	no_argument, 		0, 'd'},
      {0,		no_argument, 		0,  0},
	};
  
  int iarg = 0;
  opterr=1;
  int index;
  
  while(iarg != -1)
    {
      iarg = getopt_long(argc,argv, "x:s:b:t:c:v:p:d", longopts, &index);
      
      switch(iarg)
	{
	case 'x':
	  xml = optarg;
	  break;
	case 's':
	  set = atoi(optarg);
	  break;
	case 'b':
	  set = atof(optarg);
	  break;
	case 't':
	  tag = optarg;
	  break;
	case 'c':
	  channel = optarg;
	  break;
	case 'v':
	  var = optarg;
	  break;
	case 'p':
	  print_mode=true;
	  break;
	case 'd':
	  debug=true;
	  break;
	case '?':
	case 'h':
	  std::cout<<"Allowed arguments:"<<std::endl;
	  std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
	  std::cout<<"\t-t\t--tag\t\tA unique tag to identify the outputs."<<std::endl;
	  std::cout<<"\t-c\t--channel\tName of channel to plot the signal strength. (default: 'combined'; other options: 'np','zp' )" <<std::endl;
	  std::cout<<"\t-s\t--set\t\tA number to identify the fake data set."<<std::endl;
	  std::cout<<"\t-v\t--var\t\tparameter to plot wrt signal scaling (default: 'chi2'; other options: 'pvals' for p-values,'sign' for significance)."<<std::endl;
	  std::cout<<"\t-p\t--print\t\tRuns in print mode, making a lot more plots and Variations. (warning can take a while!) "<<std::endl;
	  std::cout<<"\t-d\t--debug\t\tRuns in debug "<<std::endl;
	  return 0;
	}
    }
  
  //POT for each datasets
  std::vector<std::string> POT = {"5.01e20 POT","7.97e20 POT","7.74e20 POT","7.75e20 POT","9e20 POT"};
  
  double  W = 800;
  double  H  = 600;
  double  T = 0.12*H;
  double  B = 0.12*H;
  double  L = 0.12*W;
  double  R = 0.04*W;
  auto c41 = new TCanvas("c","c",100,100,W,H);
  c41->SetFillColor(0);
  c41->SetBorderMode(0);
  c41->SetFrameFillStyle(0);
  c41->SetFrameBorderMode(0);
  c41->SetLeftMargin( L/W );
  c41->SetRightMargin( R/W );
  c41->SetTopMargin( T/H );
  c41->SetBottomMargin( B/H );
  c41->SetTickx(0);
  c41->SetTicky(0);
  c41->SetGrid();
  
  std::vector< std::vector< std::vector<double> > > results;
  std::vector< std::vector<double> > result;
  std::vector<double> a_result;
  std::vector<std::string> vectorname = {"median","sig1up","sig1dn","sig2up","sig2dn"};
  std::vector<std::string> hypothesisname = {"H1","data"};
  std::cout << vectorname.size() << std::endl;
  
  //get the label based on channel
  std::string label;
  if(channel=="combined") label = "nue_1e0p_numu"; 
  else if(channel=="np") label = "nue_numu"; 
  else if(channel=="zp") label = "1e0p_numu"; 
  //get the x params vector
  double min_res = 10e4;
  double max_res = -10e4;
  for(int j=0; j<hypothesisname.size(); j++){
    for(int i=0; i<vectorname.size(); i++){
       double a;
       std::string name = Form("%s_%s_%s%s_fakedata%d_strength_%s.txt",var.c_str(), vectorname[i].c_str(),tag.c_str(),label.c_str(),set,hypothesisname[j].c_str());
       std::cout << name <<std::endl;
       std::fstream myfile(name.c_str());
       if(myfile.is_open())
	 {
	   while ( myfile >> a )
	     {
	       if(debug) std::cout << "a = " << a << std::endl; 
	       a_result.push_back(a);
               if(j==1) continue;
	       if(a < min_res ) min_res = a;
	       if(a > max_res ) max_res = a;
	     } 
	   myfile.close();
	 }
       result.push_back(a_result);
       a_result.clear();
    }
    results.push_back(result);
    result.clear();
  }
  
  std::vector< std::vector<double> > sig1up, sig1dn, sig2up,  sig2dn;
  
  for(int i=0; i < 1; i++){
    sig1up.push_back(SubtractVectors(results[i][1],results[i][0]));
    sig1dn.push_back(SubtractVectors(results[i][0],results[i][2]));
    sig2up.push_back(SubtractVectors(results[i][3],results[i][0]));
    sig2dn.push_back(SubtractVectors(results[i][0],results[i][4]));
  }
  
  if(debug) std::cout << "results[0] = " << results[0][0].size() << std::endl;
  
  double isig=0.;
  int points=results[0][0].size();
  double start_bin = isig;
  double end_bin = isig + points*bin - bin;
  double ymin = *min_element(results[0][0].begin(), results[0][0].end()); 
  double ymax = *max_element(results[0][0].begin(), results[0][0].end()); 
  
  //set the y-axis title
  std::string ytitle;
  if(var=="chi2") ytitle="#Delta#chi^{2}";
  else if(var=="pvals") ytitle="p-values";
  else if(var=="sign") ytitle="#sigma";

  TH2F frame("frame",Form(";LEE signal scaling; %s;",ytitle.c_str()),points,start_bin,end_bin,1,ymin,ymax);
  TGraphAsymmErrors med(points); TGraphAsymmErrors onesig(points); TGraphAsymmErrors twosig(points); TGraphAsymmErrors fd(points);
  
  for(int k=0; k < points; k++){
    if(debug) std::cout << "point, signal scaling = " << k << ", " << isig << std::endl;
    if(debug) std::cout << "median = " << results[1][0][k] << std::endl;
    fd.SetPoint(k+1, isig, results[1][0][k]);
    if(debug) std::cout << "data = " << results[0][0][k] << std::endl;
    med.SetPoint(k+1, isig, results[0][0][k]);
    onesig.SetPoint(k+1, isig, results[0][0][k]);
    if(debug) std::cout << "-1 std dev, +1 std dev = " << sig1dn[0][k] << ", " << sig1up[0][k] <<std::endl;
    onesig.SetPointError(k+1, 0., 0., sig1dn[0][k], sig1up[0][k]);
    if(debug) std::cout << "-2 std dev, +2 std dev = " << sig2dn[0][k] << ", " << sig2up[0][k] <<std::endl;
    twosig.SetPoint(k+1, isig, results[0][0][k]);
    twosig.SetPointError(k+1, 0., 0., sig2dn[0][k], sig2up[0][k]);
    isig+=bin;
  }
  
  fd.SetLineColor(kRed);
  fd.SetLineWidth(2);
  
  med.SetLineColor(kBlack);
  med.SetLineStyle(2);
  med.SetLineWidth(2);
  
  onesig.SetFillColor(kYellow);
  onesig.SetLineColor(kYellow);
  
  twosig.SetFillColor(kGreen);
  twosig.SetLineColor(kGreen);
  
  frame.GetXaxis()->SetTitleSize(0.05);
  frame.GetYaxis()->SetTitleSize(0.05);
  frame.GetXaxis()->CenterTitle();
  frame.GetYaxis()->CenterTitle();
  frame.GetXaxis()->SetLabelSize(0.04);
  frame.GetYaxis()->SetLabelSize(0.04);
  frame.Draw(); gStyle->SetOptStat(0);
  
  twosig.Draw("same3");
  onesig.Draw("same3");
  med.Draw("same");
  fd.Draw("same");
  
  std::string legtitle;
  if( channel == "combined" ) legtitle = "1e0p0#pi + 1eNp0#pi";
  if( channel == "zp" ) legtitle = "1e0p0#pi";
  if( channel == "np" ) legtitle = "1eNp0#pi";
  
  double x1 = 0.15;
  double x2 = x1 + 0.3;
  double y2 = 0.86;
  double y1 = 0.62;

  auto legend = new TLegend(x1,y1,x2,y2);
  //legend->SetFillStyle(0);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.041);
  legend->SetTextFont(62);
  legend->SetHeader(legtitle.c_str());
  legend->SetTextSize(0.041);
  legend->SetTextFont(42);
  legend->AddEntry(&fd, "obs.","l");
  legend->AddEntry(&med, "expected median","l");
  legend->AddEntry(&onesig, "#pm 1 std. deviation","f");
  legend->AddEntry(&twosig, "#pm 2 std. deviation","f");
  legend->Draw("same");
  
  //title
  TLatex latex;
  latex.SetTextSize(0.05);
  latex.SetTextAlign(13); //align at top
  latex.SetTextFont(61);  //align at top
  latex.DrawLatex(start_bin,ymax+0.08*ymax,"MicroBooNE");
  latex.SetTextFont(52);  //align at top
  latex.SetTextSize(0.045);
  latex.DrawLatex(start_bin+1.5*bin, ymax+0.065*ymax,"Preliminary");
  latex.SetTextFont(62);  //align at top
  latex.SetTextSize(0.04);
  latex.DrawLatex(end_bin-1.1*bin, ymax+0.065*ymax, POT[set-1].c_str());
  
  TString cname = Form("signal_strength_%s%sset%d.pdf",tag.c_str(),label.c_str(),set);
  c41->Print(cname,"pdf");
  
}
