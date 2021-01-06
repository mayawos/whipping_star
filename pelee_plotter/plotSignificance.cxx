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
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
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

std::vector<double> difference, diff_fabs;

for(int i=0; i<vector1.size(); i++)
  difference.push_back( fabs(vector1[i]-vector2[i]));

return difference;

}

std::string SetPrecision(double a, bool pval){

std::string number=Form("%.1f",a);

if(!pval) return number;

if(int(a*10) >= 1) number = Form("%.1f",a);
else if(int(a*100) >= 1) number = Form("%.2f",a);
else if(int(a*1000) >= 1) number = Form("%.3f",a);
else if(int(a*10000) > 1) number = Form("%.4f",a);
else if(int(a*100000) > 1) number = Form("%.5f",a);

return number;

}

int main(int argc, char* argv[])
{

  std::string xml = "example.xml";
  std::string tag = "";
  std::string var = "chi2";
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
      {"tag", 		required_argument,	0, 't'},
      {"set", 		required_argument,	0, 's'},
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
      iarg = getopt_long(argc,argv, "x:s:t:v:p:d", longopts, &index);
      
      switch(iarg)
	{
	case 'x':
	  xml = optarg;
	  break;
	case 's':
	  set = atoi(optarg);
	  break;
	case 't':
	  tag = optarg;
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
	  std::cout<<"\t-s\t--set\t\tA number to identify the fake data set."<<std::endl;
	  std::cout<<"\t-v\t--var\t\tparameter to plot wrt signal scaling (default: 'chi2'; other options: 'pvals' for p-values,'sign' for significance)."<<std::endl;
	  std::cout<<"\t-p\t--print\t\tRuns in print mode, making a lot more plots and Variations. (warning can take a while!) "<<std::endl;
	  std::cout<<"\t-d\t--debug\t\tRuns in debug "<<std::endl;
	  return 0;
	}
    }

   std::vector<double> combine_H0, combine_H1, combine_data, np_H0, np_H1, np_data, zp_H0, zp_H1, zp_data;
   std::vector<std::string> POT = {"5.01e20 POT","7.97e20 POT","7.74e20 POT","7.75e20 POT","9e20 POT"};

   double  W = 900;
   double  H  = 800;
   double  T = 0.1*H;
   double  B = 0.12*H;
   double  L = 0.12*W;
   double  R = 0.04*W;
   auto c1 = new TCanvas("c","c",100,100,W,H);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetFrameFillStyle(0);
   c1->SetFrameBorderMode(0);
   c1->SetLeftMargin( 0.17 );
   c1->SetRightMargin( R/W );
   c1->SetTopMargin( T/H );
   c1->SetBottomMargin( B/H );
   c1->SetTickx(0);
   c1->SetTicky(0);
   c1->SetGridx(0);

   // create the arrays for the points
   Int_t n = 3;

   //get the p-vals vector
   std::vector< std::vector< std::vector<double> > > pval_results, chi2_results;
   std::vector< std::vector<double> > pval_result, chi2_result;
   std::vector<double> a_result;
   std::vector<std::string> vectorname = {"median","sig1up","sig1dn"};
   std::vector<std::string> hypothesisname = {"H0","H1","data"};
   std::cout << vectorname.size() << std::endl;

   bool pval=false;
   std::string var1, var2;
   if(var=="sign"){
     var1="sign";
     var2="sign";
   }else{
     pval=true;
     var1="pvals";
     var2="chi2";
   } 

   for(int j=0; j<hypothesisname.size(); j++){
     for(int i=0; i<vectorname.size(); i++){
       double a;
       std::string name1 = Form("%s_%s%s_fakedata%d_%s.txt",var1.c_str(), vectorname[i].c_str(),tag.c_str(),set,hypothesisname[j].c_str());
       std::cout << name1 <<std::endl;
       std::fstream myfile(name1.c_str());
       if(myfile.is_open())
       {
         while ( myfile >> a )
         {
           std::cout << var << " = " << a << std::endl; 
           a_result.push_back(a);
         } 
         myfile.close();
       }
       std::reverse(a_result.begin(),a_result.end());
       pval_result.push_back(a_result);
       a_result.clear();
     }
     pval_results.push_back(pval_result);
     pval_result.clear();
   }

   //get the chi2 vector
   double min_chi2 = 10e4;
   double max_chi2 = -10e4;
   for(int j=0; j<hypothesisname.size(); j++){
     for(int i=0; i<vectorname.size(); i++){
       double a;
       std::string name2 = Form("%s_%s%s_fakedata%d_%s.txt",var2.c_str(), vectorname[i].c_str(),tag.c_str(),set,hypothesisname[j].c_str());
       std::cout << name2 <<std::endl;
       std::fstream myfile(name2.c_str());
       if(myfile.is_open())
       {
         while ( myfile >> a )
         {
           std::cout << var << " = " << a << std::endl; 
           a_result.push_back(a);
           if(a < min_chi2 ) min_chi2 = a;
           if(a > max_chi2 ) max_chi2 = a;
         } 
         myfile.close();
       }
       std::reverse(a_result.begin(),a_result.end());
       chi2_result.push_back(a_result);
       a_result.clear();
     }
     chi2_results.push_back(chi2_result);
     chi2_result.clear();
   }

   //get the sigma -- only for H0 and H1
   std::vector< std::vector<double> > chi2_sig1up, chi2_sig1dn, pval_sig1up,  pval_sig1dn;
   for(int i=0; i < 2; i++){
     chi2_sig1up.push_back(SubtractVectors(chi2_results[i][1],chi2_results[i][0]));
     chi2_sig1dn.push_back(SubtractVectors(chi2_results[i][0],chi2_results[i][2]));
     pval_sig1up.push_back(SubtractVectors(pval_results[i][1],pval_results[i][0]));
     pval_sig1dn.push_back(SubtractVectors(pval_results[i][0],pval_results[i][2]));
   }

   std::vector<std::string> channames = {"Combined","1e0p0#pi", "1eNp0#pi"};

   double step = (max_chi2-min_chi2)/2; 
   double start_bin = 0.8*min_chi2;
   if(min_chi2<0) start_bin = 1.2*min_chi2;
   double end_bin = 1.8*max_chi2+(3*step);
   std::string xtitle = "#Delta#chi^{2}_{LEE}"; 
   if(var=="pvals") xtitle = "#Delta#chi^{2}_{LEE}"; 
   //TH2F frame("frame","; #Delta#chi^{2}_{LEE};",1,start_bin,end_bin,n+1,0,n+1);
   TH2F frame("frame","",1,start_bin,end_bin,n+1,0,n+1);
   int iChann = 0; TGraphAsymmErrors points0(n); TGraphAsymmErrors points1(n); TGraphAsymmErrors points2(n);
   for(int k=0; k < n; k++){
     TString channel = channames[k];
     points0.SetPoint(iChann, chi2_results[0][0][k], iChann+0.5);
     points1.SetPoint(iChann, chi2_results[1][0][k], iChann+0.5);
     points2.SetPoint(iChann, chi2_results[2][0][k], iChann+0.5);
     points0.SetPointError(iChann, chi2_sig1dn[0][k], chi2_sig1up[0][k], 0, 0);
     points1.SetPointError(iChann, chi2_sig1dn[1][k], chi2_sig1up[1][k], 0, 0);
     iChann++;
     frame.GetYaxis()->SetBinLabel(iChann, channel);
     frame.GetYaxis()->SetTitleFont(62);
     frame.GetXaxis()->SetTitleFont(62);
     frame.GetXaxis()->SetTitle(xtitle.c_str());
   }
   points0.SetLineColor(kGreen+2);
   points0.SetLineWidth(3);
   points0.SetMarkerStyle(20);
   points0.SetMarkerColor(kGreen+2);
   points0.SetMarkerSize(2);
   points1.SetLineColor(kBlack);
   points1.SetLineWidth(3);
   points1.SetMarkerStyle(20);
   points1.SetMarkerColor(kBlack);
   points1.SetMarkerSize(2);
   points2.SetLineColor(kRed);
   points2.SetLineWidth(3);
   points2.SetMarkerSize(3);
   points2.SetMarkerStyle(29);
   points2.SetMarkerColor(kRed);
   frame.GetXaxis()->SetTitleSize(0.05);
   frame.GetXaxis()->CenterTitle();
   frame.GetXaxis()->SetLabelSize(0.04);
   frame.GetYaxis()->SetLabelSize(0.055);
   frame.Draw(); gStyle->SetOptStat(0);
   points0.Draw("PE SAME");
   points1.Draw("PE SAME");
   points2.Draw("PE SAME");

   double x1 = 0.3;
   double x2 = x1 + 0.65;
   double y2 = 0.85;
   double y1 = 0.78;

   std::string legtitle;

   auto legend = new TLegend(x1,y1,x2,y2);
   legend->SetNColumns(3);
   //legend->SetFillStyle(0);
   legend->SetFillColor(0);
   legend->SetBorderSize(0);
   legend->SetTextSize(0.03);
   legend->SetTextFont(42);
   legend->AddEntry(&points0, "median H0","lpe");
   legend->AddEntry(&points1, "median H1","lpe");
   legend->AddEntry(&points2, "fakedata","p");
   legend->Draw("same");

   //title
   TLatex latex;
   latex.SetTextSize(0.05);
   latex.SetTextAlign(13);  //align at top
   latex.SetTextFont(61);  //align at top
   latex.DrawLatex(start_bin,4.25,"MicroBooNE");
   latex.SetTextFont(52);  //align at top
   latex.SetTextSize(0.035);
   latex.DrawLatex(start_bin+2.*step,4.18," Preliminary");
   latex.SetTextFont(62);  //align at top
   latex.SetTextSize(0.035);
   latex.DrawLatex(end_bin-1.3*step,4.2,POT[set-1].c_str());

   //calculate the p-vals text position
   double start = max_chi2; 
   std::string number;
   std::vector<int> colors={kGreen+2,kBlack,kRed};
   for(int i=0; i<3;i++ ){
   for(int k=0; k<3;k++ ){
	 //errors
	 latex.SetTextFont(62);  //align at top
	 latex.SetTextSize(0.035);
	 latex.SetTextColor(colors[i]);
	 double offset = (i*1.3*step)+step;
         number = SetPrecision(pval_results[i][0][k],pval);
         if(debug) std::cout << var << " median = " << k << ", " << pval_results[i][0][k] << std::endl;
	 latex.DrawLatex(start+offset, 0.6+k,number.c_str());
	 if(i!=2){
	   number = "+"+SetPrecision(pval_sig1up[i][k],pval);
           if(debug) std::cout << var << "+1sigma = " << pval_sig1up[i][k] << std::endl;
	   latex.SetTextSize(0.025);
	   latex.DrawLatex(start+(step/2)+offset, 0.7+k,number.c_str());
	   number = "-"+SetPrecision(pval_sig1dn[i][k],pval);
           if(debug) std::cout << var << "-1sigma = " << pval_sig1dn[i][k] << std::endl;
	   latex.SetTextSize(0.025);
	   latex.DrawLatex(start+(step/2)+offset, 0.5+k,number.c_str());
	 }
       }
   }
   
   TLine *line = new TLine(start_bin,1,end_bin,1);
   line->SetLineWidth(3);
   line->SetLineStyle(2);
   line->Draw("same");

   std::string title = Form("plotSig_%s%s_set%d.pdf",tag.c_str(),var.c_str(),set);
   c1->Print(title.c_str());

   return 0;

}
