#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <cstring>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TNtuple.h"
#include "TLine.h"
#include "TChain.h"
#include "TMath.h"
#include "TSystem.h"
#include "TMatrixT.h"
#include "TRandom.h"
#include "TError.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TMultiGraph.h"

#include "params.h"
#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBNfit.h"
#include "SBNfit3pN.h"
#include "SBNcovariance.h"
#include "SBNfeld.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;


double Median(const TH1D * h1) { 

       int n = h1->GetXaxis()->GetNbins();  
          std::vector<double>  x(n);
             h1->GetXaxis()->GetCenter( &x[0] );
                const double * y = h1->GetArray(); 
                   // exclude underflow/overflows from bin content array yG
                       return TMath::Median(n, &x[0], &y[1]); 
                       }
                   

double quick_median(std::vector<double> &v)
{
        size_t n = v.size() / 2;
        std::nth_element(v.begin(), v.begin()+n, v.end());
                return v[n];
}

double lin_interp(double x0, double x1, double y0, double y1, double x){
    return (y0*(x1-x)+y1*(x-x0))/(x1-x0);
}
/*************************************************************
 *************************************************************
 *		BEGIN sbnfit_make_covariance.cxx
 ************************************************************
 ************************************************************/
void runHelp(){
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"Modified single subchannel scaling feldman_cousins confidence belt constructor"<<std::endl;
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"--- Required arguments: ---"<<std::endl;
                std::cout<<"\t-x\t--xml\t\tInput configuration .xml file for SBNconfig"<<std::endl;
                std::cout<<"\t-t\t--tag\t\tA unique tag to identify the inputs/outputs [Default to TEST]"<<std::endl;
                std::cout<<"\t-i\t--input\t\tInput subchannel to scale (no default, required argument)"<<std::endl;
                std::cout<<"\t-g\t--grid\t\tGrid to scan, in the form 'min max num_steps' (default '1e-4 10.0 20')"<<std::endl;
                std::cout<<"\t-m\t--mode\t\tWhat mode you want to run in. Arguments are:"<<std::endl;
                std::cout<<"\t\t\t--\t feldman : Perform the pseudo universe grid scan (run first)"<<std::endl;  
                std::cout<<"\t\t\t--\t belt: Constructs the confidence belts, must be run after'feldman'"<<std::endl;
                std::cout<<"\t\t\t--\t data: Pass in an optional datafile, that will be compared to the grid"<<std::endl;
                std::cout<<"\t\t\t--\t beltdata: Same as belt, and plots line at BF value"<<std::endl;
                std::cout<<"--- Optional arguments: ---"<<std::endl;
                std::cout<<"\t-s\t--stat\t\tStatistical error only mode, will ignore any covariance matrix passed in"<<std::endl;
                std::cout<<"\t-n\t--number\t\tNumber of pseudo-experiments to simulate (default 2500)"<<std::endl; 
                std::cout<<"\t-r\t--randomseed\t\tRandomNumber Seed (default from machine)"<<std::endl; 
                std::cout<<"\t-c\t--cnp\t\tuse a Combined Newman Pearson chi2 (default false)"<<std::endl;
                std::cout<<"\t-d\t--data\t\ta data SBNspec file to input, use with mode data"<<std::endl;
                std::cout<<"\t-f\t--frequentist\t\t, perform an equivalent of lee frequentist study, grid dimension is fixed to '1e-6 2 2'"<<std::endl;
                std::cout<<"\t-y\t--detsys\t\t, add detector systematics matrix to the default covariance matrix"<<std::endl;
                std::cout<<"\t-h\t--help\t\tThis help menu."<<std::endl;
                std::cout<<"---------------------------------------------------"<<std::endl;
    return;
}

double Median(std::vector<double> vec){
  size_t vsize = vec.size();
  double median = 0.;
  if (vsize == 0){
     return median; 
  }
  
  sort(vec.begin(), vec.end()); 
  if (vsize%2 == 0){
    median = (vec[vsize / 2 - 1] + vec[vsize / 2])/2.;
  }else if(vsize%2 != 0){
    median = vec[vsize / 2];
  }else{
    median = 0.;
  }
  return median;
}


int main(int argc, char* argv[])
{

    std::string xml = "oscillate_example.xml";

    /*************************************************************
     *************************************************************
     *		Command Line Argument Reading
     ************************************************************
     ************************************************************/
    const struct option longopts[] =
    {
        {"xml", 		required_argument, 	0, 'x'},
        {"stat", 		no_argument, 		0, 's'},
        {"number", 		required_argument,	0,'n'},
        {"cnp", 		no_argument,	0,'c'},
        {"llr",         no_argument,    0,'l'},
        {"grid", 		required_argument,	0,'g'},
        {"tag", 		required_argument,	0, 't'},
        {"mode",        required_argument, 0 ,'m'},
        {"data",        required_argument, 0 ,'d'},
        {"input",       required_argument, 0 ,'i'},
        {"randomseed",        required_argument, 0 ,'r'},
        {"frequentist",        no_argument, 0 ,'f'},
        {"detsys",        no_argument, 0 ,'y'},
        {"help", 		no_argument,	0, 'h'},
        {0,			    no_argument, 		0,  0},
    };

    int iarg = 0;
    opterr=1;
    int index;

    //a tag to identify outputs and this specific run. defaults to EXAMPLE1
    std::string tag = "TEST";
    std::string mode_option;
    bool bool_stat_only = false;
    int number = 2500;
    double random_number_seed = 10;
    bool use_cnp = false;
    bool use_llr = false;
    bool do_simple_hypothesis = false;
    bool detsys = false;
        
    std::string grid_string = "1e-4 8.0 33";
    std::string input_scale_subchannel = "unset";
    std::string data_file_input = "null";

    while(iarg != -1)
    {
        iarg = getopt_long(argc,argv, "x:t:m:n:r:d:p:i:g:schfy", longopts, &index);

        switch(iarg)
        {
            case 'x':
                xml = optarg;
                break;
            case 'n':
                number = (int)strtod(optarg,NULL);
                break;
            case 't':
                tag = optarg;
                break;
            case 'g':
                grid_string = optarg; 
                break;
            case 'd':
                data_file_input = optarg; 
                break;

            case 'i':
                input_scale_subchannel = optarg;
                break;
            case 'm':
                mode_option = optarg;
                break;
            case 'c':
                use_cnp = true;
                break;
            case 'l':
                use_llr = true;
                break;
            case 'r':
                random_number_seed = (double)strtod(optarg,NULL);
                std::cout<<"Reading in random seed argument: "<<random_number_seed<<std::endl;
                break;
            case 's':
                bool_stat_only = true;
                break;
            case 'f':
                do_simple_hypothesis = true;
                break;
            case 'y':
                detsys = true;
                break;
            case '?':
            case 'h':
                runHelp(); 
                return 0;
        }
    }



    /*************************************************************
     *************************************************************
     *			Main Program Flow
     ************************************************************
     ************************************************************/
    time_t start_time = time(0);

    std::cout<<"Begining SBNfit uboone subchannel scaling Feldman Cousins confidence belt constructor for tag: "<<tag<<std::endl;

    if(input_scale_subchannel=="unset"){
        std::cout<<"Error! you must set a value for which input subchannel to scale, e.g using --input/ -i 'nu_uBooNE_1g1p_ncdelta'"<<std::endl;
        std::cout<<"Please see...."<<std::endl;
        runHelp();
        return 0;
    }
    
    NGrid mygrid;

    if( do_simple_hypothesis ) mygrid.AddDimension(input_scale_subchannel.c_str(),grid_string);
    else mygrid.AddDimension(input_scale_subchannel.c_str(),grid_string);
  
    mygrid.Print();

    SBNfeld myfeld(mygrid,tag,xml);

        std::size_t npzp_tag = tag.find("np_zp_numu_reco");
        std::size_t np_tag = tag.find("np_numu_reco");
        std::size_t zp_tag = tag.find("constrained_zp_numu_reco");

            /////////LABELS/////////////
        TString SampleName;
        if (npzp_tag < 100){
            SampleName = "1eNp + 1e0p";
        }else if(np_tag < 100){
            SampleName = "1eNp";
        }else if(zp_tag < 100){
            SampleName = "1e0p";
        }else{
            SampleName = "";
        }
        TString FDSName;
        for (int i = 1; i <= 5 ; ++i){
            std::size_t find_fds = tag.find(Form("fakedata%d",i));
            if (find_fds < 100){
                FDSName = Form("FDS%d", i);
                break;
            }else{
               FDSName=""; 
            }
        }
        TString TestStatName;
        if(use_llr){
            TestStatName = "LLR";
        }else if(use_cnp){
            TestStatName = "CNP";
        }else{
            TestStatName = "";
        }



//check if running constraint?
bool use_constraint = false;
if(tag.find("constrained") != std::string::npos ) use_constraint = true;

    if(mode_option == "feldman"){

        //setup covariance matrix 
        TFile * fsys;
        TFile * fdetsys;
        TMatrixD * cov;
        TMatrixD * covdetsys;

        //allows for stats-only test
        std::cout<<"Begininning a full Feldman-Cousins analysis for tag : "<<tag<<std::endl;

        if(bool_stat_only){
            myfeld.SetEmptyFractionalCovarianceMatrix();
            myfeld.SetStatOnly();
            std::cout<<"RUNNING Statistics uncertainty only!"<<std::endl;
        }else{
            std::cout<<"RUNNING Systematics Covariance Matrix!"<<std::endl;
            std::string tagcovar = tag;
            if(detsys && use_constraint) tagcovar = tagcovar + "_detsys";
            fsys = new TFile(Form("%s.SBNcovar.root",tagcovar.c_str()),"read");
            cov = (TMatrixD*)fsys->Get("frac_covariance");
            if( detsys && !use_constraint){ 
               std::cout<<"RUNNING Systematics with Maya's Covariance Matrix!"<<std::endl;
	           fdetsys = new TFile(Form("/uboone/data/users/wospakrk/PeLEE/SBNfit/Fakedata/%s_detsys.root",tag.c_str()),"read");
	           covdetsys = (TMatrixD*)fdetsys->Get("full_frac_covariance_matrix");
	           *cov = *cov + *covdetsys;
	        }

            myfeld.SetFractionalCovarianceMatrix(cov);
        }

        if(use_cnp) myfeld.UseCNP();
        if(use_llr) myfeld.UseLLR();
        if(do_simple_hypothesis) myfeld.doSimpleHypothesis();


        myfeld.m_subchannel_to_scale = input_scale_subchannel;
        
        myfeld.SetCoreSpectrum(tag+"_CV.SBNspec.root");
        myfeld.SetBackgroundSpectrum(tag+"_CV.SBNspec.root",input_scale_subchannel,0.0);
        myfeld.GenerateScaledSpectra();
        
        std::cout<<"Setting random seed "<<random_number_seed<<std::endl;
        myfeld.SetRandomSeed(random_number_seed);
        myfeld.SetNumUniverses(number);

        std::cout<<"Calculating the necessary SBNchi objects"<<std::endl;
        myfeld.CalcSBNchis();
        std::cout <<"DONE calculating the necessary SBNchi objects at : " << difftime(time(0), start_time)/60.0 << " Minutes.\n";

        std::cout<<"Beginning to peform FullFeldmanCousins analysis"<<std::endl;
        myfeld.FullFeldmanCousins();


    }else if(mode_option=="data"){

        std::cout<<"Begininning a real data analysis for tag : "<<tag<<std::endl;

        if(detsys) {
            myfeld.SetFractionalCovarianceMatrix(tag+"_detsys.SBNcovar.root","frac_covariance");
        }else{
            myfeld.SetFractionalCovarianceMatrix(tag+".SBNcovar.root","frac_covariance");
        }
        myfeld.m_subchannel_to_scale = input_scale_subchannel;

        if(use_cnp) myfeld.UseCNP();

        myfeld.SetCoreSpectrum(tag+"_CV.SBNspec.root");

        myfeld.SetBackgroundSpectrum(tag+"_CV.SBNspec.root",input_scale_subchannel,1.0);
        myfeld.GenerateScaledSpectra();

        std::cout<<"Setting random seed "<<random_number_seed<<std::endl;
        myfeld.SetRandomSeed(random_number_seed);
        myfeld.SetNumUniverses(number);

        std::cout<<"Calculating the necessary SBNchi objects"<<std::endl;
        myfeld.CalcSBNchis();
        std::cout <<"DONE calculating the necessary SBNchi objects at : " << difftime(time(0), start_time)/60.0 << " Minutes.\n";

        std::cout << data_file_input.c_str() << std::endl;
        SBNspec * datain = new SBNspec(data_file_input.c_str(),xml);

        myfeld.CompareToData(datain);


    }else if(mode_option == "belt"){

        if(bool_stat_only){
            myfeld.SetEmptyFractionalCovarianceMatrix();
            myfeld.SetStatOnly();
            std::cout<<"RUNNING Statistics uncertainty only!"<<std::endl;
        }else{
            if(detsys) {
                myfeld.SetFractionalCovarianceMatrix(tag+"_detsys.SBNcovar.root","frac_covariance");
            }else{
                myfeld.SetFractionalCovarianceMatrix(tag+".SBNcovar.root","frac_covariance");
            }
        }
        myfeld.m_subchannel_to_scale = input_scale_subchannel;

        if(use_cnp) myfeld.UseCNP();
        if(use_llr) myfeld.UseLLR();
/*
        myfeld.SetCoreSpectrum(tag+"_CV.SBNspec.root");
        myfeld.SetBackgroundSpectrum(tag+"_CV.SBNspec.root",input_scale_subchannel,1.0);
        myfeld.GenerateScaledSpectra();

        std::cout<<"Calculating the necessary SBNchi objects"<<std::endl;
        myfeld.CalcSBNchis();
        std::cout <<"DONE calculating the necessary SBNchi objects at : " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
*/
        TFile *f = new TFile("scan.root","recreate");
        f->cd();
        std::cout<<"Starting to peform a globalScan analysis"<<std::endl;
        std::vector<std::vector<double>> vec_grid = mygrid.GetGrid();

        TFile *fin = new TFile(("SBNfeld_output_"+tag+".root").c_str(),"read");
        //TFile *fin = new TFile("/uboone/data/users/wospakrk/FCplots/SBNfeld_output_nue_1e0p_numu_reco_e_H1_mc_fakedata1.root","read");

        //Some Manual Color Changing and such
        int plotting_true_gridpoint = 5;
        //std::vector<double> plotting_pvals = {0.68, 0.90, 0.95, 0.99};

        //std::vector<double> plotting_pvals = {0.68, 0.95};
        //std::vector<std::string> plotting_strs = {"68%","90%","95%","99%"};
        //std::vector<std::string> plotting_strs = {"68%","95%"};

        //std::vector<int> gcols = {kGreen+3,kGreen+2,kGreen-3,kGreen-9};
        //std::vector<int> gcols = {kRed-9,kBlue-9,kGreen-9};

        std::vector<double> plotting_pvals = {0.682689492, 0.90, 0.997300};
        std::vector<std::string> plotting_strs = {"1#sigma (~68%)","90%","3#sigma (~99.7%)"};
        std::vector<int> gcols = {kRed-9,kBlue-9, kGreen-9};

        std::vector<double> v_median;
        std::vector<double> v_true;
        std::vector<double> v_1sigma_p;
        std::vector<double> v_1sigma_m;
        std::vector<std::vector<double>> v_min;v_min.resize(plotting_pvals.size());
        std::vector<std::vector<double>> v_max;v_max.resize(plotting_pvals.size());

        std::cout<<"MPrinting stuff"<<std::endl;
        TH2D * f_FC = new TH2D("f_FC","f_FC",vec_grid.size(),vec_grid.front()[0],vec_grid.back()[0],vec_grid.size(),vec_grid.front()[0],vec_grid.back()[0]);

        double bfval_v[vec_grid.size()];

        for(int i=0; i< vec_grid.size(); i++){

            v_true.push_back(vec_grid[i][0]);

            //Whats the critical value?
            TTree *t =  (TTree*)fin->Get(("ttree_"+std::to_string(i)).c_str());
            TH1D * cumul = (TH1D*)fin->Get(("delta_chi2_"+std::to_string(i)+"_cumulative").c_str());

            TH1D * h_bfval = (TH1D*)fin->Get(("bf_gridvalue_"+std::to_string(i)).c_str());//Added by Ivan


            for(int p =0; p< plotting_pvals.size(); ++p){
                double plotting_pval = plotting_pvals[p];

                //First lets find a critical chi^2 for this confidence level
                double critical_delta_chi2 = 0;
                double critical_delta_chi2_mid = 0;
                for(int c = cumul->GetNbinsX()-1;  c>0 ; --c){
                    if(cumul->GetBinContent(c+1) >= plotting_pval && cumul->GetBinContent(c)< plotting_pval){
                        critical_delta_chi2_mid = cumul->GetBinCenter(c); //lin_interp(double x0, double x1, double y0, double y1, double x) = (y0*(x1-x)+y1*(x-x0))/(x1-x0);
                        critical_delta_chi2 = lin_interp(cumul->GetBinContent(c+1), cumul->GetBinContent(c), cumul->GetBinLowEdge(c+1), cumul->GetBinLowEdge(c)+cumul->GetBinWidth(c), plotting_pval);
                        break;
                   }
                }

                std::cout<<"Grid point "<<i<<" has a critical delta chi of "<<critical_delta_chi2<<"("<<critical_delta_chi2_mid<<") for a pval of "<<plotting_pval<<std::endl; 

                std::string nam = std::to_string(i)+"bfhist";
                int Nentries = t->GetEntries(); 

                const unsigned nentries = t->Draw("bf_gridvalue", ("delta_chi2<="+std::to_string(critical_delta_chi2)).c_str());
                if (nentries) {
                    double* x = t->GetV1();
                    //std::cout << "min :" << *(std::min_element(x, x+nentries)) << std::endl;
                    //std::cout << "max :" << *(std::max_element(x, x+nentries)) << std::endl;
                    v_min[p].push_back( *(std::min_element(x, x+nentries)) );
                    v_max[p].push_back( *(std::max_element(x, x+nentries)) );
                }
            }//end pval loop

            //adding median best-fit values to confidence belt plot 
            h_bfval->ComputeIntegral();
            std::vector<double> pvalues = { 0.5 };
            std::vector<double> bf_val_quantiles(pvalues.size());
            h_bfval->GetQuantiles(pvalues.size(),&bf_val_quantiles[0], &pvalues[0]);
            double bfval = bf_val_quantiles[0];

            //double bfval = h_bfval->GetMean(); //Added by Ivan
            bfval_v[i] = bfval;
            delete cumul;
        }


        TCanvas *c3 = new TCanvas("h_ono_r3");
        c3->SetFillStyle(0);

        TPad *pad = new TPad("pad", "pad", 0, 0, 0.8, 1.0);
        pad->SetRightMargin(0); // Upper and lower plot are joined
        pad->Draw();             // Draw the upper pad: pad
        pad->cd();               // pad becomes the current pad

        std::vector<TGraph*> gmaxs;
        std::vector<TGraph*> gmins;
        std::vector<TGraph*> grshades;

        TLegend * l_probs = new TLegend(0.11,0.52,0.89,0.89);//69 was 29

        TMultiGraph *mg = new TMultiGraph();
        mg->SetTitle("Feldman Cousins Confidence Belt [ " + SampleName + " - " + TestStatName + " ]");

        for(int p=plotting_pvals.size()-1; p>=0;--p){
            pad->cd();
            
            gmins.push_back(new TGraph(v_true.size(),&(v_min[p])[0], &v_true[0]));
            gmaxs.push_back(new TGraph(v_true.size(),&(v_max[p])[0], &v_true[0]));
            grshades.push_back( new TGraph(2*v_true.size()));

            for (int i=0;i<v_true.size();i++) {
                int n = v_true.size();
                grshades.back()->SetPoint(i,v_min[p][i],v_true[i]);
                grshades.back()->SetPoint(n+i,v_max[p][n-i-1],v_true[n-i-1]);
            }
            grshades.back()->SetFillColor(gcols[p]);
            gmins.back()->SetLineWidth(2);
            gmins.back()->SetLineColor(kBlack);
            gmaxs.back()->SetLineWidth(2);
            gmaxs.back()->SetLineColor(kBlack);

            l_probs->AddEntry(grshades.back(), plotting_strs[p].c_str() ,"f");
        }

        l_probs->SetHeader("#splitline{#splitline{Classical}{Confidence}}{#splitline{Level of Interval}{}}");

        for(auto &g:grshades)mg->Add(g);
        pad->cd();
        //mg->GetXaxis()->SetRangeUser(0,3);

        mg->Draw("ALF");
        
        mg->GetXaxis()->SetTitle("Measured Signal Strength (x)");
        mg->GetYaxis()->SetTitle("True Signal Strength (#mu)");
        mg->SetMinimum(v_true.front());
        mg->SetMinimum(v_true.front());

        //calculate the 2d points based on grid dimension and resolution
        double npoints = vec_grid.size();
        double maxpt = vec_grid.back()[0];
        double scale = (npoints-1)/maxpt;

        //mg->GetYaxis()->SetLimits(v_true.front(),4);
        //mg->GetXaxis()->SetLimits(v_true.front(),4);
        
        //mg->GetXaxis()->SetLimits(v_true.front(),maxpt-1); //Original
        mg->GetXaxis()->SetLimits(v_true.front(),maxpt); //org maxpt
        mg->GetYaxis()->SetLimits(v_true.front(),maxpt); 

        //mg->GetHistogram()->SetMaximum(maxpt-1);//v_true.back());  //Original
        mg->GetHistogram()->SetMaximum(maxpt);        
        mg->GetHistogram()->SetMinimum(v_true.front());     
 
        for(auto &g:gmins)g->Draw("l same");
        for(auto &g:gmaxs)g->Draw("l same");
        double u_measured[vec_grid.size()];
        double uexp = 0;

        std::cout << "//////////HERE!!!!!!!! vec_grid.size() = " << vec_grid.size() << std::endl;
        for(int k = 0; k < npoints; k++){
            uexp = ((double) k)/scale;
            //std::cout << "//////////HERE!!!!!!!! k, uexp, median = " << k << ", " << uexp << ", " << bfval_v[k] << std::endl;
            u_measured[k] = uexp;
        }

        //TGraph *bf_vals_g = new TGraph(vec_grid.size(),u_measured,bfval_v);
        TGraph *bf_vals_g = new TGraph(vec_grid.size(), bfval_v, &v_true[0]); //Given a true value, i.e Y axis point, whats the median reconstructed best fit point
        //bf_vals_g->GetXaxis()->SetRangeUser(0,3);
        bf_vals_g->Draw("same *");
        
        TLine lcross(v_true.front(),v_true.front(),npoints, npoints);

        lcross.SetLineStyle(9);
        lcross.SetLineWidth(1);
        lcross.SetLineColor(kBlack);
        lcross.Draw("same");


        TLine lv1(1.0,0.0,1.0,maxpt);
        lv1.SetLineStyle(2);
        lv1.SetLineWidth(1);
        lv1.SetLineColor(kBlack);
        //lv1.Draw("same");

        TLine lh3(0.0,3.1, 1.0,3.1);
        lh3.SetLineStyle(2);
        lh3.SetLineWidth(1);
        lh3.SetLineColor(kBlack);
        //lh3.Draw("same");

        pad->Update();
        pad->RedrawAxis();
        // TLine l;
        //  l.DrawLine(gPad->GetUxmin(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymax());
        //   l.DrawLine(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax());
        c3->cd();
        TPad *padl = new TPad("padl", "padl", 0.8, 0, 1, 1);
        padl->SetBottomMargin(0.2);
        padl->Draw();
        padl->cd();       // padl becomes the current pad
        l_probs->SetTextSize(0.0775);
        l_probs->Draw();
        l_probs->SetLineColor(kWhite);
        l_probs->SetLineWidth(0);
       

        c3->SaveAs(("FC_confidence_belt_"+tag+".pdf").c_str(),"pdf");

        std::cout<<"**************** Feldman Cousins 1D Confidence Intervals  **********************"<<std::endl;
        for(int i=0; i<v_true.size(); i++){
            std::cout<<"Grid Pt: "<<i<<", ScaleFactor: "<<vec_grid[i][0]<<std::endl;
            std::cout<<"- Median: "<<bfval_v[i]<<std::endl;
            for(int p=0; p< plotting_pvals.size();++p){
                double measured_val = vec_grid[i][0];
                std::vector<double> reg = myfeld.getConfidenceRegion(gmins[plotting_pvals.size()-p-1],gmaxs[plotting_pvals.size()-p-1],measured_val);
                std::cout<<"-- CL: "<<plotting_pvals[p]<<"  Sigma: "<<sqrt(2)*TMath::ErfInverse(plotting_pvals[p])<<"  ConfidenceInterval: ["<<reg[0]<<" -> "<<reg[1]<<"]"<<std::endl;
            }
        }

    }else if(mode_option == "beltdata"){ // Ivan Added
        
        std::cout<<"Begininning a real data analysis for tag : "<<tag<<std::endl;

        if(detsys) {
            myfeld.SetFractionalCovarianceMatrix(tag+"_detsys.SBNcovar.root","frac_covariance");
        }else{
            myfeld.SetFractionalCovarianceMatrix(tag+".SBNcovar.root","frac_covariance");
        }

        if(use_cnp) myfeld.UseCNP();
        if(use_llr) myfeld.UseLLR();
        //myfeld.SetFractionalCovarianceMatrix(tag+".SBNcovar.root","frac_covariance");
        myfeld.m_subchannel_to_scale = input_scale_subchannel;

        myfeld.SetCoreSpectrum(tag+"_CV.SBNspec.root");
        myfeld.SetBackgroundSpectrum(tag+"_CV.SBNspec.root",input_scale_subchannel,1.0);
        myfeld.GenerateScaledSpectra();

        std::cout<<"Setting random seed "<<random_number_seed<<std::endl;
        myfeld.SetRandomSeed(random_number_seed);
        myfeld.SetNumUniverses(number);

        std::cout<<"Calculating the necessary SBNchi objects"<<std::endl;
        myfeld.CalcSBNchis();
        std::cout <<"DONE calculating the necessary SBNchi objects at : " << difftime(time(0), start_time)/60.0 << " Minutes.\n";

        SBNspec * datain = new SBNspec(data_file_input.c_str(),xml);
        //myfeld.CompareToData(datain);
        std::vector<double> result = myfeld.CompareToData(datain);
        double BF_val = result.at(2); ///BF wrt mu=0

        std::vector<std::vector<double> > results = myfeld.CompareDataToAllScaleFactors(datain);

           ////Using Deltachi2 to data///////////////////////
           int num_of_gridpts = results.size();
           
           double min_dchi2 = 1000000; double BF_AllMu = 0.;
           double dchi2_data[num_of_gridpts]; double gptval_data[num_of_gridpts];
           for (int i = 0; i < num_of_gridpts; ++i){
                dchi2_data[i] = results[i][0]; gptval_data[i] = results[i][4];
                //std::cout <<"i = " << i << ", Delta Chi2 = " << dchi2_data[i] << ", " << gptval_data[i] << std::endl;
                if(dchi2_data[i] < min_dchi2) {
                    min_dchi2 = dchi2_data[i];
                    BF_AllMu = gptval_data[i]; 
                }
            }

        if(bool_stat_only){
            myfeld.SetEmptyFractionalCovarianceMatrix();
            myfeld.SetStatOnly();
            std::cout<<"RUNNING Statistics uncertainty only!"<<std::endl;
        }else{
            if(detsys) {
                myfeld.SetFractionalCovarianceMatrix(tag+"_detsys.SBNcovar.root","frac_covariance");
            }else{
                myfeld.SetFractionalCovarianceMatrix(tag+".SBNcovar.root","frac_covariance");
            }
        }
        myfeld.m_subchannel_to_scale = input_scale_subchannel;



        
        ///Same as belt option from here on out//////
        TFile *f = new TFile("scan.root","recreate");
        f->cd();
        std::cout<<"Starting to peform a globalScan analysis"<<std::endl;
        std::vector<std::vector<double>> vec_grid = mygrid.GetGrid();

        TFile *fin = new TFile(("SBNfeld_output_"+tag+".root").c_str(),"read");
        //TFile *fin = new TFile("/uboone/data/users/wospakrk/FCplots/SBNfeld_output_nue_1e0p_numu_reco_e_H1_mc_fakedata1.root","read");

        //Some Manual Color Changing and such
        int plotting_true_gridpoint = 5;
        //std::vector<double> plotting_pvals = {0.68, 0.90, 0.95, 0.99};
        //double test = 87.64;
        //string test_s = Form("%f/%",test);
        std::vector<double> plotting_pvals = {0.68 , 0.95};
        //std::vector<std::string> plotting_strs = {"68%","90%","95%","99%"};

        std::vector<std::string> plotting_strs = {"68%","95%"};
        //std::vector<int> gcols = {kGreen+3,kGreen+2,kGreen-3,kGreen-9};
        std::vector<int> gcols = {kRed-9,kBlue-9,kGreen-9};


        std::vector<double> v_median;
        std::vector<double> v_true;
        std::vector<double> v_1sigma_p;
        std::vector<double> v_1sigma_m;
        std::vector<std::vector<double>> v_min;v_min.resize(plotting_pvals.size());
        std::vector<std::vector<double>> v_max;v_max.resize(plotting_pvals.size());

        std::cout<<"MPrinting stuff"<<std::endl;
        TH2D * f_FC = new TH2D("f_FC","f_FC",vec_grid.size(),vec_grid.front()[0],vec_grid.back()[0],vec_grid.size(),vec_grid.front()[0],vec_grid.back()[0]);
        double bfval_v[vec_grid.size()];


        for(int i=0; i< vec_grid.size(); i++){

            v_true.push_back(vec_grid[i][0]);

            //Whats the critical value?
            TTree *t =  (TTree*)fin->Get(("ttree_"+std::to_string(i)).c_str());
            TH1D * cumul = (TH1D*)fin->Get(("delta_chi2_"+std::to_string(i)+"_cumulative").c_str());
            TH1D * h_bfval = (TH1D*)fin->Get(("bf_gridvalue_"+std::to_string(i)).c_str());//Added by Ivan

            for(int p =0; p< plotting_pvals.size(); ++p){
                double plotting_pval = plotting_pvals[p];

                //First lets find a critical chi^2 for this confidence level
                double critical_delta_chi2 = 0;
                double critical_delta_chi2_mid = 0;
                for(int c = cumul->GetNbinsX()-1;  c>0 ; --c){
                    if(cumul->GetBinContent(c+1) >= plotting_pval && cumul->GetBinContent(c)< plotting_pval){
                        critical_delta_chi2_mid = cumul->GetBinCenter(c); //lin_interp(double x0, double x1, double y0, double y1, double x) = (y0*(x1-x)+y1*(x-x0))/(x1-x0);
                        critical_delta_chi2 = lin_interp(cumul->GetBinContent(c+1), cumul->GetBinContent(c), cumul->GetBinLowEdge(c+1), cumul->GetBinLowEdge(c)+cumul->GetBinWidth(c), plotting_pval);
                        break;
                   }
                }

                std::cout<<"Grid point "<<i<<" has a critical delta chi of "<<critical_delta_chi2<<"("<<critical_delta_chi2_mid<<") for a pval of "<<plotting_pval<<std::endl; 

                std::string nam = std::to_string(i)+"bfhist";
                int Nentries = t->GetEntries(); 

                const unsigned nentries = t->Draw("bf_gridvalue", ("delta_chi2<="+std::to_string(critical_delta_chi2)).c_str());
                if (nentries) {
                    double* x = t->GetV1();
                    //std::cout << "min :" << *(std::min_element(x, x+nentries)) << std::endl;
                    //std::cout << "max :" << *(std::max_element(x, x+nentries)) << std::endl;
                    v_min[p].push_back( *(std::min_element(x, x+nentries)) );
                    v_max[p].push_back( *(std::max_element(x, x+nentries)) );
                }
            }//end pval loop

            //adding median best-fit values to confidence belt plot 
            h_bfval->ComputeIntegral();
            std::vector<double> pvalues = { 0.5 };
            std::vector<double> bf_val_quantiles(pvalues.size());
            h_bfval->GetQuantiles(pvalues.size(),&bf_val_quantiles[0], &pvalues[0]);
            double bfval = bf_val_quantiles[0];
            //double bfval = h_bfval->GetMean(); //Added by Ivan
            bfval_v[i] = bfval;

            delete cumul;
        }


        TCanvas *c3 = new TCanvas("h_ono_r3");
        c3->SetFillStyle(0);

        TPad *pad = new TPad("pad", "pad", 0, 0, 0.8, 1.0);
        pad->SetRightMargin(0); // Upper and lower plot are joined
        pad->Draw();             // Draw the upper pad: pad
        pad->cd();               // pad becomes the current pad

        std::vector<TGraph*> gmaxs;
        std::vector<TGraph*> gmins;
        std::vector<TGraph*> grshades;

        TLegend * l_probs = new TLegend(0.11,0.45,0.89,0.89);//69 was 29

        TMultiGraph *mg = new TMultiGraph();
        mg->SetTitle("Feldman Cousins Confidence Belt [ " + SampleName + " - " + FDSName + " - " + TestStatName + " ]");

        for(int p=plotting_pvals.size()-1; p>=0;--p){
            pad->cd();
            
            gmins.push_back(new TGraph(v_true.size(),&(v_min[p])[0], &v_true[0]));
            gmaxs.push_back(new TGraph(v_true.size(),&(v_max[p])[0], &v_true[0]));
            grshades.push_back( new TGraph(2*v_true.size()));

            for (int i=0;i<v_true.size();i++) {
                int n = v_true.size();
                grshades.back()->SetPoint(i,v_min[p][i],v_true[i]);
                grshades.back()->SetPoint(n+i,v_max[p][n-i-1],v_true[n-i-1]);
            }
            grshades.back()->SetFillColor(gcols[p]);
            gmins.back()->SetLineWidth(2);
            gmins.back()->SetLineColor(kBlack);
            gmaxs.back()->SetLineWidth(2);
            gmaxs.back()->SetLineColor(kBlack);

            l_probs->AddEntry(grshades.back(), plotting_strs[p].c_str() ,"f");
            
        }

        l_probs->SetHeader("#splitline{#splitline{Classical}{Confidence}}{#splitline{Level of}{Interval}}");

        for(auto &g:grshades)mg->Add(g);
        pad->cd();
        mg->Draw("ALF");
        
        mg->GetXaxis()->SetTitle("Measured Signal Strength (x)");
        mg->GetYaxis()->SetTitle("True Signal Strength (#mu)");
        mg->SetMinimum(v_true.front());
        mg->SetMinimum(v_true.front());

        //calculate the 2d points based on grid dimension and resolution
        double npoints = vec_grid.size();
        double maxpt = vec_grid.back()[0];
        double scale = (npoints-1)/maxpt;

        //mg->GetXaxis()->SetLimits(v_true.front(),maxpt);

        mg->GetHistogram()->SetMaximum(maxpt);//v_true.back());          
        mg->GetHistogram()->SetMinimum(v_true.front());     
 
        mg->GetXaxis()->SetLimits(v_true.front(),maxpt); 
        mg->GetYaxis()->SetLimits(v_true.front(),maxpt);
        for(auto &g:gmins)g->Draw("l same");
        for(auto &g:gmaxs)g->Draw("l same");
        double u_measured[vec_grid.size()];
        double uexp = 0;

        //std::cout << "//////////HERE!!!!!!!! vec_grid.size() = " << vec_grid.size() << std::endl;
        for(int k = 0; k < npoints; k++){
            uexp = ((double) k)/scale;
            //std::cout << "//////////HERE!!!!!!!! k, uexp, median = " << k << ", " << uexp << ", " << bfval_v[k] << std::endl;
            u_measured[k] = uexp;
        }

        //TGraph *bf_vals_g = new TGraph(vec_grid.size(),u_measured,bfval_v);
        TGraph *bf_vals_g = new TGraph(vec_grid.size(), bfval_v, &v_true[0]); //Given a true value, i.e Y axis point, whats the median reconstructed best fit point
        
        bf_vals_g->Draw("same *");
        
        TLine lcross(v_true.front(),v_true.front(),npoints, npoints);
        lcross.SetLineStyle(9);
        lcross.SetLineWidth(1);
        lcross.SetLineColor(kBlack);
        lcross.Draw("same");

        TLine* l_BF = new TLine(BF_AllMu, 0, BF_AllMu, v_true.back());  //changed to BF_AllMu 
        l_BF->SetLineStyle(2);
        l_BF->SetLineWidth(2);
        l_BF->SetLineColor(kRed);
        l_BF->Draw("same");
        l_probs->AddEntry(l_BF, Form("#splitline{Data Best}{Fit = %1.2f}",BF_AllMu),"lf");
        
        double MB_val = 1;
        TLine l_MB(MB_val, 0, MB_val, v_true.back());
        //l_MB.SetLineStyle(0);
        l_MB.SetLineWidth(3);
        l_MB.SetLineColor(kGreen+3);
        //l_MB.Draw("same");

        TBox MBband(MB_val - 0.09, 0, MB_val + 0.09, v_true.back());
        MBband.SetFillColorAlpha(kGreen-3, 0.3);
        //MBband.Draw("SAME");
        
        
        TText *text_BF = new TText(BF_val,-0.27,Form("%2.2f",BF_val));
        text_BF->SetTextSize(0.03);
        //text_BF->Draw("same");

        pad->Update();
        pad->RedrawAxis();
        // TLine l;
        //  l.DrawLine(gPad->GetUxmin(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymax());
        //   l.DrawLine(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax());

        c3->cd();
        TPad *padl = new TPad("padl", "padl", 0.8, 0, 1, 1);
        padl->SetBottomMargin(0.2);
        padl->Draw();
        padl->cd();       // padl becomes the current pad
        l_probs->SetTextSize(0.0775);
        l_probs->Draw();
        l_probs->SetLineColor(kWhite);
        l_probs->SetLineWidth(0);
       


        c3->SaveAs(("DATACOMP_FC_confidence_belt_"+tag+".pdf").c_str(),"pdf");
        std::cout << std::endl;

        Int_t bf_step; double bf68_lo = 0; double bf68_up = 0; double bf95_lo = 0; double bf95_up = 0; double bf99_lo = 0;
        std::cout << tag << std::endl;
        //std::cout<<"**************** Feldman Cousins 1D Confidence Intervals  **********************"<<std::endl;
        for(int i=0; i<v_true.size(); i++){
            if (BF_val != vec_grid[i][0]) continue;
            //std::cout<<"Grid Pt: "<<i<<", ScaleFactor: "<<vec_grid[i][0]<<std::endl;
            bf_step = i;
            for(int p=0; p< plotting_pvals.size();++p){
                double measured_val = vec_grid[i][0];
                std::vector<double> reg = myfeld.getConfidenceRegion(gmins[plotting_pvals.size()-p-1],gmaxs[plotting_pvals.size()-p-1],measured_val);
                //std::cout<<"-- CL: "<<plotting_pvals[p]<<"  Sigma: "<<sqrt(2)*TMath::ErfInverse(plotting_pvals[p])<<Form("  ConfidenceInterval: [%2.2f, %2.2f]",reg[0],reg[1])<<std::endl;
                if(p == 0){ 
                    bf68_lo = reg[0];
                    bf68_up = reg[1];
                }
                if(p == 1){ 
                    bf95_lo = reg[0];
                    bf95_up = reg[1];
                }
                if(p == 2){
                    bf99_lo = reg[0];
                }

            }
            break;
        }
        //std::cout << "bf68_lo = " << bf68_lo << ", bf68_up = " << bf68_up << ", bf95_lo = " << bf95_lo << ", bf95_up = " << bf95_up << std::endl; 
        Double_t delta_chi2; Double_t bf_gridvalue;
        TTree* ttree_bf = (TTree*)fin->Get(Form("ttree_%d",bf_step));
        
        ttree_bf->SetBranchAddress("delta_chi2",&delta_chi2);
        ttree_bf->SetBranchAddress("bf_gridvalue",&bf_gridvalue);
        
        TH2D *h2_DeltaChi2 = new TH2D("DeltaChi2","DeltaChi2",50,0,4,50,0,50);
        int nentries = ttree_bf->GetEntries();



        ///Delta Chi2 Plots
        TGraph *g = new TGraph(nentries);
        std::vector< std::pair<double,double> > bfval_deltachi2_v;
        for (int j = 0; j < nentries; ++j){
            std::pair<double,double> temp;
            ttree_bf->GetEntry(j);
            h2_DeltaChi2->Fill(bf_gridvalue,delta_chi2);
            g->SetPoint(j,bf_gridvalue,delta_chi2);
            //pDeltaChi2[j] = delta_chi2; pbf_gridvalue[j]=bf_gridvalue;
            temp.first = bf_gridvalue; temp.second = delta_chi2; bfval_deltachi2_v.push_back(temp);
        }

          sort(bfval_deltachi2_v.begin(), bfval_deltachi2_v.end());
          std::vector<double> delta_chi2_v; std::vector<double> bf_gridvalue_v;
          double old_val = bfval_deltachi2_v[0].first;
          bf_gridvalue_v.push_back(bfval_deltachi2_v[0].first);
          
         
          std::vector<double> median_dchi2_v;
          for (size_t i = 0; i < bfval_deltachi2_v.size(); ++i){
            if (bfval_deltachi2_v[i].first != old_val){
              old_val = bfval_deltachi2_v[i].first;
              delta_chi2_v.push_back(Median(median_dchi2_v));
              bf_gridvalue_v.push_back(bfval_deltachi2_v[i].first);
              median_dchi2_v.clear();
              median_dchi2_v.push_back(bfval_deltachi2_v[i].second);
            }else{
              median_dchi2_v.push_back(bfval_deltachi2_v[i].second);
            }
          }
          //delta_chi2_v.push_back(Median(median_dchi2_v));
          int vsize = delta_chi2_v.size();

  

  

          double pDeltaChi2[vsize]; double pbf_gridvalue[vsize];
          for (int i = 0; i < vsize; ++i){
            pDeltaChi2[i] = delta_chi2_v[i]; pbf_gridvalue[i] = bf_gridvalue_v[i];
            
          }

          double yMax = 10;
          TGraph *pDeltaChi2Plot = new TGraph(vsize,pbf_gridvalue,pDeltaChi2);
   
           TCanvas * cDeltaChi2Plot = new TCanvas("DeltaChi2Plot","DeltaChi2Plot",1000,800);
           //pDeltaChi2Plot->GetXaxis()->SetRange(0,100); // IVAN CHANGE
           pDeltaChi2Plot->GetXaxis()->SetLimits(0.,maxpt);
           //pDeltaChi2Plot->GetXaxis()->SetLimits(2.,4.);
           pDeltaChi2Plot->SetTitle("#Delta#chi^{2} vs. True LEE Signal Scale Factor");
           //h2_DeltaChi2->SetStats(0);
           pDeltaChi2Plot->GetXaxis()->SetTitle("True Signal Strength (#mu)");
           pDeltaChi2Plot->GetYaxis()->SetTitle("#Delta#chi^{2}");
           pDeltaChi2Plot->GetYaxis()->SetRangeUser(0,yMax);
           
           //std::cout << "IVAN ---->>> " << pbf_gridvalue[vsize-1] << std::endl;
           
           //h2_DeltaChi2->Draw("COLZ");
           pDeltaChi2Plot->Draw("A*");

           //pDeltaChi2Plot->GetXaxis()->SetRangeUser(0,pbf_gridvalue[vsize]);
           //pDeltaChi2Plot->Draw("A*");
           //cDeltaChi2Plot->Update();
           //TLine l_bf68_lo(bf68_lo, 0, bf68_lo, 50);
           //TLine l_bf68_up(bf68_up, 0, bf68_up, 50);
           

           TBox bf68(bf68_lo, 0, bf68_up, yMax);
           //TBox bf68(bf68_lo, 0, 4, yMax); // Hacky way
           bf68.SetFillColorAlpha(kRed-9, 0.2);
           bf68.Draw("SAME");

           TBox bf95_0(bf95_lo, 0, bf68_lo, yMax);
           bf95_0.SetFillColorAlpha(kBlue-9, 0.2);
           bf95_0.Draw("SAME");

           TBox bf95_1(bf68_up, 0, bf95_up, yMax);
           bf95_1.SetFillColorAlpha(kBlue-9, 0.2);
           bf95_1.Draw("SAME");
          
           TBox bf99_0(bf99_lo, 0, bf95_lo, yMax);
           bf99_0.SetFillColorAlpha(kGreen-9, 0.2);
           //bf99_0.Draw("SAME");

            
            TLine l_MB_par(MB_val, 0, MB_val, yMax);
            l_MB_par.SetLineStyle(0);
            l_MB_par.SetLineWidth(3);
            l_MB_par.SetLineColor(kGreen+3);
            l_MB_par.Draw("same");

            TBox MBband_par(MB_val - 0.09, 0, MB_val + 0.09, yMax);
            MBband_par.SetFillColorAlpha(kGreen-3, 0.3);
            MBband_par.Draw("SAME");

            TLine l_BFval(BF_val, 0, BF_val, yMax);
            l_BFval.SetLineStyle(2);
            l_BFval.SetLineWidth(2);
            l_BFval.SetLineColor(kRed);

            l_BFval.Draw("same");

           //l_bf68_lo.Draw("SAME");
           //l_bf68_up.Draw("SAME");

           cDeltaChi2Plot->SaveAs(Form("DeltaChi2Plot_%s.pdf",tag.c_str()));
        
           
           
           double bfallmu_68_lo = 0; double bfallmu_68_up = 0; double bfallmu_95_lo = 0; double bfallmu_95_up = 0;// double bf99_lo = 0;
           std::cout<<"**************** Feldman Cousins 1D Confidence Intervals  **********************"<<std::endl;
           for(int i=0; i<v_true.size(); i++){
                if (BF_AllMu != vec_grid[i][0]) continue;
                std::cout<<"Grid Pt: "<<i<<", ScaleFactor: "<<vec_grid[i][0]<<std::endl;
                bf_step = i;
                for(int p=0; p< plotting_pvals.size();++p){
                    double measured_val = vec_grid[i][0];
                    std::vector<double> reg = myfeld.getConfidenceRegion(gmins[plotting_pvals.size()-p-1],gmaxs[plotting_pvals.size()-p-1],measured_val);
                    std::cout<<"-- CL: "<<plotting_pvals[p]<<"  Sigma: "<<sqrt(2)*TMath::ErfInverse(plotting_pvals[p])<<Form("  ConfidenceInterval: [%2.2f, %2.2f]",reg[0],reg[1])<<std::endl;
                    if(p == 0){ 
                        bfallmu_68_lo = reg[0];
                        bfallmu_68_up = reg[1];
                    }
                    if(p == 1){ 
                        bfallmu_95_lo = reg[0];
                        bfallmu_95_up = reg[1];
                    }

                }
                break;
            }
           if(bfallmu_95_up < bfallmu_68_up) bfallmu_95_up = maxpt;
           if(bfallmu_95_lo > bfallmu_68_lo) bfallmu_95_lo = 0.;
           TGraph *pDeltaChi2Plot_data = new TGraph(num_of_gridpts,gptval_data,dchi2_data);
    
           TCanvas * cDeltaChi2Plot_data = new TCanvas("DeltaChi2Plot_DATA_Comp","_DATA_Comp",1000,800);
           cDeltaChi2Plot_data->cd();
           TPad *pad_data_dchi2 = new TPad("pad_data_dchi2", "pad_data_dchi2", 0, 0, 0.85, 1.0);
           
           pad_data_dchi2->Draw();
           pad_data_dchi2->cd();
           //pDeltaChi2Plot->GetXaxis()->SetRange(0,100); // IVAN CHANGE
           //pDeltaChi2Plot_data->GetXaxis()->SetLimits(0.,maxpt);
           pDeltaChi2Plot_data->GetXaxis()->SetLimits(0.,maxpt);
           pDeltaChi2Plot_data->SetTitle("#Delta#chi^{2} vs. Signal Scale Factor [ " + SampleName + " - " + FDSName + " - " + TestStatName + " ]");
           //h2_DeltaChi2->SetStats(0);
           pDeltaChi2Plot_data->GetXaxis()->SetTitle("Signal Strength (#mu)");
           pDeltaChi2Plot_data->GetYaxis()->SetTitle("#Delta#chi^{2}");//delta_chi2 = this_chi2(mu(t)) - chi2_min(mu=mu^)
           pDeltaChi2Plot_data->GetYaxis()->SetRangeUser(0,yMax);
           
           //std::cout << "IVAN ---->>> " << pbf_gridvalue[vsize-1] << std::endl;
           
           //h2_DeltaChi2->Draw("COLZ");
           pad_data_dchi2->cd();
           pDeltaChi2Plot_data->Draw("A*");

           //pDeltaChi2Plot->GetXaxis()->SetRangeUser(0,pbf_gridvalue[vsize]);
           //pDeltaChi2Plot->Draw("A*");
           //cDeltaChi2Plot->Update();
           //TLine l_bf68_lo(bf68_lo, 0, bf68_lo, 50);
           //TLine l_bf68_up(bf68_up, 0, bf68_up, 50);
           

           //TBox bfallmu_68(bfallmu_68_lo, 0, bfallmu_68_up, yMax);
           TBox* bfallmu_68 = new TBox(bfallmu_68_lo, 0, bfallmu_68_up, yMax);
           //TBox bf68(bf68_lo, 0, 4, yMax); // Hacky way
           bfallmu_68->SetFillColorAlpha(kRed-9, 0.2);
           bfallmu_68->Draw("SAME");

           //TBox bfallmu_95_0(bfallmu_95_lo, 0, bfallmu_68_lo, yMax);
           TBox* bfallmu_95_0 = new TBox(bfallmu_95_lo, 0, bfallmu_68_lo, yMax);
           bfallmu_95_0->SetFillColorAlpha(kBlue-9, 0.2);
           bfallmu_95_0->Draw("SAME");

           TBox* bfallmu_95_1 = new TBox(bfallmu_68_up, 0, bfallmu_95_up, yMax);
           bfallmu_95_1->SetFillColorAlpha(kBlue-9, 0.2);
           bfallmu_95_1->Draw("SAME");
          
//           TBox bf99_0(bf99_lo, 0, bf95_lo, yMax);
           //bf99_0.SetFillColorAlpha(kGreen-9, 0.2);
           //bf99_0.Draw("SAME");

            
 //               TLine l_MB_par(MB_val, 0, MB_val, yMax);
 //               l_MB_par.SetLineStyle(0);
 //               l_MB_par.SetLineWidth(3);
 //               l_MB_par.SetLineColor(kGreen+3);
            l_MB_par.Draw("same");

  //          TBox MBband_par(MB_val - 0.09, 0, MB_val + 0.09, yMax);
  //          MBband_par.SetFillColorAlpha(kGreen-3, 0.3);
            MBband_par.Draw("SAME");

  //          TLine l_BFval(BF_val, 0, BF_val, yMax);
  //          l_BFval.SetLineStyle(2);
  //          l_BFval.SetLineWidth(2);
  //          l_BFval.SetLineColor(kRed);
            TLine* l_BF_AllMu = new TLine(BF_AllMu, 0, BF_AllMu, yMax);
            l_BF_AllMu->SetLineStyle(2);
            l_BF_AllMu->SetLineWidth(2);
            l_BF_AllMu->SetLineColor(kRed);
        
            l_BF_AllMu->Draw("same");
            //l_probs->Draw("same");
           //l_bf68_lo.Draw("SAME");
           //l_bf68_up.Draw("SAME");
            
            //pad_data_dchi2->Update();
            cDeltaChi2Plot_data->cd();
            TPad *pad_data_dchi2_leg = new TPad("pad_data_dchi2_leg", "pad_data_dchi2_leg", 0.8, 0, 1, 1);
            pad_data_dchi2_leg->Draw("SAME");
            pad_data_dchi2_leg->cd();
            TLegend * l_data_dchi2 = new TLegend(0.11,0.59,0.89,0.89);//69 was 29
            l_data_dchi2->SetHeader("#splitline{Confidence}{Level}");
            l_data_dchi2->AddEntry(bfallmu_68, "68%", "f");
            l_data_dchi2->AddEntry(bfallmu_95_0, "95%","f");
            l_data_dchi2->AddEntry(l_BF_AllMu, Form("#splitline{Best}{Fit = %1.2f}",BF_AllMu),"lf");
            l_data_dchi2->SetTextSize(0.15);
            l_data_dchi2->SetLineColor(kWhite);
            l_data_dchi2->SetLineWidth(0);
            l_data_dchi2->Draw();
            pad_data_dchi2_leg->SetBottomMargin(0.2);
            
                  // padl becomes the current pad
            

           cDeltaChi2Plot_data->SaveAs(Form("DeltaChi2Plot_DATA_Comp_%s.pdf",tag.c_str()));
        

           

           

        

    }else{
        std::cout<<"The mode you asked for ("<<mode_option<<") is not available, the available options are.."<<std::endl;
        runHelp();

    }
    std::cout << "Fin. Total wall time: " << difftime(time(0), start_time)/60.0 << " Minutes.\n";

    return 0;

}


