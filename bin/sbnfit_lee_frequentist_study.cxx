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
#include "TStyle.h"
#include "TError.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TGraph.h"

#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBNcls.h"
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
    int iarg = 0;
    opterr=1;
    int index;
    bool sample_from_covariance = true;
    bool sample_from_collapsed = false;
    bool remove_correlations = false;
    bool stats_only = false;
    int num_MC_events = 100000;
    bool use_cnp = false;
    int which_mode = 1;

    std::string tag = "TEST";

    std::string signal_file = "EMPTY";
    std::string background_file = "EMPTY";
    std::string covariance_file = "EMPTY";
    std::string fakedata_file = "EMPTY";
    std::string ext_err_file = "EMPTY";

    bool bool_flat_det_sys = false;
    bool bool_fill_det_sys = false; //PeLEE updates
    bool bool_ext_err = false; //PeLEE updates
    double flat_det_sys_percent = 0.0;

    bool zero_off_diag = false;
    bool tester=false;
    double epsilon = 1e-12;

    const struct option longopts[] =
    {
        {"xml", 		required_argument, 	0, 'x'},
        {"covariance", 	required_argument,0,'c'},
        {"collapse", 	required_argument,	0,'j'},
        {"signal", 		required_argument,	0,'s'},
        {"mode", 	    	required_argument,	0,'m'},
        {"background", 	required_argument,	0,'b'},
        {"fakedata", 	required_argument,	0,'d'},
        {"exterr", 	required_argument,	0,'r'},
        {"tag", 	    required_argument,	0,'t'},
        {"epsilon", required_argument,0,'e'},
        {"cnp",no_argument,0,'a'},
        {"zero",no_argument,0,'z'},
        {"tester",no_argument,0,'k'},
        {"poisson", no_argument,0,'p'},
        {"flat", required_argument,0,'f'},
        {"filldetsys", required_argument,0,'d'},
        {"help",no_argument,0,'h'},
        {0,			no_argument, 		0,  0},
    };

    while(iarg != -1)
    {
        iarg = getopt_long(argc,argv, "m:a:x:n:s:e:b:d:c:f:r:t:pjkzh", longopts, &index);

        switch(iarg)
        {
            case 'k':
                tester = true; 
                break;
            case 'z':
                remove_correlations = true; 
                break;
            case 'j':
                sample_from_collapsed = true;
                break;
            case 'm':
                which_mode = (int)strtod(optarg,NULL);
                break;
            case 'x':
                xml = optarg;
                break;
            case 's':
                signal_file = optarg;
                break;
            case 'b':
                background_file = optarg;
                break;
            case 'd':
                which_mode = 2;         //PeLEE specific hacks for fakedata
                fakedata_file = optarg;
                break;
            case 'r':
                ext_err_file = optarg;
                break;
            case 'f':
                bool_flat_det_sys = true;
                flat_det_sys_percent = (double)strtod(optarg,NULL);
                break;
            case 'y':
                bool_fill_det_sys = true;
                break;
            case 'e':
                epsilon = (double)strtod(optarg,NULL);
                break;
            case 't':
                tag = optarg;
                break;
            case 'a':
                use_cnp = true;
                break;
            case 'c':
                covariance_file = optarg;
                break;
            case 'n':
                num_MC_events = (int)strtod(optarg,NULL);
                break;
            case 'p':
                sample_from_covariance = false;
                break;
            case '?':
            case 'h':
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"sbnfit_lee_frequentist_study allows for the simple hypothesis testing."<<std::endl;
                std::cout<<"---------------------------------------------------"<<std::endl;
                std::cout<<"--- Required arguments: ---"<<std::endl;
                std::cout<<"\t-x\t--xml\t\t\tInput configuration .xml file for SBNconfig"<<std::endl;
                std::cout<<"\t-t\t--tag\t\t\tA unique tag to identify the outputs [Default to TEST]"<<std::endl;
                std::cout<<"\t-s\t--signal\t\tInput signal SBNspec.root file"<<std::endl;
                std::cout<<"\t-b\t--background\t\tInput background only SBNspec.root file"<<std::endl;
                std::cout<<"\t-d\t--fakedata\t\tInput fakedata SBNspec.root file"<<std::endl;
                std::cout<<"\t-c\t--covariance\t\tInput Fractional Covariance Matrix SBNcovar.root file. If not passed, defaults to stats only!"<<std::endl;
                std::cout<<"--- Optional arguments: ---"<<std::endl;
                std::cout<<"\t-j\t--collapse\t\tSample from collapsed rather than full covariance matrix (default false, experimental!)"<<std::endl;
                std::cout<<"\t-f\t--flat\t\tAdd a flat percent systematic to fractional covariance matrix (all channels) (default false, pass in percent, i.e 5.0 for 5\% experimental)"<<std::endl;
                std::cout<<"\t-r\t--exterr\t\tAdd zero ext bin error to the stat errors (all channels) (default no additional errors)"<<std::endl;
                std::cout<<"\t-z\t--zero\t\tZero out all off diagonal elements of the systematics covariance matrix (default false, experimental!)"<<std::endl;
                std::cout<<"\t-e\t--epsilon\t\tEpsilon tolerance by which to add back to diagonal of covariance matrix if determinant is 0 (default 1e-12)"<<std::endl;
                std::cout<<"\t-n\t--number\t\tNumber of MC events for frequentist studies (default 100k)"<<std::endl;
                std::cout<<"\t-m\t--mode\t\tMode for test statistics 0: absolute chi^2, 1: delta chi^2 (default Delta Chi| obsolete, runs all concurrently)"<<std::endl;
                std::cout<<"\t-p\t--poisson\t\tUse Poissonian draws for pseudo experiments instead of from covariance matrix"<<std::endl;
                std::cout<<"\t-h\t--help\t\t\tThis help menu."<<std::endl;
                std::cout<<"---------------------------------------------------"<<std::endl;
                return 0;	
        }
    }
    if(signal_file =="EMPTY"){
        std::cout<<"Error! You must enter a signal root file with the  `--signal  XX.SBNspec.root` or `-s XX.SBNspec.root` flags "<<std::endl;
        std::cout<<"Error! Run `--help` or `-h`  for more details."<<std::endl;
        return 1;
    }
    if(background_file =="EMPTY"){
        std::cout<<"Error! You must enter a background root file with the  `--background  XX.SBNspec.root` or `-b XX.SBNspec.root`  flags "<<std::endl;
        std::cout<<"Error! Run `--help` or `-h`  for more details."<<std::endl;
        return 1;
    }
    if(covariance_file =="EMPTY"){
        std::cout<<"Note! No covariance root file with the  `--covariance  XX.SBNcovar.root` or `-c XX.SBNcovar.root`. was passed, running in stats only mode. "<<std::endl;
        stats_only = true;
        sample_from_covariance = false;
    }
    if(ext_err_file =="EMPTY"){
        std::cout<<"Note! No ext bnb error covariance matrix root file with the  `--exterr  XX.SBNcovar.root` or `-r XX.SBNcovar.root`. was passed, running without stats error. "<<std::endl;
        bool_ext_err = false;
    }
    else bool_ext_err = true;

    std::cout<<"Loading signal file : "<<signal_file<<" with xml "<<xml<<std::endl;
    SBNspec sig(signal_file,xml);

    std::cout<<"Loading background file : "<<background_file<<" with xml "<<xml<<std::endl;
    SBNspec bkg(background_file,xml);

    std::cout<<"Loading fractional covariance matrix from "<<covariance_file<<std::endl;
    std::cout << "loading frac cov matrix" << std::endl;
    TFile * fsys;
    TMatrixD * cov;
    TFile * fexterr;
    TMatrixD * exterr;
    
    if(!stats_only){
        fsys = new TFile(covariance_file.c_str(),"read");
        cov = (TMatrixD*)fsys->Get("frac_covariance");
    }
    std::cout << "found frac cov matrix: " << cov << std::endl;

    //ext error
    std::vector<double> coll_ext_err_vec, frac_coll_ext_err_vec;
    std::vector<double> ext_err_vec, frac_ext_err_vec;
    coll_ext_err_vec.resize(bkg.num_bins_total_compressed);
    frac_coll_ext_err_vec.resize(bkg.num_bins_total_compressed);
    ext_err_vec.resize(bkg.num_bins_total);
    frac_ext_err_vec.resize(bkg.num_bins_total);
    std::fill(ext_err_vec.begin(), ext_err_vec.end(), 0.0);

    bkg.CalcFullVector();
    if(bool_ext_err){
       fexterr = new TFile(ext_err_file.c_str(),"read");
       exterr = (TMatrixD*)fexterr->Get("full_covariance");
       //create SBNchi object for the ext error so that it will be consistent when collapsed
       TMatrixD collext;
       collext.ResizeTo(bkg.num_bins_total_compressed, bkg.num_bins_total_compressed);
       SBNchi chiext(bkg,exterr);
       chiext.CollapseModes(*exterr, collext);
       for( int i=0; i < bkg.num_bins_total; i++ ){ ext_err_vec[i] = sqrt((*exterr)(i,i)); frac_ext_err_vec[i] = sqrt((*exterr)(i,i))/bkg.full_vector[i]; std::cout << "err: " << ext_err_vec[i] << std::endl; }
       //if(!stats_only){for( int i=0; i < bkg.num_bins_total; i++ ){ (*cov)(i,i) += frac_ext_err_vec[i]; }}
       bkg.CollapseVector();
       for( int i=0; i< collext.GetNcols(); i++ ){ coll_ext_err_vec[i] = sqrt(collext(i,i)); frac_coll_ext_err_vec[i] = sqrt(collext(i,i)/bkg.collapsed_vector[i]); std::cout << "err coll: " << coll_ext_err_vec[i] << std::endl; }
    }
    std::cout << "found ext error matrix: " << exterr << std::endl;


    //PeLEE hacks for incorporating fake data
    TFile * fdata;
    TH1D * h_fakedata_1eNp;
    TH1D * h_fakedata_1e0p;
    TH1D * h_fakedata_numu;
    
    if(fakedata_file != "EMPTY"){
      std::cout << "Loading fakedata file: " << fakedata_file << std::endl;
      fdata = new TFile(fakedata_file.c_str(),"read");
      h_fakedata_1eNp = (TH1D*)fdata->Get("nu_uBooNE_nue_data");
      h_fakedata_1e0p = (TH1D*)fdata->Get("nu_uBooNE_1e0p_data");
      h_fakedata_numu = (TH1D*)fdata->Get("nu_uBooNE_numu_data");
    }
    std::cout << "create vector of fakedata" << std::endl; 
    //create vector of fakedata:
    std::vector<float> fakedata;
    if( which_mode==2 ){for( int k=1; k < h_fakedata_1eNp->GetNbinsX()+1; k++ ) fakedata.push_back(h_fakedata_1eNp->GetBinContent(k)); }
    if( which_mode==2 ){for( int k=1; k < h_fakedata_1e0p->GetNbinsX()+1; k++ ) fakedata.push_back(h_fakedata_1e0p->GetBinContent(k)); }
    if( which_mode==2 ){for( int k=1; k < h_fakedata_numu->GetNbinsX()+1; k++ ) fakedata.push_back(h_fakedata_numu->GetBinContent(k)); }
    if( which_mode==2 ){for( int k=1; k < h_fakedata_1eNp->GetNbinsX()+1; k++ ) std::cout << "h_fakedata_1eNp->GetBinContent k " << k << " = " <<  h_fakedata_1eNp->GetBinContent(k) << std::endl; }
    if( which_mode==2 ){for( int k=1; k < h_fakedata_1e0p->GetNbinsX()+1; k++ ) std::cout << "h_fakedata_1e0p->GetBinContent k " << k << " = " <<  h_fakedata_1e0p->GetBinContent(k) << std::endl; }
    if( which_mode==2 ){for( int k=1; k < h_fakedata_numu->GetNbinsX()+1; k++ ) std::cout << "h_fakedata_numu->GetBinContent k " << k << " = " <<  h_fakedata_numu->GetBinContent(k) << std::endl; }
    
    std::cout << "size of fakedata vector: " << fakedata.size() << std::endl;
    //End of PeLEE hacks for incorporating fake data
    
    //PELEE hack -- use the actual PELEE diagonal errors for detsys instead of flat systematics
    TMatrixD frac_flat_matrix(bkg.num_bins_total, bkg.num_bins_total);
    
    //BDT
    //std::vector<double> np_detsys = {0.2007, 0.1077, 0.0922, 0.0574, 0.0663, 0.0755, 0.0721, 0.0872, 0.0975, 0.1034, 0.2551, 0.0849, 0.1428, 0.1764, 0.1806, 0.2042, 0.1841, 0.1688, 0.1811}; //old lower stats systematics
    std::vector<double> np_detsys = {0.2454, 0.1523, 0.1546, 0.0900, 0.1500, 0.0608, 0.1530, 0.0781, 0.1114, 0.1145, 0.2005, 0.1237, 0.1871, 0.1105, 0.1349, 0.1612, 0.2478};
    //1e0p
    //std::vector<double> zp_detsys = {0.0985, 0.0985, 0.1015, 0.1015, 0.1929, 0.1929, 0.2326, 0.2326, 0.3289, 0.3289, 0.1720, 0.1720, 0.1720, 0.1720, 0.1720, 0.1720, 0.1720, 0.1720, 0.1720}; //old lower stats systematics
    std::vector<double> zp_detsys = {0.1520, 0.1520, 0.0980, 0.0980, 0.1939, 0.1939, 0.4064, 0.4064, 0.3259, 0.3259, 0.1819, 0.1819, 0.1819, 0.1819, 0.2540, 0.2540, 0.2540, 0.2540, 0.2540};
    //numu
    //std::vector<double> numu_detsys = {0.096,0.097,0.066,0.051,0.065,0.093,0.081,0.07,0.109,0.122,0.142,0.158,0.18,0.261}; //old lower stats systematics
    std::vector<double> numu_detsys = {0.1655, 0.1044, 0.1233, 0.0574, 0.0500, 0.0849, 0.0686, 0.1449, 0.1158, 0.0849, 0.0980, 0.1536, 0.1556, 0.1879, 0.2987};

    if(bool_fill_det_sys){
      std::cout << "RUNNING with PELEE systematics!" << std::endl;
      
      frac_flat_matrix.ResizeTo(bkg.num_bins_total,bkg.num_bins_total);
      frac_flat_matrix.Zero();//(bkg.num_bins_total,bkg.num_bins_total);
      int j=0;     
      
      for(auto& h: bkg.hist){
        std::string hname = h.GetName();
        for( int i=1; i < h.GetNbinsX()+1; i++ ){
          if( (hname.find("1eNp_intrinsic") != std::string::npos) || (hname.find("1eNp_lee") != std::string::npos) ){
            frac_flat_matrix(j,j) = np_detsys[i-1]*np_detsys[i-1];
          }
          else if( (hname.find("1e0p_intrinsic") != std::string::npos) || ( hname.find("1e0p_lee") != std::string::npos ) ){
            frac_flat_matrix(j,j) = zp_detsys[i-1]*zp_detsys[i-1];
          }
          else if( hname.find("numu_bnb") != std::string::npos ){
            frac_flat_matrix(j,j) = numu_detsys[i-1]*numu_detsys[i-1];
          }
	  else if( hname.find("ext") == std::string::npos ){
            frac_flat_matrix(j,j) = 0.2*0.2;
          }
          j++; 
        }
      }
      std::cout<<"Just Before"<<std::endl;
      (*cov) = (*cov)+(frac_flat_matrix);
    }
    std::cout << "Done with systematics!" << std::endl;
    
    if(bool_flat_det_sys){
      std::cout << "RUNNING with flat systematics: " << flat_det_sys_percent << "%!" << std::endl;
      for(int i=0 ; i< bkg.num_bins_total; i++){
	frac_flat_matrix(i,i)=flat_det_sys_percent*flat_det_sys_percent/10000.;
      }
      std::cout<<"Just Before"<<std::endl;
      (*cov) = (*cov)+(frac_flat_matrix);
    }

    if(remove_correlations){
      std::cout<<"WARNING! We are running in   `Remove All Off Diagional Covariances/Correlations Mode` make sure this is what you want. "<<std::endl;
      for(int i=0; i<bkg.num_bins_total;i++){ 
	for(int j=0; j<bkg.num_bins_total;j++){ 
	  if(i==j)continue;
	  (*cov)(i,j) =0.0;
	}
      }
    }
    
    if(!stats_only){
        SBNcls cls_factory(&bkg, &sig, fakedata, *cov, ext_err_vec, coll_ext_err_vec);
        cls_factory.SetTolerance(epsilon);
        //cls_factory.SetAdditionalErrors(ext_err_vec);
        if(sample_from_collapsed)  cls_factory.SetSampleFromCollapsed();
        if(sample_from_covariance) cls_factory.SetSampleCovariance();
        cls_factory.setMode(which_mode);
        if(tester){cls_factory.runConstraintTest();return 0;}
        cls_factory.CalcCLS(num_MC_events, tag);
    }else{
        SBNcls cls_factory(&bkg, &sig, fakedata, ext_err_vec, coll_ext_err_vec);
        cls_factory.SetTolerance(epsilon);
        //std::cout << "SET ADDITIONAL ERRORS" << std::endl;
        //cls_factory.SetAdditionalErrors(ext_err_vec);
        if(sample_from_collapsed)  cls_factory.SetSampleFromCollapsed();
        cls_factory.setMode(which_mode);
        if(tester){cls_factory.runConstraintTest();return 0;}
        cls_factory.CalcCLS(num_MC_events, tag);
    }

    return 0;
}
