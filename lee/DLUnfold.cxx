#include <getopt.h>

#include "TTree.h"

#include "SBNcls.h"
#include "SBNchi.h"
#include "SBNspec.h"
#include "prob.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;

int main(int argc, char* argv[])
{
  std::string xml = "";
  int iarg = 0;
  opterr=1;
  int index;
  bool gen = false;
  bool numudis = false;
  bool combined = false;
  int mass_start = -1;

  const struct option longopts[] =
    {
      {"xml", 		required_argument, 	0, 'x'},
      {"gen",	no_argument, 0, 'g'},
      {"dis",	no_argument,0,'d'},
      {"comb", no_argument,0,'c'},
      {"part", required_argument,0,'p'},
      {0,			no_argument, 		0,  0},
    };

  while(iarg != -1)
    {
      iarg = getopt_long(argc,argv, "x:bscp:g", longopts, &index);

      switch(iarg)
	{
	case 'x':
	  xml = optarg;
	  break;
	case 'g':
	  gen = true;
	  break;
	case 'd':
	  numudis = true;
	  break;
	case 'c':
	  combined = true;
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

  std::string tag = "DL";

  // If we want to do some tests with efficiency and purity
  double eff_scale = 1.;
  double bkg_scale = 1.;	


  // Let's do this quick and sloppy
  int nBins = 4;
  double binedges[] = {100,500,900,1300,1800};
  //double binedges[] = {0,200,250,300,350,400,450,500,600,800,3000};
  TH1D * sigUnfolded = new TH1D("sigUnfolded","",nBins,binedges);

  // Fill up our histogram with signal events weighted from our friend tree
  TFile *infile = new TFile("/uboone/data/users/yatesla/othersys/inputs_to_sbnfit/arxiv/input_to_sbnfit_August2018.root","READ");
  TTree *TIntrinsic = (TTree*)infile->Get("nue_intrinsic_tree");
  // Grab true energy branch
  double ereco, etrue,unfoldedWgt;

  std::map<std::string, std::vector<double> > *f_weights;
  f_weights = new std::map<std::string, std::vector<double>>;

  TIntrinsic->AddFriend("nue_intrinsic_tree","unfolder.root");
  TIntrinsic->SetBranchAddress("reco_energy",&ereco);
  TIntrinsic->SetBranchAddress("true_nu_energy",&etrue);
  TIntrinsic->SetBranchAddress("unfoldedWgt",&unfoldedWgt);
  TIntrinsic->SetBranchAddress("weights", &f_weights);
  
  std::cout << "Creating new histo for unfolded signal" << std::endl;
  
  double POTscale = .572;
  for(int i = 0; i < TIntrinsic->GetEntries(); i++){
    TIntrinsic->GetEntry(i);
    
    double weight = f_weights->at("bnbcorrection_FluxHist").at(0) * unfoldedWgt * POTscale * eff_scale;
    sigUnfolded->Fill(ereco,weight);
    std::cout << "bnbcor: " << f_weights->at("bnbcorrection_FluxHist").at(0) << ", unfoldedwgt: " << unfoldedWgt  << ", POTscale: " << POTscale  << ", eff_scale: " << eff_scale << std::endl;
  }
  
  SBNspec bkg("DL.SBNspec.root",xml);
  SBNspec sig("DL.SBNspec.root",xml);
  sig.Add("nu_uBooNE_elike_lee",sigUnfolded);
	
  // Apply scalings
  sig.Scale("intrinsic",eff_scale);	
  sig.Scale("cocktail",bkg_scale);
  sig.Scale("extbnb",bkg_scale);
	
  bkg.Scale("intrinsic",eff_scale);	
  bkg.Scale("cocktail",bkg_scale);
  bkg.Scale("extbnb",bkg_scale);


  // Stats only
  TMatrixD *cov;
  SBNchi uboone_chi_statsonly(bkg,true);

  // Stats + sys
  TFile * fsys = new TFile("DL.SBNcovar.root","read");
  cov = (TMatrixD*)fsys->Get("frac_covariance_DL");
  SBNchi uboone_chi(bkg,*cov);

  int num_MC_events = 1000000;	

  SBNcls cls_factory(&bkg, &sig,*cov);
  cls_factory.SetSampleCovariance();

  cls_factory.CalcCLS(num_MC_events, tag);

  double chi2 = uboone_chi.CalcChi(&sig);
  double chi2_statsonly = uboone_chi_statsonly.CalcChi(&sig);

  bkg.WriteOut("bkg");
  sig.WriteOut("sig");

  std::cout << "CHI2: " << chi2 << std::endl;
  std::cout << "CHI2 STATSONLY: " << chi2_statsonly << std::endl;
	
  return 0;
}
