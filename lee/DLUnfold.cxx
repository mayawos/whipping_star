#include <getopt.h>

#include "TTree.h"

#include "SBNcls.h"
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
	case '?':
	case 'h':
	  std::cout<<"Allowed arguments:"<<std::endl;
	  std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
	  return 0;
	}
    }

  std::string tag = "DL";
  bool sample_from_covariance = true;
  int num_MC_events = 2e7;

  SBNspec sig("DL.SBNspec.root",xml);
  sig.Scale("fullosc",0.0);
  sig.Scale("oscfull",0.0);
  
  SBNspec bkg("DL.SBNspec.root",xml);
  bkg.Scale("fullosc",0.0);
  bkg.Scale("oscfull",0.0);
  bkg.Scale("signal",0.0);
  
  // Stats + sys
  TFile * fsys = new TFile("DL.SBNcovar.root","read");
  TMatrixD * cov = (TMatrixD*)fsys->Get("frac_covariance_DL");
  
  SBNcls cls_factory(&bkg, &sig, *cov);
  
  if(sample_from_covariance) cls_factory.SetSampleCovariance();
  cls_factory.CalcCLS(num_MC_events, tag);
	
  return 0;
}
