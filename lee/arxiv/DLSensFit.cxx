#include <getopt.h>

#include "SBNchi.h"
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBNgenerate.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;

int main(int argc, char* argv[])
{
  std::string xml = "numu_disp.xml";
  int iarg = 0;
  opterr = 1;
  int index;
  bool gen = false;
  bool numudis = false;
  bool combined = false;
  int mass_start = -1;

  const struct option longopts[] =
    {
      {"xml" , required_argument, 0, 'x'},
      {"gen" , no_argument      , 0, 'g'},
      {"dis" , no_argument      , 0, 'd'},
      {"comb", no_argument      , 0, 'c'},
      {"part", required_argument, 0, 'p'},
      {0     , no_argument      , 0,  0 },
    };
  
  while(iarg != -1)
    {
      iarg = getopt_long(argc,argv, "x:dscp:g", longopts, &index);

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
	  std::cout << "Allowed arguments:" << std::endl;
	  std::cout << "\t-x\t--xml\t\tInput .xml file for SBNconfig" << std::endl;
	  std::cout << "\t-g\t--gen\t\t..." << std::endl;
	  std::cout << "\t-d\t--dis\t\t..." << std::endl;
	  std::cout << "\t-c\t--comb\t\t..." << std::endl;
	  std::cout << "\t-p\t--part\t\t..." << std::endl;
	  return 0;
	}
    }

  std::string tag = "DL";

  if (gen) {
    std::cout << "You cannot use GEN with this code" << std::endl;
    return 1;
  }

  //PART  2: Now that sin and sin2 libs are generated, calculate that sensitivity
  if(!gen){

    // Load up oscillation model
    NeutrinoModel nullModel(0,0,0);
		
    // Create our data with fullosc=0
    SBNspec data("DL.SBNspec.root",xml);
    data.Scale("fullosc",0.0);
    data.Scale("fullosc2",0.0);
 
    // Stats + sys
    TFile * fsys = new TFile("DL.SBNcovar.root","read");
    TMatrixD *cov = (TMatrixD*)fsys->Get("frac_covariance_DL");

    // Chi with data
    SBNchi uboone_chi(data,*cov);
    SBNchi uboone_chi_statsonly(data,true);
    
    // Background spectrum with fullosc=0 signal=0
    SBNosc osctrue(tag+"_Bkg.SBNspec.root",xml, nullModel);
    osctrue.Scale("fullosc",0.0);
    osctrue.Scale("fullosc2",0.0);
    osctrue.Scale("signal",0.0);

    // If we're doing nue appearance
    std::cout << "NUE APPEARANCE" << std::endl;
    float mnu, um, ue, sin22th;
    for(int mi = 0; mi < 50; mi++){
      for(int sin22thi = 0; sin22thi < 50; sin22thi++){

	sin22th = pow(10.,(sin22thi/float(50)*TMath::Log10(1./1e-5) + TMath::Log10(1e-5)));
	mnu = pow(10.,(mi/float(50)*TMath::Log10(10./.1) + TMath::Log10(.1)));
	ue = pow(sin22th/float(4),.5);
	um = 1.f;	// set sin22th(mu mu) = 0
	NeutrinoModel testModel(mnu,ue,um);

	std::cout << "NU MODEL: " << mnu << " " << ue << " " << um << std::endl;
					
	// osctrue = background spectrum
	SBNosc osc = osctrue;
	  
	// set numu disappearance, nue app & dis
	// but doesn't matter since um = 1, no numu disappearance
	osc.SetBothMode();
	  
	// set the oscillation model
	osc.LoadModel(testModel);
	  
	// oscillate the spectrum and add it to the bkg
	osc.OscillateThis(tag);

	// calculate chi w.r.t background uboone_chi w. stat err or w.o. stat err
	double chi2 = uboone_chi.CalcChi(&osc);
	double chi2_statsonly = uboone_chi_statsonly.CalcChi(&osc);

	std::cout << "COUNT: " << (mi*50 + sin22thi)/float(50*.50) << "%" << std::endl;
	std::cout << "ANS: " << std::setprecision(15) << pow(mnu,2) <<  " " << ue << " " << um << " " << 4*um*um*ue*ue << " " << chi2 << std::endl;
	std::cout << "ANS_STATSONLY: " << std::setprecision(15) << pow(mnu,2) <<  " " << ue << " " << um << " " << 4*um*um*ue*ue << " " << chi2_statsonly << std::endl;
      }
    }

    return 0;
  }
}