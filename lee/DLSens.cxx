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


  //PART 1: precompute all of our sin and sin2 amplitudes so we don't need to later
  if(gen){

    // write background spectrum, fullosc=0
    NeutrinoModel nullModel(0, 0, 0);
    SBNgenerate * bkgo = new SBNgenerate(xml,nullModel);
    SBNspec bkg = bkgo->spec_central_value;
    bkg.Scale("fullosc",0.0);
    bkg.Scale("signal",0.0);

    bkg.WriteOut(tag+"_Bkg");

    // write precomputed spectrum with fullosc sample
    float mnu;
    for(int mi = 0; mi < 50; mi++){
      mnu = pow(10.,(float(mi)/50.*TMath::Log10(10./.1) + TMath::Log10(.1)));

      //Model: mnu, ue4, um4
      //we're precomputing, so we don't really care about the u's
      NeutrinoModel testModel(mnu, 1, 1);

      // on construction it makes 3 SBNspecs, 1 sin amp, 1 sin2 amp, 1 CV oscilatted
      SBNgenerate * gen = new SBNgenerate(xml,testModel);
      gen->spec_central_value.Scale("signal",0.0);

      // Write them to file
      gen->WritePrecomputedOscSpecs(tag);
    }
    return 1;
  }

  //PART  2: Now that sin and sin2 libs are generated, calculate that sensitivity
  if(!gen){
	
    // Create our unoscillated background
    SBNspec bkg(tag+"_Bkg.SBNspec.root",xml);
    bkg.Scale("fullosc",0.0);
    bkg.Scale("signal",0.0);
		
    //Bring in our covariance matrix!
    // Stats only.
    TMatrixD *cov;
    SBNchi uboone_chi_statsonly(bkg,true);

    // Stats + sys
    TFile * fsys = new TFile("DL.SBNcovar.root","read");
    cov = (TMatrixD*)fsys->Get("frac_covariance_DL");
    std::cout << "Frac Covariance Matrix" << std::endl;
    cov->Print();
    SBNchi uboone_chi(bkg,*cov);
    
    // Load up oscillation model
    NeutrinoModel nullModel(0,0,0);
    SBNosc osctrue(tag+"_Bkg.SBNspec.root",xml, nullModel);
    osctrue.Scale("signal",0.0);

    if(numudis){
      // If we're doing numu disappearance:
      std::cout << "NUMU DISAPPEARANCE" <<  std::endl;
      float mnu, um, ue;
      for(int mi = 0; mi < 50; mi++){
	for(int umi = 0; umi < 50; umi++){

	  um = umi/float(100)*(.5);
	  mnu = pow(10.,(mi/float(50)*TMath::Log10(10./.1) + TMath::Log10(.1)));
	  NeutrinoModel testModel(mnu,0,um);
	  std::cout << "NU MODEL: " << mnu << " " << ue << " " << um << std::endl;

	  SBNosc osc = osctrue;
	  osc.LoadModel(testModel);
	  osc.OscillateThis(tag);

	  double chi2 = uboone_chi.CalcChi(&osc);
	  double chi2_statsonly = uboone_chi_statsonly.CalcChi(&osc);
	  std::cout << "COUNT: " << (mi*100 + umi)/float(100) << "%" << std::endl;
	  std::cout << "ANS: " << pow(mnu,2) << " " << um << " " << 4*um*um*(1-um*um) << " " << chi2 << std::endl;
	  std::cout << "ANS_STATSONLY: " << pow(mnu,2) << " " << um << " " << 4*um*um*(1-um*um) << " " << chi2_statsonly << std::endl;
	}
      }
    }
    else if(combined){
      std::cout << "COMBINED FIT" <<  std::endl;
      float mnu, um, ue;
      int mass_end;
      if(mass_start > 0){
	mass_end = mass_start + 1;
      }
      else{
	mass_start = 0;
	mass_end = 50;
      }

      for(int mi = mass_start; mi < mass_end; mi++){
	for(int umi = 0; umi < 25; umi++){
	  for(int uei = 0; uei < 25; uei++){

	    um = pow(10.,(umi/float(25)*TMath::Log10(1./1e-3) + TMath::Log10(1e-3)));
	    ue = pow(10.,(uei/float(25)*TMath::Log10(1./1e-3) + TMath::Log10(1e-3)));
	    mnu = pow(10.,(mi/float(50)*TMath::Log10(10./.1) + TMath::Log10(.1)));
	    NeutrinoModel testModel(mnu,ue,um);
	    std::cout << "NU MODEL: " << mnu << " " << ue << " " << um << std::endl;

	    SBNosc osc = osctrue;
	    osc.LoadModel(testModel);
	    osc.OscillateThis(tag);

	    if(umi == 7 && mi == 7 ){
	      osc.WriteOut("osctest");
	      bkg.WriteOut("bkgtest");
	    }

	    double chi2 = uboone_chi.CalcChi(&osc);
	    double chi2_statsonly = uboone_chi_statsonly.CalcChi(&osc);
	    std::cout << "COUNT: " << (mi*25*25 + umi*25 + uei)/float(50*25*25/100) << "%" << std::endl;
	    std::cout << "ANS: " << pow(mnu,2) << " " << um << " " << ue << " " << chi2 << std::endl;
	    std::cout << "ANS_STATSONLY: " << pow(mnu,2) << " " << um << " " << ue << " " << chi2_statsonly << std::endl;
	  }
	}
      }
    } // end combined

    //
    // Nue appearance
    //
    else{
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
    }

    return 0;
  }
}
