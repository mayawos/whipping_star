#include <getopt.h>
#include <limits>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TStopwatch.h"

#include "SBNgenerate.h"

#include "DLUtil.h"

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
  float fake_m2 = -1;
  float fake_sin22th = -1;
  bool vary_background = false;

  const struct option longopts[] =
    {
      {"xml"        , required_argument, 0, 'x'},
      {"m2"         , required_argument, 0, 'm'},
      {"sin22th"    , required_argument, 0, 's'},
      {"background" , no_argument      , 0, 'b'},
      {0            , no_argument      , 0,  0 },
    };
  
  while(iarg != -1)
    {
      iarg = getopt_long(argc,argv, "x:m:s:b", longopts, &index);

      switch(iarg)
	{
	case 'x':
	  xml = optarg;
	  break;
	case 'm':
	  fake_m2 = std::atof(optarg);
	  break;
	case 's':
	  fake_sin22th = std::atof(optarg);
	  break;
	case 'p':
	  vary_background = true;
	  break;
	case '?':
	case 'h':
	  std::cout << "{\"xml\"        , required_argument, 0, \'x\'}," << std::endl;
	  std::cout << "{\"m2\"         , required_argument, 0, \'m\'}," << std::endl;
	  std::cout << "{\"sin22th\"    , required_argument, 0, \'s\'}," << std::endl;
	  std::cout << "{\"background\" , no_argument      , 0, \'b\'}," << std::endl;
	  std::cout << "{0            , no_argument      , 0,    0 }," << std::endl;

	  return 0;
	}
    }

  std::stringstream suffix_ss;
  suffix_ss << "_";
  suffix_ss << std::fixed << std::setprecision(2) << fake_m2;
  suffix_ss << "_";
  suffix_ss << std::fixed << std::setprecision(5) << fake_sin22th;
  std::string suffix = suffix_ss.str();
  std::cout << "suffix=" << suffix << std::endl;

  std::string tag = "DL";

  // No oscilaltions model 
  NeutrinoModel nullModel(0,0,0);
  
  // Need this
  SBNgenerate * bkgo = new SBNgenerate(xml,nullModel);
  SBNspec bbkg = bkgo->spec_central_value;
  bbkg.Scale("oscfull",0.0);
  bbkg.Scale("fullosc",0.0);
  bbkg.Scale("signal",0.0);
  bbkg.WriteOut("DLSens_Bkg");

  // Background spectrum with fullosc=0 oscfull=0 signal=0
  SBNosc osctrue("DLSens_Bkg.SBNspec.root", xml, nullModel);
  osctrue.Scale("fullosc",0.0);
  osctrue.Scale("oscfull",0.0);
  osctrue.Scale("signal",0.0);
  osctrue.has_been_scaled = false;
  osctrue.is_verbose = false;

  TStopwatch twatch;
  twatch.Stop();
  twatch.Reset();
  
  // Create our background
  SBNspec bkg("DL.SBNspec.root",xml,false);
  bkg.Scale("signal",0.0);
  bkg.Scale("fullosc",0.0);
  bkg.Scale("oscfull",0.0);
  bkg.WriteOut(std::string("DLChi2Draw_bkg_spec") + suffix);

  // Create our fake data, do it this odd way w/ disk write
  float fake_mnu = std::sqrt(fake_m2);
  float fake_ue  = std::pow(fake_sin22th/float(4),.5);
  float fake_um  = 1;

  std::cout << "Generating Fake model" << std::endl;
  twatch.Start();
  NeutrinoModel fake_model(fake_mnu, 1, 1);
  SBNgenerate fake_gen(xml,fake_model);
  fake_gen.spec_osc_sinsq.Scale("signal",0.0);
  fake_gen.spec_osc_sinsq.Scale("fullosc",0.0);
  fake_gen.spec_osc_sin.Scale("signal",0.0);
  fake_gen.spec_osc_sin.Scale("fullosc",0.0);  
  fake_gen.WritePrecomputedOscSpecs(tag);
  twatch.Stop();
  std::cout << "done @ t=" << twatch.RealTime() << std::endl;
  twatch.Reset();

  std::cout << "Generating Fake signal" << std::endl;
  twatch.Start();
  SBNosc fake_signal = osctrue;
  fake_model = NeutrinoModel(fake_mnu, fake_ue, fake_um);
  fake_signal.SetBothMode();
  fake_signal.LoadModel(fake_model);
  fake_signal.OscillateThis(tag);
  SBNspec data = static_cast<SBNspec>(fake_signal);
  data.Scale("intrinsic",0.0);
  data.Scale("cocktail",0.0);
  data.Scale("extbnb",0.0);
  data.Scale("fullosc",0.0);
  data.Scale("signal",0.0);
  data.WriteOut(std::string("DLChi2Draw_data_spec") + suffix);
  twatch.Stop();
  std::cout << "done @ t=" << twatch.RealTime() << std::endl;
  twatch.Reset();
  
  // Stats + sys
  TFile * fsys = new TFile("DL.SBNcovar.root","read");
  TMatrixD *cov = (TMatrixD*)fsys->Get("frac_covariance_DL");

  // Chi2 using background only covariance
  SBNchi BF_chi(bkg,*cov);
  BF_chi.core_spectrum.Add(&data);
  BF_chi.core_spectrum.CollapseVector();

  auto bkg_and_data = BF_chi.core_spectrum;
  bkg_and_data.CalcFullVector();
  bkg_and_data.CollapseVector();
  bkg_and_data.WriteOut(std::string("DLChi2Draw_bkg_and_data_spec") + suffix);

  SBNchi true_chi(bkg_and_data,*cov);

  std::string out_file_name = "DLChi2Draw_ana";
  out_file_name += suffix;
  out_file_name += std::string(".root");
  TFile* ana_file = TFile::Open(out_file_name.c_str(),"RECREATE");
  ana_file->cd();
  
  int fexp = 0;
  double _chi2 = -1;

  TTree* outtree1 = new TTree("ana_chi","");
  outtree1->Branch("fexp" , &fexp , "fexp/I");
  outtree1->Branch("chi2" , &_chi2, "chi2/D");


  int n_fexp = 1e6;

  std::cout << "Getting spec" << std::endl;
  twatch.Start();
  auto spec_fexp_v = BF_chi.SpecReturnSCVI(bkg_and_data, n_fexp, xml);
  // auto spec_fexp_v = BF_chi.SpecReturnSPVI(bkg_and_data, n_fexp, xml);
  twatch.Stop();
  std::cout << "done @ t=" << twatch.RealTime() << std::endl;
  twatch.Reset();

  for(; fexp < n_fexp; ++fexp) {

    auto& spec_fexp = spec_fexp_v[fexp];
    // std::stringstream ss;
    // ss << "spec_output_" << fexp;
    // spec_fexp.WriteOut(ss.str());
    
    // draw fake data
    twatch.Start();
    
    BF_chi.core_spectrum = spec_fexp;
    BF_chi.core_spectrum.CollapseVector();
    
    spec_fexp.CollapseVector();
    _chi2 = true_chi.CalcChi(spec_fexp);
    
    std::cout << "_chi2=" << _chi2 << std::endl;

    outtree1->Fill();

    twatch.Stop();
    std::cout << "@ fexp=" << fexp << "/" << n_fexp << " : " << twatch.RealTime() << std::endl;
    twatch.Reset();
    
  }

  ana_file->cd();
  outtree1->Write();
  ana_file->Close();

  return 0;
}
