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
  
  size_t n_mi = 50;
  size_t n_sin22thi = 50;
  size_t n_t = n_mi * n_sin22thi;

  std::vector<SBNosc> osc_v(n_t,osctrue);
  std::vector<NeutrinoModel> testModel_v(n_t,nullModel);
  std::vector<float> mnu_v(n_t,0.0);
  std::vector<float> sin22th_v(n_t,0.0);
  std::vector<float> ue_v(n_t,0.0);
  std::vector<float> um_v(n_t,0.0);
  
  std::cout << "initialize" << std::endl;
  for(int mi = 0; mi < n_mi; mi++) {
    twatch.Start();
    for(int sin22thi = 0; sin22thi < n_sin22thi ; sin22thi++) {
      float mnu     = pow(10.,(mi/float(50)*TMath::Log10(10./.1) + TMath::Log10(.1)));
      float sin22th = pow(10.,(sin22thi/float(50)*TMath::Log10(1./1e-5) + TMath::Log10(1e-5)));
      float ue = pow(sin22th/float(4),.5);
      float um = 1.f;
      size_t idx = mi*n_mi+sin22thi;

      mnu_v[idx]     = mnu;
      sin22th_v[idx] = sin22th;
      ue_v[idx]      = ue;
      um_v[idx]      = um;

      auto& testModel = testModel_v[idx];
      testModel = NeutrinoModel(mnu,ue,um);
      
      auto& osc = osc_v[idx];
      osc.SetBothMode();
      osc.LoadModel(testModel);
      osc.OscillateThis(tag);
      
      if (mi==25 and sin22thi==25) {
	std::stringstream ss;
	ss.str("");
	ss << "DLCombinedFitDOF_25_25" << suffix;
	SBNspec spec = static_cast<SBNspec>(osc);    
	spec.WriteOut(ss.str());
      }

      // std::stringstream ss;
      // ss.str("");
      // ss << "DLCombinedFitDOF_m2_sin22th_" << mnu << "_" << sin22th << suffix;
      // SBNspec spec = static_cast<SBNspec>(osc);    
      // spec.WriteOut(ss.str());


    }
    twatch.Stop();
    std::cout << "@ mi=" << mi << "/" << n_mi << " : " << twatch.RealTime() << std::endl;
    twatch.Reset();
  }
  
  std::cout << "done" << std::endl;

  // Create our background
  SBNspec bkg("DL.SBNspec.root",xml,false);
  bkg.Scale("signal",0.0);
  bkg.Scale("fullosc",0.0);
  bkg.Scale("oscfull",0.0);
  bkg.WriteOut(std::string("DLCombinedFitDOF_bkg_spec") + suffix);

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
  data.WriteOut(std::string("DLCombinedFitDOF_fake_signal_spec") + suffix);
  data.Scale("intrinsic",0.0);
  data.Scale("cocktail",0.0);
  data.Scale("extbnb",0.0);
  data.Scale("fullosc",0.0);
  data.Scale("signal",0.0);
  data.WriteOut(std::string("DLCombinedFitDOF_data_spec") + suffix);
  twatch.Stop();
  std::cout << "done @ t=" << twatch.RealTime() << std::endl;
  twatch.Reset();
  
  // Stats + sys
  TFile * fsys = new TFile("DL.SBNcovar.root","read");
  TMatrixD *cov = (TMatrixD*)fsys->Get("frac_covariance_DL");

  // Chi2 using background only covariance
  SBNchi BF_chi(bkg,*cov);
  BF_chi.core_spectrum.Add(&data);
  BF_chi.core_spectrum.CalcFullVector();
  BF_chi.core_spectrum.CollapseVector();

  auto bkg_and_data = BF_chi.core_spectrum;
  bkg_and_data.CalcFullVector();
  bkg_and_data.CollapseVector();
  bkg_and_data.WriteOut(std::string("DLCombinedFitDOF_bkg_and_data_spec") + suffix);

  SBNchi true_chi(bkg_and_data,*cov);

  std::string out_file_name = "DLCombinedFitDOF_ana";
  out_file_name += suffix;
  out_file_name += std::string(".root");
  TFile* ana_file = TFile::Open(out_file_name.c_str(),"RECREATE");
  ana_file->cd();
  
  int fexp = 0;
  int _iter = -1;

  float _chi_sin22th = -1;
  float _chi_m2 = -1;
  float _chi  = -1;

  float _L_sin22th = -1;
  float _L_m2 = -1;
  float _L = -1;

  float _true_sin22th = -1;
  float _true_m2 = -1;

  float _chi_pt_bf = -1;
  float _chi_pt_pt = -1;

  float _L_pt_bf = -1;
  float _L_pt_pt = -1;

  TTree* outtree1 = new TTree("ana","");
  outtree1->Branch("fexp"   , &fexp     , "fexp/I");
  outtree1->Branch("iter"   , &_iter    , "iter/I");

  outtree1->Branch("chi_sin22th", &_chi_sin22th , "chi_sin22th/F");
  outtree1->Branch("chi_m2"     , &_chi_m2      , "chi_m2/F");
  outtree1->Branch("chi"        , &_chi         , "chi_chi2/F");
  outtree1->Branch("L_sin22th"  , &_L_sin22th   , "L_sin22th/F");
  outtree1->Branch("L_m2"       , &_L_m2        , "L_m2/F");
  outtree1->Branch("L"          , &_L           , "L/F");


  TTree* outtree2 = new TTree("ana_slim","");
  outtree2->Branch("fexp"       , &fexp         , "fexp/I");
  outtree2->Branch("chi_sin22th", &_chi_sin22th , "chi_sin22th/F");
  outtree2->Branch("chi_m2"     , &_chi_m2      , "chi_m2/F");
  outtree2->Branch("chi"        , &_chi         , "chi/F");
  outtree2->Branch("L_sin22th"  , &_L_sin22th   , "L_sin22th/F");
  outtree2->Branch("L_m2"       , &_L_m2        , "L_m2/F");
  outtree2->Branch("L"          , &_L           , "L/F");

  outtree2->Branch("true_sin22th", &_true_sin22th , "true_sin22th/F");
  outtree2->Branch("true_m2"     , &_true_m2      , "true_m2/F");

  outtree2->Branch("chi_pt_bf"  , &_chi_pt_bf, "chi_pt_bf/F");
  outtree2->Branch("chi_pt_pt"  , &_chi_pt_pt, "chi_pt_pt/F");

  outtree2->Branch("L_pt_bf"  , &_L_pt_bf, "L_pt_bf/F");
  outtree2->Branch("L_pt_pt"  , &_L_pt_pt, "L_pt_pt/F");


  // Fitter data holder
  CombinedFit cf;

  int n_fexp = 1000;
  int n_it = 6;

  std::cout << "Getting spec" << std::endl;
  twatch.Start();
  auto spec_fexp_v = BF_chi.SpecReturnSCVI(bkg_and_data, n_fexp, xml);
  // auto spec_fexp_v = BF_chi.SpecReturnSPVI(bkg_and_data, n_fexp, xml);
  twatch.Stop();
  std::cout << "done @ t=" << twatch.RealTime() << std::endl;
  twatch.Reset();

  for(; fexp < n_fexp; ++fexp) {

    auto BF_chi_chi_exp = BF_chi;
    auto BF_chi_L_exp   = BF_chi;

    auto& spec_fexp = spec_fexp_v[fexp];
    spec_fexp.CalcFullVector();
    spec_fexp.CollapseVector();

    // chi2_throw
    _chi_pt_pt = true_chi.CalcChi(spec_fexp);
    _L_pt_pt   = true_chi.CalcChiLog(spec_fexp);

    // std::stringstream ss;
    // ss << "spec_output_" << fexp;
    // spec_fexp.WriteOut(ss.str());
    
    // std::exit(1);
    // draw fake data
    twatch.Start();
    
    BF_chi_chi_exp.core_spectrum = spec_fexp;
    BF_chi_chi_exp.core_spectrum.CalcFullVector();
    BF_chi_chi_exp.core_spectrum.CollapseVector();
    
    BF_chi_L_exp.core_spectrum = BF_chi_chi_exp.core_spectrum;
    
    // 6 iterations fit
    for(size_t it = 0; it < n_it; ++it) {

      cf.ScanBoth(osc_v,BF_chi_chi_exp,BF_chi_L_exp, n_mi, n_sin22thi);

      const auto& lowest_chi_spec = cf.LowChiSBNspec();
      const auto& lowest_L_spec   = cf.LowLSBNspec();
      
      BF_chi_chi_exp = SBNchi(lowest_chi_spec,*cov);
      BF_chi_L_exp   = SBNchi(lowest_L_spec,*cov);

      BF_chi_chi_exp.core_spectrum = spec_fexp;
      BF_chi_chi_exp.core_spectrum.CalcFullVector();
      BF_chi_chi_exp.core_spectrum.CollapseVector();

      BF_chi_L_exp.core_spectrum = BF_chi_chi_exp.core_spectrum;
    }
    
    // save the final iteration space
    for(size_t idx=0; idx < cf.RegionChi().size(); ++idx) {
      _iter = n_it - 1;

      _chi = cf.RegionChi()[idx];
      _chi_sin22th = sin22th_v[idx];
      _chi_m2 = mnu_v[idx];
      _chi_m2 *= _chi_m2;

      _L = cf.RegionL()[idx];
      _L_sin22th = sin22th_v[idx];
      _L_m2 = mnu_v[idx];
      _L_m2 *= _L_m2;

      outtree1->Fill();
    }

    auto low_chi_idx = cf.LowChiIndex();
    _chi = cf.RegionChi()[low_chi_idx];
    _chi_sin22th = sin22th_v[low_chi_idx];
    _chi_m2 = mnu_v[low_chi_idx];
    _chi_m2 *= _chi_m2;

    auto low_L_idx = cf.LowLIndex();
    _L = cf.RegionL()[low_L_idx];
    _L_sin22th = sin22th_v[low_L_idx];
    _L_m2 = mnu_v[low_L_idx];
    _L_m2 *= _L_m2;

    _true_sin22th = fake_sin22th;
    _true_m2 = fake_m2;

    _chi_pt_bf = (float) BF_chi_chi_exp.CalcChi(fake_signal);
    _L_pt_bf   = (float) BF_chi_L_exp.CalcChiLog(fake_signal);

    outtree2->Fill();

    // take the best fit point and compute the chi2 with the true distribution

    twatch.Stop();
    std::cout << "@ fexp=" << fexp << "/" << n_fexp << " : " << twatch.RealTime() << std::endl;
    twatch.Reset();
    
  }
  ana_file->cd();
  outtree1->Write();
  outtree2->Write();
  ana_file->Close();

  return 0;
}
