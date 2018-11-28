#include <getopt.h>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TStopwatch.h"

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

  const struct option longopts[] =
    {
      {"xml" , required_argument, 0, 'x'},
      {0     , no_argument      , 0,  0 },
    };
  
  while(iarg != -1)
    {
      iarg = getopt_long(argc,argv, "x:", longopts, &index);

      switch(iarg)
	{
	case 'x':
	  xml = optarg;
	  break;
	case '?':
	case 'h':
	  std::cout << "Allowed arguments:" << std::endl;
	  std::cout << "\t-x\t--xml\t\tInput .xml file for SBNconfig" << std::endl;
	  return 0;
	}
    }

  std::string tag = "DL";

  // No oscilaltions model
  NeutrinoModel nullModel(0,0,0);

  // Background spectrum with fullosc=0 oscfull=0 signal=0
  SBNosc osctrue("DLSens_Bkg.SBNspec.root", xml, nullModel);
  osctrue.Scale("fullosc",0.0);
  osctrue.Scale("oscfull",0.0);
  osctrue.Scale("signal",0.0);
  osctrue.has_been_scaled = false;

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
	ss << "DLCombinedFit_25_25";
	SBNspec spec = static_cast<SBNspec>(osc);    
	spec.WriteOut(ss.str());
      }
      
    }
    twatch.Stop();
    std::cout << "@ mi=" << mi << "/" << n_mi << " : " << twatch.RealTime() << std::endl;
    twatch.Reset();
  }
  
  std::cout << "done" << std::endl;

  // Create our background 
  SBNspec bkg("DL.SBNspec.root",xml);
  bkg.Scale("signal",0.0);
  bkg.Scale("fullosc",0.0);
  bkg.Scale("oscfull",0.0);
  bkg.WriteOut("DLCombinedFit_bkg_spec");

  // Create our lee data
  SBNspec data("DL.SBNspec.root",xml);
  data.Scale("intrinsic",0.0); // scale intrinsic nue to zero
  data.Scale("cocktail",0.0); // scale numu and nue cocktail to zero
  data.Scale("extbnb",0.0); // scale numu and nue extbnb to zero
  data.Scale("fullosc",0.0); // scale intrinsic fullosc to zero
  data.Scale("oscfull",0.0); // scale intrinsic oscfull to zero
  data.WriteOut("DLCombinedFit_data_spec");
  
  // Stats + sys
  TFile * fsys = new TFile("DL.SBNcovar.root","read");
  TMatrixD *cov = (TMatrixD*)fsys->Get("frac_covariance_DL");

  // Chi2 using background only covariance
  SBNchi BF_chi(bkg,*cov);

  // Add data on top
  BF_chi.core_spectrum.Add(&data);
  BF_chi.core_spectrum.CollapseVector();

  // Copy data spectrum for reset each iteration
  SBNspec bkg_and_data = BF_chi.core_spectrum;

  // Output tfile + ttree
  TFile* ana_file = TFile::Open("DLCombinedFit_ana.root","RECREATE");
  ana_file->cd();

  int _iter =  -1;
  int _chi_idx = -1;
  double _sin22th = -1;
  double _m2 = -1;
  double _chi2 = -1;

  TTree* outtree = new TTree("ana","");
  outtree->Branch("iter"   , &_iter    , "iter/I");
  outtree->Branch("sin22th", &_sin22th , "sin22th/D");
  outtree->Branch("m2"     , &_m2      , "m2/D");
  outtree->Branch("chi2"   , &_chi2    , "chi2/D");
  outtree->Branch("chi_idx", &_chi_idx , "chi_idx/I");

  // Fitter data holder
  CombinedFit cf;

  // 6 iterations fit
  for(size_t it = 0; it < 6; ++it) {
    std::cout << "start @it=" << it << std::endl;
    
    std::stringstream ss;
    ss << "DLCombinedFit_BF_chi_spec_" << it;
    BF_chi.core_spectrum.WriteOut(ss.str());

    cf.ScanChi2(osc_v, BF_chi, n_mi, n_sin22thi);

    const auto& lowest_spec = cf.LowSBNspec();

    BF_chi = SBNchi(lowest_spec,*cov);
    BF_chi.core_spectrum = bkg_and_data; 
    BF_chi.core_spectrum.CollapseVector();
    
    for(size_t idx=0; idx < cf.RegionChi2().size(); ++idx) {
      _iter = it;
      _chi2 = cf.RegionChi2()[idx];
      _chi_idx = idx;
      _sin22th = sin22th_v[idx];
      _m2 = mnu_v[idx];
      _m2 *= _m2;
      outtree->Fill();
    }
    
    ss.str("");
    ss << "DLCombinedFit_lowest_spec_" << it;
    lowest_spec.WriteOut(ss.str());
    
  }
  
  // Writeout
  ana_file->cd();
  outtree->Write();
  ana_file->Close();

  return 0;
}
