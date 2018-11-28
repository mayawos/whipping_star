#ifndef _DLUTIL_H_
#define _DLUTIL_H_

#include "SBNosc.h"
#include "SBNchi.h"

namespace sbn {

  class CombinedFit {
  public:

    CombinedFit() {}
    ~CombinedFit() {}
    
    void ScanChi2(const std::vector<SBNosc>& osc_v,
		  const SBNchi& BF_chi,
		  const size_t n_mi,
		  const size_t n_sin22thi);
    
    
    const std::vector<double>& RegionChi2() { return _chi2_v; }
    double LowChi2() { return *_chi_iter; }
    size_t LowChi2Index() { return _chi_idx; }
    const SBNspec& LowSBNspec() { return _lowest_spec; }

  private:

    std::vector<double> _chi2_v;
    std::vector<double>::iterator _chi_iter;
    size_t _chi_idx;
    SBNspec _lowest_spec;
    
  };
  

}

#endif
