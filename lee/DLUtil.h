#ifndef _DLUTIL_H_
#define _DLUTIL_H_

#include "SBNosc.h"
#include "SBNchi.h"

namespace sbn {

  class CombinedFit {
  public:

    CombinedFit() {}
    ~CombinedFit() {}

    void ScanBoth(const std::vector<SBNosc>& osc_v,
		  const SBNchi& BF_chi_chi,
		  const SBNchi& BF_chi_L,
		  const size_t n_mi,
		  const size_t n_sin22thi);
    
    void ScanChi(const std::vector<SBNosc>& osc_v,
		 const SBNchi& BF_chi,
		 const size_t n_mi,
		 const size_t n_sin22thi);

    void ScanL(const std::vector<SBNosc>& osc_v,
	       const SBNchi& BF_chi,
	       const size_t n_mi,
	       const size_t n_sin22thi);
    
    
    const std::vector<double>& RegionChi() { return _chi_v; }
    const std::vector<double>& RegionL()    { return _L_v; }
    
    double LowChi() { return *_chi_iter; }
    double LowL() { return *_L_iter; }

    size_t LowChiIndex() { return _chi_idx; }
    size_t LowLIndex() { return _L_idx; }

    const SBNspec& LowChiSBNspec() { return _lowest_chi_spec; }
    const SBNspec& LowLSBNspec() { return _lowest_L_spec; }

  private:

    std::vector<double> _chi_v;
    std::vector<double> _L_v;
    std::vector<double>::iterator _chi_iter;
    std::vector<double>::iterator _L_iter;

    size_t _chi_idx;
    size_t _L_idx;

    SBNspec _lowest_chi_spec;
    SBNspec _lowest_L_spec;
    
  };
  

}

#endif
