#ifndef _DLUTIL_CXX_
#define _DLUTIL_CXX_

#include <limits>

#include "DLUtil.h"

namespace sbn {

  void CombinedFit::ScanChi2(const std::vector<SBNosc>& osc_v,
			     const SBNchi& BF_chi,
			     const size_t n_mi,
			     const size_t n_sin22thi) {
    
    size_t n_t = n_mi * n_sin22thi;

    if (_chi2_v.size() != n_t) {
      _chi2_v.clear();
      _chi2_v.resize(n_t,std::numeric_limits<double>::max());
    }
    else {
      for(auto& chi2 : _chi2_v) {
	chi2 = std::numeric_limits<double>::max();
      }
    }
  
    for(size_t mi = 0; mi < n_mi; mi++) {
      for(size_t sin22thi = 0; sin22thi < n_sin22thi; sin22thi++) {
	size_t idx = mi*n_mi+sin22thi;
	const auto& osc = osc_v[idx];
	_chi2_v[idx] = BF_chi.CalcChi(osc);
      }
    }
    
    _chi_iter = std::min_element(_chi2_v.begin(),_chi2_v.end());
    _chi_idx = std::distance(_chi2_v.begin(), _chi_iter);
    _lowest_spec = static_cast<SBNspec>(osc_v[_chi_idx]);
    
    return;
  }

}

#endif
