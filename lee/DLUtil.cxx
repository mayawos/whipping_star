#ifndef _DLUTIL_CXX_
#define _DLUTIL_CXX_

#include <limits>

#include "DLUtil.h"

namespace sbn {


  void CombinedFit::ScanBoth(const std::vector<SBNosc>& osc_v,
			     const SBNchi& BF_chi_chi,
			     const SBNchi& BF_chi_L,
			     const size_t n_mi,
			     const size_t n_sin22thi) {

    size_t n_t = n_mi * n_sin22thi;

    if (_chi_v.size() != n_t) {
      _chi_v.clear();
      _chi_v.resize(n_t,std::numeric_limits<double>::max());
    }
    else {
      for(auto& chi : _chi_v) {
	chi = std::numeric_limits<double>::max();
      }
    }

    if (_L_v.size() != n_t) {
      _L_v.clear();
      _L_v.resize(n_t,std::numeric_limits<double>::max());
    }
    else {
      for(auto& L : _L_v) {
	L = std::numeric_limits<double>::max();
      }
    }

    for(size_t mi = 0; mi < n_mi; mi++) {
      for(size_t sin22thi = 0; sin22thi < n_sin22thi; sin22thi++) {
	size_t idx = mi*n_mi+sin22thi;
	const auto& osc = osc_v[idx];
	_chi_v[idx] = BF_chi_chi.CalcChi(osc);
	_L_v[idx]   = BF_chi_L.CalcChiLog(osc);
      }
    }


    _chi_iter = std::min_element(_chi_v.begin(),_chi_v.end());
    _chi_idx = std::distance(_chi_v.begin(), _chi_iter);
    _lowest_chi_spec = static_cast<SBNspec>(osc_v[_chi_idx]);
    _lowest_chi_spec.CalcFullVector();
    _lowest_chi_spec.CollapseVector();


    _L_iter = std::min_element(_L_v.begin(),_L_v.end());
    _L_idx = std::distance(_L_v.begin(), _L_iter);
    _lowest_L_spec = static_cast<SBNspec>(osc_v[_L_idx]);
    _lowest_L_spec.CalcFullVector();
    _lowest_L_spec.CollapseVector();

  }

  void CombinedFit::ScanChi(const std::vector<SBNosc>& osc_v,
			    const SBNchi& BF_chi,
			    const size_t n_mi,
			    const size_t n_sin22thi) {
    
    size_t n_t = n_mi * n_sin22thi;
    
    if (_chi_v.size() != n_t) {
      _chi_v.clear();
      _chi_v.resize(n_t,std::numeric_limits<double>::max());
    }
    else {
      for(auto& chi : _chi_v) {
	chi = std::numeric_limits<double>::max();
      }
    }
    
    for(size_t mi = 0; mi < n_mi; mi++) {
      for(size_t sin22thi = 0; sin22thi < n_sin22thi; sin22thi++) {
	size_t idx = mi*n_mi+sin22thi;
	const auto& osc = osc_v[idx];
	_chi_v[idx] = BF_chi.CalcChi(osc);
      }
    }
    
    _chi_iter = std::min_element(_chi_v.begin(),_chi_v.end());
    _chi_idx = std::distance(_chi_v.begin(), _chi_iter);
    _lowest_chi_spec = static_cast<SBNspec>(osc_v[_chi_idx]);
    _lowest_chi_spec.CalcFullVector();
    _lowest_chi_spec.CollapseVector();
    
    return;
  }


  void CombinedFit::ScanL(const std::vector<SBNosc>& osc_v,
			  const SBNchi& BF_chi,
			  const size_t n_mi,
			  const size_t n_sin22thi) {
    
    size_t n_t = n_mi * n_sin22thi;
    
    if (_L_v.size() != n_t) {
      _L_v.clear();
      _L_v.resize(n_t,std::numeric_limits<double>::max());
    }
    else {
      for(auto& L : _L_v) {
	L = std::numeric_limits<double>::max();
      }
    }
    
    for(size_t mi = 0; mi < n_mi; mi++) {
      for(size_t sin22thi = 0; sin22thi < n_sin22thi; sin22thi++) {
	size_t idx = mi*n_mi+sin22thi;
	const auto& osc = osc_v[idx];
	_L_v[idx] = BF_chi.CalcChiLog(osc);
      }
    }
    
    _L_iter = std::min_element(_L_v.begin(),_L_v.end());
    _L_idx = std::distance(_L_v.begin(), _L_iter);
    _lowest_L_spec = static_cast<SBNspec>(osc_v[_L_idx]);
    _lowest_L_spec.CalcFullVector();
    _lowest_L_spec.CollapseVector();
    
    return;
  }

}

#endif
