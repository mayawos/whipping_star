//#define EIGEN_USE_MKL_ALL
//#pragma GCC optimize("O3","unroll-loops","inline")

//#define EIGEN_USE_BLAS
//#define EIGEN_USE_MKL_VML
//#define EIGEN_USE_LAPACKE

#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <tuple>
#include <limits>
#include <filesystem> // Require C++17
#include <regex>
#include <chrono>
#include <ctime>


#include <diy/master.hpp>
#include <diy/reduce.hpp>
#include <diy/partners/merge.hpp>
#include <diy/decomposition.hpp>
#include <diy/assigner.hpp>
#include <diy/mpi.hpp>
#include <diy/serialization.hpp>
#include <diy/partners/broadcast.hpp>
#include <diy/reduce-operations.hpp>

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
 
#include "TMatrixT.h"
#include "TH1D.h"

#include "params.h"
#include "SBNconfig.h"
#include "SBNchi.h"
#include <Eigen/Dense>
#include <Eigen/SVD>
#include "tools.h"
#include "prob.h"
#include "ngrid.h"

//#undef basic_string_view
using namespace std;
namespace fs = std::filesystem;
using namespace std::chrono;


#include "opts.h"

typedef diy::DiscreteBounds Bounds;

struct FitResult {
   size_t n_iter;
   int best_grid_point;
   double last_chi_min, delta_chi;
};

struct Block
{
    static void*    create()            { return new Block; }
    static void     destroy(void* b)    { delete static_cast<Block*>(b); }

    std::vector<double> last_chi_min, delta_chi;
    std::vector<int> best_grid_point, n_iter, i_grid, i_univ;
};

typedef std::array<double, 3> GridPoint;

// TODO what if fixed param is exactly 0 --- pow(10, x) is dangerous then
class GridPoints {
   public:
      GridPoints() {}
      GridPoints(std::vector<std::vector<double>> const & m_vec_grid) {
         for (auto gp : m_vec_grid)
            _gridpoints.push_back({pow(10, gp[0]), pow(10, gp[1]), pow(10, gp[2])});
      }

      size_t NPoints()         { return _gridpoints.size(); }
      GridPoint Get(size_t index) { return _gridpoints[index]; }

   private:
      std::vector<GridPoint> _gridpoints;

};
 
std::ostream& operator << (std::ostream& os, GridPoint const & p) {
   os << "\t0: " << p[0]  << " 1: " << p[1]  << " 2: " << p[2] <<  "\n";
   return os;
}

inline Eigen::VectorXd collapseVectorEigen(Eigen::VectorXd  const & vin, sbn::SBNconfig const & conf){
   // All we want is a representation with the subchannels added together
   Eigen::VectorXd cvec(conf.num_bins_total_compressed);
   cvec.setZero();
   for (int d=0; d<conf.num_detectors;++d) {
      size_t offset_in(0), offset_out(0);
      for (int i=0; i<conf.num_channels; i++) {
          size_t nbins_chan = conf.num_bins[i];
          for (int j=0; j<conf.num_subchannels[i]; j++) {
             size_t first_in   = d*conf.num_bins_detector_block            + offset_in;
             size_t first_out  = d*conf.num_bins_detector_block_compressed + offset_out;
             cvec.segment(first_out, nbins_chan).noalias() += vin.segment(first_in, nbins_chan);
             offset_in +=nbins_chan;
          }
          offset_out += nbins_chan;
      }
   }
   return cvec;
}

class SignalGenerator {
   public:
      SignalGenerator(
            sbn::SBNconfig const & conf, std::vector<std::vector<double>> const & vec_grid, size_t dim2,
            std::vector<Eigen::VectorXd>  const & sinsq,
            std::vector<Eigen::VectorXd>  const & sin,
            Eigen::VectorXd const & core, int oscmode
            ) 
      {
         m_conf = conf;
         m_dim2 = dim2;
         m_gridpoints=GridPoints(vec_grid);
         m_sinsq = sinsq;
         m_sin = sin;
         m_core = core;
         m_oscmode = oscmode;
         retVec = Eigen::VectorXd(m_conf.num_bins_total);
      }

   Eigen::VectorXd predict(size_t i_grid, bool compressed) {
      auto const & gp = m_gridpoints.Get(i_grid);
      //sbn::NeutrinoModel this_model(gp[0]*gp[0], gp[1], gp[2], false);
      sbn::NeutrinoModel this_model(gp[0]*gp[0], 0.0, gp[2], false);   // FIXME deal with ==0
      int m_idx = massindex(i_grid);
      Oscillate(m_sinsq[m_idx], m_sin[m_idx], this_model);

      if (compressed) return collapseVectorEigen(m_core+retVec, m_conf);
      else return m_core+retVec;
   }
   
   int massindex(size_t igrd) {return int(floor( (igrd) / m_dim2 ));}
   size_t gridsize() {return m_gridpoints.NPoints();}

   void Oscillate(Eigen::VectorXd const & sf_sinsq, Eigen::VectorXd const & sf_sin, 
         sbn::NeutrinoModel & working_model) 
   {
       retVec.setZero(); // !!!
       int which_dm = 41; // FIXME This must not be hardcoded!!! 41 is the value to be used in the 1 sterile case

       double prob_mumu(0), prob_ee(0), prob_mue(0), prob_mue_sq(0), prob_muebar(0), prob_muebar_sq(0);

       // TODO Make this logic nice
       if (m_oscmode==0) {
          prob_mue       = working_model.oscAmp( 2,  1, which_dm, 1);
          prob_mue_sq    = working_model.oscAmp( 2,  1, which_dm, 2);
          prob_muebar    = working_model.oscAmp(-2, -1, which_dm, 1);
          prob_muebar_sq = working_model.oscAmp(-2, -1, which_dm, 2);
       }

       else if (m_oscmode==1) prob_mumu = working_model.oscAmp(2,2,which_dm,2);
       else {
          std::cerr << "oscillation mode has to be either 0 or 1: " << m_oscmode << "\n";
          exit(1);
       }



       double osc_amp(0), osc_amp_sq(0);
       int osc_pattern(0);
           // Iterate over channels
        size_t offset(0);
        for (int i=0; i<m_conf.num_channels; i++) {
            size_t nbins_chan = m_conf.num_bins[i];
            auto const & thisPattern = m_conf.subchannel_osc_patterns[i];
            for (int j=0; j<m_conf.num_subchannels[i]; j++){
                osc_pattern = thisPattern[j];
                switch (osc_pattern){
                    case 11:
                        osc_amp_sq = prob_ee;
                        osc_amp = 0;
                        break;
                    case -11:
                        osc_amp_sq = prob_ee;
                        osc_amp = 0;
                        break;
                    case 22:
                        osc_amp_sq = prob_mumu;
                        osc_amp = 0;
                        break;
                    case -22:
                        osc_amp_sq = prob_mumu;
                        osc_amp = 0;
                        break;
                    case 21:
                        osc_amp    = prob_mue;
                        osc_amp_sq = prob_mue_sq;
                        break;
                    case -21:
                        osc_amp    = prob_muebar;
                        osc_amp_sq = prob_muebar_sq;
                        break;
                    case 0:  // TODO ask Mark about actual logic
                        osc_amp = 0;
                        osc_amp_sq = 0;
                    default:
                        break;
                }

                // Iterate over detectors
                for (int d=0; d<m_conf.num_detectors;++d) {
                  size_t first  = d*m_conf.num_bins_detector_block + offset;
                  retVec.segment(first, nbins_chan).noalias() += osc_amp   *  sf_sin.segment(first, nbins_chan);
                  retVec.segment(first, nbins_chan).noalias() += osc_amp_sq*sf_sinsq.segment(first, nbins_chan);
                }
                offset +=nbins_chan;
            }
        }
   }

   private:
      sbn::SBNconfig m_conf;
      GridPoints m_gridpoints;
      size_t m_dim2;
      std::vector<Eigen::VectorXd> m_sinsq, m_sin;
      Eigen::VectorXd m_core, retVec;
      int m_oscmode;
};


void createDataSets(HighFive::File* file, size_t nPoints, size_t nUniverses) {
    file->createDataSet<double>("last_chi_min", HighFive::DataSpace( { nPoints*nUniverses,       1} ));
    file->createDataSet<double>("delta_chi",    HighFive::DataSpace( { nPoints*nUniverses,       1} ));
    file->createDataSet<int>("best_grid_point", HighFive::DataSpace( { nPoints*nUniverses,       1} ));
    file->createDataSet<int>("n_iter",          HighFive::DataSpace( { nPoints*nUniverses,       1} ));
    file->createDataSet<double>("n_events",          HighFive::DataSpace( { nPoints*nUniverses,       1} ));
    // Some bookkeeping why not
    file->createDataSet<int>("i_grid",          HighFive::DataSpace( {nPoints*nUniverses,        1} ));
    file->createDataSet<int>("i_univ",          HighFive::DataSpace( {nPoints*nUniverses,        1} ));
    file->createDataSet<double>("gridx",        HighFive::DataSpace( {nPoints,                   1} ));
    file->createDataSet<double>("gridy",        HighFive::DataSpace( {nPoints,                   1} ));
}

void writeGrid(HighFive::File* file, std::vector<std::vector<double> > const & coords, int mode) {
    std::vector<double> xcoord;
    std::vector<double> ycoord;

    for (size_t i=0; i< coords.size(); i++) {
       xcoord.push_back(coords[i][0]);
       if (mode==0) ycoord.push_back(coords[i][1]);
       else if (mode==1) ycoord.push_back(coords[i][2]);
       else {
          std::cerr << "Error, the mode must be either 0 or 1 a the moment: " << mode << "\n";
          exit(1);
       }
    }
    HighFive::DataSet d_gridx          = file->getDataSet("gridx");
    HighFive::DataSet d_gridy          = file->getDataSet("gridy");
    d_gridx.select(   {0, 0}, {xcoord.size(), 1}).write(xcoord);
    d_gridy.select(   {0, 0}, {ycoord.size(), 1}).write(ycoord);
}

TMatrixD readFracCovMat(std::string const & rootfile){
    TFile  fsys(rootfile.c_str(),"read");
    TMatrixD cov =  *(TMatrixD*)fsys.Get("frac_covariance");
    fsys.Close();

    for(int i =0; i<cov.GetNcols(); i++) {
        for(int j =0; j<cov.GetNrows(); j++) {
            if ( std::isnan( cov(i,j) ) )  cov(i,j) = 0.0;
        }
    }
    return cov;
}

std::vector<TH1D> readHistos(std::string const & rootfile, std::vector<string> const & fullnames) {
    std::vector<TH1D > hist;
    TFile f(rootfile.c_str(),"read");
    for (auto fn: fullnames) hist.push_back(*((TH1D*)f.Get(fn.c_str())));
    f.Close();
    return hist;
}

std::vector<double> flattenHistos(std::vector<TH1D> const & v_hist) {
   std::vector<double> ret; 
   for (auto h : v_hist) {
      for (int i=1; i<(h.GetSize()-1); ++i) ret.push_back(h.GetBinContent(i));
   }
   return ret;
}

std::vector<std::tuple<std::string, float> > getFilesAndDm(std::string const & inDir, std::string const & tag, std::string const & subthing, double const & m_min, double const & m_max ) {
    const std::regex re("[-+]?([0-9]*\\.[0-9]+|[0-9]+)");
    std::smatch match;
    std::string result;

    std::vector<std::tuple<std::string, float> > ret;

    for (const auto & entry : fs::directory_iterator(inDir)) {
      if (std::string(fs::path(entry).stem()).rfind(tag + subthing, 0) == 0) {
        std::string test     = std::string(fs::path(entry));
        std::string teststem = std::string(fs::path(entry).stem());
        
        std::size_t loc = teststem.find(subthing);
        std::string _test = teststem.substr(loc);

        if (std::regex_search(_test, match, re) && match.size() > 1) {
           float lgmsq = std::stof(match.str(0));
           if (lgmsq > m_max || lgmsq < m_min) {
              std::cerr << "\t NOT using file " << test << " with " << match.str(0) << " " << lgmsq << "\n";
              continue;
           }
           ret.push_back({test, lgmsq});
           std::cerr << "\t Using file " << test << " with " << match.str(0) << " " << lgmsq << "\n";
        }
      }
    }
    return ret;
}

std::tuple< std::vector<std::vector<double>>, std::vector<float>> mkHistoVecStd(std::string const & inDir, std::string const & tag, std::string const & objname, std::vector<string> const & fullnames, double const & m_min, double const & m_max ) {

   // The input files come unordered
   auto const & inputs = getFilesAndDm(inDir, tag, objname, m_min, m_max);
   std::vector<std::vector<double> > temp(inputs.size());

   if (inputs.size() ==0) {
      std::cerr << "Error, no valid input files, exiting. Maybe check --xmin --xmax and location -i\n";
      exit(1);
   }
   // Sort by mass, ascending
   std::vector<float> masses;
   for (auto in : inputs) masses.push_back(std::get<1>(in));
   std::sort(masses.begin(), masses.end());

   std::cerr << "Summary of mass dimension --- we have " << inputs.size() << " inputs corresponding to these mass-squared splittings: \n";
   for (auto m : masses) std::cerr << m << " ";
   std::cerr << "\n";
   for (auto in : inputs) {
      std::vector<float>::iterator it = std::find(masses.begin(), masses.end(), std::get<1>(in));
      int mass_index = std::distance(masses.begin(), it);
      temp[mass_index] = flattenHistos(readHistos(std::get<0>(in), fullnames));
   }
   
   return {temp, masses};
}


inline double calcChi(Eigen::VectorXd const & data, Eigen::VectorXd const & prediction, Eigen::MatrixXd const & C_inv ) {
   auto const & diff = data-prediction;
   return diff.transpose() * C_inv * diff;
}
inline double calcChi(Eigen::VectorXd const & diff, Eigen::MatrixXd const & C_inv ) {
   return diff.transpose() * C_inv * diff;
}

inline std::tuple<double, int> universeChi2(Eigen::VectorXd const & data, Eigen::MatrixXd const & C_inv,
   SignalGenerator signal)
{
   double chimin=std::numeric_limits<double>::infinity();
   Eigen::VectorXd diff(data.rows());
   int bestP(0);
   for (size_t i=0; i<signal.gridsize(); ++i) {
      diff = data - signal.predict(i, true);
      double chi = calcChi(diff, C_inv);
       if (chi<chimin) {
          chimin = chi;
          bestP=i;
       }
   }
   return {chimin, bestP};
}

inline std::tuple<double, int> universeChi2(Eigen::VectorXd const & data, Eigen::MatrixXd const & C_inv,
   SignalGenerator signal, std::vector<size_t> const & goodpoints)
{
   double chimin=std::numeric_limits<double>::infinity();
   int bestP(0);
   for (auto  i : goodpoints) {
       double chi = calcChi(data - signal.predict(i, true), C_inv);
       if (chi<chimin) {
          chimin = chi;
          bestP=i;
       }
   }
   return {chimin, bestP};
}

std::vector<size_t> initialScan(Eigen::VectorXd const & data, Eigen::MatrixXd const & C_inv,
   SignalGenerator signal, double maxchi2)
{
   std::vector<size_t> goodpoints;
   for (size_t i=0; i<signal.gridsize(); ++i) {
       double chi = calcChi(data - signal.predict(i, true), C_inv);
       if (chi<maxchi2) goodpoints.push_back(i);
   }
   return goodpoints;
}

TMatrixT<double> calcCovarianceMatrix(TMatrixT<double> const & M, std::vector<double> const & spec){
    TMatrixT<double> Mout( M.GetNcols(), M.GetNcols() );
    Mout.Zero();
    // systematics per scaled event
    for(int i =0; i<M.GetNcols(); i++) {
        for(int j =0; j<M.GetNrows(); j++) {
            Mout(i,j) = M(i,j)*spec[i]*spec[j];
            if (i==j) Mout(i,i) += spec[i];
        }
    }
    return Mout;
}

// Can we optimise this?
inline Eigen::MatrixXd calcCovarianceMatrix(Eigen::MatrixXd const & M, Eigen::VectorXd const & spec){
   Eigen::MatrixXd test(M.cols(), M.cols());
    for(int i =0; i<M.cols(); i++) {
        for(int j =i; j<M.rows(); j++) {
            test(i,j) = M(i,j)*spec[i]*spec[j];
            if (i==j) test(i,i) += spec[i];
            test(j,i)=test(i,j);
        }
    }
    return test;
}
inline Eigen::MatrixXd calcCovarianceMatrixFast(Eigen::MatrixXd const & M, Eigen::VectorXd const & spec) {
   Eigen::MatrixXd ret(M.cols(), M.cols());
   ret.array()    = M.array()*(spec*spec.transpose()).array();
   ret.diagonal() += spec;
   return ret;
}
Eigen::MatrixXd calcMatrix(Eigen::MatrixXd const & M, Eigen::VectorXd const & spec){
   Eigen::MatrixXd ret(M.cols(), M.cols());
   ret.array()    = M.array()*(spec*spec.transpose()).array();
   return ret;
   }

Eigen::MatrixXd cholD(Eigen::MatrixXd const & M, Eigen::VectorXd const & spec, double tol=1e-7) {
   auto in = calcMatrix(M, spec);
   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(in);
   auto const & EV = eigensolver.eigenvalues();

   for (int i=0;i<EV.rows();++i) {
       if (EV[i]<=0) {
         if (fabs(EV[i]) < tol) for (int a=0; a<in.cols(); ++a) in(a,a) += EV[i];
       }
       if (fabs(EV[i])< tol) for (int a =0; a<in.cols(); a++) in(a,a) += tol;
   }
   Eigen::LLT<Eigen::MatrixXd> llt(in);
   return llt.matrixL();
}


Eigen::VectorXd sample(Eigen::VectorXd const & spec, Eigen::MatrixXd const & LMAT, std::mt19937 & rng) {

   std::normal_distribution<double> dist_normal(0,1);
   Eigen::VectorXd RDM(spec.rows());
   for (int i=0;i<spec.rows();++i) RDM[i] = dist_normal(rng);

   return LMAT*RDM + spec;
}

Eigen::VectorXd poisson_fluctuate(Eigen::VectorXd const & spec, std::mt19937 & rng) {
   Eigen::VectorXd RDM(spec.rows());
   for (int i=0;i<spec.rows();++i) {
      std::poisson_distribution<int> dist_pois(spec[i]);
      RDM[i] = double(dist_pois(rng));
   }
   return RDM;
}

// Cholesky decomposition and solve for inverted matrix --- fastest
inline Eigen::MatrixXd invertMatrixEigen3(Eigen::MatrixXd const & M){
    return M.llt().solve(Eigen::MatrixXd::Identity(M.rows(), M.rows()));
}

inline Eigen::MatrixXd collapseSubchannels(Eigen::MatrixXd const & EE, sbn::SBNconfig const & conf){
    Eigen::MatrixXd  retMat = Eigen::MatrixXd::Zero(conf.num_bins_detector_block_compressed, conf.num_bins_detector_block_compressed);

    int mrow(0), mcol(0), mrow_out(0), mcol_out(0);
    for(int ic = 0; ic < conf.num_channels; ic++) {
        for(int jc =0; jc < conf.num_channels; jc++) {
            for(int m=0; m < conf.num_subchannels[ic]; m++) {
                for(int n=0; n< conf.num_subchannels[jc]; n++) {
                   int a, c;
                   a=mrow + n*conf.num_bins[jc];
                   c=mcol + m*conf.num_bins[ic];
                   retMat.block(mrow_out, mcol_out, conf.num_bins[jc], conf.num_bins[ic]).noalias() += EE.block(a, c, conf.num_bins[jc], conf.num_bins[ic]);
                }
            }
            mrow     += conf.num_subchannels[jc]*conf.num_bins[jc];
            mrow_out += conf.num_bins[jc];
        } // end of column loop
        mrow      = 0; // as we end this row, reSet row count, but jump down 1 column
        mrow_out  = 0;
        mcol     += conf.num_subchannels[ic]*conf.num_bins[ic];
        mcol_out += conf.num_bins[ic];
    } // end of row loop
    return retMat;
}

inline Eigen::MatrixXd collapseDetectors(Eigen::MatrixXd const & M, sbn::SBNconfig const & conf){
    Eigen::MatrixXd  retMat = Eigen::MatrixXd::Zero(conf.num_bins_mode_block_compressed, conf.num_bins_mode_block_compressed);
    auto const & nrow = conf.num_bins_detector_block;
    auto const & crow = conf.num_bins_detector_block_compressed;
    for (int m=0; m<conf.num_detectors; m++) {
        for (int n=0; n<conf.num_detectors; n++) {
            retMat.block(n*crow, m*crow, crow, crow).noalias() = collapseSubchannels(M.block(n*nrow, m*nrow, nrow, nrow), conf);
        }
    }
    return retMat;
}


inline Eigen::MatrixXd updateInvCov(Eigen::MatrixXd const & covmat, Eigen::VectorXd const & spec_full, sbn::SBNconfig const & conf) {
    auto const & cov = calcCovarianceMatrixFast(covmat, spec_full);
    auto const & out = collapseDetectors(cov, conf);
    return invertMatrixEigen3(out);
}

inline FitResult coreFC(Eigen::VectorXd const & fake_data, Eigen::VectorXd const & v_coll,
      SignalGenerator signal,
      Eigen::MatrixXd const & INVCOV,
      Eigen::MatrixXd const & covmat,
      sbn::SBNconfig const & myconf,
      double chi_min_convergance_tolerance = 0.001,
      size_t max_number_iterations = 5
      )
{
   float last_chi_min = FLT_MAX;
   int best_grid_point = -99;
   size_t n_iter = 0;

   Eigen::MatrixXd invcov = INVCOV;//std::vector<double> temp;
   
   for(n_iter = 0; n_iter < max_number_iterations; n_iter++){
       if(n_iter!=0){
           //Calculate current full covariance matrix, collapse it, then Invert.
           auto const & temp  = signal.predict(best_grid_point, false);
           invcov = updateInvCov(covmat, temp, myconf);
       }
       //Step 2.0 Find the global_minimum_for this universe. Integrate in SBNfit minimizer here, a grid scan for now.
       float chi_min = FLT_MAX;
       auto const & resuni  = universeChi2(fake_data, invcov, signal);//, goodpoints);
       chi_min = std::get<0>(resuni);
       best_grid_point = std::get<1>(resuni);
       if(n_iter!=0){
           //Step 3.0 Check to see if min_chi for this particular fake_data  has converged sufficiently
           if(fabs(chi_min-last_chi_min)< chi_min_convergance_tolerance){
               last_chi_min = chi_min;
               break;
           }
       }
       last_chi_min = chi_min;
   } // End loop over iterations

   //Now use the curent_iteration_covariance matrix to also calc this_chi here for the delta.
   float this_chi = calcChi(fake_data, v_coll, invcov);

   FitResult fr = {n_iter, best_grid_point, last_chi_min, this_chi-last_chi_min}; 
   return fr;
}

void doScan(Block* b, diy::Master::ProxyWithLink const& cp, int rank,
      sbn::SBNconfig const & myconf,
      Eigen::MatrixXd const & ECOV, Eigen::MatrixXd const & INVCOVBG,
      Eigen::VectorXd const & ecore, SignalGenerator signal,
      HighFive::File* file, std::vector<size_t> const & rankwork, 
      double tol, size_t iter, bool debug)
{

    double starttime, endtime;
    std::vector<FitResult> results;
    std::vector<int> v_grid, v_univ, v_iter, v_best;
    std::vector<double> v_last, v_dchi;
    
    results.reserve(rankwork.size());
    v_grid.reserve(rankwork.size());
    v_univ.reserve(rankwork.size());
    v_iter.reserve(rankwork.size());
    v_best.reserve(rankwork.size());
    v_last.reserve(rankwork.size());
    v_dchi.reserve(rankwork.size());
    
    system_clock::time_point t_init = system_clock::now();
    for (int i_grid : rankwork) {

       if (debug && i_grid!=0) return;
       auto const & specfull_e = signal.predict(i_grid, false);
       auto const & speccoll   = collapseVectorEigen(specfull_e, myconf);
       auto const & corecoll   = collapseVectorEigen(ecore, myconf);

       starttime = MPI_Wtime();

       results.push_back(coreFC(corecoll, speccoll,
                signal, INVCOVBG, ECOV, myconf, tol, iter));

       v_univ.push_back(0);
       v_grid.push_back(i_grid);
       
       endtime   = MPI_Wtime();
       system_clock::time_point now = system_clock::now();

       auto t_elapsed = now - t_init;
       auto t_togo = t_elapsed * (int(rankwork.size()) - i_grid)/(i_grid+1);
       auto t_eta = now + t_togo;
       std::time_t t = system_clock::to_time_t(t_eta);

       if (rank==0 && i_grid%100==0) fmt::print(stderr, "[{}] gridp {}/{} took {} seconds. ETA: {}",cp.gid(), i_grid, rankwork.size(), endtime-starttime, std::ctime(&t));
    }

    // Write to HDF5
    starttime   = MPI_Wtime();
    HighFive::DataSet d_last_chi_min    = file->getDataSet("last_chi_min"   );
    HighFive::DataSet d_delta_chi       = file->getDataSet("delta_chi"      );
    HighFive::DataSet d_best_grid_point = file->getDataSet("best_grid_point");
    HighFive::DataSet d_n_iter          = file->getDataSet("n_iter"         );
    // write out this grid and universe
    HighFive::DataSet d_i_grid          = file->getDataSet("i_grid");
    HighFive::DataSet d_i_univ          = file->getDataSet("i_univ");
    // This is for the fake data dump

    size_t d_bgn = rankwork[0];
    for (auto res : results) {
       v_iter.push_back(res.n_iter);
       v_best.push_back(res.best_grid_point);
       v_last.push_back(res.last_chi_min);
       v_dchi.push_back(res.delta_chi);
    }

    d_last_chi_min.select(     {d_bgn, 0}, {size_t(v_last.size()), 1}).write(v_last);
    d_delta_chi.select(        {d_bgn, 0}, {size_t(v_dchi.size()), 1}).write(v_dchi);
    d_best_grid_point.select(  {d_bgn, 0}, {size_t(v_best.size()), 1}).write(v_best);
    d_n_iter.select(           {d_bgn, 0}, {size_t(v_iter.size()), 1}).write(v_iter);
    d_i_grid.select(           {d_bgn, 0}, {size_t(v_grid.size()), 1}).write(v_grid);
    d_i_univ.select(           {d_bgn, 0}, {size_t(v_univ.size()), 1}).write(v_univ);
    endtime   = MPI_Wtime();
    if (cp.gid()==0) fmt::print(stderr, "[{}] Write out took {} seconds\n", cp.gid(), endtime-starttime);
}

// TODO add size_t writeEvery to prevent memory overload
void doFC(Block* b, diy::Master::ProxyWithLink const& cp, int rank,
      const char * xmldata, sbn::SBNconfig const & myconf,
      TMatrixD const & covmat, Eigen::MatrixXd const & ECOV, Eigen::MatrixXd const & INVCOVBG,
      SignalGenerator signal,
      HighFive::File* file, std::vector<int> const & rankwork, int nUniverses, 
      double tol, size_t iter, bool debug, bool noWrite=false, int msg_every=100)
{

    double starttime, endtime;
    std::vector<FitResult> results;
    std::vector<int> v_grid, v_univ, v_iter, v_best;
    std::vector<double> v_last, v_dchi, v_nevents;
    
    if (!noWrite) {
       results.reserve(rankwork.size()*nUniverses);
       v_grid.reserve(rankwork.size()*nUniverses);
       v_univ.reserve(rankwork.size()*nUniverses);
       v_iter.reserve(rankwork.size()*nUniverses);
       v_best.reserve(rankwork.size()*nUniverses);
       v_last.reserve(rankwork.size()*nUniverses);
       v_dchi.reserve(rankwork.size()*nUniverses);
       v_nevents.reserve(rankwork.size()*nUniverses);
    }
    
    //Eigen::Map<const Eigen::MatrixXd > ECOV(covmat.GetMatrixArray(), covmat.GetNrows(), covmat.GetNrows());

    system_clock::time_point t_init = system_clock::now();
    for (int i_grid : rankwork) {

       if (debug && i_grid!=0) return;
       auto const & specfull_e = signal.predict(i_grid, false);
       auto const & speccoll = collapseVectorEigen(specfull_e, myconf);
       std::mt19937 rng(cp.gid()); // Mersenne twister
       Eigen::MatrixXd const & LMAT = cholD(ECOV, specfull_e);
       

       starttime = MPI_Wtime();
       for (int uu=0; uu<nUniverses;++uu) {
          auto const & fake_data = poisson_fluctuate(sample(specfull_e, LMAT, rng), rng);//
          auto const & fake_dataC = collapseVectorEigen(fake_data, myconf); 
          results.push_back(coreFC(fake_dataC, speccoll,
                   signal, INVCOVBG, ECOV, myconf, tol, iter));
          v_univ.push_back(uu);
          v_grid.push_back(i_grid);
          v_nevents.push_back(fake_data.sum());
       }
       endtime   = MPI_Wtime();
       system_clock::time_point now = system_clock::now();

       auto t_elapsed = now - t_init;
       auto t_togo = t_elapsed * (int(rankwork.size()) - i_grid)/(i_grid+1);
       auto t_eta = now + t_togo;
       std::time_t t = system_clock::to_time_t(t_eta);

       if (rank==0 && i_grid%msg_every==0) fmt::print(stderr, "[{}] gridp {}/{} ({} universes) took {} seconds. ETA: {}",cp.gid(), i_grid, rankwork.size(), nUniverses, endtime-starttime, std::ctime(&t));
    }

    if (!noWrite) {

       // Write to HDF5
       starttime   = MPI_Wtime();
       HighFive::DataSet d_last_chi_min    = file->getDataSet("last_chi_min"   );
       HighFive::DataSet d_delta_chi       = file->getDataSet("delta_chi"      );
       HighFive::DataSet d_best_grid_point = file->getDataSet("best_grid_point");
       HighFive::DataSet d_n_iter          = file->getDataSet("n_iter"         );
       HighFive::DataSet d_n_events        = file->getDataSet("n_events"       );
       // write out this grid and universe
       HighFive::DataSet d_i_grid          = file->getDataSet("i_grid");
       HighFive::DataSet d_i_univ          = file->getDataSet("i_univ");
       // This is for the fake data dump

       size_t d_bgn = rankwork[0]*nUniverses;
       for (auto res : results) {
          v_iter.push_back(res.n_iter);
          v_best.push_back(res.best_grid_point);
          v_last.push_back(res.last_chi_min);
          v_dchi.push_back(res.delta_chi);
       }

       d_last_chi_min.select(     {d_bgn, 0}, {size_t(v_last.size()), 1}).write(v_last);
       d_delta_chi.select(        {d_bgn, 0}, {size_t(v_dchi.size()), 1}).write(v_dchi);
       d_best_grid_point.select(  {d_bgn, 0}, {size_t(v_best.size()), 1}).write(v_best);
       d_n_iter.select(           {d_bgn, 0}, {size_t(v_iter.size()), 1}).write(v_iter);
       d_n_events.select(         {d_bgn, 0}, {size_t(v_nevents.size()), 1}).write(v_nevents);
       d_i_grid.select(           {d_bgn, 0}, {size_t(v_grid.size()), 1}).write(v_grid);
       d_i_univ.select(           {d_bgn, 0}, {size_t(v_univ.size()), 1}).write(v_univ);
       endtime   = MPI_Wtime();
       if (cp.gid()==0) fmt::print(stderr, "[{}] Write out took {} seconds\n", cp.gid(), endtime-starttime);
    }
}


void doFCLoadBalanced(Block* b, diy::Master::ProxyWithLink const& cp, int rank,
      const char * xmldata, sbn::SBNconfig const & myconf,
      TMatrixD const & covmat, Eigen::MatrixXd const & ECOV, Eigen::MatrixXd const & INVCOVBG,
      SignalGenerator signal,
      HighFive::File* file, std::vector<size_t> const & rankwork, int nUniverses, 
      double tol, size_t iter, bool debug, bool noWrite=false, int msg_every=100)
{

    double starttime, endtime;
    std::vector<FitResult> results;
    std::vector<int> v_grid, v_univ, v_iter, v_best;
    std::vector<double> v_last, v_dchi, v_nevents;

    size_t pStart = rankwork[0];
    size_t uStart = rankwork[1];
    size_t pLast  = rankwork[2];
    size_t uLast  = rankwork[3];

    size_t i_begin = pStart * nUniverses + uStart;
    size_t i_end   = pLast  * nUniverses + uLast;

    //fmt::print(stderr, "[{}] a,b,c,d: {} {} {} {} start at {} end at {}  lends {}\n", rank, pStart, uStart, pLast, uLast, i_begin, i_end, i_end-i_begin);
    size_t lenDS = i_end - i_begin;

    if (!noWrite) {
       results.reserve(lenDS);
       v_grid.reserve(lenDS);
       v_univ.reserve(lenDS);
       v_iter.reserve(lenDS);
       v_best.reserve(lenDS);
       v_last.reserve(lenDS);
       v_dchi.reserve(lenDS);
       v_nevents.reserve(lenDS);
    }

    std::vector<std::vector<size_t> > RW;
    if (pStart == pLast) RW.push_back({pStart, uStart, uLast});
    else {
      RW.push_back({pStart, uStart, nUniverses});
       for (size_t _p = pStart+1; _p<pLast;++_p) {
          RW.push_back({_p, 0, nUniverses});
       }
       if (uLast>0) RW.push_back({pLast, 0, uLast});

    }

    //for (auto r: RW) {
       //fmt::print(stderr, "[{}]  {}: {} -- {} \n", rank, r[0], r[1], r[2]);
    //}

    
    ////Eigen::Map<const Eigen::MatrixXd > ECOV(covmat.GetMatrixArray(), covmat.GetNrows(), covmat.GetNrows());

    system_clock::time_point t_init = system_clock::now();
    std::mt19937 rng(cp.gid()); // Mersenne twister
    for (auto r : RW) {

       size_t i_grid = r[0];

       if (debug && i_grid!=0) return;
       auto const & specfull_e = signal.predict(i_grid, false);
       auto const & speccoll = collapseVectorEigen(specfull_e, myconf);
       Eigen::MatrixXd const & LMAT = cholD(ECOV, specfull_e);

       starttime = MPI_Wtime();
       for (size_t uu=r[1]; uu<r[2];++uu) {
          auto const & fake_data = poisson_fluctuate(sample(specfull_e, LMAT, rng), rng);//
          auto const & fake_dataC = collapseVectorEigen(fake_data, myconf); 
          results.push_back(coreFC(fake_dataC, speccoll,
                   signal, INVCOVBG, ECOV, myconf, tol, iter));
          v_univ.push_back(uu);
          v_grid.push_back(i_grid);
          v_nevents.push_back(fake_data.sum());
       }
       endtime   = MPI_Wtime();
       system_clock::time_point now = system_clock::now();

       auto t_elapsed = now - t_init;
       auto t_togo = t_elapsed * (lenDS - results.size())/(results.size()+1);
       auto t_eta = now + t_togo;
       std::time_t t = system_clock::to_time_t(t_eta);

       if (rank==0 && results.size()%msg_every==0) fmt::print(stderr, "[{}] gridp {}/{} took {} seconds. ETA: {}",cp.gid(), results.size(), lenDS, endtime-starttime, std::ctime(&t));
    }

    if (!noWrite) {

       // Write to HDF5
       starttime   = MPI_Wtime();
       HighFive::DataSet d_last_chi_min    = file->getDataSet("last_chi_min"   );
       HighFive::DataSet d_delta_chi       = file->getDataSet("delta_chi"      );
       HighFive::DataSet d_best_grid_point = file->getDataSet("best_grid_point");
       HighFive::DataSet d_n_iter          = file->getDataSet("n_iter"         );
       HighFive::DataSet d_n_events        = file->getDataSet("n_events"       );
       // write out this grid and universe
       HighFive::DataSet d_i_grid          = file->getDataSet("i_grid");
       HighFive::DataSet d_i_univ          = file->getDataSet("i_univ");
       // This is for the fake data dump

       size_t d_bgn = i_begin;//rankwork[0]*nUniverses;
       for (auto res : results) {
          v_iter.push_back(res.n_iter);
          v_best.push_back(res.best_grid_point);
          v_last.push_back(res.last_chi_min);
          v_dchi.push_back(res.delta_chi);
       }

       d_last_chi_min.select(     {d_bgn, 0}, {size_t(v_last.size()), 1}).write(v_last);
       d_delta_chi.select(        {d_bgn, 0}, {size_t(v_dchi.size()), 1}).write(v_dchi);
       d_best_grid_point.select(  {d_bgn, 0}, {size_t(v_best.size()), 1}).write(v_best);
       d_n_iter.select(           {d_bgn, 0}, {size_t(v_iter.size()), 1}).write(v_iter);
       d_n_events.select(         {d_bgn, 0}, {size_t(v_nevents.size()), 1}).write(v_nevents);
       d_i_grid.select(           {d_bgn, 0}, {size_t(v_grid.size()), 1}).write(v_grid);
       d_i_univ.select(           {d_bgn, 0}, {size_t(v_univ.size()), 1}).write(v_univ);
       endtime   = MPI_Wtime();
       if (cp.gid()==0) fmt::print(stderr, "[{}] Write out took {} seconds\n", cp.gid(), endtime-starttime);
    }
}


inline bool file_exists (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

// --- main program ---//
int main(int argc, char* argv[]) {
    diy::mpi::environment env(argc, argv);
    diy::mpi::communicator world;

    size_t nPoints=-1;
    int mode=0;
    int msg_every=100;
    size_t nUniverses=1;
    int NTEST(0);
    std::string out_file="test.hdf5";
    std::string f_BG="BOTHv2_BKG_ONLY.SBNspec.root";
    std::string f_CV="BOTHv2_CV.SBNspec.root";
    std::string f_COV="BOTHv2.SBNcovar.root";
    std::string tag="";
    std::string d_in="";
    std::string xml="";
    double xmin(-1.0);
    double xmax(1.1);
    double xwidth(0.1);
    double ymin(-2.3);
    double ymax(0.1);
    double ywidth(0.05);
    double tol(0.001);
    size_t iter(5);
    // get command line arguments
    using namespace opts;
    Options ops(argc, argv);
    ops >> Option('o', "output",     out_file,   "Output filename.");
    ops >> Option('u', "nuniverse",  nUniverses, "Number of universes");
    ops >> Option("ntest", NTEST , "Number of universes");
    ops >> Option('x', "xml",        xml,        "XML config.");
    ops >> Option("tol",             tol,        "Minimiser tolerance");
    ops >> Option("iter",            iter,       "Max number of iterations.");
    ops >> Option('t', "tag",        tag,        "Tag.");
    ops >> Option('i', "indir",      d_in,       "Input file directory.");
    ops >> Option("core",            f_CV,       "Central values filename.");
    ops >> Option('b', "background", f_BG,       "Backgrounds filename.");
    ops >> Option('c', "covmat",     f_COV,      "Covariance matrix filename.");
    ops >> Option("xmin",            xmin,       "xmin");
    ops >> Option("xmax",            xmax,       "xmax");
    ops >> Option("xwidth",          xwidth,     "xwidth");
    ops >> Option("ymin",            ymin,       "ymin");
    ops >> Option("ymax",            ymax,       "ymax");
    ops >> Option("ywidth",          ywidth,     "ywidth");
    ops >> Option("msg",             msg_every,  "Prin a progress message every m gridpoints on rank 0 to stderr.");
    ops >> Option("mode",            mode, "Mode 0 is default --- dimension 2 is electron, mode 1 is muon");
    bool debug       = ops >> Present('d', "debug", "Operate on single gridpoint only");
    bool statonly    = ops >> Present("stat", "Statistical errors only");
    bool nowrite    = ops >> Present("nowrite", "Don't write output --- for performance estimates only");
    bool simplescan    = ops >> Present("scan", "Simple scan, no FC");
    if (ops >> Present('h', "help", "Show help")) {
        std::cout << "Usage:  [OPTIONS]\n";
        std::cout << ops;
        return 1;
    }
    
    if( world.rank()==0 ) {
      fmt::print(stderr, "\n*** This is diy running SBN Feldman Cousins ***\n");
    }

    // Whole bunch of tests
    if ( world.rank() == 0 ) {
       if (int(world.size()) > nPoints) {
          std::cerr << "Impossible to run on more ranks than grid points, exiting.\n";
          exit(1);
       }
       std::vector<std::string> infiles = {f_BG, f_COV, f_CV, xml};
       for (auto f : infiles) {
          if (!file_exists(f)) { 
             std::cerr << "Specified input file " << f <<" does not exist, exiting\n";
             exit(1);
          }
       }
       if (tag=="") {
          std::cerr << "tag (-t, --tag) cannot be undefined, exiting\n";
          exit(1);
       }
       if (d_in=="") {
          std::cerr << "Input dir (-i, --indor) cannot be undefined, exiting\n";
          exit(1);
       }
    }
    
    double time0 = MPI_Wtime();

    // Read the xml file on rank 0
    std::string line, text;
    if ( world.rank() == 0 ) {
       std::ifstream in(xml);
       while(std::getline(in, line))  text += line + "\n";
    }
    // YUCK, is this really the most elegant way to broadcast a simple string???
    int textsize = text.size();
    MPI_Bcast(&textsize, 1, MPI_INT, 0, world);
    if ( world.rank() != 0 ) text.resize(textsize);
    MPI_Bcast(const_cast<char*>(text.data()), textsize, MPI_CHAR, 0, world);

    // Central configuration object
    const char* xmldata = text.c_str();
    sbn::SBNconfig myconf(xmldata, false);

    // Pre-oscillated spectra
    std::vector<double> sinsqvec, sinvec;
    std::vector<float> msqsplittings;
    std::vector<Eigen::VectorXd > sinsqvec_eig, sinvec_eig;
    int nFilesIn(0);
    if (world.rank()==0) {
       auto temp = mkHistoVecStd(d_in, tag, "_SINSQ_", myconf.fullnames, xmin, xmax);
       sinsqvec = asVector(std::get<0>(temp));
       msqsplittings = std::get<1>(temp);
       
       auto temp2 = mkHistoVecStd(d_in, tag, "_SIN_", myconf.fullnames, xmin, xmax);
       sinvec   = asVector(std::get<0>(temp2));
       if (sinsqvec.size() != sinvec.size()) {
          std::cerr << "Error, number of input files for _SINSQ_ (" << sinsqvec.size() << ") differs from _SIN_ (" << sinvec.size() << ") exiting.\n";
          exit(1);
       }
       nFilesIn = msqsplittings.size();

    }
    diy::mpi::broadcast(world, sinsqvec, 0);
    diy::mpi::broadcast(world, sinvec,   0);
    diy::mpi::broadcast(world, msqsplittings,   0);
    diy::mpi::broadcast(world, nFilesIn, 0);
    for (auto v : splitVector(sinsqvec, nFilesIn)) sinsqvec_eig.push_back(Eigen::Map<Eigen::VectorXd> (v.data(), v.size(), 1) );
    for (auto v : splitVector(sinvec  , nFilesIn))   sinvec_eig.push_back(Eigen::Map<Eigen::VectorXd> (v.data(), v.size(), 1) );

    // Core spectrum
    std::vector<double> core;
    if (world.rank()==0) {
       auto const & cvhist = readHistos(f_CV, myconf.fullnames);
       core = flattenHistos(cvhist);
    }
    diy::mpi::broadcast(world, core, 0);

    Eigen::Map<Eigen::VectorXd> ecore(core.data(), core.size(), 1);

    // Background
    std::vector<double> bgvec;
    if (world.rank() == 0) {
       std::vector<TH1D> bghist  = readHistos(f_BG, myconf.fullnames);
       bgvec = flattenHistos(bghist);
       bghist.clear();
       bghist.shrink_to_fit();
    }
    diy::mpi::broadcast(world, bgvec, 0);

    // Read the covariance matrix on rank 0 --- broadcast and subsequently buid from array
    TMatrixD covmat;

    size_t nBins(0);
    if (!statonly) {
       std::vector<double> v_covmat;

       if ( world.rank() == 0 ) {
          TMatrixD temp = readFracCovMat(f_COV);
          nBins = temp.GetNcols();
          const double *pData = temp.GetMatrixArray();
          v_covmat.assign(pData, pData + temp.GetNoElements());
       }
       // broadcast
       diy::mpi::broadcast(world, v_covmat, 0);
       diy::mpi::broadcast(world, nBins,    0);
       // Set data of TMatrix
       covmat.ResizeTo(nBins, nBins);
       covmat.SetMatrixArray(v_covmat.data());
       releaseVec(v_covmat);
    }
    else {
       nBins=bgvec.size();
       covmat.ResizeTo(nBins, nBins);
       covmat.Zero();
    }
    Eigen::Map<const Eigen::MatrixXd > ECOVMAT(covmat.GetMatrixArray(), covmat.GetNrows(), covmat.GetNrows());
    
    // Use the BG only inv cov matrix as start point
    TMatrixT<double> _cov = calcCovarianceMatrix(covmat, bgvec);
    Eigen::Map<const Eigen::MatrixXd> ecov(_cov.GetMatrixArray(), _cov.GetNcols(), _cov.GetNcols());
    auto const & _covcol = collapseDetectors(ecov, myconf);
    auto const & INVCOVBG = invertMatrixEigen3(_covcol);

    // Setup grid
    NGrid mygrid;
    size_t dim2(0);

    double mmin = msqsplittings.front()*0.5;
    double mmax = msqsplittings.back()*0.5;
    double mwidth = (msqsplittings[1] - msqsplittings[0])*0.5;

    if (world.rank()==0) std::cerr << "Mass setup for input grid: " << mmin << " < " << mmax << " width: " << mwidth << "\n";

    mygrid.AddDimension("m4", mmin, mmax, mwidth );
    //mygrid.AddDimension("m4",  xmin, xmax, xwidth); // Sterile mass ultimately limited by input files
    if (mode==0) {
       mygrid.AddDimension("ue4", ymin, ymax, ywidth);// arbirtrarily dense! mixing angle nu_e
       mygrid.AddFixedDimension("um4", 0.0);
       dim2 = mygrid.f_dimensions[1].f_N;
    }
    else if (mode==1) {
       mygrid.AddFixedDimension("ue4", 0.0);
       mygrid.AddDimension("um4", ymin, ymax, ywidth);
       dim2 = mygrid.f_dimensions[2].f_N;
    }
    else {
       std::cerr << "Error, the mode must be either 0 or 1 a the moment: " << mode << "\n";
       exit(1);
    }
    nPoints = mygrid.f_num_total_points;
    GridPoints GP(mygrid.GetGrid());

    if (world.rank()==0) mygrid.Print();

    // Finally, the signal generator
    SignalGenerator     signal(myconf, mygrid.GetGrid(), dim2, sinsqvec_eig, sinvec_eig, ecore, mode);
    
    double time1 = MPI_Wtime();
    if (world.rank()==0) fmt::print(stderr, "[{}] Input preparation took {} seconds\n",world.rank(), time1 - time0);

    if (NTEST>0) {
       auto const & sv = signal.predict(1, false);
       std::vector<double> svb(sv.data(), sv.data() + sv.rows() * sv.cols());
       
       double t0 = MPI_Wtime();
       for (int i=0;i<NTEST;++i) signal.predict(1, false);
       double t1 = MPI_Wtime();
       for (int i=0;i<NTEST;++i) collapseVectorEigen(sv, myconf);
       double t2 = MPI_Wtime();
       //for (int i=0;i<NTEST;++i) signalc.predict(1, false);
       double t3 = MPI_Wtime();
       for (int i=0;i<NTEST;++i) collapseVectorStd(svb, myconf);
       double t4 = MPI_Wtime();

       fmt::print(stderr, "\n {} calls to predict took {} seconds\n",             NTEST, t1-t0);
       fmt::print(stderr, "\n {} calls to collapseVectorEigen took {} seconds\n", NTEST, t2-t1);
       fmt::print(stderr, "\n {} calls to class predict took {} seconds\n", NTEST, t3-t2);
       fmt::print(stderr, "\n {} calls to collapseVectorStd took {} seconds\n",   NTEST, t4-t3);
       exit(1);
    }



    if( world.rank()==0 ) {
      fmt::print(stderr, "***********************************\n");
      fmt::print(stderr, "    Output will be written to {}\n", out_file);
      fmt::print(stderr, "    Points:    {}\n"               ,  nPoints);
      fmt::print(stderr, "    nBins:     {}\n"               , nBins);
      fmt::print(stderr, "    Universes: {}\n"               , nUniverses);
      fmt::print(stderr, "    Total size of dataset:  {}\n"  , nPoints*nUniverses);
      fmt::print(stderr, "    f_BG:   {}\n"                  , f_BG);
      fmt::print(stderr, "    f_CV:   {}\n"                  , f_CV);
      fmt::print(stderr, "    f_COV:  {}\n"                  , f_COV);
      fmt::print(stderr, "    iter :  {}\n"                  , iter);
      fmt::print(stderr, "    tol:    {}\n"                  , tol);
      if (statonly) fmt::print(stderr, "    S T A T  O N L Y \n");
      if (debug)    fmt::print(stderr,    "    D E B U G \n"    );
      if (nowrite)  fmt::print(stderr,  "  N O   W R I T E \n"  );
      fmt::print(stderr, "***********************************\n");
    }
    
    // Create hdf5 file structure here 
    HighFive::File* f_out  = new HighFive::File(out_file,
                        HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate,
                        HighFive::MPIOFileDriver(MPI_COMM_WORLD,MPI_INFO_NULL));

    // Create datasets needed TODO nbinsC --- can we get that from somewhere?
    createDataSets(f_out, nPoints, nUniverses);
   
    // First rank also writes the grid so we know what the poins actually are
    if (world.rank() == 0)  writeGrid(f_out, mygrid.GetGrid(), mode);

    //// Now more blocks as we have universes
    size_t blocks = world.size();//nPoints;//*nUniverses;
    if (world.rank()==0) fmt::print(stderr, "FC will be done on {} blocks, distributed over {} ranks\n", blocks, world.size());
    Bounds fc_domain(1);
    fc_domain.min[0] = 0;
    fc_domain.max[0] = blocks-1;
    diy::RoundRobinAssigner        fc_assigner(world.size(), blocks);
    diy::RegularDecomposer<Bounds> fc_decomposer(1, fc_domain, blocks);
    diy::RegularBroadcastPartners  fc_comm(    fc_decomposer, 2, true);
    diy::RegularMergePartners      fc_partners(fc_decomposer, 2, true);
    diy::Master                    fc_master(world, 1, -1, &Block::create, &Block::destroy);
    diy::decompose(1, world.rank(), fc_domain, fc_assigner, fc_master);

    //std::vector<int> work(nPoints);
    //std::iota(std::begin(work), std::end(work), 0); //0 is the starting number
    //std::vector<int> rankwork = splitVector(work, world.size())[world.rank()];
    //work.clear();
    //work.shrink_to_fit();


    // Somehow this feels overly complicatd
    size_t _S(nPoints*nUniverses);
    std::vector<size_t> _L;
    size_t maxwork = size_t(ceil(_S/world.size()));
    for (size_t r=0; r <                _S%world.size(); ++r) _L.push_back(maxwork + 1);
    for (size_t r=0; r < world.size() - _S%world.size(); ++r) _L.push_back(maxwork    );


    std::vector<size_t> _bp, _bu;
    _bp.push_back(0);
    _bu.push_back(0);

    size_t _n(0), _temp(0);
    for (size_t i=0; i<nPoints;++i) {
       for (size_t j=0; j<nUniverses;++j) {
          if (_temp == _L[_n]) {
             _bp.push_back(i);
             _bu.push_back(j);
             _temp = 0;
             _n+=1;
          }
          _temp+=1;
       }
    }
    _bp.push_back(nPoints-1);
    _bu.push_back(nUniverses);
    
    std::vector<size_t> rankworkbl;
    //fmt::print(stderr, "[{}] Got {} ranks and {} sets of work {}\n", world.rank(), world.size(), _bu.size() -1, _bp.size() -1);
    world.barrier();
    rankworkbl.push_back(_bp[world.rank()]);
    rankworkbl.push_back(_bu[world.rank()]);
    rankworkbl.push_back(_bp[world.rank()+1]);
    rankworkbl.push_back(_bu[world.rank()+1]);


    //fmt::print(stderr, "Got maxwork {}\n", maxwork);
    //for (int gp=0;gp<nPoints;++gp) {
       //for (int uv=0;uv<nUniverses;++uv) {
          //if (_n%maxwork ==0) {
             //_bp.push_back(gp);
             //_bu.push_back(uv);
          //}
          //_n++;
       //}
    //}
    //_bp.push_back(nPoints-1);
    //_bu.push_back(nUniverses-1);
    //std::vector<int> rankworkbl;

    //fmt::print(stderr, "[{}] Got {} ranks and {} sets of work {}\n", world.rank(), world.size(), _bu.size() -1, _bp.size() -1);
    //exit(0);



    double starttime = MPI_Wtime();
    if (simplescan) {
       if (world.rank()==0) fmt::print(stderr, "Start simple scan\n");
       fc_master.foreach([world, ECOVMAT, INVCOVBG, ecore, myconf,  f_out, rankworkbl, tol, iter, debug, signal](Block* b, const diy::Master::ProxyWithLink& cp)
                              { doScan(b, cp, world.rank(), myconf, ECOVMAT, INVCOVBG, ecore, signal, f_out, rankworkbl, tol, iter, debug); });
    }
    else {
       if (world.rank()==0) fmt::print(stderr, "Start FC\n");
       fc_master.foreach([world, covmat, ECOVMAT, xmldata, INVCOVBG, myconf, nUniverses, f_out, rankworkbl, tol, iter, debug, signal, nowrite, msg_every](Block* b, const diy::Master::ProxyWithLink& cp)
                              { doFCLoadBalanced(b, cp, world.rank(), xmldata, myconf, covmat, ECOVMAT, INVCOVBG, signal, f_out, rankworkbl, nUniverses, tol, iter, debug, nowrite, msg_every); });
    }
    double endtime   = MPI_Wtime(); 
    if (world.rank()==0 || world.rank()==world.size()-1) fmt::print(stderr, "[{}] that took {} seconds\n", world.rank(), endtime-starttime);
    if (world.rank()==0) fmt::print(stderr, "Output written to {}\n",out_file);

    delete f_out;
    return 0;
}
