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
#include <diy/io/block.hpp>

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Easy.hpp>

#include <cppoptlib/meta.h>
#include <cppoptlib/problem.h>
#include <cppoptlib/function.h>
#include <cppoptlib/solver/bfgssolver.h>
#include <cppoptlib/solver/lbfgssolver.h>
#include <cppoptlib/solver/lbfgsbsolver.h>
#include <cppoptlib/solver/newtondescentsolver.h>
#include <cppoptlib/solver/cmaessolver.h>
#include <cppoptlib/solver/neldermeadsolver.h>
#include <cppoptlib/solver/conjugatedgradientdescentsolver.h>
#include <cppoptlib/solver/gradientdescentsolver.h>

#include <Eigen/Dense>
#include <vector>
#include <set>
#include <iostream>

#include "TMatrixT.h"
#include "TH1D.h"

#include "params.h"
#include "SBNconfig.h"
#include "SBNchi.h"
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <unsupported/Eigen/CXX11/Tensor>
#include "tools.h"
#include "prob.h"
#include "ngrid.h"

#include <mfa/mfa.hpp>
#include <mfa/block_base.hpp>


//#undef basic_string_view
using namespace std;
namespace fs = std::filesystem;
using namespace std::chrono;

#include "opts.h"

// set input and ouptut precision here, float or double
#if 0
typedef float                          real_t;
#else
typedef double                         real_t;
#endif

using namespace cppoptlib;

class Rosenbrock : public Problem<double> {
  public:
    double value(const Eigen::VectorXd &x) {
        const double t1 = (1 - x[0]);
        const double t2 = (x[1] - x[0] * x[0]);
        return   t1 * t1 + 100 * t2 * t2;
    }
};

struct FitResult {
   size_t n_iter;
   int best_grid_point;
   double last_chi_min, delta_chi;
   std::vector<double> fakedataC, collspec;
};

// arguments to block foreach functions
struct DomainArgs : public ModelInfo
{
    DomainArgs(int dom_dim, int pt_dim) :
        ModelInfo(dom_dim, pt_dim)
    {
        starts.resize(dom_dim);
        full_dom_pts.resize(dom_dim);
        min.resize(dom_dim);
        max.resize(dom_dim);
        s.resize(pt_dim);
        f.resize(pt_dim);
    }
    vector<int>         starts;                     // starting offsets of ndom_pts (optional, usually assumed 0)
    vector<int>         full_dom_pts;               // number of points in full domain in case a subset is taken
    vector<real_t>      min;                        // minimum corner of domain
    vector<real_t>      max;                        // maximum corner of domain
    vector<real_t>      s;                          // scaling factor for each variable or any other usage
    real_t              r;                          // x-y rotation of domain or any other usage
    vector<real_t>      f;                          // frequency multiplier for each variable or any other usage
    real_t              t;                          // waviness of domain edges or any other usage
    real_t              n;                          // noise factor [0.0 - 1.0]
    string              infile;                     // input filename
    bool                multiblock;                 // multiblock domain, get bounds from block
};


// block
template <typename T>
struct Block : public BlockBase<T>
{
    static
        void* create()              { return mfa::create<Block>(); }

    static
        void destroy(void* b)       { mfa::destroy<Block>(b); }

    static
        void add(                                   // add the block to the decomposition
            int                 gid,                // block global id
            const Bounds<T>&    core,               // block bounds without any ghost added
            const Bounds<T>&    bounds,             // block bounds including any ghost region added
            const Bounds<T>&    domain,             // global data bounds
            const RCLink<T>&    link,               // neighborhood
            diy::Master&        master,             // diy master
            int                 dom_dim,            // domain dimensionality
            int                 pt_dim,             // point dimensionality
            T                   ghost_factor = 0.0) // amount of ghost zone overlap as a factor of block size (0.0 - 1.0)
    {
        mfa::add<Block, T>(gid, core, bounds, domain, link, master, dom_dim, pt_dim, ghost_factor);
    }

    static
        void save(const void* b_, diy::BinaryBuffer& bb)    { mfa::save<Block, T>(b_, bb); }
    static
        void load(void* b_, diy::BinaryBuffer& bb)          { mfa::load<Block, T>(b_, bb); }


    // read a floating point 2d scalar dataset from HDF5
    // reads masses for geometry dimension 0 from same HDF5 file
    // assigns integer values for the geometry dimension 1 from 0 to n_pts - 1
    // f = (mass, y, value)
    //template <typename V>               // type of science value being read
    void read_2d_data(
            const       diy::Master::ProxyWithLink& cp,
            DomainArgs& args,
            Eigen::Tensor<double, 3> vals,
            bool  rescale)            // rescale science values
    {

        DomainArgs* a = &args;
        int tot_ndom_pts = 1;
        this->geometry.min_dim = 0;
        this->geometry.max_dim = this->dom_dim - 1;
        int nvars = 57;
        this->vars.resize(nvars);
        this->max_errs.resize(nvars);
        this->sum_sq_errs.resize(nvars);
	for (int i = 1; i < nvars; i++) {
          this->vars[i].min_dim = this->vars[i-1].max_dim+1;
          this->vars[i].max_dim = this->vars[i].min_dim;
        }
        VectorXi ndom_pts(this->dom_dim);
        this->bounds_mins.resize(this->pt_dim);
        this->bounds_maxs.resize(this->pt_dim);
	//std::cout << "this->dom_dim, a->ndom_pts = " << this->dom_dim << ", " << a->ndom_pts.size() << std::endl;
        for (int i = 0; i < this->dom_dim ; i++)
        { 
            ndom_pts(i)     =  a->ndom_pts[i];
            tot_ndom_pts    *= ndom_pts(i);
	    //std::cout << "tot_ndom_pts, ndom_pts(" << i << ") = " << tot_ndom_pts << ", " << ndom_pts(i) << std::endl; 
	    //std::cout << "ndom_pts(" << i << "), a->ndom_pts[" << i << "] = " << ndom_pts(i) << ", " << a->ndom_pts[i] << std::endl; 
        }
        this->domain.resize(tot_ndom_pts, this->pt_dim);
	
        assert(vals(0) == ndom_pts(0));
        assert(vals(1) == ndom_pts(1));
        // set geometry values
        int n = 0;
        double min_vals = 99999, max_vals = -999.9;
        for (size_t j = 0; j < (size_t)(ndom_pts(1)); j++)
          {
            for (size_t i = 0; i < (size_t)(ndom_pts(0)); i++)
              {
		for (int m = 1; m < nvars; m++) {
                	this->domain(n, 0) = i;
                	this->domain(n, 1) = j;
                	this->domain(n, 2) = vals(m, i,j);

			if( vals(i,j) >= max_vals ) max_vals = vals(m, i,j);
			if( vals(i,j) <= min_vals ) min_vals = vals(m, i,j);
                	n++;
		}
              }
          }
        // find extent of masses and values
        for (size_t i = 0; i < (size_t)this->domain.rows(); i++)
          {
            if (i == 0 || this->domain(i, 0) < this->bounds_mins(0))
              this->bounds_mins(0) = this->domain(i, 0);
            if (i == 0 || this->domain(i, 0) > this->bounds_maxs(0))
              this->bounds_maxs(0) = this->domain(i, 0);
            if (i == 0 || this->domain(i, 1) < this->bounds_mins(1))
              this->bounds_mins(1) = this->domain(i, 1);
            if (i == 0 || this->domain(i, 1) > this->bounds_maxs(1))
              this->bounds_maxs(1) = this->domain(i, 1);
            if (i == 0 || this->domain(i, 2) < this->bounds_mins(2))
              this->bounds_mins(2) = this->domain(i, 2);
            if (i == 0 || this->domain(i, 2) > this->bounds_maxs(2))
              this->bounds_maxs(2) = this->domain(i, 2);
          }
	 // std::cout << "tot_ndom_pt, this->domain(tot_ndom_pts - 1, 1), this->dom_dim = " << tot_ndom_pts << ", " << this->domain(tot_ndom_pts - 1, 1) << ", " << this->dom_dim << std::endl;
        
        // extents
        this->bounds_mins(2) = min_vals;
        this->bounds_maxs(2) = max_vals;
        this->bounds_mins(1) = 0;
        this->bounds_maxs(1) = this->domain(tot_ndom_pts - 1, 1);
        this->bounds_mins(0) = 0;
        this->bounds_maxs(0) = this->domain(tot_ndom_pts - 1, 0);
        this->core_mins.resize(this->dom_dim);
        this->core_maxs.resize(this->dom_dim);
        for (int i = 0; i < this->dom_dim; i++)
          {
            this->core_mins(i) = this->bounds_mins(i);
            this->core_maxs(i) = this->bounds_maxs(i);
	    std::cout << "this->core_mins(" << i << "), this->core_maxs(" << i << ") = " << this->core_mins(i) << ", " << this->core_maxs(i) << std::endl;
          }
        
        this->mfa = new mfa::MFA<T>(this->dom_dim, ndom_pts, this->domain);
        
        // normalize the science variable to the same extent as max of geometry
        if (rescale)
          {
            double extent[3];
            extent[0] = this->bounds_maxs(0) - this->bounds_mins(0);
            extent[1] = this->bounds_maxs(1) - this->bounds_mins(1);
            extent[2] = this->bounds_maxs(2) - this->bounds_mins(2);
            double scale = extent[0] > extent[1] ? extent[0] : extent[1];
            cerr << "\n * rescaling values dividing by " << extent[2] * scale << " *\n" << endl;
            for (size_t i = 0; i < (size_t)this->domain.rows(); i++)
              this->domain(i, 2) = this->domain(i, 2) / extent[2] * scale;
          }
       
        // debug
        cerr << "domain extent:\n min\n" << this->bounds_mins << "\nmax\n" << this->bounds_maxs << endl;
    }
  
  
};

typedef std::array<double, 3> GridPoint;

class GridPoints {
public:
  GridPoints() {}
  GridPoints(std::vector<std::vector<double>> const & m_vec_grid) {
    for (auto gp : m_vec_grid){
      //std::cerr << "gp:  " << gp[0] << ", " << gp[1] << ", " << gp[2] << std::endl;
      _gridpoints.push_back({pow(10, gp[0]), pow(10, gp[1]), pow(10, gp[2])});}
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
  SignalGenerator( sbn::SBNconfig const & conf, std::vector<std::vector<double>> const & vec_grid, size_t dim2,
                   Eigen::MatrixXd const & sinsq,
                   Eigen::MatrixXd const & sin,
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
    auto const & gp = m_gridpoints.Get(i_grid); //grid point after translation: 10^gp[i]
    sbn::NeutrinoModel this_model(gp[0]*gp[0], gp[1], gp[2], false);
    int m_idx = massindex(i_grid);
    Oscillate(m_sinsq.row(m_idx), m_sin.row(m_idx), this_model);
    
    if (compressed) return collapseVectorEigen(m_core+retVec, m_conf);
    else return m_core+retVec;
  }

  Eigen::VectorXd predict2D(size_t i_grid, bool compressed, int &gridx, int &gridy) {
    auto const & gp = m_gridpoints.Get(i_grid); //grid point after translation: 10^gp[i]
    sbn::NeutrinoModel this_model(gp[0]*gp[0], gp[1], gp[2], false);
    int m_idx = massindex(i_grid);
    gridx = gp[2];
    gridy = gp[0];
    //std::cerr << "i_grid: " << i_grid << " mass index " << m_idx << " dim2: " << m_dim2 << " GP:" << gp <<  "\n";
    Oscillate(m_sinsq.row(m_idx), m_sin.row(m_idx), this_model);
    
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
  //std::vector<Eigen::VectorXd> m_sinsq, m_sin;
  Eigen::MatrixXd m_sinsq, m_sin;
  Eigen::VectorXd m_core, retVec;
  int m_oscmode;
};
inline void makeSignalModel(diy::Master* master, SignalGenerator signal, int nbins)
{
    // default command line arguments
    int    pt_dim       = 3;                    // dimension of input points
    int    dom_dim      = 2;                    // dimension of domain (<= pt_dim)
    int    geom_degree  = 1;                    // degree for geometry (same for all dims)
    int    vars_degree  = 2;                    // degree for science variables (same for all dims)
    int    geom_nctrl   = geom_degree + 1;     // input number of control points for geometry (same for all dims)
    int    vars_nctrl   = vars_degree + 1;     // input number of control points for geometry (same for all dims)
    int    weighted     = 1;                   // input number of control points for all science variables (same for all dims)
    real_t noise        = 0.0;                  // fraction of noise
    // set default args for diy foreach callback functions
    DomainArgs d_args(dom_dim, pt_dim);
    d_args.weighted     = weighted;
    d_args.n            = noise;
    d_args.multiblock   = false;
    d_args.verbose      = 1;
 
    for (int i = 0; i < pt_dim - dom_dim; i++)
        d_args.f[i] = 1.0;
    for (int i = 0; i < dom_dim; i++)
    {
        d_args.geom_p[i]            = geom_degree;
        d_args.vars_p[0][i]         = vars_degree;  // assuming one science variable, vars_p[0]
        d_args.geom_nctrl_pts[i]    = geom_nctrl;
    }
    d_args.ndom_pts[1]          = 26;
    d_args.ndom_pts[0]          = 26;
    d_args.vars_nctrl_pts[0][1] = 26;       // assuming one science variable, vars_nctrl_pts[0]
    d_args.vars_nctrl_pts[0][0] = 26;       // assuming one science variable, vars_nctrl_pts[0]

   Eigen::VectorXd vec_gridx(26);
   Eigen::VectorXd vec_gridy(26);
   Eigen::Tensor<double, 3> map_bin_to_grid(57,26,26);
   Eigen::ArrayXXd values(signal.gridsize(),nbins);	
   int gridx = -1;
   int gridy = -1;
   for (size_t i=0; i<signal.gridsize(); ++i) {
      values.row(i) = signal.predict2D(i, true, gridx, gridy);
      int gridy_index = i/26;
      int gridx_index = i%26;
      vec_gridy(gridy_index) = gridy;
      vec_gridx(gridx_index) = gridx;
      for( int bin=0; bin < nbins; bin++ ) map_bin_to_grid(bin,gridx_index,gridy_index) = values(i,bin); 
   }

   master->foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp){ blocks[bin]->read_2d_data(cp,d_args,map_bin_to_grid,true); });
   master->foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp){ blocks[bin]->fixed_encode_block(cp,d_args); });

}

void loadData(const char* fname, std::string what, std::vector<double> & v_buffer, int & n_rows, int & n_cols) {
  H5Easy::File file(fname, H5Easy::File::ReadOnly);
  Eigen::MatrixXd _mat = H5Easy::load<Eigen::MatrixXd>(file, what);
  n_rows = _mat.rows();
  n_cols = _mat.cols();
  v_buffer = std::vector<double>(_mat.data(), _mat.data() + _mat.rows() * _mat.cols());
}

void loadData(const char* fname, std::string what, std::vector<int> & v_buffer, int & n_rows, int & n_cols) {
  H5Easy::File file(fname, H5Easy::File::ReadOnly);
  Eigen::MatrixXi _mat      = H5Easy::load<Eigen::MatrixXi>(file, what);
  n_rows = _mat.rows();
  n_cols = _mat.cols();
  v_buffer = std::vector<int>(_mat.data(), _mat.data() + _mat.rows() * _mat.cols());
}

Eigen::MatrixXd bcMatrixXd(diy::mpi::communicator world, std::vector<double>  v_buffer, int  n_rows, int  n_cols) {
  diy::mpi::broadcast(world, v_buffer, 0);
  diy::mpi::broadcast(world, n_rows,   0);
  diy::mpi::broadcast(world, n_cols,   0);

  Eigen::Map<Eigen::MatrixXd> mat(v_buffer.data(), n_rows, n_cols);
  return mat;
}

Eigen::MatrixXi bcMatrixXi(diy::mpi::communicator world, std::vector<int>  v_buffer, int  n_rows, int  n_cols) {
  diy::mpi::broadcast(world, v_buffer, 0);
  diy::mpi::broadcast(world, n_rows,   0);
  diy::mpi::broadcast(world, n_cols,   0);
  
  Eigen::Map<Eigen::MatrixXi> mat(v_buffer.data(), n_rows, n_cols);
  return mat;
}

void createDataSets(HighFive::File* file, size_t nPoints, size_t nUniverses) {
  file->createDataSet<double>("last_chi_min", HighFive::DataSpace( { nPoints*nUniverses,       1} ));
  file->createDataSet<double>("delta_chi",    HighFive::DataSpace( { nPoints*nUniverses,       1} ));
  file->createDataSet<int>("best_grid_point", HighFive::DataSpace( { nPoints*nUniverses,       1} ));
  file->createDataSet<int>("n_iter",          HighFive::DataSpace( { nPoints*nUniverses,       1} ));
  file->createDataSet<double>("fakedataC",    HighFive::DataSpace( { nPoints*nUniverses,       57} ));
  file->createDataSet<double>("collspec",     HighFive::DataSpace( { nPoints*nUniverses,      57} ));
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
    for (auto fn: fullnames) {
       hist.push_back(*((TH1D*)f.Get(fn.c_str())));
    }
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

//need to figure out the proper way to propagate the MFA model
namespace cppoptlib {
template<typename T>
class LLR : public Problem<T> {
    private:
        Block<real_t>&     b;            //block encoded with MFA model
        diy::Master::ProxyWithLink const& cp;
        VectorXd           data;
	const MatrixX<T>   M;

  public:
    //need mfa model, fake_data, and inverse covariance matrix
    LLR( Block<real_t>& b_,  
		    diy::Master::ProxyWithLink const& cp_,
		    VectorXd data_, 
		    MatrixX<T> M_) : b(b_), 
				     cp(cp_),	
				     data(data_),	
				     M(M_) 	
	{}
    //TO DO: Figure out the correct way to pass the MFA 
    using typename Problem<T>::TVector;
    T value(const TVector &x) {
	//to calculate the chi2 we need the to calculate the difference between the fake data and MFA model in each bin of data/model (point at dom_dim == 0)
    	// evaluate point
	TVector diff(data.size()); //data size == nvars
	VectorX<T> param(b.dom_dim);
	for(int p=0; p<data.size(); p++){
 		param(p) = x[1];
 	    	int dom_dim = b->dom_dim;
     		int pt_dim  = b->pt_dim;
		//already scaled in read_2d_data
		//scale to 1
                // normalize the science variable to the same extent as max of geometry
                //double extent[3];
                //extent[0] = b->bounds_maxs(0) - b->bounds_mins(0);
                //extent[1] = b->bounds_maxs(1) - b->bounds_mins(1);
                //extent[2] = b->bounds_maxs(2) - b->bounds_mins(2);
                //double scale = extent[0] > extent[1] ? extent[0] : extent[1];
                //std::cerr << "\n * rescaling values dividing by " << extent[2] * scale << " *\n" << std::endl;
                for (size_t i = 0; i < (size_t)blocks[p]->domain.rows(); i++)
                   blocks[p]->domain(i, 2) = blocks[p]->domain(i, 2) / extent[2] * scale;
	    	int dom_dim = blocks[p].dom_dim;
    		int pt_dim  = blocks[p].pt_dim;
    		VectorX<real_t> out_pt(pt_dim);
   		// parameters of input point to evaluate
    		VectorX<real_t> in_param(dom_dim);
    		for (auto i = 0; i < dom_dim; i++){
        		in_param(i) = param[i];
    		}
		b.decode_point(cp, in_param, out_pt);
    		diff[p] = data[p]-out_pt(pt_dim-1);
		std::cerr << "p, diff, data[p], out_pt(pt_dim-1) = " << p << ", " << diff[p] << ", " << data[p] << ", " << out_pt(pt_dim-1) << std::endl;
	}
	return (diff.transpose())*M*(diff);
    }
/*
    void gradient(const TVector &x, TVector &grad) {
	for(int p=0; p<data.size(); p++){
		//scale to 1
                // normalize the science variable to the same extent as max of geometry
                double extent[3];
                extent[0] = blocks[p]->bounds_maxs(0) - blocks[p]->bounds_mins(0);
                extent[1] = blocks[p]->bounds_maxs(1) - blocks[p]->bounds_mins(1);
                extent[2] = blocks[p]->bounds_maxs(2) - blocks[p]->bounds_mins(2);
                double scale = extent[0] > extent[1] ? extent[0] : extent[1];
                std::cerr << "\n * rescaling values dividing by " << extent[2] * scale << " *\n" << std::endl;
                for (size_t i = 0; i < (size_t)blocks[p]->domain.rows(); i++)
                   blocks[p]->domain(i, 2) = blocks[p]->domain(i, 2) / extent[2] * scale;
	    	int dom_dim = blocks[p].dom_dim;
    		int pt_dim  = blocks[p].pt_dim;
    		VectorX<real_t> out_pt(pt_dim);
   		// parameters of input point to evaluate
    		VectorX<real_t> in_param(dom_dim);
    		for (auto i = 0; i < dom_dim; i++){
        		in_param(i) = param[i];
    		}
		// evaluate point
		VectorX<real_t> out_pt(pt_dim);
		VectorX<real_t> out_pt_deriv(pt_dim);
		// parameters of input point to evaluate
		VectorX<real_t> in_param(dom_dim);
    		for (auto i = 0; i < dom_dim; i++){
       			in_param(i) = param[i];
		}
    		for (auto i = 0; i < dom_dim; i++){
			b.differentiate_point(cp, in_param, 1, i, -1, out_pt_deriv);
		}
		b.decode_point(cp, in_param, out_pt);
	        //grad
		//scale the derivation u,v back to x,y
	        grad[p] = 2*scale*out_pt_deriv(pt_dim - 1)*M*(out_pt(pt_dim-1)-data[p]);
	}
    }*/

};
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
   //auto const & goodpoints  = initialScan(fake_data, Eigen::MatrixXd::Identity(INVCOV.rows(), INVCOV.rows()), signal, 1e4);
   //auto const & goodpoints  = initialScan(fake_data, invcov, signal, 1e6);
   
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
   //convert eigen vector to regular vector since I don't know how to converge it with create_datasets
   //assert that fakedataC should have the same dimension as collspec
   if(fake_data.size() != v_coll.size() ) std::cout << "check the collapsing method!" << std::endl;
   std::vector<double> fakedataC, collspec;
   for(uint i=0; i < fake_data.size(); i++) fakedataC.push_back(fake_data(i));
   for(uint i=0; i < v_coll.size(); i++) collspec.push_back(v_coll(i));

   FitResult fr = {n_iter, best_grid_point, last_chi_min, this_chi-last_chi_min, fakedataC, collspec}; 
   return fr;
}

//new implementation using MFA 
inline FitResult coreFC(Eigen::VectorXd const & fake_data, Eigen::VectorXd const & v_coll,
      Block<real_t>* b, //replace with block that has encoded MFA model
      diy::Master::ProxyWithLink const& cp,
      int i_grid,
      Eigen::MatrixXd const & INVCOV,
      Eigen::MatrixXd const & covmat,
      sbn::SBNconfig const & myconf,
      double chi_min_convergance_tolerance = 0.001,
      size_t max_number_iterations = 5
      )
{
   float global_chi_min = FLT_MAX;
   int best_grid_point = -99;
   size_t n_iter = 0;

   Eigen::MatrixXd invcov = INVCOV;//std::vector<double> temp;

   float chi_min = FLT_MAX;
   //Step 2.0 Find the global_minimum_for this universe. Integrate in SBNfit minimizer here, a grid scan previously
   //implement optimizer
   typedef double T;
   //using typename cppoptlib::Problem<T>::TVector;
   //want to pass the mfa model here which is already encode in the block
   typedef LLR<T> llr;
   llr f(*b, cp, fake_data, INVCOV);
   cppoptlib::LbfgsSolver<llr> solver;
   //minimize the function
   Eigen::VectorXd x(fake_data.size()); 
   for( int p=0; p < fake_data.size(); p++) x << i_grid;
   //solver.minimize(f,x);
  
   //store the result here
   //global_chi_min = f(x); //the global minimum chi2
   //best_grid_point = x; //point in grid which gives the minimum chi2

   //Now use the curent_iteration_covariance matrix to also calc this_chi here for the delta.
   float this_chi = calcChi(fake_data, v_coll, invcov);
   //assert that fakedataC should have the same dimension as collspec
   if(fake_data.size() != v_coll.size() ) std::cout << "check the collapsing method!" << std::endl;
   std::vector<double> fakedataC, collspec;
   for(uint i=0; i < fake_data.size(); i++) fakedataC.push_back(fake_data(i));
   for(uint i=0; i < v_coll.size(); i++) collspec.push_back(v_coll(i));

   FitResult fr = {n_iter, best_grid_point, global_chi_min, this_chi-global_chi_min, fakedataC, collspec}; 
   return fr;
}

void doScan(Block<real_t>* b, diy::Master::ProxyWithLink const& cp, int rank,
      sbn::SBNconfig const & myconf,
      Eigen::MatrixXd const & ECOV, Eigen::MatrixXd const & INVCOVBG,
      Eigen::VectorXd const & ecore, SignalGenerator signal,
      HighFive::File* file, std::vector<int> const & rankwork, 
      double tol, size_t iter, bool debug)
{

    double starttime, endtime;
    std::vector<FitResult> results;
    std::vector<int> v_grid, v_univ, v_iter, v_best;
    std::vector<double> v_last, v_dchi;
    std::vector<std::vector<double> > v_fakedataC, v_collspec;
    
    results.reserve(rankwork.size());
    v_grid.reserve(rankwork.size());
    v_univ.reserve(rankwork.size());
    v_iter.reserve(rankwork.size());
    v_best.reserve(rankwork.size());
    v_last.reserve(rankwork.size());
    v_dchi.reserve(rankwork.size());
    v_fakedataC.reserve(rankwork.size());
    v_collspec.reserve(rankwork.size());
   
    //validation
    Eigen::MatrixXd specfull_e_mat(rankwork.size(),57); 
    Eigen::MatrixXd speccoll_mat(rankwork.size(),57); 
    Eigen::MatrixXd corecoll_mat(rankwork.size(),57);

    system_clock::time_point t_init = system_clock::now();

    //remove this loop
    for (int i_grid : rankwork) {

       if (debug && i_grid!=0) return;
       auto const & specfull_e = signal.predict(i_grid, false);
       auto const & speccoll   = collapseVectorEigen(specfull_e, myconf);
       auto const & corecoll   = collapseVectorEigen(ecore, myconf);
       
       //Eigen
       specfull_e_mat.row(i_grid) = Eigen::VectorXd::Map(&specfull_e[i_grid], specfull_e_mat.size());
       speccoll_mat.row(i_grid) = Eigen::VectorXd::Map(&speccoll[0], speccoll.size());
       corecoll_mat.row(i_grid) = Eigen::VectorXd::Map(&corecoll[0], corecoll.size());

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
       v_fakedataC.push_back(res.fakedataC);
       v_collspec.push_back(res.collspec);
    }

    d_last_chi_min.select(     {d_bgn, 0}, {size_t(v_last.size()), 1}).write(v_last);
    d_delta_chi.select(        {d_bgn, 0}, {size_t(v_dchi.size()), 1}).write(v_dchi);
    d_best_grid_point.select(  {d_bgn, 0}, {size_t(v_best.size()), 1}).write(v_best);
    d_n_iter.select(           {d_bgn, 0}, {size_t(v_iter.size()), 1}).write(v_iter);
    d_i_grid.select(           {d_bgn, 0}, {size_t(v_grid.size()), 1}).write(v_grid);
    d_i_univ.select(           {d_bgn, 0}, {size_t(v_univ.size()), 1}).write(v_univ);
    endtime   = MPI_Wtime();
    //if (world.rank()==0){
	H5Easy::File file1("/code/src/comparespectrum_mpi.h5", H5Easy::File::Overwrite);
  	H5Easy::dump(file1, "specfull", specfull_e_mat);
  	H5Easy::dump(file1, "colspec", speccoll_mat);
  	H5Easy::dump(file1, "colcore", corecoll_mat); 
    //}
    if (cp.gid()==0) fmt::print(stderr, "[{}] Write out took {} seconds\n", cp.gid(), endtime-starttime);
}

void doMin(Block<real_t>* b, diy::Master::ProxyWithLink const& cp, int rank,
      const char * xmldata, sbn::SBNconfig const & myconf,
      TMatrixD const & covmat, Eigen::MatrixXd const & ECOV, Eigen::MatrixXd const & INVCOVBG,
      SignalGenerator signal,
      HighFive::File* file, std::vector<int> const & rankwork, int nUniverses, 
      double tol, size_t iter, bool debug, bool noWrite=false, int msg_every=100)
{ 

    //return 0;
}

// TODO add size_t writeEvery to prevent memory overload
void doFC(Block<real_t>* b, diy::Master::ProxyWithLink const& cp, int rank,
      const char * xmldata, sbn::SBNconfig const & myconf,
      TMatrixD const & covmat, Eigen::MatrixXd const & ECOV, Eigen::MatrixXd const & INVCOVBG,
      SignalGenerator signal,
      HighFive::File* file, std::vector<int> const & rankwork, int nUniverses, 
      double tol, size_t iter, bool debug, bool noWrite=false, int msg_every=100)
{
    std::cout << "*** doFC ***" << std::endl;
    double starttime, endtime;
    std::vector<FitResult> results;
    std::vector<int> v_grid, v_univ, v_iter, v_best;
    std::vector<double> v_last, v_dchi;
    std::vector<std::vector<double> > v_fakedataC, v_collspec;
   
    if (!noWrite) {
       results.reserve(rankwork.size()*nUniverses);
       v_grid.reserve(rankwork.size()*nUniverses);
       v_univ.reserve(rankwork.size()*nUniverses);
       v_iter.reserve(rankwork.size()*nUniverses);
       v_best.reserve(rankwork.size()*nUniverses);
       v_last.reserve(rankwork.size()*nUniverses);
       v_dchi.reserve(rankwork.size()*nUniverses);
       v_fakedataC.reserve(rankwork.size()*nUniverses);
       v_collspec.reserve(rankwork.size()*nUniverses);
    }
     
     
    //validation
    Eigen::MatrixXd specfull_e_mat(rankwork.size(),57); 
    Eigen::MatrixXd speccoll_mat(rankwork.size(),57); 
    Eigen::MatrixXd corecoll_mat(rankwork.size(),57);

    system_clock::time_point t_init = system_clock::now();
    for (int i_grid : rankwork) {

       if (debug && i_grid!=0) return;
       auto const & specfull_e = signal.predict(i_grid, false);
       auto const & speccoll = collapseVectorEigen(specfull_e, myconf);
       std::mt19937 rng(cp.gid()); // Mersenne twister
       Eigen::MatrixXd const & LMAT = cholD(ECOV, specfull_e);
     
         
       //Eigen
       specfull_e_mat.row(i_grid) = Eigen::VectorXd::Map(&specfull_e[i_grid], specfull_e_mat.size());
       speccoll_mat.row(i_grid) = Eigen::VectorXd::Map(&speccoll[i_grid], speccoll.size());

       std::vector<Block<real_t>* > b;
       //creates a linear chain of blocks like in https://github.com/diatomic/diy/blob/master/examples/simple/iexchange-particles.cpp?
       for (int uu=0; uu<nUniverses;++uu) {
          auto const & fake_data = sample(specfull_e, LMAT, rng);//
          auto const & fake_dataC = collapseVectorEigen(fake_data, myconf);
          
          results.push_back(coreFC(fake_dataC, speccoll,
                   //signal, INVCOVBG, ECOV, myconf, tol, iter)); //old implementation
                   b, cp, i_grid, INVCOVBG, ECOV, myconf, tol, iter));

          v_univ.push_back(uu);
          v_grid.push_back(i_grid);
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
       // write out this grid and universe
       HighFive::DataSet d_i_grid          = file->getDataSet("i_grid");
       HighFive::DataSet d_i_univ          = file->getDataSet("i_univ");
       // This is for the fake data dump
       HighFive::DataSet d_fakedataC      = file->getDataSet("fakedataC");
       HighFive::DataSet d_collspec        = file->getDataSet("collspec");

       size_t d_bgn = rankwork[0]*nUniverses;
       for (auto res : results) {
          v_iter.push_back(res.n_iter);
          v_best.push_back(res.best_grid_point);
          v_last.push_back(res.last_chi_min);
          v_dchi.push_back(res.delta_chi);
	  v_fakedataC.push_back(res.fakedataC);
	  v_collspec.push_back(res.collspec);
       }

       d_last_chi_min.select(     {d_bgn, 0}, {size_t(v_last.size()), 1}).write(v_last);
       d_delta_chi.select(        {d_bgn, 0}, {size_t(v_dchi.size()), 1}).write(v_dchi);
       d_best_grid_point.select(  {d_bgn, 0}, {size_t(v_best.size()), 1}).write(v_best);
       d_n_iter.select(           {d_bgn, 0}, {size_t(v_iter.size()), 1}).write(v_iter);
       d_i_grid.select(           {d_bgn, 0}, {size_t(v_grid.size()), 1}).write(v_grid);
       d_i_univ.select(           {d_bgn, 0}, {size_t(v_univ.size()), 1}).write(v_univ);
       d_fakedataC.select(        {d_bgn, 0}, {size_t(v_fakedataC.size()), 57}).write(v_fakedataC);
       d_collspec.select(         {d_bgn, 0}, {size_t(v_collspec.size()), 57}).write(v_collspec);
       endtime   = MPI_Wtime();
       if (cp.gid()==0) fmt::print(stderr, "[{}] Write out took {} seconds\n", cp.gid(), endtime-starttime);
	H5Easy::File file1("/code/src/comparespectrum_mpi.h5", H5Easy::File::Overwrite);
  	H5Easy::dump(file1, "specfull", specfull_e_mat);
  	H5Easy::dump(file1, "colspec", speccoll_mat);
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
    std::string out_file="test.h5";
    std::string f_BG="BOTHv2_BKG_ONLY.SBNspec.root";
    std::string f_CV="BOTHv2_CV.SBNspec.root";
    std::string f_COV="BOTHv2.SBNcovar.root";
    std::string tag="";
    std::string d_in="";
    std::string xml="";
    std::string infile0  = "corespectrum.h5";                 // h5 input file
    std::string infile1  = "SIN_mfa_output_deg2.h5";                 // MFA input file
    std::string infile2  = "SINSQ_mfa_output_deg2.h5";                 // MFA input file

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
    ops >> Option("msg",             msg_every,  "Print a progress message every m gridpoints on rank 0 to stderr.");
    ops >> Option("mode",            mode, "Mode 0 is default --- dimension 2 is electron, mode 1 is muon");
    bool debug       = ops >> Present('d', "debug", "Operate on single gridpoint only");
    bool statonly    = ops >> Present("stat", "Statistical errors only");
    bool nowrite    = ops >> Present("nowrite", "Don't write output --- for performance estimates only");
    bool simplescan    = ops >> Present("scan", "Simple scan, no FC");
    bool mfa    = ops >> Present("mfa", "Use MFA input file");
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


    std::vector<double> v_buff;
    int ncols, nrows;
    if (world.rank()==0) loadData(Form("%s/%s",d_in.c_str(),infile0.c_str()), "core",        v_buff,   nrows, ncols);
    Eigen::VectorXd _core = bcMatrixXd(world, v_buff, nrows, ncols);
  
    if (world.rank()==0) loadData(Form("%s/%s",d_in.c_str(),infile0.c_str()), "covmat",        v_buff,   nrows, ncols);
    Eigen::MatrixXd _covmat = bcMatrixXd(world, v_buff, nrows, ncols);
  
    if (world.rank()==0) loadData(Form("%s/%s",d_in.c_str(),infile0.c_str()), "invcovbg",        v_buff,   nrows, ncols);
    Eigen::MatrixXd _invcovbg = bcMatrixXd(world, v_buff, nrows, ncols);
  
    if (world.rank()==0) loadData(Form("%s/%s",d_in.c_str(),infile1.c_str()), "_SIN_",        v_buff,   nrows, ncols);
    std::cout << Form("%s/%s",d_in.c_str(),infile1.c_str()) << std::endl;
    Eigen::MatrixXd _sin = bcMatrixXd(world, v_buff, nrows, ncols);
  
    if (world.rank()==0) loadData(Form("%s/%s",d_in.c_str(),infile2.c_str()), "_SINSQ_",        v_buff,   nrows, ncols);
    Eigen::MatrixXd _sinsq = bcMatrixXd(world, v_buff, nrows, ncols);
   
    //std::cout << "_sinsq nrows, ncols = " << nrows << ", " << ncols << std::endl;

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

    std::cout << "sinsqvec size, sinsqvec_eig size = " << sinsqvec.size() << ", " << sinsqvec_eig.size() << std::endl;
    std::cout << "nrows, ncols = " << nrows << ", " << ncols << std::endl;
    Eigen::MatrixXd sinsqeig(nrows-1,ncols);
    Eigen::MatrixXd sineig(nrows-1,ncols);

    for (int i = 0; i < sinsqvec_eig.size(); ++i){
        sinsqeig.row(i) = Eigen::VectorXd::Map(&sinsqvec_eig[i][0], sinsqvec_eig[0].size());
        sineig.row(i) = Eigen::VectorXd::Map(&sinvec_eig[i][0], sinvec_eig[0].size());
    }

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
       //std::cout << "mode 2, ymin, ymax = " << ymin << ", " << ymax << std::endl;
    }
    else if (mode==1) {
       mygrid.AddFixedDimension("ue4", 0.0);
       mygrid.AddDimension("um4", ymin, ymax, ywidth);
       dim2 = mygrid.f_dimensions[2].f_N;
       //std::cout << "mode 1, ymin, ymax = " << ymin << ", " << ymax << std::endl;
    }
    else {
       std::cerr << "Error, the mode must be either 0 or 1 a the moment: " << mode << "\n";
       exit(1);
    }
    nPoints = mygrid.f_num_total_points;
    GridPoints GP(mygrid.GetGrid());

    if (world.rank()==0) mygrid.Print();

    // Finally, the signal generator
    Eigen::MatrixXd _sinsq_(nrows,ncols);
    Eigen::MatrixXd _sin_(nrows,ncols);

    if(!mfa){
	_sinsq_ = sinsqeig;  
	_sin_ = sineig;
    } else {
	_sinsq_ = _sinsq; 
	_sin_ = _sin;
    }

    SignalGenerator  signal(myconf, mygrid.GetGrid(), dim2, _sinsq_, _sin_, ecore, mode);
   

    //write the mass index and gridpoint
    //Eigen::MatrixXd masses;
    //Eigen::MatrixXd masses;
    /*std::vector<Eigen::MatrixXd> map_grid_speccoll;  // 3 dimensions (4x5x6)
    Eigen::MatrixXd gridpt_spec_bin;
    std::vector<std::vector<double> > m_vec_grid = mygrid.GetGrid();
    for( int bin = 0; bin < speccol.size(); bin++ ) map_grid_speccoll.push_back(speccol(bin));
    for( uint igrid = 0; igrid < nPoints; igrid++ ){
       auto const & speccol = signal.predict(igrid, true);
       int row = igrid/26;
       int col = igrid%26;
    }*/
    /*for(int i=0; i < m_vec_grid.size(); i++){
    for(int j=0; j < m_vec_grid[i].size(); j++){
	    std::cout<< " m_vec_grid.size(), m_vec_grid[i].size(), m_vec_grid[i][j] = " <<  m_vec_grid.size() << ", " << m_vec_grid[i].size() << " = " << m_vec_grid[i][j] << std::endl;
    }
    }*/

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

     std::cout << "before world.size() = " << world.size() << std::endl;
    // Create datasets needed TODO nbinsC --- can we get that from somewhere?
    createDataSets(f_out, nPoints, nUniverses);
   
     std::cout << "after world.size() = " << world.size() << std::endl;
    // First rank also writes the grid so we know what the poins actually are
    if (world.rank() == 0)  writeGrid(f_out, mygrid.GetGrid(), mode);


    //// Now more blocks as we have universes
    size_t blocks = world.size();//nPoints;//*nUniverses;
    std::cout << "blocks = " << blocks << std::endl;
    if (world.rank()==0) fmt::print(stderr, "FC will be done on {} blocks, distributed over {} ranks\n", blocks, world.size());
    Bounds<real_t> fc_domain(2);
    fc_domain.min[0] = 0;
    fc_domain.max[0] = blocks-1;
    fc_domain.min[1] = 0;
    fc_domain.max[1] = blocks-1;
    diy::FileStorage               storage("./DIY.XXXXXX");
    diy::RoundRobinAssigner        fc_assigner(world.size(), blocks);
    diy::RegularDecomposer<Bounds<real_t>> fc_decomposer(2, fc_domain, blocks);
    diy::RegularBroadcastPartners  fc_comm(    fc_decomposer, 2, true);
    diy::RegularMergePartners      fc_partners(fc_decomposer, 2, true);
    diy::Master                    fc_master(world, 1, -1, &Block<real_t>::create, &Block<real_t>::destroy, &storage, &Block<real_t>::save, &Block<real_t>::load);
    diy::ContiguousAssigner   assigner(world.size(), blocks);
    fc_decomposer.decompose(world.rank(),
                         assigner,
                         [&](int gid, const Bounds<real_t>& core, const Bounds<real_t>& bounds, const Bounds<real_t>& domain, const RCLink<real_t>& link)
                         { Block<real_t>::add(gid, core, bounds, domain, link, fc_master, 2, 3, 0.0); });

     
    //maybe set up a different MFA blocks and decomposer here?
    //std::vector< Block<real_t>* > blks;
    makeSignalModel(&fc_master,signal,57);


    std::vector<int> work(nPoints);
    std::iota(std::begin(work), std::end(work), 0); //0 is the starting number
    std::vector<int> rankwork = splitVector(work, world.size())[world.rank()];
    work.clear();
    work.shrink_to_fit();

    double starttime = MPI_Wtime();
    if (simplescan) {
       if (world.rank()==0) fmt::print(stderr, "Start simple scan\n");
       if ( !mfa ) { fc_master.foreach([world, ECOVMAT, INVCOVBG, ecore, myconf,  f_out, rankwork, tol, iter, debug, signal](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
                              { doScan(b, cp, world.rank(), myconf, ECOVMAT, INVCOVBG, ecore, signal, f_out, rankwork, tol, iter, debug); });}
       else { fc_master.foreach([world, _covmat, _invcovbg, _core, myconf,  f_out, rankwork, tol, iter, debug, signal](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
                              { doScan(b, cp, world.rank(), myconf, _covmat, _invcovbg, _core, signal, f_out, rankwork, tol, iter, debug); });}
    }
    else {
       if (world.rank()==0) fmt::print(stderr, "Start FC\n");
       if( !mfa ){ fc_master.foreach([world, covmat, ECOVMAT, xmldata, INVCOVBG, myconf, nUniverses, f_out, rankwork, tol, iter, debug, signal, nowrite, msg_every](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
                              { doFC(b, cp, world.rank(), xmldata, myconf, covmat, ECOVMAT, INVCOVBG, signal, f_out, rankwork, nUniverses, tol, iter, debug, nowrite, msg_every); });}
       else { fc_master.foreach([world, covmat, _covmat, xmldata, _invcovbg, myconf, nUniverses, f_out, rankwork, tol, iter, debug, signal, nowrite, msg_every](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
                              { doFC(b, cp, world.rank(), xmldata, myconf, covmat, _covmat, _invcovbg, signal, f_out, rankwork, nUniverses, tol, iter, debug, nowrite, msg_every); });}
    }

    double endtime   = MPI_Wtime();
    if(mfa) std::cout << "This is the MFA implementation!" << std::endl; 
    if (world.rank()==0) fmt::print(stderr, "[{}] that took {} seconds\n",world.rank(), endtime-starttime);
    if (world.rank()==0) fmt::print(stderr, "Output written to {}\n",out_file);
    

    //Print comparisons before FC procedure 
    std::string filename;
    if(!mfa) filename="/code/src/comparesinsinsq.h5";
    else filename="/code/src/comparesinsinsq_mfa.h5";
    if(world.rank()==0){
	H5Easy::File file(filename.c_str(), H5Easy::File::Overwrite);
  	H5Easy::dump(file, "_sinsq_", _sinsq_);
  	H5Easy::dump(file, "_sin_", _sin_);
  	H5Easy::dump(file, "_sinsq", _sinsq);
  	H5Easy::dump(file, "_sin", _sin);
  	H5Easy::dump(file, "sinsqeig", sinsqeig); 
  	H5Easy::dump(file, "sineig", sineig); 
  	H5Easy::dump(file, "ECOVMAT", ECOVMAT); 
  	H5Easy::dump(file, "_covmat", _covmat); 
  	H5Easy::dump(file, "ecore", ecore); 
  	H5Easy::dump(file, "_core", _core); 
    }
    if(world.rank()==0){
	H5Easy::File file("_SINSQ_.h5", H5Easy::File::Overwrite);
  	H5Easy::dump(file, "_SINSQ_", _sinsq_);
  	H5Easy::dump(file, "masses", msqsplittings);
    }
    if(world.rank()==0){
	H5Easy::File file("_SIN_.h5", H5Easy::File::Overwrite);
  	H5Easy::dump(file, "_SIN_", _sin_);
  	H5Easy::dump(file, "masses", msqsplittings);
    }

    delete f_out;
    return 0;
}
