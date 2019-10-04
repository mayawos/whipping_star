//#pragma GCC optimize("O3","unroll-loops","inline")

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


using namespace std;

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

class GridPoints {
   public:
      GridPoints(std::vector<std::vector<double>> const & m_vec_grid) {
         for (auto gp : m_vec_grid)
            _gridpoints.push_back({pow(10, gp[0]), pow(10, gp[1]), pow(10, gp[2])});
      }

      size_t NPoints()         { return _gridpoints.size(); }
      GridPoint Get(size_t index) { return _gridpoints[index]; }

   private:
      std::vector<GridPoint> _gridpoints;

};

struct SignalGenerator {
   //sbn::SBNosc osc; sbn::SBNconfig const & conf;
   sbn::SBNconfig const & conf;
   GridPoints m_gridpoints; size_t dim2;
   std::vector<Eigen::VectorXd>  const & sinsq;
   std::vector<Eigen::VectorXd>  const & sin;
   Eigen::VectorXd const & core;

   Eigen::VectorXd predict(size_t i_grid, bool compressed) {
      auto const & gp = m_gridpoints.Get(i_grid);
      sbn::NeutrinoModel this_model(gp[0]*gp[0], gp[1], gp[2], false);
      int m_idx = massindex(i_grid);
      auto ans = Oscillate(sinsq[m_idx], sin[m_idx], this_model, conf);

      ans+=core;
      if (compressed) return collapseVectorEigen(ans, conf);
      else return ans;
   }
   
   int massindex(size_t igrd) {return int(floor( (igrd) / dim2 ));}
   size_t gridsize() {return m_gridpoints.NPoints();}


   Eigen::VectorXd Oscillate(Eigen::VectorXd const & sf_sinsq, Eigen::VectorXd const & sf_sin, 
         sbn::NeutrinoModel & working_model, sbn::SBNconfig const & conf) 
   {
       Eigen::VectorXd retVec(conf.num_bins_total);
       retVec.setZero(); // !!!
       int which_dm = 41;

       double prob_mumu(0), prob_ee(0), prob_mue(0), prob_mue_sq(0), prob_muebar(0), prob_muebar_sq(0);

       prob_mue       = working_model.oscAmp( 2,  1, which_dm, 1);
       prob_mue_sq    = working_model.oscAmp( 2,  1, which_dm, 2);
       prob_muebar    = working_model.oscAmp(-2, -1, which_dm, 1);
       prob_muebar_sq = working_model.oscAmp(-2, -1, which_dm, 2);
       double osc_amp(0), osc_amp_sq(0);
       int osc_pattern(0);
           // Iterate over channels
        size_t offset(0);
        for (int i=0; i<conf.num_channels; i++) {
            size_t nbins_chan = conf.num_bins[i];
            auto const & thisPattern = conf.subchannel_osc_patterns[i];
            for (int j=0; j<conf.num_subchannels[i]; j++){
                osc_pattern = thisPattern[j];
                switch (osc_pattern){
                    case 11:
                        osc_amp_sq = prob_ee;
                        break;
                    case -11:
                        osc_amp_sq = prob_ee;
                        break;
                    case 22:
                        osc_amp_sq = prob_mumu;
                        break;
                    case -22:
                        osc_amp_sq = prob_mumu;
                        break;
                    case 21:
                        osc_amp    = prob_mue;
                        osc_amp_sq = prob_mue_sq;
                        break;
                    case -21:
                        osc_amp    = prob_muebar;
                        osc_amp_sq = prob_muebar_sq;
                        break;
                    case 0:
                    default:
                        break;
                }

                // Iterate over detectors
                for (int d=0; d<conf.num_detectors;++d) {
                  size_t first  = d*conf.num_bins_detector_block + offset;
                  retVec.segment(first, nbins_chan) += osc_amp   *  sf_sin.segment(first, nbins_chan);
                  retVec.segment(first, nbins_chan) += osc_amp_sq*sf_sinsq.segment(first, nbins_chan);
                }
                offset +=nbins_chan;
            }
        }

    return retVec;
   }
};

void createDataSets(HighFive::File* file, size_t nPoints, size_t nUniverses) {
    file->createDataSet<double>("last_chi_min", HighFive::DataSpace( { nPoints*nUniverses,       1} ));
    file->createDataSet<double>("delta_chi",    HighFive::DataSpace( { nPoints*nUniverses,       1} ));
    file->createDataSet<int>("best_grid_point", HighFive::DataSpace( { nPoints*nUniverses,       1} ));
    file->createDataSet<int>("n_iter",          HighFive::DataSpace( { nPoints*nUniverses,       1} ));
    // Some bookkeeping why not
    file->createDataSet<int>("i_grid",          HighFive::DataSpace( {nPoints*nUniverses,        1} ));
    file->createDataSet<int>("i_univ",          HighFive::DataSpace( {nPoints*nUniverses,        1} ));
    file->createDataSet<double>("gridx",        HighFive::DataSpace( {nPoints,                   1} ));
    file->createDataSet<double>("gridy",        HighFive::DataSpace( {nPoints,                   1} ));
}

void writeGrid(HighFive::File* file, std::vector<std::vector<double> > const & coords) {
    std::vector<double> xcoord;
    std::vector<double> ycoord;

    for (size_t i=0; i< coords.size(); i++) {
       xcoord.push_back(coords[i][0]);
       ycoord.push_back(coords[i][1]);
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
std::vector<std::vector<double>> mkHistoVecStd(std::string const & prefix, std::vector<string> const & fullnames) {
   // TODO can we have all spectra in a single file or something?
   std::string mass_tag = "";
   std::vector<std::vector<double> > temp;
   float dm(-2.0);
   while( dm < 2.3) {
      std::ostringstream out;
      out <<std::fixed<< std::setprecision(4) << dm;
      mass_tag = out.str();
      dm+=0.2;
      std::string fname = prefix+mass_tag+".SBNspec.root";
      temp.push_back(flattenHistos(readHistos(fname, fullnames)));
   }
   return temp;
}

std::vector<Eigen::VectorXd> mkHistoVec(std::string const & prefix, std::vector<string> const & fullnames) {
   // TODO can we have all spectra in a single file or something?
   std::string mass_tag = "";
   std::vector<std::vector<double> > temp;
   float dm(-2.0);
   while( dm < 2.3) {
      std::ostringstream out;
      out <<std::fixed<< std::setprecision(4) << dm;
      mass_tag = out.str();
      dm+=0.2;
      std::string fname = prefix+mass_tag+".SBNspec.root";
      temp.push_back(flattenHistos(readHistos(fname, fullnames)));
   }
   std::vector<Eigen::VectorXd> ret;
   for (size_t i=0;i<temp.size();++i) ret.push_back(Eigen::Map<Eigen::VectorXd> (temp[i].data(), temp[i].size(), 1) );
   return ret;
}

double calcChi(Eigen::VectorXd const & data, Eigen::VectorXd const & prediction, Eigen::MatrixXd const & C_inv ) {
   auto const & diff = data-prediction;
   return diff.transpose() * C_inv * diff;
}
double calcChi(Eigen::VectorXd const & diff, Eigen::MatrixXd const & C_inv ) {
   return diff.transpose() * C_inv * diff;
}

std::tuple<double, int> universeChi2(Eigen::VectorXd const & data, Eigen::MatrixXd const & C_inv,
   SignalGenerator signal)
{
   double chimin=std::numeric_limits<double>::infinity();
   int bestP(0);
   for (size_t i=0; i<signal.gridsize(); ++i) {
       double chi = calcChi(data - signal.predict(i, true), C_inv);
       if (chi<chimin) {
          chimin = chi;
          bestP=i;
       }
   }
   return {chimin, bestP};
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
Eigen::MatrixXd calcCovarianceMatrix(Eigen::MatrixXd const & M, Eigen::VectorXd const & spec){
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

// Cholesky decomposition and solve for inverted matrix --- fastest
Eigen::MatrixXd invertMatrixEigen3(Eigen::MatrixXd const & M){
    return M.llt().solve(Eigen::MatrixXd::Identity(M.rows(), M.rows()));
}

Eigen::MatrixXd collapseSubchannels(Eigen::MatrixXd const & EE, sbn::SBNconfig const & conf){
    Eigen::MatrixXd  retMat = Eigen::MatrixXd::Zero(conf.num_bins_detector_block_compressed, conf.num_bins_detector_block_compressed);

    int mrow(0), mcol(0), mrow_out(0), mcol_out(0);
    for(int ic = 0; ic < conf.num_channels; ic++) {
        for(int jc =0; jc < conf.num_channels; jc++) {
            for(int m=0; m < conf.num_subchannels[ic]; m++) {
                for(int n=0; n< conf.num_subchannels[jc]; n++) {
                   int a, c;
                   a=mrow + n*conf.num_bins[jc];
                   c=mcol + m*conf.num_bins[ic];
                   retMat.block(mrow_out, mcol_out, conf.num_bins[jc], conf.num_bins[ic]) += EE.block(a, c, conf.num_bins[jc], conf.num_bins[ic]);
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

Eigen::MatrixXd collapseDetectors(Eigen::MatrixXd const & M, sbn::SBNconfig const & conf){
    Eigen::MatrixXd  retMat = Eigen::MatrixXd::Zero(conf.num_bins_mode_block_compressed, conf.num_bins_mode_block_compressed);
    auto const & nrow = conf.num_bins_detector_block;
    auto const & crow = conf.num_bins_detector_block_compressed;
    for (int m=0; m<conf.num_detectors; m++) {
        for (int n=0; n<conf.num_detectors; n++) {
            retMat.block(n*crow, m*crow, crow, crow) = collapseSubchannels(M.block(n*nrow, m*nrow, nrow, nrow), conf);
        }
    }
    return retMat;
}


Eigen::MatrixXd updateInvCov(Eigen::MatrixXd const & covmat, Eigen::VectorXd const & spec_full, sbn::SBNconfig const & conf) {
    auto const & cov = calcCovarianceMatrix(covmat, spec_full);
    auto const & out = collapseDetectors(cov, conf);
    return invertMatrixEigen3(out);
}

FitResult coreFC(Eigen::VectorXd const & fake_data, Eigen::VectorXd const & v_coll,
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
       auto resuni  = universeChi2(fake_data, invcov, signal);
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

// TODO add size_t writeEvery to prevent memory overload
void doFC(Block* b, diy::Master::ProxyWithLink const& cp, int rank,
      const char * xmldata, sbn::SBNconfig const & myconf,
      TMatrixD const & covmat, Eigen::MatrixXd const & INVCOVBG,
      SignalGenerator signal,
      HighFive::File* file, std::vector<int> const & rankwork, int nUniverses, 
      double tol, size_t iter, bool debug, bool noWrite=false)
{

    double starttime, endtime;
    std::vector<FitResult> results;
    std::vector<int> v_grid, v_univ, v_iter, v_best;
    std::vector<double> v_last, v_dchi;
    
    if (!noWrite) {
       results.reserve(rankwork.size()*nUniverses);
       v_grid.reserve(rankwork.size()*nUniverses);
       v_univ.reserve(rankwork.size()*nUniverses);
       v_iter.reserve(rankwork.size()*nUniverses);
       v_best.reserve(rankwork.size()*nUniverses);
       v_last.reserve(rankwork.size()*nUniverses);
       v_dchi.reserve(rankwork.size()*nUniverses);
    }
    
    Eigen::Map<const Eigen::MatrixXd > ECOV(covmat.GetMatrixArray(), covmat.GetNrows(), covmat.GetNrows());

    for (int i_grid : rankwork) {

       if (debug && i_grid!=0) return;
       auto const & specfull_e = signal.predict(i_grid, false);
       //std::vector<double> specfull(specfull_e.data(), specfull_e.data() + specfull_e.rows() * specfull_e.cols());
       auto const & speccoll = collapseVectorEigen(specfull_e, myconf);
       //std::vector<double> speccoll = collapseVectorStd(specfull, myconf);
       sbn::SBNchi mychi(covmat, xmldata, false);
       mychi.InitRandomNumberSeeds(double(cp.gid()));

       starttime = MPI_Wtime();
       for (int uu=0; uu<nUniverses;++uu) {

          auto const & fake_data = mychi.SampleCovariance(specfull_e);
          auto const & fake_dataC = collapseVectorEigen(fake_data, myconf); 
          results.push_back(coreFC(fake_dataC, speccoll,
                   signal, INVCOVBG, ECOV, myconf, tol, iter));

          v_univ.push_back(uu);
          v_grid.push_back(i_grid);
       }
       endtime   = MPI_Wtime(); 
       if (rank==0) fmt::print(stderr, "[{}] gridp {}/{} ({} universes) took {} seconds\n",cp.gid(), i_grid, rankwork.size(), nUniverses, endtime-starttime);
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
    int mockFactor=1;
    size_t nUniverses=1;
    int NTEST(0);
    std::string out_file="test.hdf5";
    std::string f_BG="BOTHv2_BKG_ONLY.SBNspec.root";
    std::string f_CV="BOTHv2_CV.SBNspec.root";
    std::string f_COV="BOTHv2.SBNcovar.root";
    std::string tag="";
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
    ops >> Option("core",            f_CV,       "Central values filename.");
    ops >> Option('b', "background", f_BG,       "Backgrounds filename.");
    ops >> Option('c', "covmat",     f_COV,      "Covariance matrix filename.");
    ops >> Option("xmin",            xmin,       "xmin");
    ops >> Option("xmax",            xmax,       "xmax");
    ops >> Option("xwidth",          xwidth,     "xwidth");
    ops >> Option("ymin",            ymin,       "ymin");
    ops >> Option("ymax",            ymax,       "ymax");
    ops >> Option("ywidth",          ywidth,     "ywidth");
    ops >> Option("mock",            mockFactor, "Mockfactor");
    bool debug       = ops >> Present('d', "debug", "Operate on single gridpoint only");
    bool statonly    = ops >> Present("stat", "Statistical errors only");
    bool nowrite    = ops >> Present("nowrite", "Don't write output --- for performance estimates only");
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
       if (world.size() > nPoints) {
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
    std::vector<Eigen::VectorXd > sinsqvec_eig, sinvec_eig;
    int nFilesIn(0);
    if (world.rank()==0) {
       auto const & temp =mkHistoVecStd(tag+"_SINSQ_dm_", myconf.fullnames);
       nFilesIn = temp.size();
       sinsqvec = asVector(temp);
       sinvec   = asVector(mkHistoVecStd(tag+"_SIN_dm_"  , myconf.fullnames));
    }
    diy::mpi::broadcast(world, sinsqvec, 0);
    diy::mpi::broadcast(world, sinvec,   0);
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
    
    // Use the BG only inv cov matrix as start point
    TMatrixT<double> _cov = calcCovarianceMatrix(covmat, bgvec);
    Eigen::Map<const Eigen::MatrixXd> ecov(_cov.GetMatrixArray(), _cov.GetNcols(), _cov.GetNcols());
    auto const & _covcol = collapseDetectors(ecov, myconf);
    auto const & INVCOVBG = invertMatrixEigen3(_covcol);

    // Setup grid
    NGrid mygrid;
    mygrid.AddDimension("m4",  xmin, xmax, xwidth);//0.1 --- can't change easily --- sterile mass
    mygrid.AddDimension("ue4", ymin, ymax, ywidth);//0.1 --- arbirtrarily dense! mixing angle nu_e --- can change ywidth. Try 20^6
    mygrid.AddFixedDimension("um4", 0.0);
    nPoints = mygrid.f_num_total_points;
    GridPoints GP(mygrid.GetGrid());

    // Finally, the signal generator
    SignalGenerator    signal    = {myconf, GP, mygrid.f_dimensions[1].f_N, sinsqvec_eig, sinvec_eig, ecore};
    
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
       double t3 = MPI_Wtime();
       for (int i=0;i<NTEST;++i) collapseVectorStd(svb, myconf);
       double t4 = MPI_Wtime();

       fmt::print(stderr, "\n {} calls to predict took {} seconds\n",             NTEST, t1-t0);
       fmt::print(stderr, "\n {} calls to collapseVectorEigen took {} seconds\n", NTEST, t2-t1);
       fmt::print(stderr, "\n {} calls to collapseVectorStd took {} seconds\n",   NTEST, t4-t3);
       exit(1);
    }



    if( world.rank()==0 ) {
      fmt::print(stderr, "***********************************\n");
      fmt::print(stderr, "    Output will be written to {}\n", out_file);
      fmt::print(stderr, "    Points:  {}\n"                 ,  nPoints);
      fmt::print(stderr, "    nBins:  {}\n"                  , nBins);
      fmt::print(stderr, "    Universes: {}\n"               , nUniverses);
      fmt::print(stderr, "    Total size of dataset:  {}\n"  , nPoints*nUniverses);
      fmt::print(stderr, "    f_BG:  {}\n"                   , f_BG);
      fmt::print(stderr, "    f_CV:  {}\n"                   , f_CV);
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
    if (world.rank() == 0)  writeGrid(f_out, mygrid.GetGrid());

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

    std::vector<int> work(nPoints);
    std::iota(std::begin(work), std::end(work), 0); //0 is the starting number
    std::vector<int> rankwork = splitVector(work, world.size())[world.rank()];
    work.clear();
    work.shrink_to_fit();

    double starttime = MPI_Wtime();
    if (world.rank()==0) fmt::print(stderr, "Start FC\n");
    fc_master.foreach([world, covmat, xmldata, INVCOVBG, myconf, nUniverses, f_out, rankwork, tol, iter, debug, signal, nowrite](Block* b, const diy::Master::ProxyWithLink& cp)
                           { doFC(b, cp, world.rank(), xmldata, myconf, covmat, INVCOVBG, signal, f_out, rankwork, nUniverses, tol, iter, debug, nowrite); });
    double endtime   = MPI_Wtime(); 
    if (world.rank()==0) fmt::print(stderr, "[{}] FC took {} seconds\n",world.rank(), endtime-starttime);
    if (world.rank()==0) fmt::print(stderr, "Output written to {}\n",out_file);

    delete f_out;
    return 0;
}
