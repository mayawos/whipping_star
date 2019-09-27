#pragma GCC optimize("O3","unroll-loops","inline", "peel-loops")

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
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBNcovariance.h"
#include "SBNfeld.h"
#include <Eigen/Dense>
#include <Eigen/SVD>
//#include "nomad.hpp"


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

sbn::NeutrinoModel convert3p1(std::vector<double> const & ingrid){
    if(ingrid.size()!=3){
        std::cout<<"ERROR: Currently must have a 3-grid, for Dm ue4 and um4. grid.size(): "<<ingrid.size()<<std::endl;
        exit(EXIT_FAILURE);
    }
    sbn::NeutrinoModel signalModel(pow(10,ingrid[0]), pow(10,ingrid[1]), pow(10,ingrid[2]));
    return signalModel;
}

// TODO can this be simplified?
std::vector<double> collapseVector(std::vector<double> const & vin, sbn::SBNconfig const & conf){
   std::vector<double> cvec;
   cvec.reserve(conf.num_bins_total_compressed);
   //std::vector<double> cvec2(conf.num_bins_total_compressed, 0.0);

   for(int im = 0; im < conf.num_modes; im++){
       for(int id =0; id < conf.num_detectors; id++){
           int edge = id*conf.num_bins_detector_block + conf.num_bins_mode_block*im;
           for(int ic = 0; ic < conf.num_channels; ic++){
               int corner=edge;
               for(int j=0; j< conf.num_bins.at(ic); j++){
                   double tempval=0;
                   for(int sc = 0; sc < conf.num_subchannels.at(ic); sc++){
                        tempval += vin.at(j+sc*conf.num_bins.at(ic)+corner);
                        edge +=1;
                   }
                   cvec.push_back(tempval);
               }
           }
       }
   }
   //// All we want is a representation with the subchannels added together

   //std::cerr << " target size: " << cvec.size() << "\n";
   //for (int d=0; d<conf.num_detectors;++d) {
          //size_t offset(0);
          //size_t offset2(0);
      //for (int i=0; i<conf.num_channels; i++) {
          //size_t nbins_chan = conf.num_bins.at(i);
          //for (int j=0; j<conf.num_subchannels.at(i); j++) {
                //size_t first      = d*conf.num_bins_detector_block + offset;
                //size_t last       = d*conf.num_bins_detector_block + offset + nbins_chan;
                //size_t first_out  = d*conf.num_bins_detector_block_compressed + offset2;
                //size_t last_out   = d*conf.num_bins_detector_block_compressed + offset2 + nbins_chan;
                //std::cerr << "in : " << first << " : " << last << " out --- " << first_out << " : " << last_out << "\n";
                //std::transform (cvec2.begin() + first_out, cvec2.begin() + last_out, vin.begin() +first, cvec2.begin() + first_out, std::plus<double>());
                //offset +=nbins_chan;
          //}
                //offset2 +=nbins_chan;
       //}
   //}

   return cvec;
}

// TODO add mode logic if necessary
std::vector<double> collapseVectorStd(std::vector<double> const & vin, sbn::SBNconfig const & conf){
   // All we want is a representation with the subchannels added together
   std::vector<double> cvec(conf.num_bins_total_compressed, 0.0);
   for (int d=0; d<conf.num_detectors;++d) {
      size_t offset_in(0), offset_out(0);
      for (int i=0; i<conf.num_channels; i++) {
          size_t nbins_chan = conf.num_bins.at(i);
          for (int j=0; j<conf.num_subchannels.at(i); j++) {
             size_t first_in   = d*conf.num_bins_detector_block            + offset_in;
             size_t first_out  = d*conf.num_bins_detector_block_compressed + offset_out;
             std::transform (
                   cvec.begin() + first_out, cvec.begin() + first_out + nbins_chan, 
                   vin.begin()  + first_in,  cvec.begin() + first_out,
                   std::plus<double>());
             offset_in +=nbins_chan;
          }
          offset_out += nbins_chan;
      }
   }
   return cvec;
}


struct SignalGenerator {
   sbn::SBNosc osc; std::vector<std::vector<double>> const & m_vec_grid; std::string tag;
   std::unordered_map <std::string, std::vector<TH1D> > const & sinsqmap;
   std::unordered_map <std::string, std::vector<TH1D> > const & sinmap;
   const char * xmldata;
   
   std::vector<double> predict(size_t i_grid, bool compressed) {
      sbn::NeutrinoModel this_model = convert3p1(m_vec_grid[i_grid]);
      osc.LoadModel(this_model);
      osc.SetAppMode();            
      std::vector<double> ans = osc.Oscillate(tag, compressed, xmldata, sinsqmap, sinmap);
      return ans;
   }
   size_t gridsize() {return m_vec_grid.size();}

};

struct SignalGeneratorStd {
   sbn::SBNosc osc; sbn::SBNconfig const & conf; std::vector<std::vector<double>> const & m_vec_grid;
   std::unordered_map <std::string, std::vector<double> > const & sinsqmap;
   std::unordered_map <std::string, std::vector<double> > const & sinmap;
   std::vector<double> const core;
   
   std::vector<double> predict(size_t i_grid, bool compressed) {
      sbn::NeutrinoModel this_model = convert3p1(m_vec_grid[i_grid]);
      osc.LoadModel(this_model);
      osc.SetAppMode();            
      std::vector<double> ans = osc.Oscillate(sinsqmap, sinmap);
      std::transform (ans.begin(), ans.end(), core.begin(), ans.begin(), std::plus<double>());
      if (compressed) {
         return collapseVectorStd(ans, conf);
      }
      else {
         return ans;
      }
   }
   size_t gridsize() {return m_vec_grid.size();}

};

struct SignalGeneratorEigen {
   sbn::SBNosc osc; std::vector<std::vector<double>> const & m_vec_grid;
   std::unordered_map <std::string, Eigen::VectorXd > const & sinsqmap;
   std::unordered_map <std::string, Eigen::VectorXd > const & sinmap;
   
   std::vector<double> predict(size_t i_grid, bool compressed) {
      sbn::NeutrinoModel this_model = convert3p1(m_vec_grid[i_grid]);
      osc.LoadModel(this_model);
      osc.SetAppMode();            
      std::vector<double> ans = osc.Oscillate(sinsqmap, sinmap);
      return ans;
   }
   size_t gridsize() {return m_vec_grid.size();}

};


// Split input vector into pieces of similar size
std::vector<std::vector<int>> splitVector(const std::vector<int>& vec, size_t n) {
    std::vector<std::vector<int>> outVec;
    size_t length = vec.size() / n;
    size_t remain = vec.size() % n;
    size_t begin = 0;
    size_t end = 0;

    for (size_t i = 0; i < std::min(n, vec.size()); ++i) {
        end += (remain > 0) ? (length + !!(remain--)) : length;
        outVec.push_back(std::vector<int>(vec.begin() + begin, vec.begin() + end));
        begin = end;
    }
    return outVec;
}

// Concatenate vector of vectors as vector
std::vector<double> asVector(std::vector<std::vector<double> > v_in) {
    std::vector<double> allSpectra;
    allSpectra.reserve(v_in.size()*v_in[0].size());
    for (auto temp : v_in) allSpectra.insert(allSpectra.end(), temp.begin(), temp.end());
    return allSpectra;
}

std::vector<double> asVector(TMatrixT<double> const & M) {
    const double *pData = M.GetMatrixArray();
    std::vector<double> vData;
    vData.assign(pData, pData + M.GetNoElements());
    return vData;
}

template<typename T>
std::vector<T> myslice(std::vector<T> const &v, int m, int n)
{
    auto first = v.cbegin() + m;
    auto last = v.cbegin() + n + 1;
    std::vector<T> vec(first, last);
    return vec;
}

std::vector< std::vector<double> > sliceVector(std::vector<double> const & input, int nPoints) {
   std::vector< std::vector<double> > test;
   test.reserve(nPoints);
   int nBins = input.size()/nPoints;

   std::vector<double> work;
   work.reserve(nBins);

   for (int i=0;i<nPoints;++i) {
      work=myslice(input, i*nBins, (i+1)*nBins-1);
      test.push_back(work);
   }
   return test;
}

// Properly release memory
void releaseVec(std::vector<double> & vec) {
    vec.clear();
    vec.shrink_to_fit();
}

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
    return cov;
}


std::vector<TH1D> readHistos(std::string const & rootfile, std::vector<string> const & fullnames) {
    std::vector<TH1D > hist;
    TFile f(rootfile.c_str(),"read");
    for (auto fn: fullnames) {
       //std::cerr << "reading " <<fn.c_str() << "\n";
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

std::unordered_map <std::string, std::vector<TH1D> > mkHistoMap(std::string const & prefix, std::vector<string> const & fullnames) {
   // TODO can we have all spectra in a single file or something?
   std::string mass_tag = "";
   std::unordered_map <std::string, std::vector<TH1D> > hmap;
   float dm(-2.0);
   while( dm < 2.3) {
      std::ostringstream out;
      out <<std::fixed<< std::setprecision(4) << dm;
      mass_tag = out.str();
      dm+=0.2;
      std::string fname = prefix+mass_tag+".SBNspec.root";
      hmap[mass_tag] = readHistos(fname, fullnames);
   }
   return hmap;
}

std::unordered_map <std::string, std::vector<double> > mkHistoMapStd(std::string const & prefix, std::vector<string> const & fullnames) {
   // TODO can we have all spectra in a single file or something?
   std::string mass_tag = "";
   std::unordered_map <std::string, std::vector<double> > hmap;
   float dm(-2.0);
   while( dm < 2.3) {
      std::ostringstream out;
      out <<std::fixed<< std::setprecision(4) << dm;
      mass_tag = out.str();
      dm+=0.2;
      std::string fname = prefix+mass_tag+".SBNspec.root";
      std::vector<TH1D> histos = readHistos(fname, fullnames);
      hmap[mass_tag] = flattenHistos(histos);
   }
   return hmap;
}

std::unordered_map <std::string, Eigen::VectorXd > mkHistoMapEig(std::string const & prefix, std::vector<string> const & fullnames) {
   // TODO can we have all spectra in a single file or something?
   std::string mass_tag = "";
   std::unordered_map <std::string, Eigen::VectorXd > hmap;
   float dm(-2.0);
   while( dm < 2.3) {
      std::ostringstream out;
      out <<std::fixed<< std::setprecision(4) << dm;
      mass_tag = out.str();
      dm+=0.2;
      std::string fname = prefix+mass_tag+".SBNspec.root";
      std::vector<TH1D> histos = readHistos(fname, fullnames);
      std::vector<double> ret = flattenHistos(histos);
      hmap[mass_tag] = Eigen::Map<Eigen::VectorXd>(ret.data(), ret.size(), 1);
   }
   return hmap;
}


float calcChi(std::vector<double> const & data, std::vector<double> const & prediction, TMatrixT<double> const & C_inv ) {
    int n = data.size();
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > CINV(C_inv.GetMatrixArray(), n, n);
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 1> > PRED(prediction.data(), n, 1);
    // Massage data right now to convert from float to double.
    std::vector<double> temp;
    temp.assign(data.begin(), data.end());
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 1> > DATA(temp.data(), n, 1);
    Eigen::VectorXd DIFF = DATA-PRED;
    return DIFF.transpose() * CINV * DIFF; 
}

inline size_t grid_to_index(size_t ix, size_t iy) {
   // Total grid size static for now ( number of y-values in this 2D case)
   size_t NY(48);
   return ix*NY + iy;
}

// Calculate all chi2 for this universe
std::tuple<double, int> universeChi2(std::vector<double> const & data, TMatrixT<double> const & C_inv,
   SignalGeneratorStd signal)
{
   std::vector<double> temp;
   std::vector<double> result;
   result.reserve(signal.gridsize());
   double chimin=std::numeric_limits<double>::infinity();
   int bestP(0);

   for (size_t i=0; i<signal.gridsize(); ++i) {
       temp = signal.predict(i, true);
       double chi = calcChi(data, temp, C_inv);
       if (chi<chimin) {
          chimin = chi;
          bestP=i;
       }
   }
   return {chimin, bestP};
}

TMatrixT<double> calcCovarianceMatrix(TMatrixT<double> const & M, std::vector<double> const & spec){
    TMatrixT<double> Mout( M.GetNcols(), M.GetNcols() );
    // systematics per scaled event
    for(int i =0; i<M.GetNcols(); i++) {
        for(int j =0; j<M.GetNrows(); j++) {
            if ( std::isnan( M(i,j) ) )  Mout(i,j) = 0.0;
            else                         Mout(i,j) = M(i,j)*spec[i]*spec[j];

            if (i==j) Mout(i,i) += spec[i];
        }
    }
    return Mout;
}

// Only used for covariance matrix a very beginning, just being lazy here
TMatrixT<double> invertMatrix(TMatrixT<double> &M){
    TMatrixT<double> McI(M.GetNrows(),M.GetNrows());
    McI.Zero();
    TDecompSVD svd(M);
    if (!svd.Decompose()) {
        std::cerr<<" (InvertMatrix) Decomposition FAILED, matrix not symettric?, has nans?" << std::endl;
        std::cerr<<"ERROR: The matrix to invert failed a SVD decomp!"<<std::endl;
        for(int i=0; i< M.GetNrows(); i++){
            for(int j=0; j< M.GetNrows(); j++) std::cerr<<i<<" "<<j<<" "<<M(i,j)<<std::endl;
        }
    }
    else McI = svd.Invert();
    if (!McI.IsValid()) std::cout << "ERROR: The inverted matrix isnt valid! Something went wrong.." << std::endl;
    return McI;
}


// Cholesky decomposition and solve for inverted matrix --- fastest
Eigen::MatrixXd invertMatrixEigen3(TMatrixT<double> &M){
    int n = M.GetNrows();
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > COV(M.GetMatrixArray(), n, n);
    return COV.llt().solve(Eigen::MatrixXd::Identity(n,n));
}


// This is the powerhouse, takes each detector matrix filled with num_channels channels of num_subchannels[i] subchannels, and collapses it.
void collapseSubchannels(TMatrixT <double> & M, TMatrixT <double> & Mc, sbn::SBNconfig const & conf){

    std::vector< std::vector< TMatrixT<double> > > Summed(conf.num_channels, std::vector<TMatrixT<double>>(conf.num_channels) );	//Initialise a matrix of matricies, to ZERO.
    for(int ic = 0; ic < conf.num_channels; ic++){
        for(int jc =0; jc < conf.num_channels; jc++){
            Summed[ic][jc].ResizeTo(conf.num_bins[jc], conf.num_bins[ic]) ;// This is CORRECT, do not switch (ie Summed[0][1] = size (num_bins[1], num_bins[0])
            Summed[ic][jc] = 0.0;
        }
    }

    int mrow = 0.0;
    int mcol = 0.0;

    for(int ic = 0; ic < conf.num_channels; ic++){ 	 // Loop over all rows
        for(int jc =0; jc < conf.num_channels; jc++){    // Loop over all columns
            for(int m=0; m < conf.num_subchannels[ic]; m++){
                for(int n=0; n< conf.num_subchannels[jc]; n++){ //For each big block, loop over all subchannels summing together
                   int a, b, c, d;
                   a=mrow + n*conf.num_bins[jc];
                   b=mrow + n*conf.num_bins[jc] + conf.num_bins[jc]-1;
                   c=mcol + m*conf.num_bins[ic];
                   d=mcol + m*conf.num_bins[ic] + conf.num_bins[ic]-1;
                   Summed[ic][jc] +=  M.GetSub(a,b,c,d);
                }
            }
            mrow += conf.num_subchannels[jc]*conf.num_bins[jc]; // As we work our way left in columns, add on that many bins
        } // end of column loop
        mrow = 0; // as we end this row, reSet row count, but jump down 1 column
        mcol += conf.num_subchannels[ic]*conf.num_bins[ic];
    } // end of row loop

    ///********************************* And put them back toGether! ************************//
    Mc.Zero();
    mrow = 0;
    mcol = 0;

    //Repeat again for Contracted matrix
    for(int ic = 0; ic < conf.num_channels; ic++){
        for(int jc =0; jc < conf.num_channels; jc++){
            Mc.SetSub(mrow, mcol, Summed[ic][jc]);
            mrow += conf.num_bins[jc];
        }
        mrow = 0;
        mcol +=conf.num_bins[ic];
    }
}

void collapseDetectors(TMatrixT <double> & M, TMatrixT <double> & Mc, sbn::SBNconfig const & conf){
    Mc.Zero();
    int nrow = conf.num_bins_detector_block;// N_e_bins*N_e_spectra+N_m_bins*N_m_spectra;
    int crow = conf.num_bins_detector_block_compressed; //N_e_bins+N_m_bins;
    for (int m=0; m<conf.num_detectors; m++) {
        for (int n=0; n<conf.num_detectors; n++) {
            TMatrixT<double> imat(nrow,nrow);
            TMatrixT<double> imatc(crow,crow);
            imat = M.GetSub(n*nrow, n*nrow+nrow-1, m*nrow, m*nrow+nrow-1);
            collapseSubchannels(imat, imatc, conf);
            Mc.SetSub(n*crow, m*crow, imatc);
        }
    }
}

void collapseModes(TMatrixT <double> & M, TMatrixT <double> & Mc, sbn::SBNconfig const & conf){
    Mc.Zero();
    int nrow = conf.num_bins_mode_block;// (N_e_bins*N_e_spectra+N_m_bins*N_m_spectra)*N_dets;
    int crow = conf.num_bins_mode_block_compressed;// (N_e_bins+N_m_bins)*N_dets;

    for (int m=0; m<conf.num_modes ; m++) {
        for (int n=0; n<conf.num_modes; n++) {
            TMatrixT<double> imat(nrow, nrow);
            TMatrixT<double> imatc(crow, crow);
            imat = M.GetSub(n*nrow, n*nrow+nrow-1, m*nrow, m*nrow+nrow-1);
            collapseDetectors(imat, imatc, conf);
            Mc.SetSub(n*crow, m*crow, imatc);
        }
    }
}

void updateInvCov(TMatrixT<double> & invcov, TMatrixD const & covmat, std::vector<double> const & spec_full, sbn::SBNconfig const & conf, int nBinsC) {
    TMatrixT<double> _cov = calcCovarianceMatrix(covmat, spec_full);
    TMatrixT<double> _covcol;
    _covcol.ResizeTo(nBinsC, nBinsC);
    collapseModes(_cov, _covcol, conf);
    Eigen::MatrixXd temp = invertMatrixEigen3(_covcol);
    invcov.SetMatrixArray(temp.data());
}


FitResult coreFC(std::vector<double> const & fake_data, std::vector<double> const & v_coll,
      SignalGeneratorStd signal,
      TMatrixT<double> const & INVCOV, 
      TMatrixT<double> const & covmat, 
      sbn::SBNconfig const & myconf,
      double chi_min_convergance_tolerance = 0.001,
      size_t max_number_iterations = 5)
{
   float last_chi_min = FLT_MAX;
   int best_grid_point = -99;
   size_t n_iter = 0;

   TMatrixT<double> invcov = INVCOV;
   std::vector<double> temp;
   
   for(n_iter = 0; n_iter < max_number_iterations; n_iter++){
       if(n_iter!=0){
           //Calculate current full covariance matrix, collapse it, then Invert.
           temp  = signal.predict(best_grid_point, false);
           updateInvCov(invcov, covmat, temp, myconf, 90);
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

void doFC(Block* b, diy::Master::ProxyWithLink const& cp, int rank,
      NGrid mygrid, const char * xmldata, sbn::SBNconfig const & myconf,
      TMatrixD const & covmat, TMatrixT<double> const & INVCOV,
      SignalGeneratorStd signal,
      HighFive::File* file, std::vector<int> const & rankwork, int nUniverses, 
      double tol, size_t iter, bool debug)
{

    double starttime, endtime;
    std::vector<double> fake_data;
    std::vector<double> fake_dataC;
    std::vector<FitResult> results;
    results.reserve(rankwork.size()*nUniverses);
   
    std::vector<int> v_grid, v_univ, v_iter, v_best;
    v_grid.reserve(rankwork.size()*nUniverses);
    v_univ.reserve(rankwork.size()*nUniverses);
    v_iter.reserve(rankwork.size()*nUniverses);
    v_best.reserve(rankwork.size()*nUniverses);
    std::vector<double> v_last, v_dchi;
    v_last.reserve(rankwork.size()*nUniverses);
    v_dchi.reserve(rankwork.size()*nUniverses);
    

    for (int i_grid : rankwork) {

       if (debug && i_grid!=0) return;
       std::vector<double> specfull = signal.predict(i_grid, false);
       std::vector<double> speccoll = collapseVectorStd(specfull, myconf); //signal.predict(i_grid, true);
       sbn::SBNchi mychi(covmat, xmldata, false);
       mychi.InitRandomNumberSeeds(double(cp.gid()));

       starttime = MPI_Wtime();
       for (int uu=0; uu<nUniverses;++uu) {
          fake_data = mychi.SampleCovariance(specfull);
          fake_dataC = collapseVectorStd(fake_data, myconf);
          results.push_back(coreFC(fake_dataC, speccoll, signal, INVCOV, covmat, myconf, tol, iter));
          v_univ.push_back(uu);
          v_grid.push_back(i_grid);
       }
       endtime   = MPI_Wtime(); 
       if (rank==0) fmt::print(stderr, "[{}] gridp {}/{} ({} universes) took {} seconds\n",cp.gid(), i_grid, rankwork.size(), nUniverses, endtime-starttime);
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

inline bool file_exists (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

// --- main program ---//
int main(int argc, char* argv[]) {
    diy::mpi::environment env(argc, argv);
    diy::mpi::communicator world;

    size_t nBinsC=90;
    int nPoints=-1;
    int mockFactor=1;
    int nUniverses=1;
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
    ops >> Option('x', "xml",        xml,        "XML config.");
    ops >> Option("tol",             tol,        "Minimiser tolerance");
    ops >> Option("iter",            iter,       "Max number of iterations.");
    ops >> Option('t', "tag",        tag,        "Tag.");
    ops >> Option("nbinsC",          nBinsC,     "Number of collapsed bins in 2d dataset");
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
    //bool storeINVCOV = ops >> Present("storeinvcov", "Write inverted covariance matrices to hdf5");
    if (ops >> Present('h', "help", "Show help")) {
        std::cout << "Usage:  [OPTIONS]\n";
        std::cout << ops;
        return 1;
    }
    
    if( world.rank()==0 ) {
      fmt::print(stderr, "\n*** This is diy running highfivewrite ***\n");
    }

    NGrid mygrid;
    mygrid.AddDimension("m4",  xmin, xmax, xwidth);//0.1 --- can't change easily --- sterile mass
    mygrid.AddDimension("ue4", ymin, ymax, ywidth);//0.1 --- arbirtrarily dense! mixing angle nu_e --- can change ywidth. Try 20^6
    mygrid.AddFixedDimension("um4", 0.0); //0.05

    nPoints = mygrid.f_num_total_points;


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

    const char* xmldata = text.c_str();
    sbn::SBNconfig myconf(xmldata, false);

    // Central values FIXME don't read on every rank!
    std::vector<TH1D> cvhist = readHistos(f_CV, myconf.fullnames);
    sbn::SBNosc myosc(cvhist, xmldata);
    
    // TODO find a way to serialize and reconstruct TH1D that is not too painful
    double time0 = MPI_Wtime();
    std::unordered_map <std::string, std::vector<TH1D> >   sinsqmap     = mkHistoMap(   tag+"_SINSQ_dm_", myconf.fullnames);
    std::unordered_map <std::string, std::vector<TH1D> >   sinmap       = mkHistoMap(   tag+"_SIN_dm_"  , myconf.fullnames);
    std::unordered_map <std::string, std::vector<double> > sinsqmap_std = mkHistoMapStd(tag+"_SINSQ_dm_", myconf.fullnames);
    std::unordered_map <std::string, std::vector<double> > sinmap_std   = mkHistoMapStd(tag+"_SIN_dm_"  , myconf.fullnames);
    //std::unordered_map <std::string, Eigen::VectorXd >     sinsqmap_eig = mkHistoMapEig(tag+"_SINSQ_dm_", myconf.fullnames);
    //std::unordered_map <std::string, Eigen::VectorXd >     sinmap_eig   = mkHistoMapEig(tag+"_SIN_dm_"  , myconf.fullnames);
    double time1 = MPI_Wtime();
    if (world.rank()==0) fmt::print(stderr, "[{}] loading spectra the new way took {} seconds\n",world.rank(), time1 - time0);

    std::vector<double> const core  = flattenHistos(cvhist);

    SignalGenerator signal          = {myosc, mygrid.GetGrid(), tag, sinsqmap,     sinmap,     xmldata};
    SignalGeneratorStd signal_std   = {myosc, myconf, mygrid.GetGrid(), sinsqmap_std, sinmap_std, core};
    //SignalGeneratorEigen signal_eig = {myosc, mygrid.GetGrid(), sinsqmap_eig, sinmap_eig};

    double time3 = MPI_Wtime();
    for (int i =0;i<1000;++i) signal.predict(1,false);
    double time4 = MPI_Wtime();
    for (int i =0;i<1000;++i) signal_std.predict(1,false);
    double time5 = MPI_Wtime();

    
    for (int i =0;i<1000;++i) signal.predict(1,true);
    double time6 = MPI_Wtime();
    for (int i =0;i<1000;++i) signal_std.predict(1,true);
    double time7 = MPI_Wtime();

    if (world.rank()==0) fmt::print(stderr, "[{}] TH1D way took {} seconds\n",world.rank(), time4 - time3);
    if (world.rank()==0) fmt::print(stderr, "[{}] std  way took {} seconds\n",world.rank(), time5 - time4);
    if (world.rank()==0) fmt::print(stderr, "[{}] std speed-up {} \n",world.rank(), (time4 - time3)/(time5 - time4));
  
    if (world.rank()==0) fmt::print(stderr, "[{}] C TH1D way took {} seconds\n",world.rank(), time6 - time5);
    if (world.rank()==0) fmt::print(stderr, "[{}] C std  way took {} seconds\n",world.rank(), time7 - time6);
    if (world.rank()==0) fmt::print(stderr, "[{}] C std speed-up {} \n",world.rank(), (time6 - time5)/(time7 - time6));


    std::vector<double> test =signal_std.predict(1,false);
    time3 = MPI_Wtime();
    for (int i =0;i<10000;++i) collapseVector(test, myconf);
    time4 = MPI_Wtime();
    for (int i =0;i<10000;++i) collapseVectorStd(test, myconf);
    time5 = MPI_Wtime();
    
    if (world.rank()==0) fmt::print(stderr, "[{}] old way took {} seconds\n",world.rank(), time4 - time3);
    if (world.rank()==0) fmt::print(stderr, "[{}] std  way took {} seconds\n",world.rank(), time5 - time4);
    if (world.rank()==0) fmt::print(stderr, "[{}] std speed-up {} \n",world.rank(), (time4 - time3)/(time5 - time4));

    //std::vector<double> a =collapseVector(   test, myconf);
    //std::vector<double> b =collapseVectorStd(test, myconf);



    //for (int i =0;i<1;++i) signal_eig.predict(1,false);
    //double time6 = MPI_Wtime();
    //if (world.rank()==0) fmt::print(stderr, "[{}] TH1D way took {} seconds\n",world.rank(), time4 - time3);
    //if (world.rank()==0) fmt::print(stderr, "[{}] std  way took {} seconds\n",world.rank(), time5 - time4);
    //if (world.rank()==0) fmt::print(stderr, "[{}] eig  way took {} seconds\n",world.rank(), time6 - time5);
    //if (world.rank()==0) fmt::print(stderr, "[{}] std speed-up {} \n",world.rank(), (time4 - time3)/(time5 - time4));
    //if (world.rank()==0) fmt::print(stderr, "[{}] eig speed-up {} \n",world.rank(), (time4 - time3)/(time6 - time5));
    
    
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

    if( world.rank()==0 ) {
      fmt::print(stderr, "\n    Output will be written to {}\n", out_file);
      fmt::print(stderr, "\n    Points:  {}\n", nPoints);
      fmt::print(stderr, "\n    nBins:  {}\n", nBins);
      fmt::print(stderr, "\n    Universes: {}\n", nUniverses);
      fmt::print(stderr, "\n    f_BG:  {}\n", f_BG);
      fmt::print(stderr, "\n    f_CV:  {}\n", f_CV);
      fmt::print(stderr, "\n    f_COV:  {}\n", f_COV);
      fmt::print(stderr, "\n    iter :  {}\n", iter);
      fmt::print(stderr, "\n    tol:    {}\n", tol);
      if (statonly) fmt::print(stderr, "\n    S T A T  O N L Y \n");
      if (debug) fmt::print(stderr, "\n    D E B U G \n");
      fmt::print(stderr, "\n    Total size of dataset:  {}\n", nPoints*nUniverses);
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

    // First diy setup to load spectra
    size_t blocks = nPoints;
    Bounds domain(1);
    domain.min[0] = 0;
    domain.max[0] = blocks-1;
    diy::RoundRobinAssigner        assigner(world.size(), blocks);
    diy::RegularDecomposer<Bounds> decomposer(1, domain, blocks);
    diy::RegularBroadcastPartners  comm(    decomposer, 2, true);
    diy::RegularMergePartners      partners(decomposer, 2, true);
    diy::Master master(world, 1, -1, &Block::create, &Block::destroy);
    diy::decompose(1, world.rank(), domain, assigner, master);
    

    // Add the BG
    TMatrixT<double> _cov = calcCovarianceMatrix(covmat, bgvec);
    TMatrixT<double> _covcol;
    _covcol.ResizeTo(nBinsC, nBinsC);
    collapseModes(_cov, _covcol, myconf);
    TMatrixT<double> INVCOV = invertMatrix(_covcol);


    //// Now more blocks as we have universes
    blocks = world.size();//nPoints;//*nUniverses;
    if (world.rank()==0) fmt::print(stderr, "FC will be done on {} blocks, distributed over {} ranks\n", blocks, world.size());
    Bounds fc_domain(1);
    fc_domain.min[0] = 0;
    fc_domain.max[0] = blocks-1;
    diy::RoundRobinAssigner        fc_assigner(world.size(), blocks);
    diy::RegularDecomposer<Bounds> fc_decomposer(1, fc_domain, blocks);
    diy::RegularBroadcastPartners  fc_comm(    fc_decomposer, 2, true);
    diy::RegularMergePartners      fc_partners(fc_decomposer, 2, true);
    diy::Master                    fc_master(world, 1, -1, &Block::create, &Block::destroy);
    diy::decompose(1, world.rank(), fc_domain, fc_assigner, fc_master);//, share_face, wrap, ghosts);


    std::vector<int> work(nPoints);
    std::iota(std::begin(work), std::end(work), 0); //0 is the starting number
    std::vector<int> rankwork = splitVector(work, world.size())[world.rank()];

    double starttime = MPI_Wtime();
    fc_master.foreach([world, mygrid, covmat, xmldata, INVCOV, myconf, nUniverses, f_out, rankwork, tol, iter, debug, signal_std](Block* b, const diy::Master::ProxyWithLink& cp)
                           { doFC(b, cp, world.rank(), mygrid, xmldata, myconf, covmat, INVCOV, signal_std, f_out, rankwork, nUniverses, tol, iter, debug); });
    double endtime   = MPI_Wtime(); 
    if (world.rank()==0) fmt::print(stderr, "[{}] FC took {} seconds\n",world.rank(), endtime-starttime);


    delete f_out;
    return 0;
}
