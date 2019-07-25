
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include <functional>
#include <numeric>
#include <stdexcept>

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
#include "SBNfit.h"
#include "SBNcovariance.h"
#include "SBNfeld.h"

using namespace std;

#include "opts.h"

typedef diy::DiscreteBounds Bounds;

struct Block
{
    static void*    create()            { return new Block; }
    static void     destroy(void* b)    { delete static_cast<Block*>(b); }

    void show_link(const diy::Master::ProxyWithLink& cp)
    {
      diy::RegularLink<Bounds>* link = static_cast<diy::RegularLink<Bounds>*>(cp.link());
      std::cout << "Block (" << cp.gid() << "): "
                << link->core().min[0]   << ' ' << link->core().min[1]   << ' ' << link->core().min[2] << " - "
                << link->core().max[0]   << ' ' << link->core().max[1]   << ' ' << link->core().max[2] << " : "
                << link->bounds().min[0] << ' ' << link->bounds().min[1] << ' ' << link->bounds().min[2] << " - "
                << link->bounds().max[0] << ' ' << link->bounds().max[1] << ' ' << link->bounds().max[2] << " : "
                << link->size()   << ' ' //<< std::endl
                << std::dec
                << std::endl;
    }
};

// TODO make this work with std::vector
std::vector<TH1D> readHistos(std::string const & rootfile, std::vector<string> const & fullnames) {
     std::vector<TH1D > hist;
     //Contruct from a prexisting histograms that exist in a rootfile
     TFile f(rootfile.c_str(),"read");

     //Loop over all filenames that should be there, and load up the histograms.
     for (auto fn: fullnames) hist.push_back(*((TH1D*)f.Get(fn.c_str())));
     f.Close();

     return hist;
}

// Assume order is correct
std::vector<double> flattenHistos(std::vector<TH1D> const & v_hist) {
   std::vector<double> temp;
   for (auto hist : v_hist) {
      for (size_t i=1; i< (hist.GetNbinsX()+1); ++i) {
         temp.push_back(hist.GetBinContent(i));
      }
   }
   return temp;
}

// That's the chi2 calculation
float calcChi(std::vector<float> const & data, std::vector<double> const & prediction, TMatrixT<double> const & C_inv ){
    float tchi = 0;

    for(int i =0; i<data.size(); i++){
        for(int j =0; j<data.size(); j++)  tchi += (data[i]-prediction[i]) * C_inv[i][j] * (data[j]-prediction[j]);
    }
    return tchi;
}


TMatrixD readFracCovMat(std::string const & rootfile){
    TFile  fsys(rootfile.c_str(),"read");
    TMatrixD cov =  *(TMatrixD*)fsys.Get("frac_covariance");
    fsys.Close();
    return cov;
}


unique_ptr<sbn::SBNspec> loadPreOscillatedSpectrum(int i_uni, string const & tag, string const & xml, NGrid const & mygrid, const char * xmldata, std::vector<TH1D> const & bghist,  std::vector<TH1D> const & cvhist, TMatrixD const & covmat) {

    sbn::SBNfeld myfeld(mygrid, tag, xmldata);

    myfeld.SetCoreSpectrum(cvhist, xmldata);

    double random_number_seed = -1;
    //std::cout<<"Setting random seed "<<random_number_seed<<std::endl;
    myfeld.SetRandomSeed(random_number_seed);


    // FIXME the dtor of SBNfeld only works when these functions are called :(
    myfeld.SetFractionalCovarianceMatrix(covmat);
    myfeld.LoadBackgroundSpectrum(xmldata, bghist);
    return myfeld.LoadPreOscillatedSpectrum(i_uni, xmldata);
};

void loadSpectrum(Block* b, diy::Master::ProxyWithLink const& cp, int size, int rank, string tag, string xml, NGrid const & mygrid, HighFive::File* f_out, const char * xmldata, std::vector<TH1D> const & bghist,  std::vector<TH1D> const & cvhist,TMatrixD const & covmat) {
    auto myspec = loadPreOscillatedSpectrum(cp.gid(), tag, xml, mygrid, xmldata, bghist, cvhist, covmat);
    std::vector<double> full = myspec->full_vector;
    std::vector<double> coll = myspec->collapsed_vector;

    HighFive::DataSet spec = f_out->getDataSet("specfull");
    spec.select(   {std::size_t(cp.gid()), 0}, {1, std::size_t(full.size())}).write(full);

    HighFive::DataSet speccoll = f_out->getDataSet("speccoll");
    speccoll.select(   {std::size_t(cp.gid()), 0}, {1, std::size_t(coll.size())}).write(coll);
 
}


// This is doing one univere for one gridpoint
void doFCsmart(Block* b, diy::Master::ProxyWithLink const& cp, int size, int rank, string tag, string xml, NGrid mygrid, HighFive::File* f_out, int num_universes, bool dry, const char * xmldata, std::vector<TH1D> const & bghist, std::vector<TH1D> const & cvhist, TMatrixD const & covmat, std::vector<double> const & bgvec) {
    // Some arithmetic to figure out the gridpoint and universe from cp.gid
    int i_grid = cp.gid() % mygrid.f_num_total_points;
    int i_univ = floor(cp.gid()/mygrid.f_num_total_points);


    int nBinsFull =covmat.GetNcols();
    // These HDF5 dataset have all the full/collapsed vector spectra as loaded in loadSpectrum
    HighFive::DataSet g_specfull = f_out->getDataSet("specfull");
    HighFive::DataSet g_speccoll = f_out->getDataSet("speccoll");
    std::vector<double> specfull = {nBinsFull};
    g_specfull.select(   {i_grid, 0}, {1, nBinsFull}).read(specfull);


    double sumf = std::accumulate(specfull.begin(), specfull.end(), 0.0);
    if (sumf >1e10) {
       std::cerr << "sum of specfull: " << sumf << "\n";
       for (auto v : specfull ) std::cerr << " " << v;
       std::cerr << "\n";
       //abort();
    }

    sbn::SBNspec myspec(specfull, xmldata, i_grid, false);
    
    double starttime, endtime;
    starttime = MPI_Wtime();
    //sbn::SBNfeld myfeld(mygrid, tag, xmldata);
    //myfeld.SetCoreSpectrum(cvhist, xmldata);
    //myfeld.SetFractionalCovarianceMatrix(covmat);

    //// TODO: Setting the seed to cp.gid()???
    //double random_number_seed = -1;
    //myfeld.SetRandomSeed(random_number_seed);
    
    //auto myspec = myfeld.LoadPreOscillatedSpectrum(i_grid, xmldata); // FIXME We need to get around the second loading here.

    //myfeld.LoadBackgroundSpectrum(xmldata, bghist);
    TMatrixT<double> stat_only_matrix(nBinsFull, nBinsFull);
    stat_only_matrix.Zero();

    sbn::SBNchi* mychi = nullptr; 

    // TODO make statonly a CL option and pass
    //if(myfeld.statOnly()){
         //mychi = new sbn::SBNchi(myspec, stat_only_matrix, xmldata, false);//, myfeld.seed());
     //}else{
         mychi = new sbn::SBNchi(myspec, covmat, xmldata, false);//, myfeld.seed()); FIXME pass covmat directly???
     //}
    int nBinsColl = myspec.collapsed_vector.size();

    // Single grid point multiuniverse FC
    int max_number_iterations = 5;
    double chi_min_convergance_tolerance = 0.001;

    // FIXME * CalcCovMat --- make free function
    //       * Understand CollapseModes
    //       * ReloadCoreSpectrum --- why?
    //       * InvertMatrix --- make free function and DONT EXIT IN LIBRARY!!!
    //Ok take the background only spectrum and form a background only covariance matrix. CalcCovarianceMatrix includes stats
    //std::vector<double> test = flattenHistos(bghist);
    //std::cerr << "bgvec size: " << bgvec.size() << "\n";
    //for (auto b : bgvec ) {
       //std::cerr << "  " << b;
    //}
    //std::cerr << "\n";
    TMatrixT<double> background_full_covariance_matrix = mychi->CalcCovarianceMatrix(covmat, bgvec);
    TMatrixT<double> background_collapsed_covariance_matrix(nBinsColl, nBinsColl);
    mychi->CollapseModes(background_full_covariance_matrix, background_collapsed_covariance_matrix);    
    TMatrixT<double> inverse_background_collapsed_covariance_matrix = mychi->InvertMatrix(background_collapsed_covariance_matrix);


    // These are for write out of FC results
    HighFive::DataSet d_last_chi_min    = f_out->getDataSet("last_chi_min"   );
    HighFive::DataSet d_delta_chi       = f_out->getDataSet("delta_chi"      );
    HighFive::DataSet d_best_grid_point = f_out->getDataSet("best_grid_point");
    HighFive::DataSet d_n_iter          = f_out->getDataSet("n_iter"         );
    
    // write out this grid and universe
    HighFive::DataSet d_i_grid          = f_out->getDataSet("i_grid");
    HighFive::DataSet d_i_univ          = f_out->getDataSet("i_univ");

    // This is for the fake data dump
    HighFive::DataSet d_fakedata          = f_out->getDataSet("fakedata");
    
    // These are our work vectors, we will fill them with numbers from hdf5 as needed  
    //std::vector<double> specfull = {nBinsFull};
    std::vector<double> speccoll = {nBinsColl};
    std::vector< std::vector<double> > allColl;


    for(size_t r =0; r < mygrid.f_num_total_points; r++){
       // Load spectrum number r from HDF into speccoll
        g_speccoll.select(   {r, 0}, {1, nBinsColl}).read(speccoll);
        allColl.push_back(speccoll);
    }
    
    //sbn::SBNchi * mychi  = mychi; 


    //step 0. Make a fake-data-experiment for this point, drawn from covariance
    std::vector<float> fake_data= mychi->SampleCovariance(&myspec);
    double sum = std::accumulate(fake_data.begin(), fake_data.end(), 0.0);
    if (sum <10) {
       std::cerr << "sum of fakedata: " << sum << "\n";
       myspec.PrintFullVector();
       std::cerr << "\n";
       myspec.PrintCollapsedVector();
       abort();

    }

    
    d_fakedata.select(   {size_t(cp.gid()), 0}, {1, size_t(nBinsColl)}).write(fake_data);
    float last_chi_min = FLT_MAX;
    int best_grid_point = -99;
    
    TMatrixT<double> inverse_current_collapsed_covariance_matrix = inverse_background_collapsed_covariance_matrix;  
    size_t n_iter = 0;
    
    if (!dry) {
       for(n_iter = 0; n_iter < max_number_iterations; n_iter++){
    
           //Step 1. What covariance matrix do we use?
           //For first iteration, use the precalculated background only inverse covariance matrix.
           //For all subsequent iterations what is the full covariance matrix? Use the last best grid point.
           if(n_iter!=0){
               //Calculate current full covariance matrix, collapse it, then Invert.
               // Fill vector specfull with the full vector from hdf5 at position best_grid_point
               g_specfull.select(   {best_grid_point, 0}, {1, nBinsFull}).read(specfull);
               TMatrixT<double> current_full_covariance_matrix = mychi->CalcCovarianceMatrix(covmat, specfull);
               TMatrixT<double> current_collapsed_covariance_matrix(nBinsColl, nBinsColl);
               mychi->CollapseModes(current_full_covariance_matrix, current_collapsed_covariance_matrix);

               inverse_current_collapsed_covariance_matrix = mychi->InvertMatrix(current_collapsed_covariance_matrix);
           }
    
           //Step 2.0 Find the global_minimum_for this universe. Integrate in SBNfit minimizer here, a grid scan for now.
           float chi_min = FLT_MAX;
    
           for(size_t r =0; r < mygrid.f_num_total_points; r++){
              // Load spectrum number r from HDF into speccoll
               //g_speccoll.select(   {r, 0}, {1, nBinsColl}).read(speccoll);
               //float chi_tmp = myfeld.CalcChi(fake_data, speccoll,  inverse_current_collapsed_covariance_matrix);
               float chi_tmp = calcChi(fake_data, allColl[r],  inverse_current_collapsed_covariance_matrix);
               if(chi_tmp < chi_min){
                   best_grid_point = r;
                   chi_min = chi_tmp;
               }
           }
    
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
       float this_chi = calcChi(fake_data, myspec.collapsed_vector,inverse_current_collapsed_covariance_matrix);
    
       //step 4 calculate the delta_chi for this universe
       // Write out numbers of interest
       std::vector<double> v_last_chi_min    = { last_chi_min };
       std::vector<double> v_delta_chi       = { this_chi-last_chi_min };
       std::vector<int>    v_best_grid_point = { best_grid_point };
       std::vector<int>    v_n_iter          = { n_iter };
       std::vector<int>    v_i_grid          = { i_grid };
       std::vector<int>    v_i_univ          = { i_univ };
       d_last_chi_min.select(   {size_t(cp.gid()), 0}, {1,1}).write( v_last_chi_min   );
       d_delta_chi.select(      {size_t(cp.gid()), 0}, {1,1}).write( v_delta_chi      );
       d_best_grid_point.select({size_t(cp.gid()), 0}, {1,1}).write( v_best_grid_point);
       d_n_iter.select(         {size_t(cp.gid()), 0}, {1,1}).write( v_n_iter         );
       d_i_grid.select(         {size_t(cp.gid()), 0}, {1,1}).write( v_i_grid         );
       d_i_univ.select(         {size_t(cp.gid()), 0}, {1,1}).write( v_i_univ         );
    }

    delete mychi;
    //delete true_chi; 
    endtime   = MPI_Wtime(); 
    if (rank==0) fmt::print(stderr, "[{}] gridp {} univ {} iteration {}  took {} seconds, chi2min: {}\n",cp.gid(), i_grid, i_univ, n_iter, endtime-starttime, last_chi_min);
}



// --- main program ---//
int main(int argc, char* argv[])
{
    diy::mpi::environment env(argc, argv);
    diy::mpi::communicator world;

    size_t nBlocks = 0;
    size_t nBins=246;
    size_t nBinsC=90;
    int nPoints=1000;
    int nUniverses=1;
    std::string out_file="test.hdf5";
    std::string in_file="";
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
    // get command line arguments
    using namespace opts;
    Options ops(argc, argv);
    //ops >> Option('p', "npoints",   nPoints,   "Number of parameter points");
    ops >> Option('u', "nuniverse", nUniverses, "Number of universes");
    ops >> Option('b', "nblocks",   nBlocks,    "Number of blocks");
    //ops >> Option('n', "nbins",     nBins,      "Number of bins in 2d dataset");
    ops >> Option("nbinsC",    nBinsC,     "Number of collapsed bins in 2d dataset");
    ops >> Option('o', "output",    out_file,  "Output filename.");
    ops >> Option('g', "background",    f_BG,  "Backgrounds filename.");
    ops >> Option("core",    f_CV,  "Central values filename.");
    ops >> Option('c', "covmat",    f_COV,  "Covariance matrix filename.");
    ops >> Option('f', "fin",    in_file,  "Output filename.");
    ops >> Option('t', "tag",    tag,  "Tag.");
    ops >> Option('x', "xml",    xml,  "XML config.");
    ops >> Option("xmin",    xmin,   "xmin");
    ops >> Option("xmax",    xmax,   "xmax");
    ops >> Option("xwidth",  xwidth, "xwidth");
    ops >> Option("ymin",    ymin,   "ymin");
    ops >> Option("ymax",    ymax,   "ymax");
    ops >> Option("ywidth",  ywidth, "ywidth");
    bool verbose     = ops >> Present('v', "verbose", "verbose output");
    bool dryrun     = ops >> Present("dry", "dry run");
    if (ops >> Present('h', "help", "Show help"))
    {
        std::cout << "Usage:  [OPTIONS]\n";
        std::cout << ops;
        return 1;
    }

    
    NGrid mygrid;

    mygrid.AddDimension("m4",  xmin, xmax, xwidth);//0.1
    mygrid.AddDimension("ue4", ymin, ymax, ywidth);//0.1
    mygrid.AddFixedDimension("um4",0.0); //0.05

    
    nPoints = mygrid.f_num_total_points;

    std::string line,text;
    std::ifstream in(xml);
    while(std::getline(in, line))  text += line + "\n";
    const char* xmldata = text.c_str();

    sbn::SBNconfig myconf(xmldata, false);

    // Background
    std::vector<TH1D> bghist = readHistos(f_BG, myconf.fullnames);
    std::vector<double> bgvec = flattenHistos(bghist);

    // Central values
    std::vector<TH1D> cvhist = readHistos(f_CV, myconf.fullnames);

    // Read the covariance matrix on every rank --- we'll think about broadcast later. Maybe Eigen?
    TMatrixD covmat = readFracCovMat(f_COV);


    nBins =  covmat.GetNcols();

    
    if( world.rank()==0 ) {
      fmt::print(stderr, "\n*** This is diy running highfivewrite ***\n");
      fmt::print(stderr, "\n    Output will be written to {}\n", out_file);
      fmt::print(stderr, "\n    Points:  {}\n", nPoints);
      fmt::print(stderr, "\n    nBins:  {}\n", nBins);
      fmt::print(stderr, "\n    Universes: {}\n", nUniverses);
      fmt::print(stderr, "\n    f_BG:  {}\n", f_BG);
      fmt::print(stderr, "\n    f_CV:  {}\n", f_CV);
      fmt::print(stderr, "\n    f_COV:  {}\n", f_COV);
      fmt::print(stderr, "\n    Total size of dataset:  {}\n", nPoints*nUniverses);
      fmt::print(stderr, "***********************************\n");
    }
    
    // Create hdf5 file structure here 
    HighFive::File* f_out  = new HighFive::File(out_file,
                        HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate,
                        HighFive::MPIOFileDriver(MPI_COMM_WORLD,MPI_INFO_NULL));

    f_out->createDataSet<double>("specfull", HighFive::DataSpace( { nPoints, nBins  }));
    f_out->createDataSet<double>("speccoll", HighFive::DataSpace( { nPoints, nBinsC }));

    f_out->createDataSet<double>("fakedata",     HighFive::DataSpace( { nPoints*nUniverses, nBinsC  }));
   
    f_out->createDataSet<double>("last_chi_min", HighFive::DataSpace( {nPoints*nUniverses, 1} ));
    f_out->createDataSet<double>("delta_chi",    HighFive::DataSpace( {nPoints*nUniverses, 1} ));
    f_out->createDataSet<int>("best_grid_point", HighFive::DataSpace( {nPoints*nUniverses, 1} ));
    f_out->createDataSet<int>("n_iter",          HighFive::DataSpace( {nPoints*nUniverses, 1} ));
   
    // Some bookkeeping why not 
    f_out->createDataSet<int>("i_grid",          HighFive::DataSpace( {nPoints*nUniverses, 1} ));
    f_out->createDataSet<int>("i_univ",      HighFive::DataSpace( {nPoints*nUniverses, 1} ));
    f_out->createDataSet<double>("gridx",      HighFive::DataSpace( {nPoints, 1} ));
    f_out->createDataSet<double>("gridy",      HighFive::DataSpace( {nPoints, 1} ));
   
    // First rank also writes the grid so we know what the poins actually are
    if (world.rank() == 0) {
       std::vector<std::vector<double>> coords = mygrid.GetGrid();
       std::vector<double> xcoord;
       std::vector<double> ycoord;

       for (int i=0; i< nPoints; i++) {
          xcoord.push_back(coords[i][0]);
          ycoord.push_back(coords[i][1]);
       }
       HighFive::DataSet d_gridx          = f_out->getDataSet("gridx");
       HighFive::DataSet d_gridy          = f_out->getDataSet("gridy");
       d_gridx.select(   {0, 0}, {xcoord.size(), 1}).write(xcoord);
       d_gridy.select(   {0, 0}, {ycoord.size(), 1}).write(ycoord);
    }



    // First diy setup to load specra 
    size_t blocks = nPoints;
    Bounds domain;
    domain.min[0] = 0;
    domain.max[0] = blocks-1;
    diy::RoundRobinAssigner        assigner(world.size(), blocks);
    diy::RegularDecomposer<Bounds> decomposer(1, domain, blocks);
    diy::RegularBroadcastPartners  comm(    decomposer, 2, true);
    diy::RegularMergePartners      partners(decomposer, 2, true);
    diy::Master master(world, 1, -1, &Block::create, &Block::destroy);
    diy::decompose(1, world.rank(), domain, assigner, master);//, share_face, wrap, ghosts);

    double starttime, endtime;
    starttime = MPI_Wtime();
    master.foreach([world, tag, xml, mygrid, f_out, xmldata, bghist, cvhist, covmat](Block* b, const diy::Master::ProxyWithLink& cp)
                           {loadSpectrum(b, cp, world.size(), world.rank(), tag, xml, mygrid, f_out, xmldata, bghist, cvhist, covmat); });
    endtime   = MPI_Wtime(); 
    if (world.rank()==0) fmt::print(stderr, "[{}] loading spectra took {} seconds\n",world.rank(), endtime-starttime);


    //// Now more blocks as we have universes
    blocks = nPoints*nUniverses;
    if (world.rank()==0) fmt::print(stderr, "FC will be done on {} blocks, distributed over {} ranks\n", blocks, world.size());
    Bounds fc_domain;
    fc_domain.min[0] = 0;
    fc_domain.max[0] = blocks-1;
    diy::RoundRobinAssigner        fc_assigner(world.size(), blocks);
    diy::RegularDecomposer<Bounds> fc_decomposer(1, fc_domain, blocks);
    diy::RegularBroadcastPartners  fc_comm(    fc_decomposer, 2, true);
    diy::RegularMergePartners      fc_partners(fc_decomposer, 2, true);
    diy::Master                    fc_master(world, 1, -1, &Block::create, &Block::destroy);
    diy::decompose(1, world.rank(), fc_domain, fc_assigner, fc_master);//, share_face, wrap, ghosts);

    starttime = MPI_Wtime();
    fc_master.foreach([world, tag, xml, mygrid, f_out, nUniverses, dryrun, xmldata, bghist, cvhist, covmat, bgvec](Block* b, const diy::Master::ProxyWithLink& cp)
                           {doFCsmart(b, cp, world.size(), world.rank(), tag, xml, mygrid, f_out, nUniverses, dryrun, xmldata, bghist, cvhist, covmat, bgvec); });
    endtime   = MPI_Wtime(); 
    if (world.rank()==0) fmt::print(stderr, "[{}] FC took {} seconds\n",world.rank(), endtime-starttime);

    delete f_out;
    return 0;
}

