
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
    //int nUni;
    //int nPts;
    ////std::vector<int> myUnis;
    ////std::vector<int> myPts;
    //int dsIndex(int uni, int pt) {
};

void singlePoint(Block* b, diy::Master::ProxyWithLink const& cp, int size, int rank, string tag, string xml, int numUniverses, NGrid mygrid) {
    
    sbn::SBNfeld myfeld(mygrid,tag,xml);
    std::cout<<"Begininning a full Feldman-Cousins analysis for tag : "<<tag<<std::endl;

    myfeld.SetCoreSpectrum(tag+"_BKG_ONLY.SBNspec.root");
    myfeld.SetFractionalCovarianceMatrix(tag+".SBNcovar.root","frac_covariance");

    double random_number_seed = -1;
    std::cout<<"Setting random seed "<<random_number_seed<<std::endl;
    myfeld.SetRandomSeed(random_number_seed);
    std::cout<<"Loading precomputed spectra"<<std::endl;
    myfeld.LoadPreOscillatedSpectra();
    //std::cout <<"DONE loading precomputed spectra at : " << difftime(time(0), start_time)/60.0 << " Minutes.\n";
    myfeld.LoadBackgroundSpectrum();

    myfeld.SetNumUniverses(numUniverses);

    std::cout<<"Calculating the necessary SBNchi objects"<<std::endl;
    myfeld.CalcSBNchis();

   std::cout<<"Beginning to peform Single Grid PointFeldmanCousins analysis on pt: "<<cp.gid()+5<<std::endl;
   myfeld.PointFeldmanCousins((size_t) cp.gid()+5);
   
   // HERE we need something that purges myfeld of all root stuff.
}


void process_block(Block* b, diy::Master::ProxyWithLink const& cp, int size, int rank, size_t nbins, bool verbose, HighFive::File* f_out, const vector<double> data, const vector<double> specdata)
{
   HighFive::DataSet test = f_out->getDataSet("dataset");
   test.select({std::size_t(cp.gid())}, {1}).write(data);
   

   HighFive::DataSet spec = f_out->getDataSet("dataset2d");
   spec.select(   {std::size_t(    cp.gid()), 0}, {1, nbins}).write(specdata);
}

//
//
void loadSpectrum(Block* b, diy::Master::ProxyWithLink const& cp, int size, int rank, string tag, string xml, NGrid mygrid, HighFive::File* f_out) {

    sbn::SBNfeld myfeld(mygrid, tag, xml);

    myfeld.SetCoreSpectrum(tag+"_BKG_ONLY.SBNspec.root");

    double random_number_seed = -1;
    myfeld.SetRandomSeed(random_number_seed);

    sbn::SBNspec* myspec = myfeld.LoadPreOscillatedSpectrum(cp.gid());
    std::vector<double> full = myspec->full_vector;
    std::vector<double> coll = myspec->collapsed_vector;

    HighFive::DataSet spec = f_out->getDataSet("specfull");
    spec.select(   {std::size_t(cp.gid()), 0}, {1, std::size_t(full.size())}).write(full);

    HighFive::DataSet speccoll = f_out->getDataSet("speccoll");
    speccoll.select(   {std::size_t(cp.gid()), 0}, {1, std::size_t(coll.size())}).write(coll);
    
    delete myspec;

}

//void loadSpectrum2(Block* b, diy::Master::ProxyWithLink const& cp, int size, int rank, string tag, string xml, NGrid mygrid, HighFive::File* f_out) {

    //sbn::SBNfeld myfeld(mygrid, tag, xml);

    //myfeld.SetCoreSpectrum(tag+"_BKG_ONLY.SBNspec.root");
    //myfeld.SetFractionalCovarianceMatrix(tag+".SBNcovar.root","frac_covariance");

    //double random_number_seed = -1;
    //myfeld.SetRandomSeed(random_number_seed);

    //myfeld.LoadBackgroundSpectrum();
    //sbn::SBNspec* myspec = myfeld.LoadPreOscillatedSpectrum(cp.gid());
    //std::vector<double> full = myspec->full_vector;
    //std::vector<double> coll = myspec->collapsed_vector;

    //HighFive::DataSet spec = f_out->getDataSet("specfull");
    //spec.select(   {std::size_t(cp.gid()), 0}, {1, std::size_t(full.size())}).write(full);

    //HighFive::DataSet speccoll = f_out->getDataSet("speccoll");
    //speccoll.select(   {std::size_t(cp.gid()), 0}, {1, std::size_t(coll.size())}).write(coll);
    
    //sbn::SBNchi* mychi = new sbn::SBNchi(*myspec, *myfeld.fullFracCovMat(), myfeld.xmlname, false, myfeld.seed());
    ////Ok take the background only spectrum and form a background only covariance matrix. CalcCovarianceMatrix includes stats
    //TMatrixT<double> C_full = mychi->CalcCovarianceMatrix(myfeld.fullFracCovMat(), full);
    //TMatrixT<double> C_coll(myfeld.bgBinsCompressed(), myfeld.bgBinsCompressed());
    //mychi->CollapseModes(C_full, C_coll);
    //TMatrixT<double> C_coll_inv(myfeld.bgBinsCompressed(), myfeld.bgBinsCompressed()); 
    //C_coll_inv = mychi->InvertMatrix(C_coll);

    //std::vector<double> c_flat;//(myfeld.bgBinsCompressed()*myfeld.bgBinsCompressed());
    //for (int i=0; i<myfeld.bgBinsCompressed();i++) {
       //for (int j=0; j<myfeld.bgBinsCompressed();j++) {
          //c_flat.push_back(C_coll_inv[i][j]);
       //}
    //}
    //HighFive::DataSet d_invcollcov = f_out->getDataSet("invcollcov");
    //d_invcollcov.select(   {std::size_t(cp.gid()), 0}, {1, c_flat.size()}).write(c_flat);
    //delete myspec;
    //delete mychi;

//}

// FC for single point
// TODO: * seed, max #it, chi_min conv ---> block members?
//       * some sort of progress report
void doFC(Block* b, diy::Master::ProxyWithLink const& cp, int size, int rank, string tag, string xml, NGrid mygrid, HighFive::File* f_out, int num_universes) {
    
    sbn::SBNfeld myfeld(mygrid, tag, xml);

    myfeld.SetCoreSpectrum(tag+"_BKG_ONLY.SBNspec.root");
    myfeld.SetFractionalCovarianceMatrix(tag+".SBNcovar.root","frac_covariance");

    double random_number_seed = -1;
    myfeld.SetRandomSeed(random_number_seed);
    
    sbn::SBNspec* myspec = myfeld.LoadPreOscillatedSpectrum(cp.gid());

    myfeld.LoadBackgroundSpectrum();
    TMatrixT<double> stat_only_matrix(myfeld.num_bins_total, myfeld.num_bins_total);
    stat_only_matrix.Zero();

    sbn::SBNchi* mychi=NULL; 

    if(myfeld.statOnly()){
         mychi = new sbn::SBNchi(*myspec, stat_only_matrix, myfeld.xmlname, false, myfeld.seed());
     }else{
         mychi = new sbn::SBNchi(*myspec, *myfeld.fullFracCovMat(), myfeld.xmlname, false, myfeld.seed());
     }
    int nBinsFull = myspec->full_vector.size();
    int nBinsColl = myspec->collapsed_vector.size();

    // Single grid point multiuniverse FC
    int max_number_iterations = 5;
    double chi_min_convergance_tolerance = 0.001;

    //Ok take the background only spectrum and form a background only covariance matrix. CalcCovarianceMatrix includes stats
    TMatrixT<double> background_full_covariance_matrix = mychi->CalcCovarianceMatrix(myfeld.fullFracCovMat(), *myfeld.bgSpectrum());
    TMatrixT<double> background_collapsed_covariance_matrix(myfeld.bgBinsCompressed(), myfeld.bgBinsCompressed());
    mychi->CollapseModes(background_full_covariance_matrix, background_collapsed_covariance_matrix);    
    TMatrixT<double> inverse_background_collapsed_covariance_matrix = mychi->InvertMatrix(background_collapsed_covariance_matrix);   

    //fmt::print(stderr, "[{}] start FC\n", cp.gid());

    // These HDF5 dataset have all the full/collapsed vector spectra as loaded in loadSpectrum
    HighFive::DataSet g_specfull = f_out->getDataSet("specfull");
    HighFive::DataSet g_speccoll = f_out->getDataSet("speccoll");

    // These are for write out of FC results
    HighFive::DataSet d_last_chi_min    = f_out->getDataSet("last_chi_min"   );
    HighFive::DataSet d_delta_chi       = f_out->getDataSet("delta_chi"      );
    HighFive::DataSet d_best_grid_point = f_out->getDataSet("best_grid_point");
    HighFive::DataSet d_n_iter          = f_out->getDataSet("n_iter"         );

    // This is for the fake data dump
    HighFive::DataSet d_fakedata          = f_out->getDataSet("fakedata");
      
    // These are our work vectors, we will fill them with numbers from hdf5 as needed  
    std::vector<double> specfull = {nBinsFull};
    std::vector<double> speccoll = {nBinsColl};

    sbn::SBNspec* true_spec = myspec; 
    sbn::SBNchi * true_chi  = mychi; 

    size_t offset = 0;
    for(size_t i=0; i< num_universes; i++){
        offset = i*mygrid.f_num_total_points;

        //step 0. Make a fake-data-experimet for this point, drawn from covariance
        std::vector<float> fake_data= true_chi->SampleCovariance(true_spec);
        d_fakedata.select(   {offset + size_t(cp.gid()), 0}, {1, size_t(nBinsColl)}).write(fake_data);
        float last_chi_min = FLT_MAX;
        int best_grid_point = -99;

        TMatrixT<double> inverse_current_collapsed_covariance_matrix = inverse_background_collapsed_covariance_matrix;  
        size_t n_iter = 0;
        for(n_iter = 0; n_iter < max_number_iterations; n_iter++){

            //Step 1. What covariance matrix do we use?
            //For first iteration, use the precalculated background only inverse covariance matrix.
            //For all subsequent iterations what is the full covariance matrix? Use the last best grid point.
            if(n_iter!=0){
                //Calculate current full covariance matrix, collapse it, then Invert.
                // Fill vector specfull with the full vector from hdf5 at position best_grid_point
                g_specfull.select(   {best_grid_point, 0}, {1, nBinsFull}).read(specfull);
                TMatrixT<double> current_full_covariance_matrix = true_chi->CalcCovarianceMatrix(myfeld.fullFracCovMat(), specfull);
                TMatrixT<double> current_collapsed_covariance_matrix(myfeld.bgBinsCompressed(), myfeld.bgBinsCompressed());
                true_chi->CollapseModes(current_full_covariance_matrix, current_collapsed_covariance_matrix);    
                inverse_current_collapsed_covariance_matrix = true_chi->InvertMatrix(current_collapsed_covariance_matrix);   
            }

            //Step 2.0 Find the global_minimum_for this universe. Integrate in SBNfit minimizer here, a grid scan for now.
            float chi_min = FLT_MAX;

            for(size_t r =0; r < mygrid.f_num_total_points; r++){
               // Load spectrum number r from HDF into speccoll
                g_speccoll.select(   {r, 0}, {1, nBinsColl}).read(speccoll);
                float chi_tmp = myfeld.CalcChi(fake_data, speccoll,  inverse_current_collapsed_covariance_matrix);
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
        float this_chi = myfeld.CalcChi(fake_data, true_spec->collapsed_vector,inverse_current_collapsed_covariance_matrix);

        //step 4 calculate the delta_chi for this universe
        // Write out numbers of interest
        std::vector<double> v_last_chi_min    = { last_chi_min };
        std::vector<double> v_delta_chi       = { this_chi-last_chi_min };
        std::vector<int>    v_best_grid_point = { best_grid_point };
        std::vector<int>    v_n_iter          = { n_iter };
        d_last_chi_min.select(   {offset + size_t(cp.gid()), 0}, {1,1}).write( v_last_chi_min   );
        d_delta_chi.select(      {offset + size_t(cp.gid()), 0}, {1,1}).write( v_delta_chi      );
        d_best_grid_point.select({offset + size_t(cp.gid()), 0}, {1,1}).write( v_best_grid_point);
        d_n_iter.select(         {offset + size_t(cp.gid()), 0}, {1,1}).write( v_n_iter         );
    } // End loop over universes

    delete mychi;
    delete myspec;
}

// This is doing one univere for one gridpoint
void doFCsmart(Block* b, diy::Master::ProxyWithLink const& cp, int size, int rank, string tag, string xml, NGrid mygrid, HighFive::File* f_out, int num_universes) {
    // Some arithmetic to figure out the gridpoint and universe from cp.gid
    int i_grid = cp.gid() % mygrid.f_num_total_points;
    int i_univ = floor(cp.gid()/mygrid.f_num_total_points);
     
    //fmt::print(stderr, "[{}] doing grid point {}\n", cp.gid(), i_grid);

    double starttime, endtime;
    starttime = MPI_Wtime();
    sbn::SBNfeld myfeld(mygrid, tag, xml);
    myfeld.SetCoreSpectrum(tag+"_BKG_ONLY.SBNspec.root");
    myfeld.SetFractionalCovarianceMatrix(tag+".SBNcovar.root","frac_covariance");

    // TODO: Setting the seed to cp.gid()???
    double random_number_seed = -1;
    myfeld.SetRandomSeed(random_number_seed);
    
    sbn::SBNspec* myspec = myfeld.LoadPreOscillatedSpectrum(i_grid);

    myfeld.LoadBackgroundSpectrum();
    TMatrixT<double> stat_only_matrix(myfeld.num_bins_total, myfeld.num_bins_total);
    stat_only_matrix.Zero();

    sbn::SBNchi* mychi=NULL; 

    if(myfeld.statOnly()){
         mychi = new sbn::SBNchi(*myspec, stat_only_matrix, myfeld.xmlname, false, myfeld.seed());
     }else{
         mychi = new sbn::SBNchi(*myspec, *myfeld.fullFracCovMat(), myfeld.xmlname, false, myfeld.seed());
     }
    int nBinsFull = myspec->full_vector.size();
    int nBinsColl = myspec->collapsed_vector.size();

    // Single grid point multiuniverse FC
    int max_number_iterations = 5;
    double chi_min_convergance_tolerance = 0.001;

    //Ok take the background only spectrum and form a background only covariance matrix. CalcCovarianceMatrix includes stats
    TMatrixT<double> background_full_covariance_matrix = mychi->CalcCovarianceMatrix(myfeld.fullFracCovMat(), *myfeld.bgSpectrum());
    TMatrixT<double> background_collapsed_covariance_matrix(myfeld.bgBinsCompressed(), myfeld.bgBinsCompressed());
    mychi->CollapseModes(background_full_covariance_matrix, background_collapsed_covariance_matrix);    
    TMatrixT<double> inverse_background_collapsed_covariance_matrix = mychi->InvertMatrix(background_collapsed_covariance_matrix);


    //endtime   = MPI_Wtime(); 
    //fmt::print(stderr, "[{}] preparation took {} seconds\n",cp.gid(), endtime-starttime);

    //fmt::print(stderr, "[{}] start FC\n", cp.gid());

    // These HDF5 dataset have all the full/collapsed vector spectra as loaded in loadSpectrum
    HighFive::DataSet g_specfull = f_out->getDataSet("specfull");
    HighFive::DataSet g_speccoll = f_out->getDataSet("speccoll");

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
    std::vector<double> specfull = {nBinsFull};
    std::vector<double> speccoll = {nBinsColl};
    std::vector< std::vector<double> > allColl;


   for(size_t r =0; r < mygrid.f_num_total_points; r++){
      // Load spectrum number r from HDF into speccoll
       g_speccoll.select(   {r, 0}, {1, nBinsColl}).read(speccoll);
       allColl.push_back(speccoll);
   }


    sbn::SBNspec* true_spec = myspec; 
    sbn::SBNchi * true_chi  = mychi; 

    //HighFive::DataSet d_invcollcov          = f_out->getDataSet("invcollcov");
    //std::vector<double> c_flat;
    //c_flat.reserve(myfeld.bgBinsCompressed()*myfeld.bgBinsCompressed());

     //step 0. Make a fake-data-experimet for this point, drawn from covariance
     std::vector<float> fake_data= true_chi->SampleCovariance(true_spec);
     d_fakedata.select(   {size_t(cp.gid()), 0}, {1, size_t(nBinsColl)}).write(fake_data);
     float last_chi_min = FLT_MAX;
     int best_grid_point = -99;

     TMatrixT<double> inverse_current_collapsed_covariance_matrix = inverse_background_collapsed_covariance_matrix;  
     size_t n_iter = 0;
     for(n_iter = 0; n_iter < max_number_iterations; n_iter++){

         //Step 1. What covariance matrix do we use?
         //For first iteration, use the precalculated background only inverse covariance matrix.
         //For all subsequent iterations what is the full covariance matrix? Use the last best grid point.
         if(n_iter!=0){
             //Calculate current full covariance matrix, collapse it, then Invert.
             // Fill vector specfull with the full vector from hdf5 at position best_grid_point
             g_specfull.select(   {best_grid_point, 0}, {1, nBinsFull}).read(specfull);
             TMatrixT<double> current_full_covariance_matrix = true_chi->CalcCovarianceMatrix(myfeld.fullFracCovMat(), specfull);
             TMatrixT<double> current_collapsed_covariance_matrix(myfeld.bgBinsCompressed(), myfeld.bgBinsCompressed());
             true_chi->CollapseModes(current_full_covariance_matrix, current_collapsed_covariance_matrix);    
             inverse_current_collapsed_covariance_matrix = true_chi->InvertMatrix(current_collapsed_covariance_matrix);
             //d_invcollcov.select(   {best_grid_point, 0}, {1, myfeld.bgBinsCompressed()*myfeld.bgBinsCompressed()}).read(c_flat);
             //int vpos=0;
             //for (int i=0; i<myfeld.bgBinsCompressed();i++) {
                //for (int j=0; j<myfeld.bgBinsCompressed();j++) {
                   //inverse_current_collapsed_covariance_matrix[i][j] = c_flat[vpos];
                   //vpos++;
                //}
             //}

             

         }

         //Step 2.0 Find the global_minimum_for this universe. Integrate in SBNfit minimizer here, a grid scan for now.
         float chi_min = FLT_MAX;

         for(size_t r =0; r < mygrid.f_num_total_points; r++){
            // Load spectrum number r from HDF into speccoll
             //g_speccoll.select(   {r, 0}, {1, nBinsColl}).read(speccoll);
             //float chi_tmp = myfeld.CalcChi(fake_data, speccoll,  inverse_current_collapsed_covariance_matrix);
             float chi_tmp = myfeld.CalcChi(fake_data, allColl[r],  inverse_current_collapsed_covariance_matrix);
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
     float this_chi = myfeld.CalcChi(fake_data, true_spec->collapsed_vector,inverse_current_collapsed_covariance_matrix);

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

    delete mychi;
    delete myspec;
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
    std::string tag="";
    std::string xml="";
    // get command line arguments
    using namespace opts;
    Options ops(argc, argv);
    //ops >> Option('p', "npoints",   nPoints,   "Number of parameter points");
    ops >> Option('u', "nuniverse", nUniverses, "Number of universes");
    ops >> Option('b', "nblocks",   nBlocks,   "Number of blocks");
    ops >> Option('n', "nbins",   nBins,   "Number of bins in 2d dataset");
    ops >> Option('n', "nbinsC",   nBinsC,   "Number of collapsed bins in 2d dataset");
    ops >> Option('o', "output",    out_file,  "Output filename.");
    ops >> Option('f', "fin",    in_file,  "Output filename.");
    ops >> Option('t', "tag",    tag,  "Tag.");
    ops >> Option('x', "xml",    xml,  "XML config.");
    bool verbose     = ops >> Present('v', "verbose", "verbose output");
    if (ops >> Present('h', "help", "Show help"))
    {
        std::cout << "Usage:  [OPTIONS]\n";
        std::cout << ops;
        return 1;
    }

    
    NGrid mygrid;

    mygrid.AddDimension("m4", -1.0, 1.1, 0.1);//0.1
    mygrid.AddDimension("ue4", -2.3, 0.1, 0.05);//0.1
    mygrid.AddFixedDimension("um4",0.0); //0.05

    
    nPoints = mygrid.f_num_total_points;

    
    if( world.rank()==0 ) {
      fmt::print(stderr, "\n*** This is diy running highfivewrite ***\n");
      fmt::print(stderr, "\n    Output will be written to {}\n", out_file);
      fmt::print(stderr, "\n    Points:  {}\n", nPoints);
      fmt::print(stderr, "\n    Universes: {}\n", nUniverses);
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
    master.foreach([world, tag, xml, mygrid, f_out](Block* b, const diy::Master::ProxyWithLink& cp)
                           {loadSpectrum(b, cp, world.size(), world.rank(), tag, xml, mygrid, f_out); });
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
    fc_master.foreach([world, tag, xml, mygrid, f_out, nUniverses](Block* b, const diy::Master::ProxyWithLink& cp)
                           {doFCsmart(b, cp, world.size(), world.rank(), tag, xml, mygrid, f_out, nUniverses); });
    endtime   = MPI_Wtime(); 
    if (world.rank()==0) fmt::print(stderr, "[{}] FC took {} seconds\n",world.rank(), endtime-starttime);
    return 0;
}

