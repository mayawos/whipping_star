
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

//#include <TH1D.h>
//#include <TFile.h>

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


void singlePoint(Block* b, diy::Master::ProxyWithLink const& cp, int size, int rank, string tag, string xml, int numUniverses) {
    
    NGrid mygrid;

    mygrid.AddDimension("m4", -1.0, 1.1, 0.1);//0.1
    mygrid.AddDimension("ue4", -2.3, 0.1, 0.05);//0.1
    mygrid.AddFixedDimension("um4",0.0); //0.05
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

    //if(grid_pt ==-1){
        //std::cout<<"Beginning to peform FullFeldmanCousins analysis"<<std::endl;
        //myfeld.FullFeldmanCousins();
    //}else if(grid_pt>=0){
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


// --- main program ---//
int main(int argc, char* argv[])
{
    diy::mpi::environment env(argc, argv);
    diy::mpi::communicator world;

    size_t nBlocks = 0;
    size_t nBins=10;
    int nPoints=1000;
    int nUniverses=1000;
    std::string out_file="test.hdf5";
    std::string in_file="";
    std::string tag="";
    std::string xml="";
    // get command line arguments
    using namespace opts;
    Options ops(argc, argv);
    ops >> Option('p', "npoints",   nPoints,   "Number of parameter points");
    ops >> Option('u', "nuniverse", nUniverses, "Number of universes");
    ops >> Option('b', "nblocks",   nBlocks,   "Number of blocks");
    ops >> Option('n', "nbins",   nBins,   "Number of bins in 2d dataset");
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

    
    // Create hdf5 file structure here 
    //HighFive::File* f_out  = new HighFive::File(out_file,
			//HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate,
			//HighFive::MPIOFileDriver(MPI_COMM_WORLD,MPI_INFO_NULL));


    //std::vector<size_t> ds_dims(1);
    //ds_dims[0] = nPoints*nUniverses;
    //ds_dims[1] = 1;
    //f_out->createDataSet<double>("dataset", HighFive::DataSpace(ds_dims));

    //std::vector<size_t> spec_dims(2);
    //spec_dims[0] = nPoints*nUniverses;
    //spec_dims[1] = nBins;
    //f_out->createDataSet<double>("dataset2d", HighFive::DataSpace(spec_dims));

    
    size_t blocks;
    if (nBlocks==0) blocks= nPoints;//world.size() * threads;
    else blocks=nBlocks;


    int dim =1;
    diy::FileStorage storage("./DIY.XXXXXX"); // used for blocks moved out of core
    Bounds domain;
    for (int i = 0; i < dim; ++i) {
      domain.min[i] = 0;
      domain.max[i] = blocks-1;
    }
    ////// choice of contiguous or round robin assigner
    diy::ContiguousAssigner   assigner(world.size(), blocks);
    //// decompose the domain into blocks
    //// This is a DIY regular way to assign neighbors. You can do this manually.
    diy::RegularDecomposer<Bounds> decomposer(dim, domain, blocks);

    int k = 2;       // the radix of the k-ary reduction tree
    diy::RegularBroadcastPartners comm(decomposer, k, true);

    diy::RegularMergePartners  partners(decomposer,  // domain decomposition
                                        k,           // radix of k-ary reduction
                                        true); // contiguous = true: distance doubling

    diy::Master master(world, 1, -1, &Block::create, &Block::destroy);
    diy::decompose(dim, world.rank(), domain, assigner, master);//, share_face, wrap, ghosts);


    if( world.rank()==0 ) {
      fmt::print(stderr, "\n*** This is diy running highfivewrite ***\n");
      fmt::print(stderr, "\n    Output will be written to {}\n", out_file);
      fmt::print(stderr, "\n    DIY blocks:              {}\n", blocks);
      fmt::print(stderr, "\n    Points:  {}\n", nPoints);
      fmt::print(stderr, "\n    Universes: {}\n", nUniverses);
      fmt::print(stderr, "\n    Total size of dataset:  {}\n", nPoints*nUniverses);
      fmt::print(stderr, "***********************************\n");
    }


    //master.foreach([world, verbose, nBins, f_out, data, specdata](Block* b, const diy::Master::ProxyWithLink& cp)
                           //{process_block(b, cp, world.size(), world.rank(), nBins, verbose, f_out, data, specdata); });
    master.foreach([world, tag, xml, nUniverses](Block* b, const diy::Master::ProxyWithLink& cp)
                           {singlePoint(b, cp, world.size(), world.rank(), tag, xml, nUniverses); });

    //delete f_out;
    return 0;
}

