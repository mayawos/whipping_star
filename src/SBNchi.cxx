#include "SBNchi.h"

#ifdef USE_GPU
#include "openacc.h"
#include "openacc_curand.h"
#endif
using namespace sbn;


std::vector<TMatrixT<double>> splitNormShape(TMatrixT<double> & Min,std::vector<double> & fullvec){
    int ncol = Min.GetNrows();
    std::vector<TMatrixT<double>> ans;
    for(int i=0; i<3;i++){
        ans.push_back(TMatrixT<double>(ncol,ncol));
        ans.back().Zero();
    }

    // ans[0] is shape, ans[1] is mixed, ans[2] is norm

    double Nt = std::accumulate(fullvec.begin(), fullvec.end(),0.0);

    for(int i=0; i<ncol; i++){
        for(int j=0; j<ncol; j++){

            ans[0](i,j)  = Min(i,j);
            ans[1](i,j)  = 0.0;
            ans[2](i,j)  = 0.0;
            for(int k=0; k<ncol;k++){
                ans[0](i,j)-fullvec[j]/Nt*Min(i,k) - fullvec[i]/Nt*Min(k,j);
                ans[1](i,j)+fullvec[j]/Nt*Min(i,k) + fullvec[i]/Nt*Min(k,j);

                for(int l=0; l<ncol; l++){
                    ans[0](i,j) += fullvec[i]*fullvec[j]/(Nt*Nt)*Min(k,l);
                    ans[1](i,j) += -2*fullvec[i]*fullvec[j]/(Nt*Nt)*Min(k,l);
                    ans[2](i,j) += fullvec[i]*fullvec[j]/(Nt*Nt)*Min(k,l);
                }

            };
        }
    }

    return ans;
}


/***********************************************
 *		Constructors
 * ********************************************/

SBNchi::SBNchi(std::string xml) : SBNconfig(xml,false){};
SBNchi::SBNchi(SBNspec in, TMatrixT<double> matrix_systematicsin) : SBNchi(in,matrix_systematicsin,true){}
SBNchi::SBNchi(SBNspec in, TMatrixT<double> matrix_systematicsin, bool is_verbose) : SBNchi(in,matrix_systematicsin,in.xmlname, is_verbose){}
SBNchi::SBNchi(SBNspec in, TMatrixT<double> matrix_systematicsin, std::string inxml, bool is_verbose) : SBNchi(in, matrix_systematicsin, inxml,  is_verbose,-1){}
SBNchi::SBNchi(SBNspec in, TMatrixT<double> matrix_systematicsin, std::string inxml, bool is_verbose, double random_seed) : SBNconfig(inxml, is_verbose), core_spectrum(in){

    last_calculated_chi = -9999999;
    is_stat_only= false;


    pseudo_from_collapsed = false;
    matrix_collapsed.ResizeTo(num_bins_total_compressed, num_bins_total_compressed);
    matrix_systematics.ResizeTo(num_bins_total, num_bins_total);
    matrix_fractional_covariance.ResizeTo(num_bins_total, num_bins_total);

    TMatrixD m = matrix_systematicsin;
    for (int i = 0; i < used_bins.size(); i++){
        TMatrixDColumn(m,i) = TMatrixDColumn(m,used_bins.at(i));
        m.ResizeTo(used_bins.size(),used_bins.size());
    }

    m_cmin = -999;
    m_cmax = -999;

    matrix_fractional_covariance = m;
    matrix_systematics.Zero();
    max_sample_chi_val =150.0;
    m_tolerance = 1e-8;

    this->InitRandomNumberSeeds(random_seed);
    this->ReloadCoreSpectrum(&core_spectrum);
}

//Alternative constrctors
SBNchi::SBNchi(SBNspec in, std::string newxmlname) : SBNconfig(newxmlname), core_spectrum(in){
    is_stat_only = false;

    matrix_collapsed.ResizeTo(num_bins_total_compressed, num_bins_total_compressed);
    matrix_fractional_covariance.ResizeTo(num_bins_total, num_bins_total);

    if(fullnames.size() !=in.fullnames.size()){
        std::cerr<<"ERROR: SBNchi::SBNchi | Selected covariance matrix and background spectrum are different sizes!"<<std::endl;
        exit(EXIT_FAILURE);
    }else{
        for(int i=0; i< fullnames.size(); i++){
            if(fullnames[i]!=in.fullnames[i]){
                std::cerr<<"ERROR: SBNchi::SBNchi | Spectrum and Covariance matrix have different (or different order) subchannels!"<<std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }
    m_cmin = -999;
    m_cmax = -999;


    m_tolerance = 1e-7;
    pseudo_from_collapsed = false;
    max_sample_chi_val =150.0;
    matrix_fractional_covariance = FillSystematicsFromXML();
    last_calculated_chi = -9999999;
    core_spectrum.CollapseVector();
    matrix_systematics.ResizeTo(num_bins_total,num_bins_total);
    matrix_systematics.Zero();
    matrix_systematics=matrix_fractional_covariance;

    this->InitRandomNumberSeeds();

    this->ReloadCoreSpectrum(&in);
}

SBNchi::SBNchi(SBNspec in) : SBNchi(in,true){}


SBNchi::SBNchi(SBNspec in, bool is_is_stat_only): SBNconfig(in.xmlname), core_spectrum(in), is_stat_only(is_is_stat_only){

    last_calculated_chi = -9999999;

    matrix_collapsed.ResizeTo(num_bins_total_compressed, num_bins_total_compressed);
    matrix_systematics.ResizeTo(num_bins_total, num_bins_total);
    matrix_fractional_covariance.ResizeTo(num_bins_total, num_bins_total);
m_cmin = -999;
    m_cmax = -999;


    m_tolerance = 1e-7;
    max_sample_chi_val =150.0;
    this->InitRandomNumberSeeds();

    pseudo_from_collapsed = false;

    if(is_is_stat_only){
        matrix_fractional_covariance.Zero();
        matrix_systematics.Zero();

    }else{
        matrix_fractional_covariance = FillSystematicsFromXML();
        matrix_systematics.Zero();
    }

    this->ReloadCoreSpectrum(&core_spectrum);

}


/***********************************************
 *		Rest for now
 * ********************************************/
void SBNchi::InitRandomNumberSeeds(){
    this->InitRandomNumberSeeds(-1);
}

void SBNchi::InitRandomNumberSeeds(double seed){

    if(seed<0){
        rangen_twister = new std::mt19937(random_device_seed());
        rangen_linear = new std::minstd_rand(random_device_seed());
        rangen_carry = new std::ranlux24_base(random_device_seed());
        rangen = new TRandom3(0);
    }else{
        rangen_twister = new std::mt19937(seed);
        rangen_linear = new std::minstd_rand(seed);
        rangen_carry = new std::ranlux24_base(seed);
        rangen = new TRandom3(seed);
    }


    m_dist_normal=new std::normal_distribution<float>;
    std::normal_distribution<float> dtemp(0.0,1.0);
    m_dist_normal->param(dtemp.param());

}

int SBNchi::ReloadCoreSpectrum(SBNspec *bkgin){
    otag = "SBNchi::ReloadCoreSpectrum\t|| ";

    bool is_fractional = true;
    cholosky_performed = false;

    if(is_verbose)std::cout<<otag<<"Begininning to reload core spec! First Set new core spec"<<std::endl;
    core_spectrum = *bkgin;
    core_spectrum.CollapseVector();

    if(is_verbose)std::cout<<otag<<" || Clear all previous chi^2 data"<<std::endl;
    vec_last_calculated_chi.clear();
    vec_last_calculated_chi.resize(num_bins_total_compressed, std::vector<double>( num_bins_total_compressed,0) );

    //Reset matrix_systematics to fractional
    if(is_verbose) std::cout<<otag<<" Reseting matrix_systematics to matrix_fractional_covariance"<<std::endl;
    matrix_systematics = matrix_fractional_covariance;

    if(matrix_systematics.GetNcols()!=num_bins_total ){
        std::cout<<otag<<"ERROR: trying to pass a matrix to SBNchi that isnt the right size"<<std::endl;
        std::cout<<otag<<"ERROR: num_bins_total: "<<num_bins_total<<" and matrix is: "<<matrix_systematics.GetNcols()<<std::endl;
        exit(EXIT_FAILURE);
    }

    if(is_verbose)std::cout<<otag<<"Go from fracCovariance to fullCovariance. matrix_systematics.GetNcols(): "<<matrix_systematics.GetNcols()<<" matrix_systematics.GetNrows(): "<<matrix_systematics.GetNrows()<<" core->fullvec.size(): "<<core_spectrum.full_vector.size()<<std::endl;
    // systematics per scaled event
    for(int i =0; i<matrix_systematics.GetNcols(); i++)
    {
        //std::cout<<"KRAK: "<<core_spectrum.full_vector.at(i)<<std::endl;
        for(int j =0; j<matrix_systematics.GetNrows(); j++)
        {
            if(is_fractional){
                if(std::isnan(matrix_systematics(i,j)))
                    matrix_systematics(i,j) = 0;
                else
                    matrix_systematics(i,j) = matrix_systematics(i,j)*core_spectrum.full_vector.at(i)*core_spectrum.full_vector.at(j);
                    //if(i==j) std::cout << "1 mc err, mc stats err  =  " << matrix_systematics(i,j) << ", " << pow(core_spectrum.full_error[i],2) << std::endl;
                    //if(i==j) matrix_systematics(i,j) += pow(core_spectrum.full_error[i],2); //we are including the MC intrinsic error in the constraint so commenting this out for now
                    //if(i==j) std::cout << "2 mc err, mc stats err  =  " << matrix_systematics(i,j) << ", " << pow(core_spectrum.full_error[i],2) << std::endl;
            }
        }
    }

    if(is_verbose)std::cout<<otag<<"Filling stats into cov matrix"<<std::endl;
    // Fill stats from the back ground vector
    TMatrixT <double> Mstat(num_bins_total, num_bins_total);
    Mstat.Zero();
    FillStatsMatrix(Mstat, core_spectrum.full_vector);

    if(Mstat.IsSymmetric()){
        if(is_verbose)std::cout<<otag<<"Stat matrix is symmetric (it is just diagonal core)"<<std::endl;
    }else{
        std::cout<<otag<<"ERROR: SBNchi::FormCovarianceMatrix, stats  is not symmetric!"<<std::endl;
        exit(EXIT_FAILURE);
    }


    /*
       TMatrixD Mcorr = matrix_systematics;
    //TEMP
    TFile *f = new TFile("gloop.root","recreate");
    f->cd();
    for(int i=0; i<Mcorr.GetNrows();i++){
    for(int j=0; j<Mcorr.GetNrows();j++){
    Mcorr(i,j) = matrix_systematics(i,j)/(sqrt(matrix_systematics(i,i))*sqrt(matrix_systematics(j,j)));
    }
    }
    TH2D Hcorr(Mcorr);
    TH2D Hsys(matrix_systematics);

    f->cd();
    Hcorr.Write("corr");
    Hsys.Write("total");
    f->Close();
    */



    //And then define the total covariance matrix in all its glory
    TMatrixT <double> Mtotal(num_bins_total,num_bins_total);
    Mtotal.Zero();


    if(is_stat_only){
        if(is_verbose)std::cout<<otag<<"Using stats only in covariance matrix"<<std::endl;
        Mtotal = Mstat;
    }else{
        if(is_verbose)std::cout<<otag<<" Using stats+sys in covariance matrix"<<std::endl;
        Mtotal = Mstat + matrix_systematics;
    }


    //Also going to do matrix_systematics_collapsed;

    m_matrix_systematics_collapsed.ResizeTo(num_bins_total_compressed, num_bins_total_compressed);
    CollapseModes(matrix_systematics, m_matrix_systematics_collapsed);


    if(is_verbose)std::cout<<otag<<"Mstat: "<<Mstat.GetNrows()<<" x "<<Mstat.GetNcols()<<std::endl;
    if(is_verbose)std::cout<<otag<<"matrix_systematics: "<<matrix_systematics.GetNrows()<<" x "<<matrix_systematics.GetNcols()<<std::endl;
    if(is_verbose)std::cout<<otag<<"Mtotal: "<<Mtotal.GetNrows()<<" x "<<Mtotal.GetNcols()<<std::endl;

    if(Mtotal.IsSymmetric() ){
        if(is_verbose)	std::cout<<otag<<"Total Mstat +matrix_systematics is symmetric"<<std::endl;
    }else{

        double tol = m_tolerance;
        double biggest_deviation = 0;
        int bi=0;
        int bj=0;

        if(is_verbose)std::cout<<otag<<"WARNING: Stats + sys result appears to be not symmetric!"<<std::endl;
        for(int i=0; i<Mtotal.GetNrows(); i++){
            for(int j=0; j<Mtotal.GetNcols(); j++){
                double dev = fabs(Mtotal(i,j)-Mtotal(j,i));
                if(dev>biggest_deviation){
                    biggest_deviation = 2*dev/(fabs(Mtotal(i,j))+fabs(Mtotal(j,i)));
                    bi=i;
                    bj=j;
                }
                if(Mtotal(i,j)!=Mtotal(i,j)){

                    std::cout<<"ERROR: we have NAN's  Better check your inputs."<<std::endl;
                    exit(EXIT_FAILURE);

                }
            }
        }

        if(is_verbose) std::cout<<otag<<"WARNING: Biggest Relative Deviation from symmetry is i:"<<bi<<" j: "<<bj<<" of order "<<biggest_deviation<<" M(j,i)"<<Mtotal(bj,bi)<<" M(i,j)"<<Mtotal(bi,bj)<<std::endl;

        if(biggest_deviation >tol){

            std::cout<<"ERROR: Thats too unsymettric, killing process. Better check your inputs. Tolerance is "<< tol << "and the deviation is: " << biggest_deviation << std::endl;

            exit(EXIT_FAILURE);
        }else{

            if(is_verbose)	std::cout<<otag<<"WARNING: Thats within tolderence. Continuing."<<std::endl;
        }
    }

    TMatrixT<double > Mctotal(num_bins_total_compressed,num_bins_total_compressed);

    CollapseModes(Mtotal, Mctotal);

    matrix_collapsed = Mctotal;

    vec_matrix_collapsed = TMatrixDToVector(Mctotal);
    double invdet=0;

    TMatrixD McI(num_bins_total_compressed,num_bins_total_compressed);
    McI.Zero();

    if(is_verbose) std::cout<<otag<<" About to do a SVD decomposition"<<std::endl;
    TDecompSVD svd(Mctotal);
    if (!svd.Decompose()) {

        std::cout <<otag<<"Decomposition failed, matrix not symettric?, has nans?" << std::endl;
        std::cout<<otag<<"ERROR: The matrix to invert failed a SVD decomp!"<<std::endl;

        for(int i=0; i< num_bins_total_compressed; i++){
            for(int j=0; j< num_bins_total_compressed; j++){
                std::cout<<i<<" "<<j<<" "<<Mctotal(i,j)<<std::endl;
            }
        }

        exit(EXIT_FAILURE);
        return 0;
    } else {
        McI = svd.Invert();
    }
    if( !McI.IsValid()){
        std::cout<<otag<<"ERROR: The inverted matrix isnt valid! Something went wrong.."<<std::endl;
        exit(EXIT_FAILURE);

    }

    if(is_verbose)std::cout<<otag<<"SUCCESS! Inverted."<<std::endl;
    // matrix_systematics.Print();
    // McI.Print();
    //Mtotal.Print();

    vec_matrix_inverted = TMatrixDToVector(McI);


    // test for validity
    bool is_small_negative_eigenvalue = false;
    double tolerence_positivesemi = m_tolerance;


    //if a matrix is (a) real and (b) symmetric (checked above) then to prove positive semi-definite, we just need to check eigenvalues and >=0;

    TMatrixDEigen eigen (Mctotal);
    TVectorD eigen_values = eigen.GetEigenValuesRe();


    for(int i=0; i< eigen_values.GetNoElements(); i++){
        if(eigen_values(i)<0){
            is_small_negative_eigenvalue = true;
            if(fabs(eigen_values(i))> tolerence_positivesemi ){
                std::cout<<otag<<" collapsed covariance matrix contains (at least one)  negative eigenvalue: "<<eigen_values(i)<<std::endl;
                Mctotal.Print();
                std::cout<<otag<<" full covariance "<<std::endl;
                Mtotal.Print();
                exit(EXIT_FAILURE);
            }
        }
    }


    if(is_small_negative_eigenvalue){
        if(is_verbose)	std::cout<<otag<<"Generated covariance matrix is (allmost) positive semi-definite. It did contain small negative values of absolute value <= :"<<tolerence_positivesemi<<std::endl;
    }else{
        if(is_verbose)	std::cout<<otag<<"Generated covariance matrix is also positive semi-definite."<<std::endl;
    }

    core_spectrum.CollapseVector();

    return 0;
}



/*********************************************
 *		Different Chi^2 calculations
 * ******************************************/

//Standard chi^2 calculation
double SBNchi::CalcChi(SBNspec *sigSpec){
    double tchi = 0;

    if(sigSpec->collapsed_vector.size()==0){
        //        if(is_verbose)	std::cout<<"WARNING: SBNchi::CalcChi, inputted sigSpec has un-compressed vector, I am doing it now, but this is inefficient!"<<std::endl;
        sigSpec->CollapseVector();
    }

    int k=0;

    for(int i =0; i<num_bins_total_compressed; i++){
        for(int j =0; j<num_bins_total_compressed; j++){
            k++;

            if(i==j && fabs(vec_matrix_inverted.at(i).at(j)) > 1e16 && vec_matrix_inverted[i][j] < 0){
                std::cout<<"ERROR: SBNchi::CalcChi || diagonal of inverse covariance is negative! : "<<vec_matrix_inverted[i][j]<<" @ ("<<i<<","<<j<<")"<<std::endl;
            }
            vec_last_calculated_chi.at(i).at(j) =(core_spectrum.collapsed_vector.at(i)-sigSpec->collapsed_vector.at(i))*vec_matrix_inverted.at(i).at(j)*(core_spectrum.collapsed_vector.at(j)-sigSpec->collapsed_vector.at(j) );
            tchi += vec_last_calculated_chi.at(i).at(j);
        }
    }

    last_calculated_chi = tchi;
    return tchi;
}


float SBNchi::CalcChi(std::vector<float> * sigVec){
    float tchi = 0;

#ifndef _OPENACC
    if(sigVec->size() != num_bins_total_compressed ){
        std::cerr<<"ERROR: SBNchi::CalcChi(std::vector<float>) ~ your inputed vector does not have correct dimensions"<<std::endl;
        std::cerr<<"sigVec.size(): "<<sigVec->size()<<" num_bins_total_compressed: "<<num_bins_total_compressed<<std::endl;
        exit(EXIT_FAILURE);
    }
#endif

    for(int i =0; i<num_bins_total_compressed; i++){
        for(int j =0; j<num_bins_total_compressed; j++){
            tchi += (core_spectrum.collapsed_vector[i]-(*sigVec)[i])*vec_matrix_inverted[i][j]*(core_spectrum.collapsed_vector[j]-(*sigVec)[j] );
        }
    }

    //last_calculated_chi = tchi;
    return tchi;
}


double SBNchi::CalcChi(std::vector<double> * sigVec){
    double tchi = 0;

#ifndef _OPENACC
    if(sigVec->size() != num_bins_total_compressed ){
        std::cerr<<"ERROR: SBNchi::CalcChi(std::vector<double>) ~ your inputed vector does not have correct dimensions"<<std::endl;
        std::cerr<<"sigVec.size(): "<<sigVec->size()<<" num_bins_total_compressed: "<<num_bins_total_compressed<<std::endl;
        exit(EXIT_FAILURE);
    }
#endif

    for(int i =0; i<num_bins_total_compressed; i++){
        for(int j =0; j<num_bins_total_compressed; j++){
            tchi += (core_spectrum.collapsed_vector[i]-(*sigVec)[i])*vec_matrix_inverted[i][j]*(core_spectrum.collapsed_vector[j]-(*sigVec)[j] );
        }
    }

    //last_calculated_chi = tchi;
    return tchi;
}

double SBNchi::CalcChi(double* sigVec){
    double tchi = 0;

    for(int i =0; i<num_bins_total_compressed; i++){
        for(int j =0; j<num_bins_total_compressed; j++){
            tchi += (core_spectrum.collapsed_vector[i]-(sigVec)[i])*vec_matrix_inverted[i][j]*(core_spectrum.collapsed_vector[j]-(sigVec)[j] );
        }
    }

    //last_calculated_chi = tchi;
    return tchi;
}

float SBNchi::CalcChi(float **invert_matrix, float* core, float *sig){
    float tchi = 0;

    for(int i =0; i<num_bins_total_compressed; i++){
        for(int j =0; j<num_bins_total_compressed; j++){
            tchi += (core[i]-sig[i])*invert_matrix[i][j]*(core[j]-sig[j] );
            //if(i==j)std::cout << "1. invert_matrix: " << invert_matrix[i][j] << ", " << core[i] << ", " << sig[i] << ", " << core[j] << ", " << sig[j] << ", " << tchi << std::endl;
        }
    }
   // std::cout << "tchi = " << tchi << std::endl;

   return tchi;
}

float SBNchi::PoissonLogLiklihood(float * pred, float *data){
    float ans = 0;

    int n = num_bins_total_compressed;
    for(int i =0; i< n; i++){
        ans +=   ( data[i] >0.001 ? 2.0*(pred[i]-data[i]  + data[i]*log(data[i]/pred[i]) ) : 2*pred[i]);
    }
    return ans;
}


double SBNchi::CalcChi(double **invert_matrix, double* core, double *sig){
    double tchi = 0;

    for(int i =0; i<num_bins_total_compressed; i++){
        for(int j =0; j<num_bins_total_compressed; j++){
            if(i==j)std::cout << "2.invert_matrix: " << invert_matrix[i][j] << std::endl;
            tchi += (core[i]-sig[i])*invert_matrix[i][j]*(core[j]-sig[j]);
        }
    }

    return tchi;
}

double SBNchi::CalcChi(TMatrixT<double> M_invert, std::vector<double>& spec, std::vector<double>& data){

    double tchi = 0;
    for(int i=0; i< num_bins_total_compressed; i++){
        for(int j=0; j< num_bins_total_compressed ;j++){
            if(i==j)std::cout << "invert_matrix: " << M_invert[i][j] << std::endl;
            tchi += M_invert(i,j)*(spec[i]- data[i])*(spec[j]-data[j]);
        }
    }

    return tchi;
}


//same as above but passing in a vector instead of whole SBNspec
double SBNchi::CalcChi(std::vector<double> sigVec){
    double tchi = 0;

    if(sigVec.size() != num_bins_total_compressed ){
        std::cerr<<"ERROR: SBNchi::CalcChi(std::vector<double>) ~ your inputed vector does not have correct dimensions"<<std::endl;
        std::cerr<<"sigVec.size(): "<<sigVec.size()<<" num_bins_total_compressed: "<<num_bins_total_compressed<<std::endl;
        exit(EXIT_FAILURE);
    }

    for(int i =0; i<num_bins_total_compressed; i++){
        for(int j =0; j<num_bins_total_compressed; j++){
            tchi += (core_spectrum.collapsed_vector[i]-sigVec[i])*vec_matrix_inverted[i][j]*(core_spectrum.collapsed_vector[j]-sigVec[j] );
        }
    }

    last_calculated_chi = tchi;
    return tchi;
}

//A log-lilihood based one used @ MiniBooNE
double SBNchi::CalcChiLog(SBNspec *sigSpec){
    double tchi = 0;

    if(sigSpec->collapsed_vector.size()==0){
        std::cout<<"WARNING: SBNchi::CalcChi, inputted sigSpec has un-compressed vector, I am doing it now, but this is inefficient!"<<std::endl;
        sigSpec->CollapseVector();
    }

    for(int i =0; i<num_bins_total_compressed; i++){
        for(int j =0; j<num_bins_total_compressed; j++){
            vec_last_calculated_chi.at(i).at(j) =(core_spectrum.collapsed_vector[i]-sigSpec->collapsed_vector[i])*vec_matrix_inverted[i][j]*(core_spectrum.collapsed_vector[j]-sigSpec->collapsed_vector[j] );
            tchi += vec_last_calculated_chi.at(i).at(j);
        }
    }

    double absDetM = log(fabs(matrix_collapsed.Determinant()));

    last_calculated_chi = tchi+absDetM;
    return tchi+absDetM;
}


double SBNchi::CalcChi(SBNspec *sigSpec, SBNspec *obsSpec){
    double tchi=0;
    if(sigSpec->collapsed_vector.size()==0){
        std::cout<<"WARNING: SBNchi::CalcChi, inputted sigSpec has un-compressed vector, I am doing it now, but this is inefficient!"<<std::endl;
        sigSpec->CollapseVector();
    }

    if(obsSpec->collapsed_vector.size()==0){
        std::cout<<"WARNING: SBNchi::CalcChi, inputted obsSpec has un-compressed vector, I am doing it now, but this is inefficient!"<<std::endl;
        obsSpec->CollapseVector();
    }

    for(int i =0; i<num_bins_total_compressed; i++){
        for(int j =0; j<num_bins_total_compressed; j++){
            tchi += (obsSpec->collapsed_vector[i]-sigSpec->collapsed_vector[i])*vec_matrix_inverted[i][j]*(obsSpec->collapsed_vector[j]-sigSpec->collapsed_vector[j] );
        }
    }

    last_calculated_chi = tchi;
    return tchi;
}



/**************************************************************************
 *			Collapsing code
 * ************************************************************************/


//This is the powerhouse, takes each detector matrix filled with num_channels channels of num_subchannels[i] subchannels, and collapses it.
void SBNchi::CollapseSubchannels(TMatrixT <double> & M, TMatrixT <double> & Mc){
    bool debug = false;
    if(debug)	std::cout<<"Starting:M "<<M.GetNcols()<<" "<<M.GetNrows()<<" "<<115<<std::endl;
    if(debug)	std::cout<<"Starting:Mc "<<Mc.GetNcols()<<" "<<Mc.GetNrows()<<" "<<30<<std::endl;

    std::vector<std::vector<TMatrixT<double>>> Summed(num_channels, std::vector<TMatrixT<double>>(num_channels) );	//Initialise a matrix of matricies, to ZERO.
    for(int ic = 0; ic < num_channels; ic++){
        for(int jc =0; jc < num_channels; jc++){
            Summed[ic][jc].ResizeTo(num_bins[jc],num_bins[ic]) ;// This is CORRECT, do not switch (ie Summed[0][1] = size (num_bins[1], num_bins[0])
            Summed[ic][jc] = 0.0;
        }
    }

    int mrow = 0.0;
    int mcol = 0.0;

    for(int ic = 0; ic < num_channels; ic++){ 	 //Loop over all rows
        for(int jc =0; jc < num_channels; jc++){ //Loop over all columns

            if(debug)std::cout<<"Diagonal! : "<<ic<<" "<<jc<<" mcol is: "<<mcol<<" mrow is: "<<mrow<<std::endl;

            for(int m=0; m < num_subchannels[ic]; m++){
                for(int n=0; n< num_subchannels[jc]; n++){ //For each big block, loop over all subchannels summing toGether
                    Summed[ic][jc] +=  M.GetSub(mrow+n*num_bins[jc] ,mrow + n*num_bins[jc]+num_bins[jc]-1, mcol + m*num_bins[ic], mcol+ m*num_bins[ic]+num_bins[ic]-1 );
                }
            }


            mrow += num_subchannels[jc]*num_bins[jc];//As we work our way left in columns, add on that many bins
        }//end of column loop

        mrow = 0; // as we end this row, reSet row count, but jump down 1 column
        mcol += num_subchannels[ic]*num_bins[ic];
    }//end of row loop

    ///********************************* And put them back toGether! ************************//
    Mc.Zero();
    mrow = 0;
    mcol = 0;

    //Repeat again for Contracted matrix
    for(int ic = 0; ic < num_channels; ic++){
        for(int jc =0; jc < num_channels; jc++){

            Mc.SetSub(mrow,mcol,Summed[ic][jc]);
            mrow += num_bins[jc];
        }

        mrow = 0;
        mcol +=num_bins[ic];
    }

    return;
}



//This is the detector layer, Take a given mode and run over each detector V detector sub matrix
void SBNchi::CollapseDetectors(TMatrixT <double> & M, TMatrixT <double> & Mc){

    Mc.Zero();
    int nrow = num_bins_detector_block;// N_e_bins*N_e_spectra+N_m_bins*N_m_spectra;
    int crow = num_bins_detector_block_compressed; //N_e_bins+N_m_bins;

    for(int m =0; m< num_detectors; m++){
        for(int n =0; n< num_detectors; n++){
            TMatrixT<double> imat(nrow,nrow);
            TMatrixT<double> imatc(crow,crow);

            imat = M.GetSub(n*nrow,n*nrow+nrow-1, m*nrow,m*nrow+nrow-1);
            CollapseSubchannels(imat,imatc);
            Mc.SetSub(n*crow,m*crow,imatc);
        }
    }

    return;
}

//This is the Mode layer, Take a given full matrix and runs over each Mode V Mode sub matrix
void SBNchi::CollapseModes(TMatrixT <double> & M, TMatrixT <double> & Mc){

    Mc.Zero();
    int nrow = num_bins_mode_block;// (N_e_bins*N_e_spectra+N_m_bins*N_m_spectra)*N_dets;
    int crow=  num_bins_mode_block_compressed;// (N_e_bins+N_m_bins)*N_dets;

    for(int m =0; m< num_modes ; m++){
        for(int n =0; n< num_modes; n++){

            TMatrixT<double> imat(nrow,nrow);
            TMatrixT<double> imatc(crow,crow);

            imat = M.GetSub(n*nrow,n*nrow+nrow-1, m*nrow,m*nrow+nrow-1);

            CollapseDetectors(imat,imatc);
            Mc.SetSub(n*crow,m*crow,imatc);

        }
    }

    return;
}

TMatrixT<double> SBNchi::CalcCovarianceMatrix(TMatrixT<double>*M, std::vector<double>& spec, bool add_stats){

    TMatrixT<double> Mout(M->GetNcols(), M->GetNcols() );

    for(int i =0; i<M->GetNcols(); i++)
    {
        for(int j =0; j<M->GetNrows(); j++)
        {
            if(  std::isnan( (*M)(i,j) )){
                Mout(i,j) = 0.0;
            }else{

                Mout(i,j) = (*M)(i,j)*spec[i]*spec[j];
            }
            if(add_stats){  if(i==j) Mout(i,i) += spec[i]; }  //stats part
        }
    }
    return Mout;
}

TMatrixT<double> SBNchi::CalcCovarianceMatrix(TMatrixT<double>*M, std::vector<double>& spec, std::vector<double> &mcerr, bool add_stats){

    TMatrixT<double> Mout(M->GetNcols(), M->GetNcols() );

    for(int i =0; i<M->GetNcols(); i++)
    {
        for(int j =0; j<M->GetNrows(); j++)
        {
            if(  std::isnan( (*M)(i,j) )){
                Mout(i,j) = 0.0;
            }else{

                Mout(i,j) = (*M)(i,j)*spec[i]*spec[j];
                if(i==j) Mout(i,j) += mcerr[i]*mcerr[i]; 
            }
	    if(add_stats) {
               if(i==j) Mout(i,i) += spec[i];   //stats part
            }
        }
    }
    return Mout;
}
TMatrixT<double> SBNchi::CalcCovarianceMatrix(TMatrixT<double>*M, std::vector<double>& spec){

    TMatrixT<double> Mout(M->GetNcols(), M->GetNcols() );

    for(int i =0; i<M->GetNcols(); i++)
    {
        for(int j =0; j<M->GetNrows(); j++)
        {
            if(  std::isnan( (*M)(i,j) )){
                Mout(i,j) = 0.0;
            }else{

                Mout(i,j) = (*M)(i,j)*spec[i]*spec[j];
            }
            if(i==j) Mout(i,i) += spec[i];   //stats part
        }
    }
    return Mout;
}
//here spec is full vector of MC, spec_collapse is collapsed vector of MC, datavec is collapsed vector of data
TMatrixT<double> SBNchi::CalcCovarianceMatrixCNP(TMatrixT<double> M, std::vector<double>& spec, std::vector<double>& spec_collapse, const std::vector<double>& datavec ){

    if(M.GetNcols() != spec.size()){
        std::cout << "ERROR: your input vector does not have the right dimenstion  " << std::endl; 
        std::cout << "Fractional Matrix size :"<< M.GetNcols() << " || Input Full Vector size "<< spec.size() << std::endl;  
        exit(EXIT_FAILURE);
    }

    TMatrixT<double> M_temp(M.GetNcols(), M.GetNcols() );
    TMatrixT<double> Mout(spec_collapse.size(), spec_collapse.size()); //collapsed covariance matrix

    //systematic apart 
    for(int i =0; i<M.GetNcols(); i++)
    {
        for(int j =0; j<M.GetNrows(); j++)
        {
            if(  std::isnan( M(i,j) )){
                M_temp(i,j) = 0.0;
            }else{

                M_temp(i,j) = M(i,j)*spec[i]*spec[j];
            }
        }
    }

    CollapseModes(M_temp, Mout);
    //add stats part	
    for(int i=0; i< spec_collapse.size(); i++){
        Mout(i,i) +=   ( datavec[i] >0.001 ? 3.0/(1.0/datavec[i] +  2.0/spec_collapse[i])  : spec_collapse[i]/2.0 ); 
        //Mout(i,i) +=   spec_collapse[i];//( datavec[i] >0.001 ? 3.0/(1.0/datavec[i] +  2.0/spec_collapse[i])  : spec_collapse[i]/2.0 ); 
    }
    return Mout;
}


//here spec is full vector of MC, spec_collapse is collapsed vector of MC, datavec is collapsed vector of data
TMatrixT<double> SBNchi::CalcCovarianceMatrixCNP(TMatrixT<double> *M, std::vector<double>& spec, std::vector<double>& spec_collapse, std::vector<double>& spec_mcerr, const std::vector<float>& datavec, bool add_stats ){

    std::cout << "am using covariance matrix with mc err" << std::endl;
    if(M->GetNcols() != spec.size()){
        std::cout << "ERROR: your input vector does not have the right dimenstion  " << std::endl; 
        std::cout << "Fractional Matrix size :"<< M->GetNcols() << " || Input Full Vector size "<< spec.size() << std::endl;  
        exit(EXIT_FAILURE);
    }

    TMatrixT<double> M_temp(M->GetNcols(), M->GetNcols() );
    TMatrixT<double> Mout(spec_collapse.size(), spec_collapse.size()); //collapsed covariance matrix

    //systematic apart 
    for(int i =0; i<M->GetNcols(); i++)
    {
        for(int j =0; j<M->GetNrows(); j++)
        {
            if(  std::isnan( (*M)(i,j) )){
                M_temp(i,j) = 0.0;
            }else{

                M_temp(i,j) = (*M)(i,j)*spec[i]*spec[j];
                if(i==j) M_temp(i,j) += spec_mcerr[i]*spec_mcerr[i];
            }
        }
    }

    CollapseModes(M_temp, Mout);
    //add stats part
    if(add_stats){	
    for(int i=0; i< spec_collapse.size(); i++){
        Mout(i,i) +=   ( datavec[i] >0.001 ? 3.0/(1.0/datavec[i] +  2.0/spec_collapse[i])  : spec_collapse[i]/2.0 ); 
        //Mout(i,i) +=   spec_collapse[i];//( datavec[i] >0.001 ? 3.0/(1.0/datavec[i] +  2.0/spec_collapse[i])  : spec_collapse[i]/2.0 ); 
    }
    }
    return Mout;
}


TMatrixT<double> SBNchi::CalcCovarianceMatrixCNP(TMatrixT<double>*M, std::vector<double>& spec, const std::vector<float>& datavec ){

    TMatrixT<double> Mout(M->GetNcols(), M->GetNcols() );

    for(int i =0; i<M->GetNcols(); i++)
    {
        for(int j =0; j<M->GetNrows(); j++)
        {
            if(std::isnan( (*M)(i,j) )){
                Mout(i,j) = 0.0;
            }else{
                Mout(i,j) = (*M)(i,j)*spec[i]*spec[j];
            }
            if(i==j) Mout(i,i) +=   ( datavec[i] >0.001 ? 3.0/(1.0/datavec[i] +  2.0/spec[i])  : spec[i]/2.0 );
        }
    }
    return Mout;
}

float SBNchi::CalcChi_CNP(float * pred, float* data){

    is_verbose = false;
    TMatrixT<double> inverse_collapsed = m_matrix_systematics_collapsed;
    //Add on CNP terms proportional to data
    for(int j =0; j<num_bins_total_compressed; j++)
    {
        inverse_collapsed(j,j) +=  ( data[j] >0.001 ? 3.0/(1.0/data[j] +  2.0/pred[j])  : pred[j]/2.0 );
    }

    inverse_collapsed = this->InvertMatrix(inverse_collapsed);   

    float tchi = 0.0;
    for(int i =0; i<num_bins_total_compressed; i++){
        for(int j =0; j<num_bins_total_compressed; j++){
            tchi += (pred[i]-data[i])*inverse_collapsed(i,j)*(pred[j]-data[j]);
        }
    }

    return tchi;
}

TMatrixT<double> SBNchi::CalcCovarianceMatrix(TMatrixT<double>*M, TVectorT<double>& spec, TVectorT<double> &err){

    TMatrixT<double> Mout( M->GetNcols(), M->GetNcols() );
    // systematics per scaled event
    for(int i =0; i<M->GetNcols(); i++)
    {
        //std::cout<<"KRAK: "<<core_spectrum.full_vector.at(i)<<std::endl;
        for(int j =0; j<M->GetNrows(); j++)
        {
            if(  std::isnan( (*M)(i,j))){
                Mout(i,j) = 0.0;
            }else{
                Mout(i,j) = (*M)(i,j)*spec(i)*spec(j);
            }
        }
    }
    return Mout;
}




TMatrixT<double> SBNchi::CalcCovarianceMatrix(TMatrixT<double>*M, TVectorT<double>& spec){

    TMatrixT<double> Mout( M->GetNcols(), M->GetNcols() );
    // systematics per scaled event
    for(int i =0; i<M->GetNcols(); i++)
    {
        //std::cout<<"KRAK: "<<core_spectrum.full_vector.at(i)<<std::endl;
        for(int j =0; j<M->GetNrows(); j++)
        {
            if(  std::isnan( (*M)(i,j))){
                Mout(i,j) = 0.0;
            }else{
                Mout(i,j) = (*M)(i,j)*spec(i)*spec(j);
            }
            if(i==j) Mout(i,i) +=spec(i);
        }
    }
    return Mout;
}



TMatrixT<double> SBNchi::CalcCovarianceMatrix(TMatrixT<double>*M, TVectorT<double>& spec,bool add_stats){

    TMatrixT<double> Mout( M->GetNcols(), M->GetNcols() );
    // systematics per scaled event
    for(int i =0; i<M->GetNcols(); i++)
    {
        //std::cout<<"KRAK: "<<core_spectrum.full_vector.at(i)<<std::endl;
        for(int j =0; j<M->GetNrows(); j++)
        {
            if(  std::isnan( (*M)(i,j))){
                Mout(i,j) = 0.0;
            }else{
                Mout(i,j) = (*M)(i,j)*spec(i)*spec(j);
            }
            if(add_stats){           if(i==j) Mout(i,i) +=spec(i);}
        }
    }
    return Mout;
}


TMatrixT<double> SBNchi::InvertMatrix(TMatrixT<double> &M){

    double invdet=0;

    TMatrixT<double> McI(M.GetNrows(),M.GetNrows());
    McI.Zero();

    if(is_verbose) std::cout<<otag<<" About to do a SVD decomposition"<<std::endl;
    TDecompSVD svd(M);

    if (!svd.Decompose()){
        std::cout<<otag<<" (InvertMatrix) Decomposition failed, matrix not symettric?, has nans?" << std::endl;
        std::cout<<otag<<"ERROR: The matrix to invert failed a SVD decomp!"<<std::endl;

        for(int i=0; i< M.GetNrows(); i++){
            for(int j=0; j< M.GetNrows(); j++){
                std::cout<<i<<" "<<j<<" "<<M(i,j)<<std::endl;
            }
        }

        exit(EXIT_FAILURE);

    } else {
        McI = svd.Invert();
    }
    if( !McI.IsValid()){
        std::cout<<otag<<"ERROR: The inverted matrix isnt valid! Something went wrong.."<<std::endl;
        exit(EXIT_FAILURE);

    }


    return McI;

}



/**************************************************************************
 *			Misc
 * ************************************************************************/

int SBNchi::FillCovarianceMatrix(TMatrixT<double>*in){
    in->ResizeTo(num_bins_total,num_bins_total);
    for(int i=0; i<num_bins_total;i++){
        for(int j=0; j<num_bins_total;j++){
            (*in)(i,j) = matrix_systematics(i,j);
        }
    }

    return 0;
}


int SBNchi::FillCollapsedCovarianceMatrix(TMatrixT<double>*in){
    in->ResizeTo(num_bins_total_compressed,num_bins_total_compressed) ;
    for(int i=0; i<num_bins_total_compressed;i++){
        for(int j=0; j<num_bins_total_compressed;j++){
            (*in)(i,j) = vec_matrix_collapsed.at(i).at(j) - (i==j? core_spectrum.collapsed_vector.at(i) : 0.0 );
        }
    }

    return 0;
}


int SBNchi::FillCorrelationMatrix(TMatrixT<double>*in){
    in->ResizeTo(num_bins_total,num_bins_total) ;
    for(int i=0; i<num_bins_total;i++){
        for(int j=0; j<num_bins_total;j++){
            double val = matrix_systematics(i,j);
            double p1 = sqrt(matrix_systematics(j,j)); 
            double p2 = sqrt(matrix_systematics(i,i)); 
            (*in)(i,j) = val/(p1*p2);
        }
    }
    return 0;
}

int SBNchi::FillCollapsedCorrelationMatrix(TMatrixT<double>*in){
    in->ResizeTo(num_bins_total_compressed,num_bins_total_compressed) ;
    for(int i=0; i<num_bins_total_compressed;i++){
        for(int j=0; j<num_bins_total_compressed;j++){
            double val = vec_matrix_collapsed.at(i).at(j) - (i==j? core_spectrum.collapsed_vector.at(i) : 0.0 ); 
            double p1 = sqrt(vec_matrix_collapsed.at(j).at(j)- core_spectrum.collapsed_vector.at(j)); 
            double p2 = sqrt(vec_matrix_collapsed.at(i).at(i)- core_spectrum.collapsed_vector.at(i)); 
            (*in)(i,j) = val/(p1*p2);
        }
    }

    return 0;
}


int SBNchi::FillFractionalMatrix(TMatrixT<double>*in){
    in->ResizeTo(num_bins_total,num_bins_total) ;
    for(int i=0; i<num_bins_total;i++){
        for(int j=0; j<num_bins_total;j++){
            (*in)(i,j) = ( matrix_systematics(i,j))/(core_spectrum.full_vector.at(i)*core_spectrum.full_vector.at(j));
        }
    }

    return 0;
}


int SBNchi::FillCollapsedFractionalMatrix(TMatrixT<double>*in){
    in->ResizeTo(num_bins_total_compressed,num_bins_total_compressed) ;
    for(int i=0; i<num_bins_total_compressed;i++){
        for(int j=0; j<num_bins_total_compressed;j++){
            (*in)(i,j) = ( vec_matrix_collapsed.at(i).at(j) - (i==j? core_spectrum.collapsed_vector.at(i) : 0.0 ) )/(core_spectrum.collapsed_vector.at(i)*core_spectrum.collapsed_vector.at(j));
        }
    }

    return 0;
}


TMatrixT<double> * SBNchi::GetCollapsedMatrix(){
    TMatrixT<double> * tmp = new TMatrixT<double>(num_bins_total_compressed,num_bins_total_compressed);
    for(int i=0; i<num_bins_total_compressed;i++){
        for(int j=0; j<num_bins_total_compressed;j++){
            (*tmp)(i,j) = vec_matrix_collapsed.at(i).at(j);
        }
    }

    return tmp;
}




void SBNchi::FakeFillMatrix(TMatrixT <double> &M){
    //Fills a square matrix of dim matrix_size with random numbers for now.
    std::uniform_real_distribution<double> dist(0,1);

    int matrix_size=M.GetNrows();
    if(M.GetNrows()!=M.GetNcols()){std::cout<<"#ERROR: not a square matrix!"<<std::endl;}
    for(int i=0; i<matrix_size; i++){
        for (int j = i;j<matrix_size;j++){
            M(i,j)= dist(*rangen_twister);
            M(j,i)=M(i,j);
        }
    }
    return ;
}


std::vector<std::vector<double >> SBNchi::TMatrixDToVector(TMatrixT <double > Min)
{
    int dimension =  Min.GetNrows();

    std::vector<std::vector<double >>  ans(dimension, std::vector<double>(dimension));

    for(int i = 0; i< dimension; i++){
        for(int k = 0; k< dimension; k++){
            ans[i][k]=Min(i,k);
            if(ans[i][k]==-0){
                ans[i][k]=0;
            }
        }
    }
    return ans;
}

/*void SBNchi::FillStatsMatrix(TMatrixT <double> &M, std::vector<double> diag, std::vector<double> add_stats_err){
    int matrix_size = M.GetNrows();

    std::cout << "matrix size, diag size = " << matrix_size << ", " << diag.size() << ", " << add_stats_err.size() << std::endl;
    if(matrix_size != diag.size()){std::cout<<"#ERROR: FillStatsMatrix, matrix not equal to diagonal"<<std::endl;}
    if(M.GetNrows()!=M.GetNcols()){std::cout<<"#ERROR: not a square matrix!"<<std::endl;}

    M.Zero();

    for(int i=0; i<matrix_size; i++)
    {
        std::cout << "add stats err = " << add_stats_err.at(i) << std::endl;
        M(i,i) = diag.at(i) + (add_stats_err.at(i)*add_stats_err.at(i));
    }

    return ;
}*/

void SBNchi::FillStatsMatrix(TMatrixT <double> &M, std::vector<double> diag ){
    int matrix_size = M.GetNrows();

    std::cout << "matrix size, diag size = " << matrix_size << ", " << diag.size() << std::endl;
    if(matrix_size != diag.size()){std::cout<<"#ERROR: FillStatsMatrix, matrix not equal to diagonal"<<std::endl;}
    if(M.GetNrows()!=M.GetNcols()){std::cout<<"#ERROR: not a square matrix!"<<std::endl;}

    M.Zero();

    for(int i=0; i<matrix_size; i++)
    {
        M(i,i) = diag.at(i);
        std::cout << "stats error = " <<  M(i,i) << std::endl;
    }

    return ;
}



TMatrixT<double> SBNchi::FillSystematicsFromXML(){
    return FillSystematicsFromXML(correlation_matrix_rootfile, correlation_matrix_name);
}


TMatrixT<double > SBNchi::FillSystematicsFromXML(std::string rootname, std::string matname){
    //Pretty much obsolete now, should fill directly really.
    std::cout<<"SBNchi::FillSystematicsFromXML || filling from "<<rootname<<std::endl;

    TMatrixT<double> temp2(num_bins_total,num_bins_total);
    TFile *fm= new TFile(rootname.c_str());

    TMatrixT<float> * temp = (TMatrixT <float>* )fm->Get(matname.c_str());
    //TMatrixT<double> * temp = (TMatrixT <double>* )fm->Get(matname.c_str());

    std::vector<std::vector<double>> mcont;

    for(int p:used_bins){
        std::vector<double> tvec;
        for(int u:used_bins){
            tvec.push_back( (*temp)(p,u) );
        }
        mcont.push_back(tvec);
    }

    for(int i =0; i<num_bins_total; i++)
    {
        for(int j =0; j<num_bins_total; j++)
        {
            temp2(i,j)=mcont[i][j];
        }
    }
    delete temp;

    std::cout<<"SBNchi::FillSystematicsFromXML || loaded with dim : "<<temp2.GetNcols()<<" "<<temp2.GetNrows()<<std::endl;

    fm->Close();
    delete fm;

    if(temp2.IsSymmetric()){
        if(is_verbose)std::cout<<"Inputted fracCov covariance matrix is symmetric"<<std::endl;
    }else{
        std::cerr<<"ERROR: SBNchi::FillSystematicsFromXML, matrix_systematics input is not symmetric!"<<std::endl;
        //exit(EXIT_FAILURE);
    }

    return temp2;

}





TH2D* SBNchi::GetChiogram(){
    TH2D *tmp = new TH2D("chi-o-gram","chi-o-gram",num_bins_total_compressed,0, num_bins_total_compressed ,num_bins_total_compressed,0, num_bins_total_compressed);

    for(int i =0; i<num_bins_total_compressed; i++){
        for(int j =0; j<num_bins_total_compressed; j++){

            tmp->SetBinContent(i+1, j+1, vec_last_calculated_chi.at(i).at(j));
        }
    }

    return tmp;
}

int SBNchi::PrintMatricies(std::string tag){
    TFile* fout = new TFile(("SBNfit_collapsed_matrix_plots_"+tag+".root").c_str(),"recreate");
    fout->cd();

    gStyle->SetOptStat(0);

    TMatrixD full, frac, corr;
    this->FillCollapsedCovarianceMatrix(&full);
    this->FillCollapsedFractionalMatrix(&frac);
    this->FillCollapsedCorrelationMatrix(&corr);

    corr.Write("collapsed_correlation");

    gStyle->SetPalette(kLightTemperature);

    frac.Write("collapsed_fractional_covariance");
    TH2D h2_frac(frac);
    //h2_frac.Write();
    h2_frac.SetName("frac");
    TCanvas *c_frac = new TCanvas("collapsed fractional covariance matrix");
    TPad* p_frac = (TPad*)c_frac->cd();
    //p_frac->SetLogz();
    c_frac->SetFixedAspectRatio();
    h2_frac.Draw("colz");
    h2_frac.SetTitle("Collapsed fractional covariance matrix");
    h2_frac.GetXaxis()->SetTitle("Reco Bin i");
    h2_frac.GetYaxis()->SetTitle("Reco Bin j");
    if(m_cmin !=-999 || m_cmax !=-999)   h2_frac.GetZaxis()->SetRangeUser(m_cmin,m_cmax);

    c_frac->SetRightMargin(0.150);

    int use_frac =0;

    for(int im =0; im<num_modes; im++){
        for(int id =0; id<num_detectors; id++){
            for(int ic = 0; ic < num_channels; ic++){
                TLine *lv = new TLine(0, num_bins.at(ic)+use_frac, num_bins_total_compressed, num_bins.at(ic)+use_frac);
                TLine *lh = new TLine(num_bins.at(ic)+use_frac,0, num_bins.at(ic)+use_frac, num_bins_total_compressed);
                lv->SetLineWidth(1.5);
                lh->SetLineWidth(1.5);
                use_frac+=num_bins.at(ic);
                lv->Draw();
                lh->Draw();
            }
        }
    }
    c_frac->Write();
    c_frac->SaveAs(("SBNfit_collapsed_fractional_covariance_"+tag+".SBNplot.pdf").c_str(),"pdf");

    for(int i=0; i<h2_frac.GetNbinsX(); i++){
        std::cout<<"Collapsed Frac "<<i<<" "<<sqrt(h2_frac.GetBinContent(i+1,i+1))*100.0<<std::endl;
    }

    for(int i=0; i<core_spectrum.collapsed_vector.size(); i++){
        std::cout<<sqrt(core_spectrum.collapsed_vector.at(i))/core_spectrum.collapsed_vector.at(i)*100.0<<" ";
    }std::cout<<std::endl;


    gStyle->SetPalette(kLightTemperature);

    full.Write("collapsed_covariance");
    TH2D h2_full(full);
    h2_full.SetName("full");
    TCanvas *c_full = new TCanvas("collapsed covariance matrix");
    c_full->cd();
    c_full->SetFixedAspectRatio();
    h2_full.Draw("colz");
    h2_full.SetTitle("Collapsed covariance matrix");
    h2_full.GetXaxis()->SetTitle("Reco Bin i");
    h2_full.GetYaxis()->SetTitle("Reco Bin j");

    c_full->SetRightMargin(0.150);

    int use_full =0;
    for(int im =0; im<num_modes; im++){
        for(int id =0; id<num_detectors; id++){
            for(int ic = 0; ic < num_channels; ic++){
                TLine *lv = new TLine(0, num_bins.at(ic)+use_full, num_bins_total_compressed, num_bins.at(ic)+use_full);
                TLine *lh = new TLine(num_bins.at(ic)+use_full,0, num_bins.at(ic)+use_full, num_bins_total_compressed);
                lv->SetLineWidth(1.5);
                lh->SetLineWidth(1.5);
                use_full+=num_bins.at(ic);
                lv->Draw();
                lh->Draw();
            }
        }
    }
    c_full->Write();
    c_full->SaveAs(("SBNfit_collapsed_covariance_"+tag+".SBNplot.pdf").c_str(),"pdf");

    gStyle->SetPalette(kLightTemperature);
    TH2D h2_corr(corr);
    h2_corr.SetName("corr");
    //h2_corr.Write();
    TCanvas *c_corr = new TCanvas("collapsed correlation matrix");
    c_corr->cd();
    c_corr->SetFixedAspectRatio();
    h2_corr.Draw("colz");
    h2_corr.SetTitle("Collapsed correlation matrix");
    h2_corr.GetXaxis()->SetTitle("Reco Bin i");
    h2_corr.GetYaxis()->SetTitle("Reco Bin j");
    h2_corr.GetZaxis()->SetRangeUser(0.0,1);
    c_corr->SetRightMargin(0.150);

    int use_corr =0;
    for(int im =0; im<num_modes; im++){
        for(int id =0; id<num_detectors; id++){
            for(int ic = 0; ic < num_channels; ic++){
                TLine *lv = new TLine(0, num_bins.at(ic)+use_corr, num_bins_total_compressed, num_bins.at(ic)+use_corr);
                TLine *lh = new TLine(num_bins.at(ic)+use_corr,0, num_bins.at(ic)+use_corr, num_bins_total_compressed);
                lv->SetLineWidth(1.5);
                lh->SetLineWidth(1.5);
                use_corr+=num_bins.at(ic);
                lv->Draw();
                lh->Draw();
            }
        }
    }
    c_corr->Write();
    c_corr->SaveAs(("SBNfit_collapsed_correlation_"+tag+".SBNplot.pdf").c_str(),"pdf");



    TCanvas *chiogram = new TCanvas("Chi-o-gram","Chi-o-gram");
    chiogram->cd();
    chiogram->SetFixedAspectRatio();

    TH2D * h_chiogram = (TH2D*)this->GetChiogram();
    h_chiogram->Draw("colz");
    h_chiogram->SetTitle("Reco Bin i");
    h_chiogram->SetTitle("Reco Bin j");

    int use_chio =0;
    for(int im =0; im<num_modes; im++){
        for(int id =0; id<num_detectors; id++){
            for(int ic = 0; ic < num_channels; ic++){
                TLine *lv = new TLine(0, num_bins.at(ic)+use_chio, num_bins_total_compressed, num_bins.at(ic)+use_chio);
                TLine *lh = new TLine(num_bins.at(ic)+use_chio,0, num_bins.at(ic)+use_chio, num_bins_total_compressed);
                lv->SetLineWidth(1.5);
                lh->SetLineWidth(1.5);
                use_chio+=num_bins.at(ic);
                lv->Draw();
                lh->Draw();
            }
        }
    }
    chiogram->Write();

    gStyle->SetPalette(kBird);

    TMatrixD ufull, ufrac, ucorr;
    this->FillCovarianceMatrix(&ufull);
    this->FillFractionalMatrix(&ufrac);
    this->FillCorrelationMatrix(&ucorr);

    plot_one(ufull,"SBNfit_uncollapsed_covariance_matrix_"+tag,fout,true,false,false);
    plot_one(ufrac,"SBNfit_uncollapsed_fractional_covariance_matrix"+tag,fout,true,false,false);
    plot_one(ucorr,"SBNfit_uncollapsed_correlation_matrix"+tag,fout,true,false,true);

    fout->Close();
    return 0;
}
int SBNchi::plot_one(TMatrixD matrix, std::string tag, TFile *fin, bool plot_pdf, bool indiv, bool is_corr){
    fin->cd();
    if(indiv){
        TDirectory *individualDir = fin->GetDirectory("individualDir"); 
        if (!individualDir) { 
            individualDir = fin->mkdir("individualDir");       
        }
        fin->cd(); 
        individualDir->cd();
    }
    if(is_corr) gStyle->SetPalette(kLightTemperature);


    if(is_corr){
        for(int i=0; i<matrix.GetNrows();i++){
            for(int j=0; j<matrix.GetNrows();j++){
                double val = matrix(i,j);
                if(val!=val || isinf(val) || std::isnan(val)){
                    matrix(i,j)=0.0;
                    matrix(i,i)=1.0;
                    matrix(j,j)=1.0;
                }
            }

        }
    }

    TH2D h2_full(matrix);
    h2_full.SetName((tag+"_th2d").c_str());
    TCanvas *c_full = new TCanvas((tag+"_canvas").c_str());
    TPad *p_full = (TPad*)c_full->cd();
    c_full->SetFixedAspectRatio();
    h2_full.Draw("colz");
    h2_full.SetTitle(tag.c_str());
    h2_full.GetXaxis()->SetTitle("Global Bin Number");
    h2_full.GetYaxis()->SetTitle(" ");
    h2_full.GetYaxis()->SetLabelSize(0);
    //p_full->SetLogz();
    if(is_corr){
        h2_full.GetZaxis()->SetRangeUser(0.4,1);
    }
    else{
        h2_full.GetZaxis()->SetRangeUser(-0.25,0.25);
    }    

    c_full->SetFrameFillColor(kWhite);
    c_full->SetFillColor(kWhite);
    p_full->SetFillColor(kWhite);


    c_full->SetRightMargin(0.150);
    c_full->SetLeftMargin(0.10);//0.250
    c_full->SetTopMargin(0.10);
    int use_full =0;

    double percent_left = 0.05;//0.15
    double nice_shift = num_bins_total*0.02;

    for(int im =0; im<num_modes; im++){
        for(int id =0; id<num_detectors; id++){
            for(int ic = 0; ic < num_channels; ic++){
                for(int isc = 0; isc < num_subchannels.at(ic); isc++){


                    std::string mode_det = mode_names[im] +" " +detector_names[id];
                    std::string chan_sub = channel_names[ic]+" "+subchannel_names[ic][isc];


                    TText * tmd = new TText(-num_bins_total*percent_left*0.15, use_full+nice_shift*0.5, (mode_det+" "+chan_sub).c_str() );

                    //TText * tmd = new TText(use_full*1.05, num_bins_total*1.015, chan_sub.c_str());
                    //TText * tcs = new TText(use_full*1.05, num_bins_total*1.055, mode_det.c_str());
                    tmd->SetTextColor(kBlack);
                    //tcs->SetTextColor(kBlack);
                    tmd->SetTextSize(0.03);
                    tmd->SetTextAlign(31);
                    //tcs->SetTextSize(0.03);

                    //dont plot names ta the moment
                    //tmd->Draw();
                    //tcs->Draw();


                    /*
                       TText * tlow_bin = new TText(-num_bins_total*percent_left, use_full+nice_shift*0.5, to_string_prec(bin_edges[ic].front(),0).c_str());
                       TText * thigh_bin = new TText(-num_bins_total*percent_left, (use_full+num_bins[ic])-nice_shift*1.4, to_string_prec(bin_edges[ic].back(),0).c_str());
                       tlow_bin->SetTextSize(0.02);
                       thigh_bin->SetTextSize(0.02);
                       tlow_bin->Draw();
                       thigh_bin->Draw();

                       TText * tunit = new TText(-num_bins_total*percent_left, use_full+0.5*num_bins[ic], channel_units[ic].c_str());
                       tunit->SetTextSize(0.03);
                       tunit->Draw();
                       */

                    if(isc<num_subchannels[ic]-1){
                        TLine *lscv = new TLine(-num_bins_total*percent_left, num_bins.at(ic)+use_full, num_bins_total, num_bins.at(ic)+use_full);
                        TLine *lsch = new TLine(num_bins.at(ic)+use_full,0, num_bins.at(ic)+use_full, num_bins_total*1.045);
                        lscv->SetLineWidth(3);
                        lsch->SetLineWidth(3);
                        lscv->SetLineColor(kRed);
                        lsch->SetLineColor(kRed);
                        lscv->SetLineStyle(9);
                        lsch->SetLineStyle(9);

                        //Going to drop the little ones for now
                        //lscv->Draw();
                        //lsch->Draw();

                        use_full+=num_bins.at(ic);

                    }
                }
                TLine *lv = new TLine(-num_bins_total*percent_left, num_bins.at(ic)+use_full, num_bins_total, num_bins.at(ic)+use_full);
                TLine *lh = new TLine(num_bins.at(ic)+use_full,0, num_bins.at(ic)+use_full, num_bins_total*1.045);
                lv->SetLineWidth(2);
                lh->SetLineWidth(2);
                lv->SetLineColor(kBlack);
                lh->SetLineColor(kBlack);
                use_full+=num_bins.at(ic);
                lv->Draw();
                lh->Draw();

            }
        }
    }


    c_full->Write();
    if(plot_pdf) c_full->SaveAs((tag+".pdf").c_str(),"pdf");




    return 0;
}


int SBNchi::PerformCholoskyDecomposition(SBNspec *specin){
    specin->CalcFullVector();
    is_verbose=false;
    double tol = m_tolerance;

    std::cout<<" Starting Cholosky Decomp, tolderance is "<<tol<<" and sample_collapse "<<pseudo_from_collapsed<<std::endl;

    TMatrixD U  = matrix_fractional_covariance;

    for(int i =0; i<U.GetNcols(); i++)
    {
        for(int j =0; j<U.GetNrows(); j++)
        {
            if(std::isnan(U(i,j)))
                U(i,j) = 0;
            else
                U(i,j)=U(i,j)*specin->full_vector.at(i)*specin->full_vector.at(j);

            //Comment this in if you want to sample from full gaussian
            //if(i==j)U(i,i)+=specin->full_vector.at(i);
        }
    }

    //New bit, do we collapse and sample from collapsed or full! Debugging April2020 Collab Meeting @ Zarkos info from MiniBooNE ERA
    int n_t = (pseudo_from_collapsed ? num_bins_total_compressed : num_bins_total);
    TMatrixT<double > U_use(n_t,n_t);
    if(pseudo_from_collapsed){
        CollapseModes(U, U_use);
    }else{
        U_use = U;
    }

    //Stats error are NOT added back in herebut treat them as Poisson later. Seems better
    //TMatrixDEigen eigen (U); // This was original, but caused a lot of "Error in <MakeSchurr>: too many iterations". Move to explicit symmetric matrix
    //Is this Really the best way to construct?!?
    // TMatrixDSym U_explicit_sym(n_t);
    // for(int i=0; i< n_t;i++){
    //   for(int j=i; j< n_t;j++){
    //     U_explicit_sym[i][j] = U(i,j);
    // }
    // }

    bool was_modified = false;

    //Get a Determinant, check things.
    double det = U_use.Determinant(); 
    std::cout<<"SBNchi::CholeskyDecomposition\t||\t Checking determinant of Covariance Matrix: U "<<det<<std::endl; 

    if(det < tol){
        std::cout<<"SBNchi::CholeskyDecomposition\t||\t This determinant is below tolerance of : "<<m_tolerance<<std::endl;
        std::cout<<"SBNchi::CholeskyDecomposition\t||\t Going to add this back to diagonal of full covariance matrix : "<<tol<<std::endl;
        for(int a =0; a<U_use.GetNcols(); a++){
            U_use(a,a) += tol;
            was_modified = true;
        }
        det = U_use.Determinant(); 
        std::cout<<"SBNchi::CholeskyDecomposition\t||\t The modified determinant is now: "<<det<<std::endl;
    }


    TMatrixDEigen eigen(U_use);
    TVectorD eigen_values = eigen.GetEigenValuesRe();
    TVectorD eigen_values_IM = eigen.GetEigenValuesIm();

    int n_zeros = 0;
    for(int i=0; i< eigen_values.GetNoElements(); i++){
        std::cout<<"SBNchi::CholeskyDecomposition\t||\t Eigenvalue "<<i<<" is Re: "<<eigen_values(i)<<" Im: "<<eigen_values_IM(i)<<std::endl;
        if(eigen_values(i)<tol){
            if(fabs(eigen_values(i)) > 0){
                if(is_verbose)std::cout<<"SBNchi::CholeskyDecomposition\t|| cov has a very small, < "<<tol<<" , negative eigenvalue. Adding it back to diagonal of : "<<eigen_values(i)<<std::endl;

                for(int a =0; a<U_use.GetNcols(); a++){
                    U_use(a,a) += eigen_values(i);
                    was_modified = true;
                }

            }else{
                std::cout<<"SBNchi::SampleCovariance\t|| 0 or negative eigenvalues! error: Value "<<eigen_values(i)<<" Tolerence "<<tol<<std::endl;
                //U_use.Print();
                //std::cout<<"Hmm"<<std::endl;
                //exit(EXIT_FAILURE);
            }
        }
    }


    //If everything is OK, lets pass this matrix back to SBNchi for use.
    if(was_modified){
        std::cout<<"We had to add on small diagonal terms to covariance matrix to Decompose it. Adding back to primary fractional covariance for consistency"<<std::endl;
        std::cout<<"This potentially causes an infinite loop. Check"<<std::endl;
        for(int i=0; i< n_t; i++){
            for(int j=0; j< n_t; j++){
                double mi = specin->full_vector.at(i)*specin->full_vector.at(j);
                matrix_fractional_covariance(i,j) = (mi==0 ? 0 : U_use(i,j)/(mi));
            }
        }
        this->ReloadCoreSpectrum(specin);
    }


    //Seconndly attempt a Cholosky Decomposition
    TDecompChol * chol = new TDecompChol(U_use,tol);
    bool worked = chol->Decompose();

    if(!worked){
        std::cout<<"SBNchi::SampleCovariance\t|| Cholosky Decomposition Failed. Tolerance is set at "<<tol<<std::endl;
        //U_use.Print();
        std::cout<<"With eigens"<<std::endl;

        for(int i=0; i< eigen_values.GetNoElements(); i++){
            std::cout<<eigen_values(i)<<std::endl;
        }

        exit(EXIT_FAILURE);
    }

    TMatrixT<float> upper_trian(n_t,n_t);
    matrix_lower_triangular.ResizeTo(n_t,n_t);
    upper_trian = chol->GetU();
    matrix_lower_triangular = upper_trian;
    matrix_lower_triangular.T();

    vec_matrix_lower_triangular.resize(n_t, std::vector<float>(n_t));
    for(int i=0; i< n_t; i++){
        for(int j=0; j< n_t; j++){
            vec_matrix_lower_triangular[i][j] = matrix_lower_triangular[i][j];
            //            std::cout<<"Flormph "<<i<<" "<<j<<" "<<vec_matrix_lower_triangular[i][j]<<" "<<U_use(i,j)<<std::endl;
        }
    }


    //New Check, rebuild matrix and compare
    TMatrixDSym rebuild = chol->GetMatrix();
    TMatrixD rebuild2 = matrix_lower_triangular*upper_trian;

    for(int i=0; i< n_t; i++){
        for(int j=0; j< n_t; j++){
            double fd1 = fabs(rebuild(i,j)-U_use(i,j))/fabs(rebuild(i,j)+U_use(i,j));
            double fd2 = fabs(rebuild2(i,j)-U_use(i,j))/fabs(rebuild2(i,j)+U_use(i,j));
            if( fd1>m_tolerance || fd2>m_tolerance){
                //std::cout<<"ERROR the rebuilt matrix after Cholesky Decomp is not the same as the original."<<std::endl;
                //std::cout<<i<<" "<<j<<" Original: "<<U_use(i,j)<<" , Rebuild: "<<rebuild(i,j)<<" , Diff "<<fd1<<std::endl;
                //std::cout<<i<<" "<<j<<" Original: "<<U_use(i,j)<<" , Rebuild2: "<<rebuild2(i,j)<<" , Diff "<<fd2<<std::endl;
                //exit(EXIT_FAILURE);
            }
        }
    }
    std::cout<<"Rebuilt matrix is the same as input after Cholesky Decomp"<<std::endl;



    cholosky_performed = true;	
    delete chol;
    return 0;
}


TH1D SBNchi::SampleCovarianceVaryInput(SBNspec *specin, int num_MC, double maxchi){ 
    max_sample_chi_val = maxchi;
    std::vector<double>  tmp;
    return SampleCovarianceVaryInput(specin,num_MC,&tmp);
}

TH1D SBNchi::SampleCovarianceVaryInput(SBNspec *specin, int num_MC, std::vector<double> * chival){
    if(!cholosky_performed) this->PerformCholoskyDecomposition(specin); 

    float** a_vec_matrix_lower_triangular = new float*[num_bins_total];
    float** a_vec_matrix_inverted = new float*[num_bins_total_compressed];

    for(int i=0; i < num_bins_total; i++){
        a_vec_matrix_lower_triangular[i] = new float[num_bins_total];
    }

    for(int i=0; i < num_bins_total_compressed; i++){
        a_vec_matrix_inverted[i] = new float[num_bins_total_compressed];
    }

    for(int i=0; i < num_bins_total; i++){
        for(int j=0; j < num_bins_total; j++){
            a_vec_matrix_lower_triangular[i][j] = vec_matrix_lower_triangular[i][j]; 
        }
    }

    for(int i=0; i< num_bins_total_compressed; i++){
        for(int j=0; j< num_bins_total_compressed; j++){
            a_vec_matrix_inverted[i][j] = vec_matrix_inverted[i][j]; 
        }
    }

    float *a_specin = new float[num_bins_total];
    float *a_corein = new float[num_bins_total_compressed];

    for(int i=0; i< num_bins_total; i++){
        a_specin[i] = specin->full_vector[i];
    }

    for(int i=0; i< num_bins_total_compressed; i++) {
        a_corein[i] = core_spectrum.collapsed_vector[i];
    }

    TH1D ans("","",std::max(200,(int)max_sample_chi_val),0,max_sample_chi_val );
    //ans.GetXaxis()->SetCanExtend(kTRUE);
    is_verbose = false;

    std::vector<float> vec_chis (num_MC, 0.0);

    float* a_vec_chis  = (float*)vec_chis.data();
    int num_chival = chival->size();
    float* a_chival = new float[num_chival];

    int *nlower = new int[num_chival];
    for(int i=0; i< num_chival; i++){
        nlower[i]=0; 
        a_chival[i] = chival->at(i);
    }

    std::vector < float > gaus_sample_v(num_bins_total), sampled_fullvector_v(num_bins_total);
    std::vector<float> collapsed_v(num_bins_total_compressed, 0.0);

    float* gaus_sample = new float[num_bins_total];
    float* sampled_fullvector = new float[num_bins_total] ;
    float* collapsed = new float[num_bins_total_compressed];

    //  float gaus_sample[54];
    // float sampled_fullvector[54];
    //  float collapsed[38];

    //We will need a uniform dist and a Gaussian
    std::uniform_int_distribution<int> dist_int(0,pow(2,32));
    std::normal_distribution<float> dist_normal(0,1);

#ifdef USE_GPU
    unsigned long long seed[num_MC];
    unsigned long long seq = 0ULL;
    unsigned long long offset = 0ULL;
    curandState_t state;

    for(int i=0; i<num_MC; ++i) {
        seed[i] = dist_int(*rangen_twister);
    }
#endif

#ifdef USE_GPU
#pragma acc parallel loop  private(gaus_sample[:54],sampled_fullvector[:54],collapsed[:36],state) \
    copyin(this[0:1],							\
            a_specin[:num_bins_total],					\
            a_vec_matrix_lower_triangular[:num_bins_total][:num_bins_total],\ 
            a_corein[:num_bins_total_compressed],				\
            a_vec_matrix_inverted[:num_bins_total_compressed][:num_bins_total_compressed],	\
            seed[0:num_MC],						\
            a_chival[:num_chival],						\
            this->a_num_bins[:num_channels],				\
            this->a_num_subchannels[:num_channels])			\    
        copyout(a_vec_chis[:num_MC]) \
        copy(nlower[:num_chival])
#endif

        for(int i=0; i < num_MC;i++){

#ifdef USE_GPU
            unsigned long long seed_sd = seed[i];
            curand_init(seed_sd, seq, offset, &state);

            for(int a=0; a<num_bins_total; a++) {
                gaus_sample[a]= curand_normal(&state);
            }
#else
            for(int a=0; a<num_bins_total; a++) {
                gaus_sample[a]= dist_normal(*rangen_twister);
            }      
#endif

            for(int j = 0; j < num_bins_total; j++){
                sampled_fullvector[j] = a_specin[j];
                for(int k = 0; k < num_bins_total; k++){
                    sampled_fullvector[j] += a_vec_matrix_lower_triangular[j][k] * gaus_sample[k];
                }

                if(sampled_fullvector[j]<0) sampled_fullvector[j]=0.0;

                sampled_fullvector[j] = rangen->Poisson(sampled_fullvector[j]);
                //sampled_fullvector[j] = rangen->Poisson(a_specin[j]);
                //std::cout<<"P: "<<a_specin[j]<<" "<<sampled_fullvector[j]<<std::endl;
            }

            this->CollapseVectorStandAlone(sampled_fullvector, collapsed);

            a_vec_chis[i] = this->CalcChi(a_vec_matrix_inverted, a_corein, collapsed);
            //Just to get some pvalues that were asked for.

            for(int j=0; j< num_chival; j++){
#pragma acc atomic update
                if(a_vec_chis[i]>=a_chival[j]) nlower[j]++;
            }

        }

    is_verbose = true;


    for(int i=0; i<num_MC; i++){
        //       if (i<(int)1e3) 
        //         std::cout << "@i=" << a_vec_chis[i] << std::endl;
        ans.Fill(a_vec_chis[i]);
    }
    for(int n =0; n< num_chival; n++){
        chival->at(n) = nlower[n]/(double)num_MC;
    }



    delete[] a_corein;
    delete[] a_specin;
    delete[] nlower;

    for(int i=0; i < num_bins_total; i++){
        delete[] a_vec_matrix_lower_triangular[i];
    }

    for(int i=0; i < num_bins_total_compressed; i++){
        delete[] a_vec_matrix_inverted[i];  
    }

    delete[] a_vec_matrix_lower_triangular;
    delete[] a_vec_matrix_inverted;

    delete[] gaus_sample;
    delete[] sampled_fullvector;
    delete[] collapsed;

    return ans;
}


int SBNchi::CollapseVectorStandAlone(std::vector<double> * full_vector, std::vector<double> *collapsed_vector){
    for(int im = 0; im < num_modes; im++){
        for(int id =0; id < num_detectors; id++){
            int edge = id*num_bins_detector_block + num_bins_mode_block*im; // This is the starting index for this detector
            int out_edge = edge;
            int tmp_chan = 0;
            for(int ic = 0; ic < num_channels; ic++){
                int corner=edge;

                for(int j=0; j< num_bins[ic]; j++){

                    double tempval=0;

                    for(int sc = 0; sc < num_subchannels[ic]; sc++){
                        tempval += (*full_vector)[j+sc*num_bins[ic]+corner];
                        edge +=1;	//when your done with a channel, add on every bin you just summed
                    }
                    //we can size this vector beforehand and get rid of all push_back()

                    int collapsed_index = tmp_chan+out_edge;
                    (*collapsed_vector)[collapsed_index] = tempval;
                    tmp_chan++;
                }
            }
        }
    }


    return 0;
}

int SBNchi::CollapseVectorStandAlone(float* full_vector, float *collapsed_vector){

    //int tmp_num_bins[3] = {25,25,6};
    //int tmp_num_subchannels[3] = {2,1,1};

    int collapsed_index = 0;
    for(int im = 0; im < num_modes; im++){
        for(int id =0; id < num_detectors; id++){
            int edge = id*num_bins_detector_block + num_bins_mode_block*im; // This is the starting index for this detector
            int out_edge = edge;
            //            int chan = 0;
            for(int ic = 0; ic < num_channels; ic++){
                int corner=edge;

                for(int j=0; j< this->a_num_bins[ic]; j++){

                    float tempval=0;

                    for(int sc = 0; sc < this->a_num_subchannels[ic]; sc++){
                        tempval += full_vector[j+sc*this->a_num_bins[ic]+corner];
                        //tempval += full_vector[j+sc*tmp_num_bins[ic]+corner];
                        edge +=1;	//when your done with a channel, add on every bin you just summed
                    }
                    //we can size this vector beforehand and get rid of all push_back()

                    // int collapsed_index = chan+out_edge;
                    collapsed_vector[collapsed_index] = tempval;
                    collapsed_index++;
                }
            }
        }
    }


    return 0;
}


int SBNchi::CollapseVectorStandAlone(double* full_vector, double *collapsed_vector){

    //int tmp_num_bins[3] = {25,25,6};
    //int tmp_num_subchannels[3] = {2,1,1};

    for(int im = 0; im < num_modes; im++){
        for(int id =0; id < num_detectors; id++){
            int edge = id*num_bins_detector_block + num_bins_mode_block*im; // This is the starting index for this detector
            int out_edge = edge;
            int chan = 0;
            for(int ic = 0; ic < num_channels; ic++){
                int corner=edge;

                for(int j=0; j< this->a_num_bins[ic]; j++){

                    double tempval=0;

                    for(int sc = 0; sc < this->a_num_subchannels[ic]; sc++){
                        tempval += full_vector[j+sc*this->a_num_bins[ic]+corner];
                        //tempval += full_vector[j+sc*tmp_num_bins[ic]+corner];
                        edge +=1;	//when your done with a channel, add on every bin you just summed
                    }
                    //we can size this vector beforehand and get rid of all push_back()

                    int collapsed_index = chan+out_edge;
                    collapsed_vector[collapsed_index] = tempval;
                    chan++;
                }
            }
        }
    }
    return 0;
}

std::vector<float> SBNchi::GeneratePseudoExperiment(){

    core_spectrum.CollapseVector();
    if(!cholosky_performed || is_stat_only) PerformCholoskyDecomposition(&core_spectrum); 

    int n_t =  (pseudo_from_collapsed ? num_bins_total_compressed : num_bins_total); 
    std::vector<float> sampled(n_t);
    is_verbose = false;

    std::vector<double> v_gaus(n_t,0.0);
    for(int i=0; i< n_t; ++i){
        //v_gaus[i] = rangen->Gaus(0,1);
        //v_gaus[i] = rangen->Uniform(-1,1); 
        v_gaus[i] = (*m_dist_normal)(*rangen_twister); 
    }

    for(int i=0; i< n_t; ++i){
        sampled[i] = (pseudo_from_collapsed ? core_spectrum.collapsed_vector[i] : core_spectrum.full_vector[i]); 
        if(!is_stat_only){
            for(int j=0; j<n_t; ++j){
                sampled[i] += vec_matrix_lower_triangular[i][j]*v_gaus[j];
            }
        }
    }
    //Now poisson fluctuate the sampled spectrum
    for(int j=0; j<n_t; ++j){
        std::poisson_distribution<int> dist_pois(sampled[j]);
        sampled[j] = float(dist_pois(*rangen_twister));
    }

    std::vector<float> collapsed(num_bins_total_compressed,0.0);
    if(pseudo_from_collapsed){
        collapsed = sampled;
    }else{
        this->CollapseVectorStandAlone(&sampled[0], &collapsed[0]);
    }
    return collapsed;
}


std::vector<float> SBNchi::SampleCovariance(SBNspec *specin){
    if(!cholosky_performed) this->PerformCholoskyDecomposition(specin); 

    int n_t = specin->full_vector.size();

    TVectorT<float> u(n_t);
    for(int i=0; i<n_t; i++){
        u(i) = specin->full_vector[i];
    }

    std::normal_distribution<float> dist_normal(0,1);
    is_verbose = false;

    TVectorT<float> gaus_sample(n_t);
    TVectorT<float> multi_sample(n_t);
    for(int a=0; a<n_t; a++){
        gaus_sample(a) = dist_normal(*rangen_twister);	
    }

    multi_sample = u + matrix_lower_triangular*gaus_sample;

    std::vector<float> sampled_fullvector(n_t,0.0);
    for(int j=0; j<n_t; j++){
        sampled_fullvector[j] = multi_sample(j);
    }

    std::vector<float> collapsed(num_bins_total_compressed,0.0);
    this->CollapseVectorStandAlone(&sampled_fullvector[0], &collapsed[0]);

    return collapsed;
}


TH1D SBNchi::SamplePoissonVaryInput(SBNspec *specin, int num_MC, double maxchi){ 
    max_sample_chi_val = maxchi;
    std::vector<double>  tmp = {};
    return SamplePoissonVaryInput(specin,num_MC,&tmp);
}


//This one varies the input comparative spectrum, and as sucn has  only to calculate the matrix_systematics once
TH1D SBNchi::SamplePoissonVaryInput(SBNspec *specin, int num_MC, std::vector<double> *chival){

    float** a_vec_matrix_inverted = new float*[num_bins_total_compressed];

    for(int i=0; i < num_bins_total_compressed; i++){
        a_vec_matrix_inverted[i] = new float[num_bins_total_compressed];
    }
    for(int i=0; i< num_bins_total_compressed; i++){
        for(int j=0; j< num_bins_total_compressed; j++){
            a_vec_matrix_inverted[i][j] = vec_matrix_inverted[i][j]; 
        }
    }

    float *a_specin = new float[num_bins_total];
    float *a_corein = new float[num_bins_total_compressed];

    for(int i=0; i< num_bins_total; i++){
        a_specin[i] = specin->full_vector[i];
    }

    for(int i=0; i< num_bins_total_compressed; i++) {
        a_corein[i] = core_spectrum.collapsed_vector[i];
    }

    std::vector<float> vec_chis (num_MC, 0.0);

    float* a_vec_chis  = (float*)vec_chis.data();
    int num_chival = chival->size();
    float* a_chival = new float[num_chival];

    int *nlower = new int[num_chival];
    for(int i=0; i< num_chival; i++){
        nlower[i]=0; 
        a_chival[i]=chival->at(i); 
    }
    float* sampled_fullvector = new float[num_bins_total] ;
    float* collapsed = new float[num_bins_total_compressed];

    //So save the core one that we will sample for
    //ans.GetXaxis()->SetCanExtend(kTRUE);
    is_verbose = false;

    std::vector< std::poisson_distribution<int>> dist_pois;
    //std::vector< std::normal_distribution<float>> dist_pois;
    for(int j = 0; j < num_bins_total; j++){
        //for tesing purposes
        dist_pois.push_back(std::poisson_distribution<int>(a_specin[j])); 
        //dist_pois.push_back(std::normal_distribution<float>(a_specin[j],sqrt(a_specin[j])));
    }

    for(int i=0; i < num_MC;i++){

        for(int j = 0; j < num_bins_total; j++){

            //float p = dist_pois[j](*rangen_twister); 
            int p = dist_pois[j](*rangen_twister); 

            sampled_fullvector[j] =  (float)p;
            //std::cout<<"P: "<<a_specin[j]<<" "<<sampled_fullvector[j]<<" "<<p<<std::endl;
        }

        this->CollapseVectorStandAlone(sampled_fullvector, collapsed);
        a_vec_chis[i] = this->CalcChi(a_vec_matrix_inverted, a_corein, collapsed);

        for(int j=0; j< num_chival; j++){
            if(a_vec_chis[i]>=a_chival[j]) nlower[j]++;
        }

    }

    TH1D ans("","",std::max(200,(int)max_sample_chi_val),0,max_sample_chi_val);


    for(int i=0; i<num_MC; i++){
        ans.Fill(a_vec_chis[i]);
    }

    for(int n =0; n< num_chival; n++){
        chival->at(n) = nlower[n]/(double)num_MC;
    }

    is_verbose = true;

    delete[] a_corein;
    delete[] a_specin;
    delete[] nlower;

    for(int i=0; i < num_bins_total_compressed; i++){
        delete[] a_vec_matrix_inverted[i];  
    }

    delete[] a_vec_matrix_inverted;

    delete[] sampled_fullvector;
    delete[] collapsed;

    return ans;


}

TH1D SBNchi::SamplePoisson_NP(SBNspec *specin, SBNchi &chi_h0, SBNchi & chi_h1, int num_MC, double maxchi,int which_sample){ 
    max_sample_chi_val = maxchi;
    std::vector<double>  tmp = {};
    return SamplePoisson_NP(specin,chi_h0,chi_h1,num_MC,&tmp, which_sample);
}


TH1D SBNchi::SamplePoisson_NP(SBNspec *specin, SBNchi &chi_h0, SBNchi & chi_h1, int num_MC, std::vector<double> *chival, int which_sample){

    float** h0_vec_matrix_inverted = new float*[num_bins_total_compressed];
    float** h1_vec_matrix_inverted = new float*[num_bins_total_compressed];

    for(int i=0; i < num_bins_total_compressed; i++){
        h0_vec_matrix_inverted[i] = new float[num_bins_total_compressed];
        h1_vec_matrix_inverted[i] = new float[num_bins_total_compressed];
    }
    for(int i=0; i< num_bins_total_compressed; i++){
        for(int j=0; j< num_bins_total_compressed; j++){
            h0_vec_matrix_inverted[i][j] = chi_h0.vec_matrix_inverted[i][j]; 
            h1_vec_matrix_inverted[i][j] = chi_h1.vec_matrix_inverted[i][j]; 
        }
    }

    float *a_specin = new float[num_bins_total];

    float *h0_corein = new float[num_bins_total_compressed];
    float *h1_corein = new float[num_bins_total_compressed];


    for(int i=0; i< num_bins_total; i++){
        a_specin[i] = specin->full_vector[i];
    }

    for(int i=0; i< num_bins_total_compressed; i++) {
        h0_corein[i] = chi_h0.core_spectrum.collapsed_vector[i];
        h1_corein[i] = chi_h1.core_spectrum.collapsed_vector[i];
    }

    std::vector<float> vec_chis (num_MC, 0.0);

    float* a_vec_chis  = (float*)vec_chis.data();
    int num_chival = chival->size();
    float* a_chival = new float[num_chival];

    int *nlower = new int[num_chival];
    for(int i=0; i< num_chival; i++){
        nlower[i]=0; 
        a_chival[i]=chival->at(i); 
    }
    float* sampled_fullvector = new float[num_bins_total] ;
    float* collapsed = new float[num_bins_total_compressed];

    float min_delta_chi = 99999;

    //So save the core one that we will sample for
    //ans.GetXaxis()->SetCanExtend(kTRUE);
    is_verbose = false;

    std::vector< std::poisson_distribution<int>> dist_pois;
    //std::vector< std::normal_distribution<float>> dist_pois;
    for(int j = 0; j < num_bins_total; j++){
        //for tesing purposes
        dist_pois.push_back(std::poisson_distribution<int>(a_specin[j])); 
        //dist_pois.push_back(std::normal_distribution<float>(a_specin[j],sqrt(a_specin[j])));
    }

    std::cout<<otag<<" Starting to generate "<<num_MC<<" pseudo universes according to poisson distribution"<<std::endl;
    for(int i=0; i < num_MC;i++){

        if(which_sample==0){//Poisson Mode
            for(int j = 0; j < num_bins_total; j++){
                //float p = dist_pois[j](*rangen_twister); 
                int p = dist_pois[j](*rangen_twister); 
                sampled_fullvector[j] =  (float)p;
                //std::cout<<"P: "<<a_specin[j]<<" "<<sampled_fullvector[j]<<" "<<p<<std::endl;
            }

            this->CollapseVectorStandAlone(sampled_fullvector, collapsed);

        }else if(which_sample==1){//Covariance Sampling

            std::vector<float> exp  = this->GeneratePseudoExperiment();

            for(int j = 0; j < num_bins_total_compressed; j++){
                collapsed[j] = exp[j];
            }
        }

        float val_chi_h0  = chi_h0.CalcChi(h0_vec_matrix_inverted, h0_corein, collapsed);
        float val_chi_h1  = chi_h1.CalcChi(h1_vec_matrix_inverted, h1_corein, collapsed);

        a_vec_chis[i] = val_chi_h0-val_chi_h1;

        //std::cout<<" "<<val_chi_h0<<" "<<val_chi_h1<<std::endl;

        /*
           std::cout<<" Pseudo Data "<<std::endl;
           for(int m=0; m<num_bins_total_compressed; m++) std::cout<<" "<<collapsed[m]<<" ";
           std::cout<<std::endl;
           std::cout<<" H0 "<<val_chi_h0<<std::endl;
           for(int m=0; m<num_bins_total_compressed; m++) std::cout<<" "<<h0_corein[m]<<" ";
           std::cout<<std::endl;
           std::cout<<" H1 "<<val_chi_h1<<std::endl;
           for(int m=0; m<num_bins_total_compressed; m++) std::cout<<" "<<h1_corein[m]<<" ";
           std::cout<<std::endl;
           std::cout<<" H0_ERR "<<std::endl;
           for(int m=0; m<num_bins_total_compressed; m++) std::cout<<" "<<h0_vec_matrix_inverted[m][m]<<" ";
           std::cout<<std::endl;
           std::cout<<" H1_ERR "<<std::endl;
           for(int m=0; m<num_bins_total_compressed; m++) std::cout<<" "<<h1_vec_matrix_inverted[m][m]<<" ";
           std::cout<<std::endl;
           std::cout<<"Delta H0-H1 "<<a_vec_chis[i]<<std::endl;
           */
        if(a_vec_chis[i] < min_delta_chi) min_delta_chi = a_vec_chis[i];

        for(int j=0; j< num_chival; j++){
            if(a_vec_chis[i]>=a_chival[j]) nlower[j]++;
        }

    }

    TH1D ans("","",std::max(200,(int)max_sample_chi_val),min_delta_chi,max_sample_chi_val);


    for(int i=0; i<num_MC; i++){
        ans.Fill(a_vec_chis[i]);
    }

    for(int n =0; n< num_chival; n++){
        chival->at(n) = nlower[n]/(double)num_MC;
    }

    is_verbose = true;

    delete[] h1_corein;
    delete[] h0_corein;
    delete[] a_specin;
    delete[] nlower;

    for(int i=0; i < num_bins_total_compressed; i++){
        delete[] h1_vec_matrix_inverted[i];  
        delete[] h0_vec_matrix_inverted[i];  
    }

    delete[] h1_vec_matrix_inverted;
    delete[] h0_vec_matrix_inverted;

    delete[] sampled_fullvector;
    delete[] collapsed;

    return ans;


}


std::vector<CLSresult> SBNchi::Mike_NP(SBNspec *specin, SBNchi &chi_h0, SBNchi & chi_h1, int num_MC, int which_sample, int id){

    std::cout << "SBNchi::Mike_NP" << std::endl;
    std::vector<CLSresult> v_results(5);

    float** h0_vec_matrix_inverted = new float*[num_bins_total_compressed];
    float** h1_vec_matrix_inverted = new float*[num_bins_total_compressed];

    for(int i=0; i < num_bins_total_compressed; i++){
        h0_vec_matrix_inverted[i] = new float[num_bins_total_compressed];
        h1_vec_matrix_inverted[i] = new float[num_bins_total_compressed];
    }

    for(int i=0; i< num_bins_total_compressed; i++){
        for(int j=0; j< num_bins_total_compressed; j++){
            h0_vec_matrix_inverted[i][j] = chi_h0.vec_matrix_inverted[i][j]; 
            h1_vec_matrix_inverted[i][j] = chi_h1.vec_matrix_inverted[i][j]; 
        }
    }

    float *a_specin = new float[num_bins_total];

    float *h0_corein = new float[num_bins_total_compressed];
    float *h1_corein = new float[num_bins_total_compressed];


    for(int i=0; i< num_bins_total; i++){
        a_specin[i] = specin->full_vector[i];
    }

    for(int i=0; i< num_bins_total_compressed; i++) {
        h0_corein[i] = chi_h0.core_spectrum.collapsed_vector[i];
        h1_corein[i] = chi_h1.core_spectrum.collapsed_vector[i];
    }

    std::vector<float> vec_chis (num_MC, 0.0);
    std::vector<float> vec_pois (num_MC, 0.0);
    std::vector<float> vec_cnp (num_MC, 0.0);
    std::vector<float> vec_h0 (num_MC, 0.0);
    std::vector<float> vec_h1 (num_MC, 0.0);

    float* a_vec_chis  = (float*)vec_chis.data();
    float* a_vec_pois  = (float*)vec_pois.data();
    float* a_vec_cnp  = (float*)vec_cnp.data();
    float* a_vec_h0  = (float*)vec_h0.data();
    float* a_vec_h1  = (float*)vec_h1.data();


    float* sampled_fullvector = new float[num_bins_total] ;
    float* collapsed = new float[num_bins_total_compressed];
    float* cv_collapsed = new float[num_bins_total_compressed];

    is_verbose = false;
    this->CollapseVectorStandAlone(a_specin, cv_collapsed);

    for(int i=0; i < num_MC;i++){

        //Generate our spectra
        if(which_sample==0){//Poisson Mode
            for(int j = 0; j < num_bins_total_compressed; j++){
                std::poisson_distribution<int> dist_pois(cv_collapsed[j]);
                collapsed[j] = float(dist_pois(*rangen_twister));
            }
        }else if(which_sample==1){//Covariance Sampling
            auto exp =  this->GeneratePseudoExperiment();
            for(int j=0; j < num_bins_total_compressed; j++){
                collapsed[j] = exp[j];
            }
        }

        //Base Default Chi
        //if(i%1000==0){ for(int b=0; b < num_bins_total_compressed; b++) std::cout << "h0_vec_matrix_inverted: " << h0_vec_matrix_inverted[b][b] << ", " << h0_corein[b] << ", " << collapsed[b] << std::endl;
        //for(int b=0; b < num_bins_total_compressed; b++) std::cout << "h1_vec_matrix_inverted: " << h0_vec_matrix_inverted[b][b] << ", " << h1_corein[b] << ", " << collapsed[b] << std::endl;}
        float val_chi_h0  = chi_h0.CalcChi(h0_vec_matrix_inverted, h0_corein, collapsed );
        float val_chi_h1  = chi_h1.CalcChi(h1_vec_matrix_inverted, h1_corein, collapsed );

        //Poisson Log Likli
        float val_pois_h0  = chi_h0.PoissonLogLiklihood(h0_corein, collapsed);

        float val_pois_h1  = chi_h1.PoissonLogLiklihood(h1_corein, collapsed);

        //CNP, going to need to recalculate and reinvert.
        float val_cnp_h0  = chi_h0.CalcChi_CNP(h0_corein, collapsed);
        float val_cnp_h1  = chi_h1.CalcChi_CNP(h1_corein, collapsed);

        a_vec_chis[i] = val_chi_h0 - val_chi_h1;
        a_vec_pois[i] = val_pois_h0 - val_pois_h1;
        a_vec_cnp[i]  = val_cnp_h0 - val_cnp_h1;
        a_vec_h0[i] = val_chi_h0;
        a_vec_h1[i] = val_chi_h1;

        if(i%1000==0) std::cout<<"Pseudo-Experiment: "<<i<<"/"<<num_MC<<" DeltaChi: "<<a_vec_chis[i]<<" PoisLogLiki: "<<a_vec_pois[i]<<" CNP_chi: "<<a_vec_cnp[i]<<std::endl;

        if(a_vec_chis[i] < v_results[0].m_min_value) v_results[0].m_min_value = a_vec_chis[i];
        if(a_vec_pois[i] < v_results[1].m_min_value) v_results[1].m_min_value = a_vec_pois[i];
        if(a_vec_cnp[i]  < v_results[2].m_min_value) v_results[2].m_min_value = a_vec_cnp[i];
        if(val_chi_h0  < v_results[3].m_min_value) v_results[3].m_min_value = val_chi_h0;
        if(val_chi_h1  < v_results[4].m_min_value) v_results[4].m_min_value = val_chi_h1;

        if(a_vec_chis[i] > v_results[0].m_max_value) v_results[0].m_max_value = a_vec_chis[i];
        if(a_vec_pois[i] > v_results[1].m_max_value) v_results[1].m_max_value = a_vec_pois[i];
        if(a_vec_cnp[i]  > v_results[2].m_max_value) v_results[2].m_max_value = a_vec_cnp[i];
        if(val_chi_h0  > v_results[3].m_max_value) v_results[3].m_max_value = val_chi_h0;
        if(val_chi_h1  > v_results[4].m_max_value) v_results[4].m_max_value = val_chi_h1;

    }

    for(int i=0; i<3;i++){
        std::cout<<"Res "<<i<<" "<<v_results[i].m_max_value<<" "<<v_results[i].m_min_value<<std::endl;
    }

    TH1D ans0(("0"+std::to_string(id)).c_str(),("0"+std::to_string(id)).c_str(),std::max(200,(int)v_results[0].m_max_value),v_results[0].m_min_value,v_results[0].m_max_value);
    TH1D ans1(("1"+std::to_string(id)).c_str(),("1"+std::to_string(id)).c_str(),std::max(200,(int)v_results[1].m_max_value),v_results[1].m_min_value,v_results[1].m_max_value);
    TH1D ans2(("2"+std::to_string(id)).c_str(),("2"+std::to_string(id)).c_str(),std::max(200,(int)v_results[2].m_max_value),v_results[2].m_min_value,v_results[2].m_max_value);
    TH1D ans3(("3"+std::to_string(id)).c_str(),("3"+std::to_string(id)).c_str(),std::max(200,(int)v_results[3].m_max_value),v_results[3].m_min_value,v_results[3].m_max_value);
    TH1D ans4(("4"+std::to_string(id)).c_str(),("4"+std::to_string(id)).c_str(),std::max(200,(int)v_results[4].m_max_value),v_results[4].m_min_value,v_results[4].m_max_value);

    for(int i=0; i<num_MC; i++){
        ans0.Fill(a_vec_chis[i]);
        ans1.Fill(a_vec_pois[i]);
        ans2.Fill(a_vec_cnp[i]);
        ans3.Fill(a_vec_h0[i]);
        ans4.Fill(a_vec_h1[i]);
    }
    v_results[0].m_pdf = ans0;
    v_results[1].m_pdf = ans1;
    v_results[2].m_pdf = ans2;
    v_results[3].m_pdf = ans3;
    v_results[4].m_pdf = ans4;

    v_results[0].m_values = vec_chis;
    v_results[1].m_values = vec_pois;
    v_results[2].m_values = vec_cnp;
    v_results[3].m_values = vec_h0;
    v_results[4].m_values = vec_h1;

    is_verbose = true;

    delete[] h1_corein;
    delete[] h0_corein;
    delete[] a_specin;

    for(int i=0; i < num_bins_total_compressed; i++){
        delete[] h1_vec_matrix_inverted[i];  
        delete[] h0_vec_matrix_inverted[i];  
    }

    delete[] h1_vec_matrix_inverted;
    delete[] h0_vec_matrix_inverted;

    delete[] sampled_fullvector;
    delete[] collapsed;

    return v_results;

}

//Do the same as above function, but taking into account the fakedata
std::vector<CLSresult> SBNchi::Mike_NP_fakedata(SBNspec *specin, std::vector<float> fakedata, std::vector<float> &chidata, SBNchi &chi_h0, SBNchi & chi_h1, int num_MC, int which_sample, int id ){

    std::vector<CLSresult> v_results(9);

    float** h0_vec_matrix_inverted = new float*[num_bins_total_compressed];
    float** h1_vec_matrix_inverted = new float*[num_bins_total_compressed];

    for(int i=0; i < num_bins_total_compressed; i++){
        h0_vec_matrix_inverted[i] = new float[num_bins_total_compressed];
        h1_vec_matrix_inverted[i] = new float[num_bins_total_compressed];
    }
    for(int i=0; i< num_bins_total_compressed; i++){
        for(int j=0; j< num_bins_total_compressed; j++){
            h0_vec_matrix_inverted[i][j] = chi_h0.vec_matrix_inverted[i][j]; 
            h1_vec_matrix_inverted[i][j] = chi_h1.vec_matrix_inverted[i][j];
        }
    }

    float *a_specin = new float[num_bins_total];
    float *a_fakedata = new float[num_bins_total_compressed];

    float *h0_corein = new float[num_bins_total_compressed];
    float *h1_corein = new float[num_bins_total_compressed];


    for(int i=0; i< num_bins_total; i++){
        a_specin[i] = specin->full_vector[i];
    }

    for(int i=0; i< num_bins_total_compressed; i++) {
        h0_corein[i] = chi_h0.core_spectrum.collapsed_vector[i];
        h1_corein[i] = chi_h1.core_spectrum.collapsed_vector[i];
        a_fakedata[i] = fakedata[i];
    }

    std::vector<float> vec_chis (num_MC, 0.0);
    std::vector<float> vec_pois (num_MC, 0.0);
    std::vector<float> vec_cnp (num_MC, 0.0);
    std::vector<float> vec_chih0 (num_MC, 0.0);
    std::vector<float> vec_chih1 (num_MC, 0.0);
    std::vector<float> vec_poish0 (num_MC, 0.0);
    std::vector<float> vec_poish1 (num_MC, 0.0);
    std::vector<float> vec_cnph0 (num_MC, 0.0);
    std::vector<float> vec_cnph1 (num_MC, 0.0);

    float* a_vec_chis  = (float*)vec_chis.data();
    float* a_vec_pois  = (float*)vec_pois.data();
    float* a_vec_cnp  = (float*)vec_cnp.data();
    float* a_vec_chih0  = (float*)vec_chih0.data();
    float* a_vec_chih1  = (float*)vec_chih1.data();
    float* a_vec_poish0  = (float*)vec_poish0.data();
    float* a_vec_poish1  = (float*)vec_poish1.data();
    float* a_vec_cnph0  = (float*)vec_cnph0.data();
    float* a_vec_cnph1  = (float*)vec_cnph1.data();


    float* sampled_fullvector = new float[num_bins_total] ;
    float* collapsed = new float[num_bins_total_compressed];
    float* cv_collapsed = new float[num_bins_total_compressed];

    is_verbose = false;
    this->CollapseVectorStandAlone(a_specin, cv_collapsed);


    std::cout<<otag<<" Starting to generate "<<num_MC<<" Pseduo-Experiments."<<std::endl;
    for(int i=0; i < num_MC;i++){

        //Generate our spectra
        if(which_sample==0){//Poisson Mode
            for(int j = 0; j < num_bins_total_compressed; j++){
                std::poisson_distribution<int> dist_pois(cv_collapsed[j]);
                collapsed[j] = float(dist_pois(*rangen_twister));
            }
        }else if(which_sample==1){//Covariance Sampling
            auto exp =  this->GeneratePseudoExperiment();
            for(int j=0; j < num_bins_total_compressed; j++){
                collapsed[j] = exp[j];
            }
        }

        //Base Default Chi
        float val_chi_h0  = chi_h0.CalcChi(h0_vec_matrix_inverted, h0_corein, collapsed);
        float val_chi_h1  = chi_h1.CalcChi(h1_vec_matrix_inverted, h1_corein, collapsed);

        //Poisson Log Likli
        float val_pois_h0  = chi_h0.PoissonLogLiklihood(h0_corein, collapsed);
        float val_pois_h1  = chi_h1.PoissonLogLiklihood(h1_corein, collapsed);

        //CNP, going to need to recalculate and reinvert.
        float val_cnp_h0  = chi_h0.CalcChi_CNP(h0_corein, collapsed);
        float val_cnp_h1  = chi_h1.CalcChi_CNP(h1_corein, collapsed);

        a_vec_chis[i] = val_chi_h0 - val_chi_h1;
        a_vec_pois[i] = val_pois_h0 - val_pois_h1;
        a_vec_cnp[i]  = val_cnp_h0 - val_cnp_h1;
        a_vec_chih0[i] = val_chi_h0;
        a_vec_chih1[i] = val_chi_h1;
        a_vec_poish0[i] = val_pois_h0;
        a_vec_poish1[i] = val_pois_h1;
        a_vec_cnph0[i] = val_cnp_h0;
        a_vec_cnph1[i] = val_cnp_h1;
/*
        a_vec_chis[i] = val_chi_h0;
        a_vec_pois[i] = val_pois_h0;
        a_vec_cnp[i]  = val_cnp_h0;
        a_vec_h0[i] = val_chi_h0;
        a_vec_h1[i] = val_chi_h1;
*/
        if(i%1000==0) std::cout<<"Pseudo-Experiment: "<<i<<"/"<<num_MC<<" DeltaChi: "<<a_vec_chis[i]<<" PoisLogLiki: "<<a_vec_pois[i]<<" CNP_chi: "<<a_vec_cnp[i]<<std::endl;

        if(a_vec_chis[i] < v_results[0].m_min_value) v_results[0].m_min_value = a_vec_chis[i];
        if(a_vec_pois[i] < v_results[1].m_min_value) v_results[1].m_min_value = a_vec_pois[i];
        if(a_vec_cnp[i]  < v_results[2].m_min_value) v_results[2].m_min_value = a_vec_cnp[i];
        if(val_chi_h0  < v_results[3].m_min_value) v_results[3].m_min_value = val_chi_h0;
        if(val_chi_h1  < v_results[4].m_min_value) v_results[4].m_min_value = val_chi_h1;
        if(val_pois_h0  < v_results[5].m_min_value) v_results[5].m_min_value = val_pois_h0;
        if(val_pois_h1  < v_results[6].m_min_value) v_results[6].m_min_value = val_pois_h1;
        if(val_cnp_h0  < v_results[7].m_min_value) v_results[7].m_min_value = val_cnp_h0;
        if(val_cnp_h1  < v_results[8].m_min_value) v_results[8].m_min_value = val_cnp_h1;

        if(a_vec_chis[i] > v_results[0].m_max_value) v_results[0].m_max_value = a_vec_chis[i];
        if(a_vec_pois[i] > v_results[1].m_max_value) v_results[1].m_max_value = a_vec_pois[i];
        if(a_vec_cnp[i]  > v_results[2].m_max_value) v_results[2].m_max_value = a_vec_cnp[i];
        if(val_chi_h0  > v_results[3].m_max_value) v_results[3].m_max_value = val_chi_h0;
        if(val_chi_h1  > v_results[4].m_max_value) v_results[4].m_max_value = val_chi_h1;
        if(val_pois_h0  < v_results[5].m_min_value) v_results[5].m_min_value = val_pois_h0;
        if(val_pois_h1  < v_results[6].m_min_value) v_results[6].m_min_value = val_pois_h1;
        if(val_cnp_h0  < v_results[7].m_min_value) v_results[7].m_min_value = val_cnp_h0;
        if(val_cnp_h1  < v_results[8].m_min_value) v_results[8].m_min_value = val_cnp_h1;

    }

    //--- fakedata ---
    float val_chi_h0fakedata  = chi_h0.CalcChi(h0_vec_matrix_inverted, h0_corein, a_fakedata);
    float val_chi_h1fakedata  = chi_h1.CalcChi(h1_vec_matrix_inverted, h1_corein, a_fakedata);
    float val_pois_h0fakedata  = chi_h0.PoissonLogLiklihood(h0_corein, a_fakedata);
    float val_pois_h1fakedata  = chi_h1.PoissonLogLiklihood(h1_corein, a_fakedata);
    float val_cnp_h0fakedata  = chi_h0.CalcChi_CNP(h0_corein, a_fakedata);
    float val_cnp_h1fakedata  = chi_h1.CalcChi_CNP(h1_corein, a_fakedata);

    float deltachi_fakedata = val_chi_h0fakedata - val_chi_h1fakedata;
    float deltapoischi_fakedata = val_pois_h0fakedata - val_pois_h1fakedata;
    float deltacnp_fakedata = val_cnp_h0fakedata - val_cnp_h1fakedata;

    std::cout << "    Fake Data delta Chi2 " << std::endl;
    std::cout << "    ==================== " << std::endl;
    std::cout << "    Pearson : " << val_chi_h0fakedata << " - " << val_chi_h1fakedata << " = " << deltachi_fakedata << std::endl;
    std::cout << "    Poisson : " << val_pois_h0fakedata << " - " << val_pois_h1fakedata << " = " << deltapoischi_fakedata << std::endl;
    std::cout << "    CNP : " << val_cnp_h0fakedata << " - " << val_cnp_h1fakedata << " = " << deltacnp_fakedata << std::endl;

    chidata.push_back(deltachi_fakedata);
    chidata.push_back(deltapoischi_fakedata);
    chidata.push_back(deltacnp_fakedata);
    chidata.push_back(val_chi_h0fakedata);
    chidata.push_back(val_chi_h1fakedata);
    chidata.push_back(val_pois_h0fakedata);
    chidata.push_back(val_pois_h1fakedata);
    chidata.push_back(val_cnp_h0fakedata);
    chidata.push_back(val_cnp_h1fakedata);
    // --------
 
    for(int i=0; i<3;i++){
        std::cout<<"Res "<<i<<" "<<v_results[i].m_max_value<<" "<<v_results[i].m_min_value<<std::endl;
    }

    TH1D ans0(("0"+std::to_string(id)).c_str(),("0"+std::to_string(id)).c_str(),std::max(200,(int)v_results[0].m_max_value),v_results[0].m_min_value,v_results[0].m_max_value);
    TH1D ans1(("1"+std::to_string(id)).c_str(),("1"+std::to_string(id)).c_str(),std::max(200,(int)v_results[1].m_max_value),v_results[1].m_min_value,v_results[1].m_max_value);
    TH1D ans2(("2"+std::to_string(id)).c_str(),("2"+std::to_string(id)).c_str(),std::max(200,(int)v_results[2].m_max_value),v_results[2].m_min_value,v_results[2].m_max_value);
    TH1D ans3(("3"+std::to_string(id)).c_str(),("3"+std::to_string(id)).c_str(),std::max(200,(int)v_results[3].m_max_value),v_results[3].m_min_value,v_results[3].m_max_value);
    TH1D ans4(("4"+std::to_string(id)).c_str(),("4"+std::to_string(id)).c_str(),std::max(200,(int)v_results[4].m_max_value),v_results[4].m_min_value,v_results[4].m_max_value);
    TH1D ans5(("5"+std::to_string(id)).c_str(),("5"+std::to_string(id)).c_str(),std::max(200,(int)v_results[5].m_max_value),v_results[5].m_min_value,v_results[5].m_max_value);
    TH1D ans6(("6"+std::to_string(id)).c_str(),("6"+std::to_string(id)).c_str(),std::max(200,(int)v_results[6].m_max_value),v_results[6].m_min_value,v_results[6].m_max_value);
    TH1D ans7(("7"+std::to_string(id)).c_str(),("7"+std::to_string(id)).c_str(),std::max(200,(int)v_results[7].m_max_value),v_results[7].m_min_value,v_results[7].m_max_value);
    TH1D ans8(("8"+std::to_string(id)).c_str(),("8"+std::to_string(id)).c_str(),std::max(200,(int)v_results[8].m_max_value),v_results[8].m_min_value,v_results[8].m_max_value);

    for(int i=0; i<num_MC; i++){
        ans0.Fill(a_vec_chis[i]);
        ans1.Fill(a_vec_pois[i]);
        ans2.Fill(a_vec_cnp[i]);
        ans3.Fill(a_vec_chih0[i]);
        ans4.Fill(a_vec_chih1[i]);
        ans5.Fill(a_vec_poish0[i]);
        ans6.Fill(a_vec_poish1[i]);
        ans7.Fill(a_vec_cnph0[i]);
        ans8.Fill(a_vec_cnph1[i]);
    }
    v_results[0].m_pdf = ans0;
    v_results[1].m_pdf = ans1;
    v_results[2].m_pdf = ans2;
    v_results[3].m_pdf = ans3;
    v_results[4].m_pdf = ans4;
    v_results[5].m_pdf = ans5;
    v_results[6].m_pdf = ans6;
    v_results[7].m_pdf = ans7;
    v_results[8].m_pdf = ans8;

    v_results[0].m_values = vec_chis;
    v_results[1].m_values = vec_pois;
    v_results[2].m_values = vec_cnp;
    v_results[3].m_values = vec_chih0;
    v_results[4].m_values = vec_chih1;
    v_results[5].m_values = vec_poish0;
    v_results[6].m_values = vec_poish1;
    v_results[7].m_values = vec_cnph0;
    v_results[8].m_values = vec_cnph1;

    is_verbose = true;

    delete[] h1_corein;
    delete[] h0_corein;
    delete[] a_specin;

    for(int i=0; i < num_bins_total_compressed; i++){
        delete[] h1_vec_matrix_inverted[i];  
        delete[] h0_vec_matrix_inverted[i];  
    }

    delete[] h1_vec_matrix_inverted;
    delete[] h0_vec_matrix_inverted;

    delete[] sampled_fullvector;
    delete[] collapsed;

    return v_results;

}



/*
   std::vector<double> SBNchi::SampleCovarianceVaryInput_getpval(SBNspec *specin, int num_MC, std::vector<double> chival){
   if(!cholosky_performed) this->PerformCholoskyDecomposition(specin); 

   int n_t = specin->full_vector.size();
   std::vector<int> nlower(chival.size(),0);

   TVectorT<double> u(n_t);
   for(int i=0; i<n_t; i++){
   u(i) = specin->full_vector.at(i);
   }

   TRandom3 * rangen = new TRandom3(0);


   TH1D ans("","",100,0,100);
   ans.GetXaxis()->SetCanExtend(kTRUE);
   is_verbose = false;
   for(int i=0; i < num_MC;i++){

   TVectorT<double> gaus_sample(n_t);
   TVectorT<double> multi_sample(n_t);
   for(int a=0; a<n_t; a++){
   gaus_sample(a) = rangen->Gaus(0,1);	
   }

   multi_sample = u + matrix_lower_triangular*gaus_sample;

   std::vector<double> sampled_fullvector(n_t,0.0);
   for(int i=0; i<n_t; i++){
   sampled_fullvector.at(i) = multi_sample(i);
   }
   SBNspec sampled_spectra(sampled_fullvector, specin->xmlname ,false);

   sampled_spectra.CollapseVector(); //this line important isnt it!

   double thischi = this->CalcChi(&sampled_spectra);
   ans.Fill(thischi);

   for(int j=0; j< chival.size(); j++){
   if(thischi>=chival.at(j)) nlower.at(j)++;
   }


   if(i%1000==0) std::cout<<"SBNchi::SampleCovarianceVaryInput(SBNspec*, int) on MC :"<<i<<"/"<<num_MC<<". Ans: "<<thischi<<std::endl;
   }
   is_verbose = true;

   std::vector<double> pval;
   for(auto n: nlower){
   pval.push_back(n/(double)num_MC);

   }

   return pval;

   }


//This one varies the input comparative spectrum, and as sucn has  only to calculate the matrix_systematics once
std::vector<double> SBNchi::SamplePoissonVaryInput_getpval(SBNspec *specin, int num_MC, std::vector<double> chival){
std::vector<int> nlower(chival.size(),0);

TRandom3 *rangen = new TRandom3(0);

TH1D ans("","",100,0,100);
//So save the core one that we will sample for
ans.GetXaxis()->SetCanExtend(kTRUE);
is_verbose = false;
for(int i=0; i < num_MC;i++){

SBNspec tmp = *specin;
tmp.ScalePoisson(rangen);
tmp.CollapseVector(); //this line important isnt it!
//tmp.PrintFullVector();

double thischi = this->CalcChi(&tmp);
ans.Fill(thischi);

for(int j=0; j< chival.size(); j++){
    if(thischi>=chival.at(j)) nlower.at(j)++;
}

if(i%1000==0) std::cout<<"SBNchi::SamplePoissonVaryInput(SBNspec*, int) on MC :"<<i<<"/"<<num_MC<<". Ans: "<<thischi<<std::endl;
}
std::vector<double> pval;
for(auto n: nlower){
    pval.push_back(n/(double)num_MC);

}

is_verbose = true;
return pval;


}
*/

//This one varies the core spectrum, and as sucn has to recalculate the matrix_systematics each stem
TH1D SBNchi::SamplePoissonVaryCore(SBNspec *specin, int num_MC){
    double center = this->CalcChi(specin);
    int nlower=0;

    TRandom3 *rangen = new TRandom3(0);

    TH1D ans("MCans","MCans",100,center-100,center+200);
    //So save the core one that we will sample for
    SBNspec core  = core_spectrum;

    is_verbose = false;
    for(int i=0; i<num_MC;i++){

        SBNspec tmp = core;
        tmp.ScalePoisson(rangen);
        this->ReloadCoreSpectrum(&tmp);
        double thischi = this->CalcChi(specin);
        ans.Fill(thischi);
        if(thischi<=center)nlower++;

        if(i%1000==0) std::cout<<"SBNchi::SamplePoissonVaryCore(SBNspec*, int) on MC :"<<i<<"/"<<num_MC<<". Ans: "<<thischi<<std::endl;
    }
    std::cout<<"pval: "<<nlower/(double)num_MC<<std::endl;

    is_verbose = true;
    return ans;
}

//Pelee specific detsys
void SBNchi::FillDetSysMatrix(TMatrixT <double> &M, SBNspec core_spectrum, bool frac){
    int matrix_size = M.GetNrows();

    if(matrix_size != (core_spectrum.full_vector).size()){std::cout<<"#ERROR: FillStatsMatrix, matrix not equal to diagonal"<<std::endl;}
    if(M.GetNrows()!=M.GetNcols()){std::cout<<"#ERROR: not a square matrix!"<<std::endl;}

    M.Zero();

    int j=0;

    //BDT
    std::vector<double> np_detsys = {0.2007, 0.1077, 0.0922, 0.0574, 0.0663, 0.0755, 0.0721, 0.0872, 0.0975, 0.1034, 0.2551, 0.0849, 0.1428, 0.1764, 0.1806, 0.2042, 0.1841, 0.1688, 0.1811}; //RealAnalysis
    //std::vector<double> np_detsys = {0.2454, 0.1523, 0.1546, 0.0900, 0.1500, 0.0608, 0.1530, 0.0781, 0.1114, 0.1145, 0.2005, 0.1237, 0.1871, 0.1105, 0.1349, 0.1612, 0.2478}; //FD
    //1e0p
    std::vector<double> zp_detsys = {0.0985, 0.0985, 0.1015, 0.1015, 0.1929, 0.1929, 0.2326, 0.2326, 0.3289, 0.3289, 0.1720, 0.1720, 0.1720, 0.1720, 0.1720, 0.1720, 0.1720, 0.1720, 0.1720}; //RealAnalysis
    //std::vector<double> zp_detsys = {0.1520, 0.1520, 0.0980, 0.0980, 0.1939, 0.1939, 0.4064, 0.4064, 0.3259, 0.3259, 0.1819, 0.1819, 0.1819, 0.1819, 0.2540, 0.2540, 0.2540, 0.2540, 0.2540}; //FD
    //numu
    std::vector<double> numu_detsys = {0.096,0.097,0.066,0.051,0.065,0.093,0.081,0.07,0.109,0.122,0.142,0.158,0.18,0.261};//RealAnalysis
    //std::vector<double> numu_detsys = {0.1655, 0.1044, 0.1233, 0.0574, 0.0500, 0.0849, 0.0686, 0.1449, 0.1158, 0.0849, 0.0980, 0.1536, 0.1556, 0.1879, 0.2987}; //FD

    //initialize vectors to store index of lee and intrinsic to make the correlation
    std::vector<int> col_intrinsic_np, col_lee_np, col_intrinsic_zp, col_lee_zp;
    for(auto& h: core_spectrum.hist){
        std::string hname = h.GetName();
        for( int i=1; i < h.GetNbinsX()+1; i++ ){
          if( hname.find("1eNp_bg_intrinsic") != std::string::npos ){
            std::cout << "Fill 1eNp signal detsys error, histo, bin number, matrix column = " << h.GetName() << ", " << i << ", " << j <<std::endl;
            if(!frac) M(j,j) = np_detsys[i-1]*np_detsys[i-1]*h.GetBinContent(i)*h.GetBinContent(i);
            else M(j,j) = np_detsys[i-1]*np_detsys[i-1];
            std::cout << ", " << np_detsys[i-1] << ", " << h.GetBinContent(i) << ", " << M(j,j) << std::endl;
            col_intrinsic_np.push_back(j);
          }
          else if( hname.find("1eNp_sig_lee") != std::string::npos ){
            std::cout << "Fill 1eNp signal detsys error, histo, bin number, matrix column = " << h.GetName() << ", " << i << ", " << j << std::endl;
            if(!frac) M(j,j) = np_detsys[i-1]*np_detsys[i-1]*h.GetBinContent(i)*h.GetBinContent(i);
            else M(j,j) = np_detsys[i-1]*np_detsys[i-1];
            std::cout << ", " << np_detsys[i-1] << ", " << h.GetBinContent(i) << ", " << M(j,j) << std::endl;
            col_lee_np.push_back(j);
          }
          else if( hname.find("1e0p_bg_intrinsic") != std::string::npos ){
            std::cout << "Fill 1e0p signal detsys error, histo, bin number, matrix column = " << h.GetName() << ", " << i << ", " << j << std::endl;
            if(!frac) M(j,j) = zp_detsys[i-1]*zp_detsys[i-1]*h.GetBinContent(i)*h.GetBinContent(i);
            else M(j,j) = zp_detsys[i-1]*zp_detsys[i-1];
            std::cout << ", " << zp_detsys[i-1] << ", " << h.GetBinContent(i) << ", " << M(j,j)*h.GetBinContent(i)*h.GetBinContent(i) << std::endl;
            col_intrinsic_zp.push_back(j);
          }
          else if( hname.find("1e0p_sig_lee") != std::string::npos){
            std::cout << "Fill 1e0p signal detsys error, histo, bin number, matrix column = " << h.GetName() << ", " << i << ", " << j << std::endl;
            if(!frac) M(j,j) = zp_detsys[i-1]*zp_detsys[i-1]*h.GetBinContent(i)*h.GetBinContent(i);
            else M(j,j) = zp_detsys[i-1]*zp_detsys[i-1];
            std::cout << ", " << zp_detsys[i-1] << ", " << h.GetBinContent(i) << ", " << M(j,j)*h.GetBinContent(i)*h.GetBinContent(i) << std::endl;
            col_lee_zp.push_back(j);
          }
          else if( hname.find("numu_bnb") != std::string::npos ){
            std::cout << "Fill numu signal detsys error, histo, bin number, matrix column = " << h.GetName() << ", " << i << ", " << j << std::endl;
            if(!frac) M(j,j) = numu_detsys[i-1]*numu_detsys[i-1]*h.GetBinContent(i)*h.GetBinContent(i);
            else M(j,j) = numu_detsys[i-1]*numu_detsys[i-1];
            std::cout << ", " << numu_detsys[i-1] << ", " << h.GetBinContent(i) << ", " << M(j,j) << std::endl;
            //col_numu_bnb.push_back(j); 
          }
          else if( hname.find("numu_dirt") != std::string::npos ){
            std::cout << "Fill numu signal detsys error, histo, bin number, matrix column = " << h.GetName() << ", " << i << ", " << j << std::endl;
            if(!frac) M(j,j) = 0.2*0.2*h.GetBinContent(i)*h.GetBinContent(i);
            else M(j,j) = 0.2*0.2;
            std::cout << ", " << numu_detsys[i-1] << ", " << h.GetBinContent(i) << ", " << M(j,j) << std::endl;
            //col_numu_bnb.push_back(j); 
          }
	  else if( hname.find("ext") == std::string::npos && hname.find("data") == std::string::npos ){
          std::cout << "Fill bg detsys error, histo, bin number, matrix column = " << h.GetName() << ", " << i << ", " << j << std::endl;
            if(!frac)M(j,j) = 0.2*0.2*h.GetBinContent(i)*h.GetBinContent(i);
            else M(j,j) = 0.2*0.2;
            std::cout << ", 0.2 , " << h.GetBinContent(i) << ", " << M(j,j)*h.GetBinContent(i)*h.GetBinContent(i) << std::endl;
          }
	  else if( hname.find("ext") != std::string::npos ){ std::cout << "@@@@@@@@ histo name, bin, number of event: " << h.GetName() << ", " << i << ", " << h.GetBinContent(i) << std::endl;}
          j++; 
        }
    }
    //fill correlation between nue and lee
    core_spectrum.CalcFullVector();
    //1eNp
    std::cout << "col_intrinsic_np.size() = " << col_intrinsic_np.size() << std::endl;
    for(int i=0; i<col_intrinsic_np.size(); i++){
      std::cout << "i, col_lee_np[i] = " << i << ", " << col_lee_np[i] << std::endl;
      M(col_intrinsic_np[i],col_lee_np[i]) = sqrt(M(col_intrinsic_np[i],col_intrinsic_np[i])*M(col_lee_np[i],col_lee_np[i])); 
      M(col_lee_np[i],col_intrinsic_np[i]) = sqrt(M(col_lee_np[i],col_lee_np[i])*M(col_intrinsic_np[i],col_intrinsic_np[i])); 
    }
    //1eNp
    std::cout << "col_intrinsic_zp.size() = " << col_intrinsic_zp.size() << std::endl;
    for(int i=0; i<col_intrinsic_zp.size(); i++){
      M(col_intrinsic_zp[i],col_lee_zp[i]) = sqrt(M(col_intrinsic_zp[i],col_intrinsic_zp[i])*M(col_lee_zp[i],col_lee_zp[i])); 
      M(col_lee_zp[i],col_intrinsic_zp[i]) = sqrt(M(col_lee_zp[i],col_lee_zp[i])*M(col_intrinsic_zp[i],col_intrinsic_zp[i])); 
    }

    //Print Detsys error
    //std::cout << "========DETSYS ERROR MATRIX========" << std::endl;
    //M.Print();

    return ;
}

