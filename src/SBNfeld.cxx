#include "SBNfeld.h"
using namespace sbn;


NeutrinoModel SBNfeld::convert3p1(std::vector<double> ingrid){
    if(ingrid.size()!=3){
        std::cout<<"ERROR: Currently must have a 3-grid, for Dm ue4 and um4. grid.size(): "<<ingrid.size()<<std::endl;
        exit(EXIT_FAILURE);
    }

    NeutrinoModel signalModel(pow(10,ingrid[0]),pow(10,ingrid[1]),pow(10,ingrid[2]));

    return signalModel;
}

int SBNfeld::GenerateOscillatedSpectra(){

    //This will loop over the predefined grid for oscillation studies, and calculate it

    int n_mass = m_grid.f_dimensions[0].GetNPoints();
    bool cv_runone = true;
    std::cout<<"Beginining to Generate all oscillated spectra for : "<<n_mass<<" Mass splittings."<<std::endl;

    for(size_t t =0; t< n_mass; t++){

        std::cout<<"SBNfeld::GenerateOcillatedSpectra()\t\t||\t\t On Mass Splitting number "<<t<<" which is "<<m_grid.f_dimensions[0].GetPoint(t)<<std::endl;
        //need to convert from this gridpoints to a neutrinoModel (Blarg, don't like this step but needed at the moment)
        NeutrinoModel this_model(pow(10,m_grid.f_dimensions[0].GetPoint(t)),1,1); // = this->convert3p1(m_vec_grid.at(t)); 
        this_model.Printall();

        // on construction it makes 3 SBNspecs, 1 sin amp, 1 sin2 amp, 1 CV oscilatted
        SBNgenerate gen(this->xmlname,this_model);

        // Write them to file
        gen.WritePrecomputedOscSpecs(this->tag);

        if(cv_runone){
            gen.spec_central_value.WriteOut(this->tag+"_CV");
            cv_runone= false;
        }
    }
}


int SBNfeld::GenerateScaledSpectra(){

    //Rather than oscillations, we just want to pick 1 subchannel, and scale it up and down. The channel to scale is defined in xml, scale_signal="true"
    // Basically, we take the core spectra, scale it by each gridpoint, and push that back

    std::cout<<"Beginining to Generate all Scaled spectra."<<std::endl;
    m_cv_spec_grid.clear();
    m_cv_spec_grid.resize(m_grid.f_num_total_points);
    m_core_spectrum->CalcFullVector();

    std::vector<std::vector<double>> vgrid = m_grid.GetGrid();

    for(size_t t =0; t < m_grid.f_num_total_points; t++){

        std::cout<<"SBNfeld::GenerateScaledSpectra()\t\t||\t\t On scaling point "<<t<<std::endl;

        for(int i=0 ; i<m_grid.f_num_dimensions; i++){
            std::cout<<"SBNfeld::GenerateScaledSpectra()\t\t||\t\t  which is "<<m_grid.f_dimensions[i].f_name<<" "<<vgrid[t][i]<<std::endl;
        }

        m_cv_spec_grid[t] = new SBNspec(m_core_spectrum->xmlname,t); //m_core_spectrum->full_vector, m_core_spectrum->xmlname, t, false);
        (*m_cv_spec_grid[t]) = (*m_core_spectrum);

        for(int i=0 ; i<m_grid.f_num_dimensions; i++){
            m_cv_spec_grid[t]->Scale(m_grid.f_dimensions[i].f_name, vgrid[t][i]);
        }

        m_cv_spec_grid[t]->CalcFullVector();
        m_cv_spec_grid[t]->CollapseVector();
    }
}



int SBNfeld::GenerateBackgroundSpectrum(){

    std::cout<<"SBNfeld::GenerateBackgroundSpectrum()\t\t||\t\t Generating a background 3+0 spectra"<<std::endl;

    NeutrinoModel this_model = this->convert3p1(m_vec_grid[0]); 
    NeutrinoModel background_only_model(this_model.mNu[0],0.0,0.0); // quirk, this works better

    m_core_spectrum->LoadModel(background_only_model);

    //Is this a troublesome line?!? Shouldn't be right?!
    //m_core_spectrum->SetAppMode();

    std::vector<double> ans = m_core_spectrum->Oscillate(this->tag, false);
    SBNspec background(ans, m_core_spectrum->xmlname,false);
    background.Scale("fullosc",0.0);
    background.Scale("ncdelta",0.0);
    background.WriteOut(this->tag+"_BKG_ONLY");

    std::cout<<"SBNfeld::GenerateBackgroundSpectrum()\t\t||\t\t Done."<<std::endl;
    return 0;

}

int SBNfeld::GenerateBackgroundScaledSpectrum(){

    std::cout<<"SBNfeld::GenerateBackgroundScaledSpectrum()\t\t||\t\t Generating a background spectra"<<std::endl;

    SBNspec background(m_core_spectrum->xmlname,false);
    background.WriteOut(this->tag+"_BKG_ONLY");

    std::cout<<"SBNfeld::GenerateBackgroundScaledSpectrum()\t\t||\t\t Done."<<std::endl;
    return 0;

}


int SBNfeld::SetCoreSpectrum(std::string file){

    std::cout<<"Set Core Spectrum1"<<std::endl;
    m_core_spectrum= new SBNosc(file,this->xmlname);
    m_bool_core_spectrum_set = true;
    std::cout<<"Set Core Spectrum2"<<std::endl;
    return 0;
}

int SBNfeld::AddFlatDetSystematic(double percent){

    std::cout<<"Adding a "<<percent<<" percent to the diagonal "<<percent*percent<<std::endl;
    for(int j=0; j< m_full_fractional_covariance_matrix->GetNrows(); j++){
        std::cout<<j<<" before "<<(*m_full_fractional_covariance_matrix)(j,j);
        (*m_full_fractional_covariance_matrix)(j,j)+=percent*percent;   
        std::cout<<" after "<<(*m_full_fractional_covariance_matrix)(j,j)<<std::endl;

    }

    return 0 ;
}

int SBNfeld::SetFractionalCovarianceMatrix(std::string filename,std::string matrixname){
    TFile * fsys = new TFile(filename.c_str(),"read");
    m_full_fractional_covariance_matrix = (TMatrixD*)fsys->Get(matrixname.c_str());
    fsys->Close();

    return 0;
}

int SBNfeld::SetFractionalCovarianceMatrix(TMatrixT<double> * in){

    m_full_fractional_covariance_matrix = in;
    return 0;
}

int SBNfeld::SetEmptyFractionalCovarianceMatrix(){

    m_full_fractional_covariance_matrix = new TMatrixT<double>(num_bins_total,num_bins_total);
    m_full_fractional_covariance_matrix->Zero();
    return 0;
}

int SBNfeld::SetRandomSeed(double sin){
    m_random_seed = sin;
    std::cout<<"SBNfeld::SetRandomSeed()\t\t||\t\tSetting random seed to: "<<sin<<std::endl;
    return 0;
}

int SBNfeld::LoadPreOscillatedSpectrum(int which_pt){

    std::cout<<"Loading Just Grid Pt: "<<which_pt<<" ";
    for(int k=0; k<m_vec_grid[which_pt].size();k++){
        std::cout<<" "<<m_vec_grid[which_pt][k];
    }
    std::cout<<std::endl;

    //need to convert from this gridpoints to a neutrinoModel (Blarg, don't like this step but needed at the moment)
    NeutrinoModel this_model = this->convert3p1(m_vec_grid[which_pt]); 

    // this_model.Printall();
    //And load thus model into our spectra. At this point its comuted all the necessary mass-splittins and which frequencies they are
    m_core_spectrum->LoadModel(this_model);
    //m_core_spectrum->SetAppMode();

    //And apply this oscillaion! Adding to it the bkgSpec that it was initilised with.
    //NOTE we want to return the FULL spectrum, not compressed so we can calculate the covariance matrix, hense the false in this Oscilate
    std::vector<double> ans = m_core_spectrum->Oscillate(this->tag, false);
    std::cout<<"Spectrum: ";
    for(int p=0; p<ans.size();p++){
        std::cout<<" "<<ans[p];
    }
    SBNspec * thispoint = new SBNspec(ans, m_core_spectrum->xmlname,which_pt, false);

    thispoint->ScaleAll(global_scale);
    thispoint->CollapseVector();

    //make a print out of this exact spectrum as compared to the "core" spectrum
    std::string tlog  = std::to_string(m_vec_grid[which_pt][0])+"_"+std::to_string(m_vec_grid[which_pt][1])+"_"+std::to_string(m_vec_grid[which_pt][2]);
    m_core_spectrum->CompareSBNspecs(thispoint,tlog); 

    return 0;
}

int SBNfeld::LoadPreOscillatedSpectra(){
    //This is where we load the precalculated spectra and assemble into a actuall oscillate spectrum
    //we have a core spectrum loaded on the start of the SBNfeld class. This is the fullosc+intrinsics...etc..
    // i.e m_cv_spec_grid;

    m_cv_spec_grid.clear();
    std::cout<<"Beginining to Grab all oscillated spectra"<<std::endl;
    m_cv_spec_grid.resize(m_num_total_gridpoints);

    for(size_t t =0; t < m_num_total_gridpoints; t++){
        std::cout<<"SBNfeld::LoadPreOcillatedSpectra()\t\t||\t\t On spectrum "<<t<<"/"<<m_num_total_gridpoints;
        for(int k=0; k<m_vec_grid[t].size();k++){
            std::cout<<" "<<m_vec_grid[t][k];
        }
        std::cout<<std::endl;

        //need to convert from this gridpoints to a neutrinoModel (Blarg, don't like this step but needed at the moment)
        NeutrinoModel this_model = this->convert3p1(m_vec_grid[t]); 

        // this_model.Printall();
        //And load thus model into our spectra. At this point its comuted all the necessary mass-splittins and which frequencies they are
        m_core_spectrum->LoadModel(this_model);
        //m_core_spectrum->SetAppMode();

        //And apply this oscillaion! Adding to it the bkgSpec that it was initilised with.
        //NOTE we want to return the FULL spectrum, not compressed so we can calculate the covariance matrix, hense the false in this Oscilate
        std::vector<double> ans = m_core_spectrum->Oscillate(this->tag, false);
        std::cout<<"Spectrum: ";
        for(int p=0; p<ans.size();p++){
            std::cout<<" "<<ans[p];
        }
        std::cout<<std::endl;
        m_cv_spec_grid[t] = new SBNspec(ans, m_core_spectrum->xmlname,t, false);
        m_cv_spec_grid[t]->ScaleAll(global_scale);
        m_cv_spec_grid[t]->CollapseVector();


        if(m_bool_print_comparasons && t ==490){// t==1668 
            //make a print out of this exact spectrum as compared to the "core" spectrum
            std::string tlog  = std::to_string(m_vec_grid[t][0])+"_"+std::to_string(m_vec_grid[t][1])+"_"+std::to_string(m_vec_grid[t][2]);
            m_core_spectrum->CompareSBNspecs(m_cv_spec_grid[t],tlog); 
        }
    }

    return 0;
}


int SBNfeld::LoadBackgroundSpectrum(){
    m_background_spectrum= new SBNosc(this->tag+"_BKG_ONLY.SBNspec.root",this->xmlname);
    m_background_spectrum->ScaleAll(global_scale);
    m_bool_background_spectrum_set = true;

    m_background_spectrum->CollapseVector();
    m_tvec_background_spectrum = new TVectorT<double>(m_background_spectrum->full_vector.size(), &(m_background_spectrum->full_vector)[0]);
    m_tvec_background_mcerr = new TVectorT<double>(m_background_spectrum->full_error.size(), &(m_background_spectrum->full_error)[0]);

    if(m_background_spectrum->full_vector.size() !=  m_full_fractional_covariance_matrix->GetNcols()){

        std::cout<<"SBNfeld::LoadBackgroundSpectrum || ERROR!! background spectrum is of length : "<<m_background_spectrum->full_vector.size()<<" and frac matrix is of size "<<m_full_fractional_covariance_matrix->GetNcols()<<std::endl;
        exit(EXIT_FAILURE);
    }


    m_background_chi = new SBNchi(*m_background_spectrum, *m_full_fractional_covariance_matrix, this->xmlname, false);
    return 0;
}

int SBNfeld::SetBackgroundSpectrum(std::string filein, std::string scale_nam, double val){

    m_background_spectrum= new SBNosc(filein,this->xmlname);
    m_background_spectrum->ScaleAll(global_scale);
    m_background_spectrum->Scale(scale_nam,val);
    m_bool_background_spectrum_set = true;

    m_background_spectrum->CollapseVector();
    m_tvec_background_spectrum = new TVectorT<double>(m_background_spectrum->full_vector.size(), &(m_background_spectrum->full_vector)[0]);
    m_tvec_background_mcerr = new TVectorT<double>(m_background_spectrum->full_error.size(), &(m_background_spectrum->full_error)[0]);

    if(m_background_spectrum->full_vector.size() !=  m_full_fractional_covariance_matrix->GetNcols()){

        std::cout<<"SBNfeld::LoadBackgroundSpectrum || ERROR!! background spectrum is of length : "<<m_background_spectrum->full_vector.size()<<" and frac matrix is of size "<<m_full_fractional_covariance_matrix->GetNcols()<<std::endl;
        exit(EXIT_FAILURE);
    }

    m_background_chi = new SBNchi(*m_background_spectrum, *m_full_fractional_covariance_matrix, this->xmlname, false);


    return 0;
}


int SBNfeld::LoadBackgroundSpectrum(std::string filein){
    m_background_spectrum= new SBNosc(filein,this->xmlname);
    m_bool_background_spectrum_set = true;

    m_background_spectrum->CollapseVector();
    m_tvec_background_spectrum = new TVectorT<double>(m_background_spectrum->full_vector.size(), &(m_background_spectrum->full_vector)[0]);
    m_tvec_background_mcerr = new TVectorT<double>(m_background_spectrum->full_error.size(), &(m_background_spectrum->full_error)[0]);

    m_background_chi = new SBNchi(*m_background_spectrum, *m_full_fractional_covariance_matrix, this->xmlname, false);
    return 0;
}


int SBNfeld::CalcSBNchis(){
    //This is where we will calculate a vector of SBNchi's
    //i.e m_sbnchi_grid;

    TMatrixT<double> stat_only_matrix(num_bins_total,num_bins_total);
    stat_only_matrix.Zero();

    for(size_t t =0; t < m_num_total_gridpoints; t++){
        //Setup a SBNchi for this true point
        //std::cout<<"Setting up grid point SBNchi "<<t<<std::endl;
        if(m_bool_stat_only){
            m_sbnchi_grid.push_back(new SBNchi(*m_cv_spec_grid.at(t), stat_only_matrix, this->xmlname, false, m_random_seed)) ;
            m_sbnchi_grid.back()->is_stat_only = true;
        }else{
            m_sbnchi_grid.push_back(new SBNchi(*m_cv_spec_grid.at(t), *m_full_fractional_covariance_matrix, this->xmlname, false, m_random_seed)) ;
        }


    }

    return 0;
}


int SBNfeld::FullFeldmanCousins(){
    int num_universes = m_num_universes;

    //Ok take the background only spectrum and form a background only covariance matrix. CalcCovarianceMatrix includes stats
    TMatrixT<double> background_full_covariance_matrix = m_sbnchi_grid[0]->CalcCovarianceMatrix(m_full_fractional_covariance_matrix, *m_tvec_background_spectrum);
    TMatrixT<double> background_collapsed_covariance_matrix(m_background_spectrum->num_bins_total_compressed, m_background_spectrum->num_bins_total_compressed);
    m_sbnchi_grid[0]->CollapseModes(background_full_covariance_matrix, background_collapsed_covariance_matrix);    
    TMatrixT<double> inverse_background_collapsed_covariance_matrix = m_sbnchi_grid[0]->InvertMatrix(background_collapsed_covariance_matrix);   

    TFile *fout =  new TFile(("SBNfeld_output_"+tag+".root").c_str(),"recreate");
    fout->cd();

    for(size_t t =0; t < m_num_total_gridpoints; t++){

        time_t start_time = time(0);

        SBNspec * true_spec = m_cv_spec_grid.at(t);
        SBNchi  * true_chi = m_sbnchi_grid.at(t);          


        if(m_bool_stat_only)true_chi->is_stat_only = true;


        double t_val = m_grid.f_dimensions[0].GetPoint(t);

        std::cout<<"Starting on point "<<t<<"/"<<m_num_total_gridpoints<<std::endl;

        double sum = std::accumulate(true_spec->full_vector.begin(),true_spec->full_vector.end(),0.0);
        std::cout<<"Sum of events "<<sum<<std::endl;

        std::vector<double> vec_delta_chi(num_universes,0);
        std::vector<double> vec_chi_min(num_universes,0);
        std::vector<double> vec_bf_val(num_universes,0);
        std::vector<int> vec_bf_pt(num_universes,0);

        double delta_chi_critical = DBL_MIN;
        //Some output informatation
        double tree_delta_chi = 0;
        double tree_chi_min = 0;
        double tree_bf_value=0;
        int tree_bf_pt=0;

        fout->cd();
        TTree t_outtree(("ttree_"+std::to_string(t)).c_str(),("ttree_"+std::to_string(t)).c_str());
        t_outtree.Branch("delta_chi2",&tree_delta_chi);
        t_outtree.Branch("chi2_min",&tree_chi_min);
        t_outtree.Branch("bf_gridvalue",&tree_bf_value);       
        t_outtree.Branch("bf_gridpoint",&tree_bf_pt);       

        for(size_t i=0; i< num_universes; i++){


            const std::vector<float>  fake_data = true_chi->GeneratePseudoExperiment();
            std::vector<double> ans = this->PerformIterativeGridFit(fake_data,t,inverse_background_collapsed_covariance_matrix);

            vec_delta_chi[i] = ans[1]-ans[2];
            vec_chi_min[i] = ans[2];
            vec_bf_val[i] = m_grid.f_dimensions[0].GetPoint((int)ans[0]);
            vec_bf_pt[i] = (int)ans[0];

            //some root output
            tree_delta_chi = vec_delta_chi[i];
            tree_chi_min = vec_chi_min[i];
            tree_bf_value= vec_bf_val[i]; 
            tree_bf_pt = (int)ans[0];

            t_outtree.Fill();
        }

        std::cout <<"Finished universe loop for grid pt "<<t<<" at  walltime time : " << difftime(time(0), start_time) << " Seconds.\n";

        double delta_chi_min =t_outtree.GetMinimum("delta_chi2");
        double delta_chi_max =t_outtree.GetMaximum("delta_chi2");

        double bf_chi2_min =t_outtree.GetMinimum("chi2_min");
        double bf_chi2_max =t_outtree.GetMaximum("chi2_min");

        double bf_value_min = t_outtree.GetMinimum("bf_gridvalue");
        double bf_value_max = t_outtree.GetMaximum("bf_gridvalue");

        double bf_pt_min = t_outtree.GetMinimum("bf_gridpoint");
        double bf_pt_max = t_outtree.GetMaximum("bf_gridpoint");

        std::cout<<"For this point, minimum delta chi value is "<<delta_chi_min<<" max is "<<delta_chi_max<<std::endl;
        std::cout<<"For this point, minimum  chi^2 value at bf is "<<bf_chi2_min<<" max is "<<bf_chi2_max<<std::endl;
        std::cout<<"For this point, minimum  bf scale value is "<<bf_value_min<<" max is "<<bf_value_max<<std::endl;
        std::cout<<"For this point, minimum  bf grid pt is "<<bf_pt_min<<" max is "<<bf_pt_max<<std::endl;

        int nbins_plot = num_universes/4;
        std::string identifier = "_"+std::to_string(t);

        TH1D h_delta_chi(("delta_chi2"+identifier).c_str(),("delta_chi2"+identifier).c_str(),nbins_plot,0.0,delta_chi_max*1.01);  // This will store all the delta_chi's from each universe for this g_true point
        for(double&v:vec_delta_chi) h_delta_chi.Fill(v);

        //Lets grab a cumulative probabilty function
        TH1D* h_cumulative_delta_chi = (TH1D*)h_delta_chi.GetCumulative();
        h_cumulative_delta_chi->Scale(1.0/h_cumulative_delta_chi->GetBinContent(h_cumulative_delta_chi->GetNbinsX()));

        TH1D h_chi_min(("bf_chi2"+identifier).c_str(),("bf_chi2"+identifier).c_str(),nbins_plot,0.0,bf_chi2_max*1.01);  // This will store all the chi_mins from each universe for this g_true point
        for(double&v:vec_chi_min) h_chi_min.Fill(v);

        TH1D h_bf_val(("bf_value"+identifier).c_str(),("bf_value"+identifier).c_str(),nbins_plot,bf_value_min*0.8,bf_value_max*1.2);  // This will store all the bf pts from each universe for this g_true point
        for(double&v:vec_bf_val)  h_bf_val.Fill(v);

        TH1D h_bf_pt(("bf_point"+identifier).c_str(),("bf_point"+identifier).c_str(),nbins_plot, bf_pt_min*0.9, bf_pt_max*1.1);  // This will store all the bf pts from each universe for this g_true point
        for(int&v:vec_bf_pt)  h_bf_pt.Fill((double)v);

        fout->cd();
        h_delta_chi.Write();
        h_cumulative_delta_chi->Write();
        h_bf_val.Write();
        h_chi_min.Write();
        t_outtree.Write();

        //Step 3: We now have a distrubution of Delta_chi's based on many pseudo-fake data experiments. This should approximate a 
        //chi^2 distrubution with 2 DOF but there can be quite large variations. 

        double sig1 = 0.5-(0.6827)/2.0;
        double sig2 = 0.5-(0.9545)/2.0;
        std::vector<double> prob_values = {0.01,sig2,0.05,0.1,sig1,0.5,1-sig1,0.9,0.95,1-sig2,0.99};
        std::vector<double> delta_chi_quantiles(prob_values.size());	
        std::vector<double> chi_min_quantiles(prob_values.size());	
        std::vector<double> bf_val_quantiles(prob_values.size());	
        std::vector<double> bf_pt_quantiles(prob_values.size());	

        h_delta_chi.ComputeIntegral(); 
        h_delta_chi.GetQuantiles(prob_values.size(), &delta_chi_quantiles[0], &prob_values[0]);

        h_chi_min.ComputeIntegral(); 
        h_chi_min.GetQuantiles(prob_values.size(), &chi_min_quantiles[0], &prob_values[0]);

        h_bf_val.ComputeIntegral(); 
        h_bf_val.GetQuantiles(prob_values.size(), &bf_val_quantiles[0], &prob_values[0]);

        h_bf_pt.ComputeIntegral(); 
        h_bf_pt.GetQuantiles(prob_values.size(), &bf_pt_quantiles[0], &prob_values[0]);

        //Silly way to save the 3 useful values.
        //Need to transform to actual integer gridpoints, alas.
        TVectorD bfsave(3); bfsave[0]=(int)std::round(bf_pt_quantiles[4]); bfsave[1]=(int)std::round(bf_pt_quantiles[5]);bfsave[2]=(int)std::round(bf_pt_quantiles[6]);
        bfsave.Write(("median_values"+identifier).c_str()); 

        std::cout<<"Grid Point: "<<t;
        for(auto &p: m_vec_grid.at(t)){
            std::cout<<" "<<p;
        }
        std::cout<<std::endl;

        for(int i=0; i< prob_values.size(); i++){
            std::cout<<"--delta_quantile "<<delta_chi_quantiles[i]<<" chimin_quantile "<<chi_min_quantiles[i]<<" bf_val_quantil "<<bf_val_quantiles[i]<<" @prob "<<prob_values[i]<<std::endl;
        }

        delete h_cumulative_delta_chi;
        std::cout<<"Finished This Point"<<std::endl;
    }

    fout->Close();
    return 0;
};


int SBNfeld::CompareToData(SBNspec *datain){

    //Ok take the background only spectrum and form a background only covariance matrix. CalcCovarianceMatrix includes stats
    TMatrixT<double> background_full_covariance_matrix = m_sbnchi_grid[0]->CalcCovarianceMatrix(m_full_fractional_covariance_matrix, *m_tvec_background_spectrum);
    TMatrixT<double> background_collapsed_covariance_matrix(m_background_spectrum->num_bins_total_compressed, m_background_spectrum->num_bins_total_compressed);
    m_sbnchi_grid[0]->CollapseModes(background_full_covariance_matrix, background_collapsed_covariance_matrix);    
    TMatrixT<double> inverse_background_collapsed_covariance_matrix = m_sbnchi_grid[0]->InvertMatrix(background_collapsed_covariance_matrix);   

    datain->CollapseVector();
    const std::vector<float>  fake_data = datain->f_collapsed_vector;
    std::vector<double> ans = this->PerformIterativeGridFit(fake_data,0,inverse_background_collapsed_covariance_matrix,true);

    TFile *fin = new TFile(("SBNfeld_output_"+tag+".root").c_str(),"read");


    std::cout<<"------ Print DeltaChi^2 ------ "<<std::endl;
    std::vector<double> rall;
    std::vector<double> vall;
    std::vector<double> fcall;
    std::vector<double> wilkscall;

    for(int k=3; k<ans.size();k++){
        int i = k-3;
        double val = m_grid.f_dimensions[0].GetPoint(i);

        TTree *t =  (TTree*)fin->Get(("ttree_"+std::to_string(i)).c_str());

        double pval_fc = (1.0-t->GetEntries(("delta_chi2 > "+std::to_string(ans[k]-ans[2])).c_str())/(double)t->GetEntries())*100.0;
        double pval_wilks = (1.0-TMath::Prob(ans[k]-ans[2],1))*100.0;

        std::cout<<i<<" val: "<<val<<" dchi: "<<ans[k]-ans[2]<<" pval_wilks: "<<pval_wilks<<" pval_fc: "<<pval_fc<<std::endl;

        fcall.push_back(pval_fc);
        wilkscall.push_back(pval_wilks);
        rall.push_back(ans[k]-ans[2]);
        vall.push_back(val);
    }

    fin->Close();

    std::cout<<"------ DONE Print DeltaChi^2 ------ "<<std::endl;
    double delta_chi = ans[1]-ans[2];
    double chi_min = ans[2];
    double bf_val = m_grid.f_dimensions[0].GetPoint((int)ans[0]);
    double bf_pt = (int)ans[0];


    TCanvas*c = new TCanvas("c","c",1600,800);
    c->Divide(2,1);
    c->cd(1);

    TGraph *g = new TGraph(rall.size(),&vall[0],&rall[0]);
    g->SetLineColor(kBlack);
    g->SetLineWidth(2);
    g->Draw("al");
    g->SetTitle("");

    g->GetXaxis()->SetTitle("x_{#Delta} (NC #Delta Radiative Scaling)");
    g->GetYaxis()->SetTitle("#Delta #chi^{2} (data | x_{#Delta})");
    g->GetXaxis()->SetRangeUser(vall.front(),vall.back());

    TLatex * qnam = new TLatex();
    qnam->SetTextSize(0.045);
    qnam->SetTextAlign(12);  //align at top
//  qnam->SetTextAngle(-0);
    qnam->DrawLatexNDC(0.25,0.8,("Best Fit: x_{#Delta}= " +  to_string_prec(bf_val,2) ).c_str());
    qnam->DrawLatexNDC(0.25,0.75,("#chi^{2}_{Min} = " +  to_string_prec(chi_min,1) +  " ( "+ std::to_string(m_background_spectrum->num_bins_total_compressed-1)+" dof)" ).c_str());

    c->Update();

    TPad*p2 = (TPad*)c->cd(2);
    //p2->SetLogy();

    TGraph *gfc = new TGraph(rall.size(),&vall[0],&fcall[0]);
    TGraph *gwilks = new TGraph(rall.size(),&vall[0],&wilkscall[0]);
    gfc->SetLineColor(kRed-7);
    gwilks->SetLineColor(kBlue-7);
    gfc->SetLineWidth(2);
    gwilks->SetLineWidth(2);
    gwilks->Draw("al");
    gfc->Draw("l same");

    TLegend *leg = new TLegend(0.52,0.11,0.92,0.3);
    leg->AddEntry(gwilks,"Wilks Theorem","l");
    leg->AddEntry(gfc,"Feldman-Cousins","l");
    leg->SetLineWidth(0);
    leg->SetLineColor(kWhite);
    leg->SetFillStyle(0);
    leg->Draw();

    gwilks->SetTitle("");
    gwilks->GetHistogram()->SetMaximum(100.0);
    gwilks->GetXaxis()->SetRangeUser(vall.front(),vall.back());
    gwilks->GetXaxis()->SetTitle("x_{#Delta} (NC #Delta Radiative Scaling)");
    gwilks->GetYaxis()->SetTitle("Confidence Level (%)");

    TLine *l68 = new TLine(vall.front(),68.0,vall.back(),68.0);
    l68->SetLineColor(kBlack);
    l68->SetLineWidth(1);
    l68->SetLineStyle(9);
    l68->Draw();

    TLine *l90 = new TLine(vall.front(),90.0,vall.back(),90.0);
    l90->SetLineColor(kBlack);
    l90->SetLineWidth(1);;
    l90->SetLineStyle(2);
    l90->Draw();

    TLine *l95 = new TLine(vall.front(),95.0,vall.back(),95.0);
    l95->SetLineColor(kBlack);
    l95->SetLineWidth(1);
    l95->SetLineStyle(3);
    l95->Draw();

    c->SaveAs(("dataFC_"+tag+".pdf").c_str(),"pdf");



    //Some BF 
    m_cv_spec_grid[bf_pt]->CompareSBNspecs(background_collapsed_covariance_matrix,datain, "Data_Comparason_Feld_"+tag);
    std::cout<<"DATA_Comparason_Point : Delta Chi "<<delta_chi<<" Chi^Min "<<chi_min<<" BF_val "<<bf_val<<" BF_PT "<<bf_pt<<std::endl;
    return 0;
};




int SBNfeld::UpdateInverseCovarianceMatrixCNP(size_t best_grid_point, const std::vector<float> &datavec, TMatrixT<double>& inverse_collapsed, SBNchi * helper){
    TMatrixT<double> coll = helper->CalcCovarianceMatrixCNP(m_full_fractional_covariance_matrix, m_cv_spec_grid[best_grid_point]->full_vector, m_cv_spec_grid[best_grid_point]->collapsed_vector, m_cv_spec_grid[best_grid_point]->full_error ,datavec);
    inverse_collapsed = helper->InvertMatrix(coll);   
    return 0;
}


int SBNfeld::UpdateInverseCovarianceMatrix(size_t best_grid_point, TMatrixT<double>& inverse_collapsed, SBNchi * helper){
    TMatrixT<double> full = helper->CalcCovarianceMatrix(m_full_fractional_covariance_matrix, m_cv_spec_grid[best_grid_point]->full_vector, m_cv_spec_grid[best_grid_point]->full_error);
    helper->CollapseModes(full, inverse_collapsed);    
    inverse_collapsed = helper->InvertMatrix(inverse_collapsed);   
    return 0;
}


std::vector<double> SBNfeld::PerformIterativeGridFit(const std::vector<float> &datavec, const size_t grid_pt, const TMatrixT<double>& inverse_background_collapsed_covariance_matrix){
    return PerformIterativeGridFit(datavec, grid_pt, inverse_background_collapsed_covariance_matrix, false);
}
std::vector<double> SBNfeld::PerformIterativeGridFit(const std::vector<float> &datavec, const size_t grid_pt, const TMatrixT<double>& inverse_background_collapsed_covariance_matrix, bool returnall){

    double last_chi_min = DBL_MAX;
    int best_grid_point = -99;

    TMatrixT<double> inverse_current_collapsed_covariance_matrix = inverse_background_collapsed_covariance_matrix;  

    SBNchi  * grid_chi = m_sbnchi_grid.at(grid_pt); 
    std::vector<double> rall;

    for(size_t n_iter = 0; n_iter < m_max_number_iterations; n_iter++){
        rall.clear();
        //Step 1. What covariance matrix do we use?
        //For first iteration, use the precalculated background only inverse covariance matrix.
        //For all subsequent iterations what is the full covariance matrix? Use the last best grid point.
        if(n_iter!=0){
            //Calculate current full covariance matrix, collase it, then Invert. 
            if(m_use_CNP){
                UpdateInverseCovarianceMatrixCNP(best_grid_point, datavec, inverse_current_collapsed_covariance_matrix ,grid_chi);
            }else{
                UpdateInverseCovarianceMatrix(best_grid_point, inverse_current_collapsed_covariance_matrix ,grid_chi);
            }
        }

        //Step 2.0 Find the global_minimum_for this universe. Integrate in SBNfit minimizer here, a grid scan for now.
        double chi_min = DBL_MAX;
        for(size_t r =0; r < m_num_total_gridpoints; r++){

            double chi_tmp = this->CalcChi(datavec, m_cv_spec_grid[r]->collapsed_vector, inverse_current_collapsed_covariance_matrix);

            //            std::cout<<"Iter: "<<n_iter<<" "<<r<<" "<<chi_tmp<<std::endl;

            if(chi_tmp <= chi_min){
                best_grid_point = r;
                chi_min = chi_tmp;
            }
            if(returnall)rall.push_back(chi_tmp);
        }

        if(n_iter!=0){
            //std::cout<<"On iter: "<<n_iter<<" w/ chi^2: "<<chi_min<<" lastchi^2: "<<last_chi_min<<" diff() "<<fabs(chi_min-last_chi_min)<<" tol: "<<m_chi_min_convergance_tolerance<<" best_grid_point: "<<best_grid_point<<std::endl;

            //Step 3.0 Check to see if mgrid_chi for this particular fake_data  has converged sufficiently

            if(fabs(chi_min-last_chi_min)< m_chi_min_convergance_tolerance){
                last_chi_min = chi_min;
                break;
            }
        }
        last_chi_min = chi_min;
    }

    //Now use the curent_iteration_covariance matrix to also calc this_chi here for the delta.
    double this_chi = this->CalcChi(datavec, m_cv_spec_grid[grid_pt]->collapsed_vector, inverse_current_collapsed_covariance_matrix);

    //returns the BF grid, the chi^2 and the minimum_chi at the BF. 
    std::vector<double> ans = {(double)best_grid_point, this_chi, last_chi_min};

    if(returnall) ans.insert(ans.end(), rall.begin(), rall.end());

    return ans;
}

int SBNfeld::PointFeldmanCousins(size_t grid_pt){

    int num_universes = m_num_universes;

    //Ok take the background only spectrum and form a background only covariance matrix. CalcCovarianceMatrix includes stats
    TMatrixT<double> background_full_covariance_matrix = m_sbnchi_grid[0]->CalcCovarianceMatrix(m_full_fractional_covariance_matrix, *m_tvec_background_spectrum);
    TMatrixT<double> background_collapsed_covariance_matrix(m_background_spectrum->num_bins_total_compressed, m_background_spectrum->num_bins_total_compressed);
    m_sbnchi_grid[0]->CollapseModes(background_full_covariance_matrix, background_collapsed_covariance_matrix);    
    TMatrixT<double> inverse_background_collapsed_covariance_matrix = m_sbnchi_grid[0]->InvertMatrix(background_collapsed_covariance_matrix);   

    if(grid_pt > m_num_total_gridpoints){
        std::cout<<"ERROR! SBNfeld::PointFeldmanCousins() || input grid_point ("<<grid_pt<<") is above num_total_gridpoints ("<<m_num_total_gridpoints<<"). Exiting."<<std::endl;
        exit(EXIT_FAILURE);
    }

    std::cout<<"Starting on grid point "<<grid_pt<<"/"<<m_num_total_gridpoints<<std::endl;

    SBNspec * true_spec = m_cv_spec_grid[grid_pt]; 
    SBNchi  * true_chi = m_sbnchi_grid[grid_pt]; 

    for(size_t i=0; i< num_universes; i++){

        //step 0. Make a fake-data-experimet for this point, drawn from covariance.
        std::vector<float> fake_data= true_chi->SampleCovariance(true_spec);

        float last_chi_min = FLT_MAX;
        int best_grid_point = -99;

        TMatrixT<double> inverse_current_collapsed_covariance_matrix = inverse_background_collapsed_covariance_matrix;  

        size_t n_iter = 0;
        for(n_iter = 0; n_iter < m_max_number_iterations; n_iter++){

            //Step 1. What covariance matrix do we use?
            //For first iteration, use the precalculated background only inverse covariance matrix.
            //For all subsequent iterations what is the full covariance matrix? Use the last best grid point.
            if(n_iter!=0){
                //Calculate current full covariance matrix, collase it, then Invert. 
                TMatrixT<double> current_full_covariance_matrix = true_chi->CalcCovarianceMatrix(m_full_fractional_covariance_matrix,m_cv_spec_grid[best_grid_point]->full_vector);
                TMatrixT<double> current_collapsed_covariance_matrix(num_bins_total_compressed, num_bins_total_compressed);
                true_chi->CollapseModes(current_full_covariance_matrix, current_collapsed_covariance_matrix);    
                inverse_current_collapsed_covariance_matrix = true_chi->InvertMatrix(current_collapsed_covariance_matrix);   
            }

            //Step 2.0 Find the global_minimum_for this universe. Integrate in SBNfit minimizer here, a grid scan for now.
            float chi_min = FLT_MAX;

            for(size_t r =0; r < m_num_total_gridpoints; r++){
                float chi_tmp = this->CalcChi(fake_data, m_cv_spec_grid[r]->collapsed_vector,  inverse_current_collapsed_covariance_matrix);
                if(chi_tmp < chi_min){
                    best_grid_point = r;
                    chi_min = chi_tmp;
                }
            }

            if(n_iter!=0){
                //std::cout<<"On iter: "<<n_iter<<" of uni "<<i<<"/"<<num_universes<<" w/ chi^2: "<<chi_min<<" lastchi^2: "<<last_chi_min<<" diff() "<<fabs(chi_min-last_chi_min)<<" tol: "<<m_chi_min_convergance_tolerance<<" best_grid_point: "<<best_grid_point<<std::endl;

                //Step 3.0 Check to see if min_chi for this particular fake_data  has converged sufficiently
                if(fabs(chi_min-last_chi_min)< m_chi_min_convergance_tolerance){
                    last_chi_min = chi_min;
                    break;
                }
            }else{
                //std::cout<<"On iter: "<<n_iter<<" chi^2: "<<chi_min<<std::endl;
            }

            last_chi_min = chi_min;
        }

        //Now use the curent_iteration_covariance matrix to also calc this_chi here for the delta.
        float this_chi = this->CalcChi(fake_data, true_spec->collapsed_vector,inverse_current_collapsed_covariance_matrix);

        //step 4 calculate the delta_chi for this universe

        std::cout<<grid_pt<<" "<<i<<" "<<last_chi_min<<" "<<this_chi-last_chi_min<<" "<<best_grid_point<<" "<<n_iter<<std::endl;
    }

    std::cout<<"Finished Grid Point: "<<grid_pt;
    for(auto &p: m_vec_grid[grid_pt]){
        std::cout<<" "<<p;
    }
    std::cout<<std::endl;


    return 0;
};



float SBNfeld::CalcChi(const std::vector<float>& data, const std::vector<double>& prediction,const TMatrixT<double> & inverse_covariance_matrix ){
    float tchi = 0;

    for(int i =0; i<num_bins_total_compressed; i++){
        for(int j =0; j<num_bins_total_compressed; j++){
            tchi += (data[i]-prediction[i])*inverse_covariance_matrix(i,j)*(data[j]-prediction[j]);
        }
    }

    return tchi;
};


std::vector<double> SBNfeld::GlobalScan(){
    return this->GlobalScan(-1);
}

std::vector<double> SBNfeld::GlobalScan(int which_pt){

    //Ok take the background only spectrum and form a background only covariance matrix. CalcCovarianceMatrix includes stats

    std::vector<double> ans;

    TMatrixT<double> background_full_covariance_matrix = m_sbnchi_grid[0]->CalcCovarianceMatrix(m_full_fractional_covariance_matrix, *m_tvec_background_spectrum);
    TMatrixT<double> background_collapsed_covariance_matrix(m_background_spectrum->num_bins_total_compressed, m_background_spectrum->num_bins_total_compressed);
    m_sbnchi_grid[0]->CollapseModes(background_full_covariance_matrix, background_collapsed_covariance_matrix);    
    TMatrixT<double> inverse_background_collapsed_covariance_matrix = m_sbnchi_grid[0]->InvertMatrix(background_collapsed_covariance_matrix);   


    SBNspec * m_observed_spectrum;
    if(which_pt <0){
        m_observed_spectrum = m_background_spectrum;
    }else{
        m_observed_spectrum = m_cv_spec_grid.at(which_pt);   
    }

    for(size_t t =0; t < m_num_total_gridpoints; t++){

        std::cout<<"Starting on point "<<t<<"/"<<m_num_total_gridpoints<<std::endl;

        SBNspec * reco_spec = m_cv_spec_grid.at(t); 
        SBNchi  * reco_chi = m_sbnchi_grid.at(t); 
        std::vector<double> reco_spec_vec =reco_spec->full_vector;

        double chiSq;
        if(which_pt<0){
            chiSq = m_background_chi->CalcChi(reco_spec);  //exclusion
        }else{
            chiSq = m_sbnchi_grid[which_pt]->CalcChi(reco_spec);  //allowed

        }

        ans.push_back(chiSq);
        std::cout<<"ANS: "<<t<<" "<<chiSq;
        for(int k=0; k<m_vec_grid[t].size();k++){
            std::cout<<" "<<m_vec_grid[t][k];
        }

        std::cout<<std::endl;
    }



    return ans;
};


std::vector<double> SBNfeld::GlobalScan2(int which){

    return this->GlobalScan(m_cv_spec_grid[which]);
}

std::vector<double> SBNfeld::GlobalScan(SBNspec * observed_spectrum){

    std::vector<double> ans_out;

    //Ok take the background only spectrum and form a background only covariance matrix. CalcCovarianceMatrix includes stats

    TMatrixT<double> background_full_covariance_matrix = m_sbnchi_grid[0]->CalcCovarianceMatrix(m_full_fractional_covariance_matrix, *m_tvec_background_spectrum);
    TMatrixT<double> background_collapsed_covariance_matrix(m_background_spectrum->num_bins_total_compressed, m_background_spectrum->num_bins_total_compressed);
    m_sbnchi_grid[0]->CollapseModes(background_full_covariance_matrix, background_collapsed_covariance_matrix);    
    TMatrixT<double> inverse_background_collapsed_covariance_matrix = m_sbnchi_grid[0]->InvertMatrix(background_collapsed_covariance_matrix);   

    //first calculate a chi^2 minimum and BF

    observed_spectrum->CollapseVector();
    std::vector<double> ans = this->PerformIterativeGridFit(observed_spectrum->f_collapsed_vector,0,inverse_background_collapsed_covariance_matrix);
    double chi_min = ans[2];
    size_t bf = ans[0];

    std::cout<<"Best fit is point "<<bf<<" chi^2 min of "<<chi_min<<"  is point ";
    for(int k=0; k<m_vec_grid[bf].size();k++){
        std::cout<<" "<<m_vec_grid[bf][k];
    }
    std::cout<<std::endl;


    SBNspec * bf_spec = m_cv_spec_grid[bf];

    bf_spec->CollapseVector();
    TVectorT<double>* m_tvec_bf_spectrum = new TVectorT<double>(bf_spec->full_vector.size(), &(bf_spec->full_vector)[0]);


    TMatrixT<double> bf_full_covariance_matrix = m_sbnchi_grid[0]->CalcCovarianceMatrix(m_full_fractional_covariance_matrix, *m_tvec_bf_spectrum);
    TMatrixT<double> bf_collapsed_covariance_matrix(bf_spec->num_bins_total_compressed, bf_spec->num_bins_total_compressed);
    m_sbnchi_grid[0]->CollapseModes(bf_full_covariance_matrix, bf_collapsed_covariance_matrix);    

    TMatrixT<double> bf_full_covariance_matrix_NOSTAT = m_sbnchi_grid[0]->CalcCovarianceMatrix(m_full_fractional_covariance_matrix, *m_tvec_bf_spectrum,false);
    TMatrixT<double> bf_collapsed_covariance_matrix_NOSTAT(bf_spec->num_bins_total_compressed, bf_spec->num_bins_total_compressed);
    m_sbnchi_grid[0]->CollapseModes(bf_full_covariance_matrix_NOSTAT, bf_collapsed_covariance_matrix_NOSTAT);    


    TMatrixT<double> inverse_bf_collapsed_covariance_matrix = m_sbnchi_grid[0]->InvertMatrix(bf_collapsed_covariance_matrix);   


    std::vector<double> chis;
    for(size_t t =0; t < m_num_total_gridpoints; t++){
        std::cout<<"Starting on point "<<t<<"/"<<m_num_total_gridpoints<<std::endl;
        SBNchi *test_chi = m_sbnchi_grid.at(t); 
        SBNspec *test_spec = m_cv_spec_grid.at(t);
        test_spec->CollapseVector();

        //double chiSq = test_chi->CalcChi(observed_spectrum); ignore this
        double chiSq   = this->CalcChi(observed_spectrum->f_collapsed_vector, test_spec->collapsed_vector, inverse_bf_collapsed_covariance_matrix);


        double deltaChi = chiSq-chi_min;
        chis.push_back(deltaChi);

        std::cout<<"ANS: "<<t<<" "<<chiSq<<" "<<deltaChi;
        for(int k=0; k<m_vec_grid[t].size();k++){
            std::cout<<" "<<m_vec_grid[t][k];
        }
        std::cout<<std::endl;

        ans_out.push_back(chiSq);

    }


    std::vector<std::vector<double>> vgrid = m_grid.GetGrid();
    TFile *f = new TFile("margin_test.root","recreate");
    f->cd();
    for(int i=0 ; i<m_grid.f_num_dimensions; i++){
        std::vector<double> margin;
        for(auto &pt: m_grid.f_dimensions[i].f_points){
            double min_val = DBL_MAX;
            for(size_t t =0; t < m_num_total_gridpoints; t++){
                if(vgrid[t][i]==pt){
                    if(chis[t]<=min_val)min_val=chis[t];
                }
            }
            margin.push_back(min_val);
        }
        TGraph g(margin.size(),&(m_grid.f_dimensions[i].f_points)[0], &margin[0]);
        g.Write();
    }
    f->Close();

    std::cout << "Check 1" << std::endl;
    m_cv_spec_grid[bf]->CompareSBNspecs(bf_collapsed_covariance_matrix_NOSTAT,observed_spectrum, "Ans_"+tag);
    std::cout << "Check 2" << std::endl;
    return ans_out;




};



int SBNfeld::SetStatOnly(){
    m_bool_stat_only = true;
    return 0;
}

int SBNfeld::RasterScan(){
    return 0;
};
int SBNfeld::SetNumUniverses(int n){
    m_num_universes = n;
    return 0;
}  


std::vector<TGraph*> SBNfeld::MakeFCMaps(std::string filein,size_t i){
    //This is a terrible hack;
    std::vector<TGraph*> ans;

    TFile* f= new TFile(filein.c_str(),"read");
    f->cd();

    {
        size_t t = i;
        std::cout<<"On Map "<<t<<" out of "<<m_num_total_gridpoints<<std::endl;
        TTree *tre = (TTree*)f->Get( ("ttree_"+std::to_string(t)).c_str() ); 
        double max = tre->GetMaximum("delta_chi2");
        double n = tre->GetEntries();
        std::vector<double>x,y;

        for(double i=0; i<max; i+=max/10000.0 ){
            double ni = tre->GetEntries(("delta_chi2<"+std::to_string(i)).c_str());
            y.push_back(1-ni/n);
            x.push_back(i);
        }
        ans.push_back(new TGraph(x.size(),&x[0],&y[0]));
        delete tre;
    }

    f->Close();
    return ans;
}

std::vector<TGraph*> SBNfeld::MakeFCMaps(std::string filein){
    //This is a terrible hack;
    std::vector<TGraph*> ans;

    TFile* f= new TFile(filein.c_str(),"read");
    f->cd();


    for(size_t t =0; t < m_num_total_gridpoints; t++){
        std::cout<<"On Map "<<t<<" out of "<<m_num_total_gridpoints<<std::endl;
        TTree *tre = (TTree*)f->Get( ("ttree_"+std::to_string(t)).c_str() ); 
        double max = tre->GetMaximum("delta_chi2");
        double n = tre->GetEntries();
        std::vector<double>x,y;

        for(double i=0; i<max; i+=max/10000.0 ){
            double ni = tre->GetEntries(("delta_chi2<"+std::to_string(i)).c_str());
            y.push_back(1-ni/n);
            x.push_back(i);
        }
        ans.push_back(new TGraph(x.size(),&x[0],&y[0]));
    }

    f->Close();
    return ans;
}


std::vector<TGraph*> SBNfeld::LoadFCMaps(std::string filein){

    std::vector<TGraph*> ans;
    TFile* f= new TFile(filein.c_str(),"read");
    f->cd();


    for(size_t t =0; t < m_num_total_gridpoints; t++){
        std::cout<<"On Map "<<t<<" out of "<<m_num_total_gridpoints<<std::endl;
        TH1D *cumul = (TH1D*)f->Get( ("delta_chi2_"+std::to_string(t)+"_cumulative").c_str() ); 
        std::vector<double>x,y;
        for(double i=1; i<=cumul->GetNbinsX(); ++i ){
            y.push_back(1.0-cumul->GetBinContent(i));
            x.push_back(cumul->GetBinCenter(i));
        }
        ans.push_back(new TGraph(x.size(),&x[0],&y[0]));
        delete cumul;
    }

    f->Close();
    return ans;

}

std::vector<double> SBNfeld::getConfidenceRegion(TGraph *gmin, TGraph *gmax, double val){

    std::vector<double> ans =  {gmin->Eval(val), gmax->Eval(val)};

    for(auto &a:ans){
        if(a< m_grid.f_dimensions[0].f_min) a = m_grid.f_dimensions[0].f_min;
        if(a> m_grid.f_dimensions[0].f_max) a = m_grid.f_dimensions[0].f_max;
    }

    return {std::min(ans[0],ans[1]), std::max(ans[0],ans[1])};
}


