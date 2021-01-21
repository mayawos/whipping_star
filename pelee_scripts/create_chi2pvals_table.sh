#!bin/bash

#How to run this script for full syst + mc stats error:
# sh <path to the the chi2summary table> constrained_ _with_mc_err 
#How to run this script for full syst only:
# sh <path to the the chi2summary table> constrained_

#loop over hypothesis
#for x in 'H1matrix'; do

    #clear tables
    rm -f pvalstable_${2}fakedatasets_${3}.tex
    
    #create latex table formatting  
    echo "\\begin{table}[ht]" >> pvalstable_${2}fakedatasets_${3}.tex;
    echo "\\begin{center}" >> pvalstable_${2}fakedatasets_${3}.tex;
    echo "\\begin{tabular}{{|c|c|c|c|c|c|c|c|c|c|}}" >> pvalstable_${2}fakedatasets_${3}.tex;
  
    #loop over fake data sets
    for n in 1 2 3 4 5; do
	
	echo "\\multicolumn{9}{c}{\textbf{Fake Data Set ${n}}} \\\\" >> pvalstable_${2}fakedatasets_${3}.tex;
	echo "\\hline" >> pvalstable_${2}fakedatasets_${3}.tex;
	echo "\\multirow{2}{2cm}{channel} & \\multicolumn{3}{c|}{\$H_0\$} & \\multicolumn{3}{c|}{\$H_1\$} & \\multicolumn{2}{c|}{\$H_0 - H_1\$}  \\\\" >> pvalstable_${2}fakedatasets_${3}.tex;
	echo "\\cline{2-9}" >> pvalstable_${2}fakedatasets_${3}.tex;
	echo " & \$\\chi^2\$ & \$\\chi^2\$/dof & p-values & \$\\chi^2\$ & \$\\chi^2\$/dof & p-values & \$\\Delta\\chi^2\$ & p-values \\\\" >> pvalstable_${2}fakedatasets_${3}.tex;
	echo "\\hline" >> pvalstable_${2}fakedatasets_${3}.tex;
	
	#find the chisq and pvals from the frequentist study
        #loop over channels
	#for t in nue_numu 1e0p_numu nue_1e0p_numu; do 
	for t in np_numu zp_numu np_zp_numu; do 
	#for t in np_zp_numu; do 
	    
	    echo "Fake data set $n";
            echo "${1}/chi2summary_${2}${t}_reco_e_H1_mc_fakedata${n}${3}_syst_detsys.txt"	    
	    #=====================================================
	    h0_sign=`head -57 ${1}/chi2summary_${2}${t}_reco_e_H1_mc_fakedata${n}${3}_syst_detsys.txt | tail -1 | awk -F, '{printf "%.4f\n", $5}'` 
	    h1_sign=`head -58 ${1}/chi2summary_${2}${t}_reco_e_H1_mc_fakedata${n}${3}_syst_detsys.txt | tail -1 | awk -F, '{printf "%.4f\n", $5}'`
	    delta_sign=`head -59 ${1}/chi2summary_${2}${t}_reco_e_H1_mc_fakedata${n}${3}_syst_detsys.txt | tail -1 | awk -F, '{printf "%.4f\n", $5}'`
	    #====================================================
	    h0_pvals=`head -57 ${1}/chi2summary_${2}${t}_reco_e_H1_mc_fakedata${n}${3}_syst_detsys.txt | tail -1 | awk -F, '{printf "%.4f\n", $4}'` 
	    h1_pvals=`head -58 ${1}/chi2summary_${2}${t}_reco_e_H1_mc_fakedata${n}${3}_syst_detsys.txt | tail -1 | awk -F, '{printf "%.4f\n", $4}'` 
	    delta_pvals=`head -59 ${1}/chi2summary_${2}${t}_reco_e_H1_mc_fakedata${n}${3}_syst_detsys.txt | tail -1 | awk -F, '{printf "%.4f\n", $4}'` 
	    #====================================================
	    h0_chi2=`head -57 ${1}/chi2summary_${2}${t}_reco_e_H1_mc_fakedata${n}${3}_syst_detsys.txt | tail -1 | awk -F: '{printf "%.2f\n", $2}'`
	    h1_chi2=`head -58 ${1}/chi2summary_${2}${t}_reco_e_H1_mc_fakedata${n}${3}_syst_detsys.txt | tail -1 | awk -F: '{printf "%.2f\n", $2}'` 
	    delta_chi2=`head -59 ${1}/chi2summary_${2}${t}_reco_e_H1_mc_fakedata${n}${3}_syst_detsys.txt | tail -1 | awk -F: '{printf "%.2f\n", $2}'` 
	    #====================================================
	    
	    if [[ "`echo ${t}`" == "np_numu" ]]; then 
		channame="1eNp";
		dof=7;
	    elif [[ "`echo ${t}`" == "zp_numu" ]]; then
		channame="1e0p";
		dof=7;
	    elif [[ "`echo ${t}`" == "np_zp_numu" ]]; then
		channame="combined";  
		dof=14;
	    fi
	    
	    ##calculate chi-sq/dof
	    h0_chisqdof=`echo ${h0_chi2} / ${dof} | bc -l | awk '{printf "%.2f\n", $1}'`   
	    h1_chisqdof=`echo ${h1_chi2} / ${dof} | bc -l | awk '{printf "%.2f\n", $1}'`   
	    
	    ##write the actual values for each channel
	    echo "${channame} & ${h0_chi2} & ${h0_chisqdof} & ${h0_pvals} & ${h1_chi2} & ${h1_chisqdof} & ${h1_pvals} & ${delta_chi2} & ${delta_pvals} \\\\">> pvalstable_${2}fakedatasets_${3}.tex;
	    echo "\\hline" >> pvalstable_${2}fakedatasets_${3}.tex;
	    
	done #end loop over channels
	
	echo "\\multicolumn{9}{c}{} \\\\" >> pvalstable_${2}fakedatasets_${3}.tex;
	echo "\\multicolumn{9}{c}{} \\\\" >> pvalstable_${2}fakedatasets_${3}.tex;
	
    done #end loop over fake data sets
    
    #end table formatting
    echo "\\end{tabular}" >> pvalstable_${2}fakedatasets_${3}.tex;
    echo "\\end{center}" >> pvalstable_${2}fakedatasets_${3}.tex;
    echo "\\label{tab:pval_fakedata_${x}}" >> pvalstable_${2}fakedatasets_${3}.tex;
    echo "\\end{table}" >> pvalstable_${2}fakedatasets_${3}.tex;
    echo "" >> pvalstable_${2}fakedatasets_${3}.tex;
    
#done #end loop over hypothesis used in cov matrix creation
