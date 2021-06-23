#!/bin/bash

usage="$(basename "$0") [-h] [-c -t -e -d -m -a -l -g -x -y -z -M -C -S] -- script to run the simple hypothesis test

where:
        -h  show this help text
        -c  input channel name of interest (options: '1eNp','1e0p','combined', default: running all channels)
    	-t  tag to use for the xml file (default: mc_collabOctFullMCunisim; multiple tags can be provided enclosed with single quotes)
    	-e  execute constraint (default: CreateConstrainedMatrix, other option: CreateConstrainedFarSBMatrix)
        -d  add detsys error and perform the full systematics frequentist test by default
        -m  add mc stats error
        -a  add zero bin ext error
        -l  constraint mode (1: numu only constraint; 2: far SB only constraint; 3: far SB + numu constraint)
        -g  perform gof study
        -x  perform the stats only frequentist test
        -y  perform the flux+xsec systematics frequentist test
        -z  make the systematics tables (note:this includes creating the flux,genie,g4 only covariance table)
        -M  DO NOT make covariance matrix and spectrum
	-C  DO NOT run constraint procedure
	-S  DO NOT run sensitivity test"


#flags default
detsys=false; mcstats=false; zerobin=false; statsonly=false; syst=false; makesystable=false; gof=false;
makecovariance=true; runconstraint=true; runsensitivity=true; constrmode=1;

while getopts ':dmaxyzgMCSh:c:t:e:l:' option; do
    case "$option" in
        h) echo "$usage"
            exit 1
            ;;
        c) channel=$OPTARG
            ;;
        t) tag=$OPTARG
           ;;
        e) exe=$OPTARG
            ;;
        l) constrmode=$OPTARG
            ;;
        d) detsys=true
            ;;
        m) mcstats=true
            ;;
        a) zerobin=true
            ;;
        g) gof=true
            ;;
        x) statsonly=true
            ;;
        y) syst=true
            ;;
        z) makesystable=true
            ;;
        M) makecovariance=false
            ;;
        C) runconstraint=false
            ;;
        S) runsensitivity=false
            ;;
        :) printf "missing argument for -%s\n" "$OPTARG" >&2
            echo "$usage" >&2
            exit 1
            ;;
        \?) printf "illegal option: -%s\n" "$OPTARG" >&2
            echo "$usage" >&2
            exit 1
            ;;
    esac
done

if [ $OPTIND -eq 1 ]; then 
    echo "NO ARGUMENTS ARE PASSED!" >&2 ; 
    echo "See the help instructions" >&2; 
    echo "" >&2
    echo "$usage" >&2;
    exit 1
fi

shift $((OPTIND-1))

#define the default samples and tag
sample="np_numu_reco_e zp_numu_reco_e np_zp_numu_reco_e"
#channel names
if [[ "`echo ${channel}`" == "1eNp" ]]
then
    sample="np_numu_reco_e"
    nm="np_numu"
elif [[ "`echo ${channel}`" == "1e0p" ]]
then
    sample="zp_numu_reco_e"
    nm="zp_numu"
elif [[ "`echo ${channel}`" == "combined" ]]
then
    sample="np_zp_numu_reco_e"
    nm="np_zp_numu"
elif [[ "`echo ${channel}`" == "numu0p" ]]
then
    sample="zp_numu0p_reco_e"
    nm="numu0p"
elif [[ "`echo ${channel}`" == "numuNp" ]]
then
    sample="np_numuNp_reco_e"
    nm="numuNp"
fi

echo "$tag"
if [[ "`echo ${tag}`" == "" ]];
then
    tag="mc_collabOctFullMCunisim"; 
fi

echo "$exe"
if [[ "`echo ${exe}`" == "" ]];
then
    exe="CreateConstrainedMatrix"; 
fi

for t in ${sample}
do
    for s in ${tag} #this tag identify the xml with the correct treatment of the unisim error
    do
        #tag used to identify different configuration for systematics tables, 
	#sensitivity plots, and logfiles of sensitivity/pvals results
        m="_with_mc_err"
	a="_with_zerobin_err"
	
	if [[ "$makecovariance" == "true" ]]; then
	    make
	    ./sbnfit_make_covariance -x ${t}_H1_${s}.xml -t ${t}_H1_${s} #-p
	    if [[ "$makesystable" == "true" ]]; then
		./sbnfit_make_covariance -x ${t}_H1_${s}_fluxonly.xml -t ${t}_H1_${s}_fluxonly #-p
		./sbnfit_make_covariance -x ${t}_H1_${s}_genieonly.xml -t ${t}_H1_${s}_genieonly #-p
		./sbnfit_make_covariance -x ${t}_H1_${s}_g4only.xml -t ${t}_H1_${s}_g4only #-p
	    fi
	    ./sbnfit_make_spec -x ${t}_H1_DATA_${s}.xml -t ${t}_H1_DATA_${s} #-p
            mv ${t}_H1_DATA_${s}_CV.SBNspec.root ${t}_H1_DATA_${s}.SBNspec.root
            ./sbnfit_make_spec -x constrained_${t}_H1_DATA_${s}.xml -t constrained_${t}_H1_DATA_${s}
            cp constrained_${t}_H1_DATA_${s}_CV.SBNspec.root /uboone/data/users/wospakrk/PeLEE/SBNfit/Unblinding/
	fi #end make covariance
	
	if [[ "$runconstraint" == "true" ]]; then
	    #go to the build/examples directory to run ${exe}.cxx
	    cd ../examples/
	    make
            
            #Constraint procedures
            #Flags explanation:
	    #  -n, -z, -c : channels: 1eNp+numu, 1e0p+numu, 1eNp+1e0p+numu
	    #  -d : add PeLEE detsys error to the diagonals of full covariance matrix, then collapse, and constrain
	    #  -m : add mc stats error to the diagonal of the  full covariance matrix, then collapse, and constrain
	    #  -a : add zero ext bin error to the diagonal of the full covariance matrix, then collapse, and constrain
	    if [[ "`echo $t`" == "np_numu_reco_e" || "`echo $t`" == "np_numuNp_reco_e" ]]; then
		if [[ "$detsys" == "true"  ]]; then
	            if [[ "$mcstats" == "true" ]]; then
	                if [[ "$zerobin" == "true" ]]; then
                            ./${exe} -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/${t}_H1_DATA_${s}.SBNspec.root -n -d -m -a -l ${constrmode};
                            mv mc_stats_err.tex ${t}_${s}${m}${a}_detsys.tex
			fi;
                    fi;
		fi;
	    elif [[ "`echo $t`" == "zp_numu_reco_e" || "`echo $t`" == "zp_numu0p_reco_e" ]]; then
	        if [[ "$detsys" == "true"  ]]; then
	            if [[ "$mcstats" == "true" ]]; then
	                if [[ "$zerobin" == "true" ]]; then
                            ./${exe} -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/${t}_H1_DATA_${s}.SBNspec.root -z -d -m -a -l ${constrmode};
                            mv mc_stats_err.tex ${t}_${s}${m}${a}_detsys.tex
			fi;
                    fi;
		fi;
	    elif [[ "`echo $t`" == "np_zp_numu_reco_e" ]]; then
		if [[ "$makesystable" == "true" ]]; then
		    ./${exe} -x ../bin/${t}_H1_${s}_fluxonly.xml -t ${t}_H1_${s}_fluxonly -f ../bin/${t}_H1_DATA_${s}.SBNspec.root -c 
		    ./${exe} -x ../bin/${t}_H1_${s}_genieonly.xml -t ${t}_H1_${s}_genieonly -f ../bin/${t}_H1_DATA_${s}.SBNspec.root -c 
		    ./${exe} -x ../bin/${t}_H1_${s}_g4only.xml -t ${t}_H1_${s}_g4only -f ../bin/${t}_H1_DATA_${s}.SBNspec.root -c
		fi 
	       if [[ "$detsys" == "true"  ]]; then
	            if [[ "$mcstats" == "true" ]]; then
			if [[ "$zerobin" == "true" ]]; then
                            ./${exe} -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/${t}_H1_DATA_${s}.SBNspec.root -c -d -m -a -l ${constrmode};
                            ls -lrt mc_stats_err.tex
                            mv mc_stats_err.tex mc_stats_${t}_${s}${m}${a}_detsys.tex
			    if [[ "$makesystable" == "true" ]]; then
				#create the systematics table
				echo "sh ../../pelee_scripts/make_beforeconstraint_table.sh 1eNp ${t} ${s}"
				sh ../../pelee_scripts/make_beforeconstraint_table.sh 1eNp ${t} ${s} 
				sh ../../pelee_scripts/make_beforeconstraint_table.sh 1e0p ${t} ${s} 
				sh ../../pelee_scripts/make_beforeconstraint_table.sh numu ${t} ${s}
			    fi
			    #create the before-after systematics table 
			    sh ../../pelee_scripts/make_beforeafter_constraint_table.sh 1eNp ${t} ${s} ${m} ${a} 
			    sh ../../pelee_scripts/make_beforeafter_constraint_table.sh 1e0p ${t} ${s} ${m} ${a}
			fi;
                    fi;
		fi;
	   fi    
	fi #end run constraint procedure
	
	if [[ "$runsensitivity" == "true" ]]; then	
            #move to bin directory to write the sensitivity results
	    cd ../bin
	    make
	    
            #copy the constrained/unconstrained Covariance Matrices and Histograms/Spectras to ../bin for sensitivity test
	    cp ../examples/collapsed_zero_ext_bin_${t}_H1_${s}${m}${a}.SBNcovar.root .
	    cp ../examples/unconstrained_${t}_H1_${s}${m}${a}_detsys.SBNspec.root unconstrained_${t}_H1_${s}${m}${a}_detsys_CV.SBNspec.root
	    cp unconstrained_${t}_H1_${s}${m}${a}_detsys_CV.SBNspec.root /uboone/data/users/wospakrk/PeLEE/SBNfit/Unblinding/
	    cp ../examples/unconstrained_${t}_H1_${s}${m}${a}_detsys.SBNcovar.root .
	    cp unconstrained_${t}_H1_${s}${m}${a}_detsys.SBNcovar.root /uboone/data/users/wospakrk/PeLEE/SBNfit/Unblinding/
	    cp ../examples/constrained_${t}_H1_${s}${m}${a}_detsys.SBNspec.root constrained_${t}_H1_${s}${m}${a}_detsys_CV.SBNspec.root
	    cp constrained_${t}_H1_${s}${m}${a}_detsys_CV.SBNspec.root /uboone/data/users/wospakrk/PeLEE/SBNfit/Unblinding/
	    cp ../examples/constrained_${t}_H1_${s}${m}${a}.SBNspec.root constrained_${t}_H1_${s}${m}${a}_CV.SBNspec.root
	    cp constrained_${t}_H1_${s}${m}${a}_detsys_CV.SBNspec.root /uboone/data/users/wospakrk/PeLEE/SBNfit/Unblinding/constrained_${t}_H1_${s}${m}${a}_CV.SBNspec.root
	    cp ../examples/constrained_${t}_H1_${s}${m}${a}.SBNcovar.root .
	    cp constrained_${t}_H1_${s}${m}${a}.SBNcovar.root /uboone/data/users/wospakrk/PeLEE/SBNfit/Unblinding/
	    cp ../examples/constrained_${t}_H1_${s}${m}${a}_detsys.SBNcovar.root .
	    cp constrained_${t}_H1_${s}${m}${a}_detsys.SBNcovar.root /uboone/data/users/wospakrk/PeLEE/SBNfit/Unblinding/
            
            #scale the lee channels to 0.0 for the different configuration
	    #./sbnfit_scale_spec -x const_${t}_H1_${s}.xml --input constrained_${t}_H1_${s}${m}${a}_CV.SBNspec.root --scalestring "lee" -v 0.0 -t constrained_${t}_H1_${s}${m}${a}_BKG
	    ./sbnfit_scale_spec -x const_${t}_H1_${s}.xml --input constrained_${t}_H1_${s}${m}${a}_detsys_CV.SBNspec.root --scalestring "lee" -v 0.0 -t constrained_${t}_H1_${s}${m}${a}_detsys_BKG
	    #./sbnfit_scale_spec -x const_${t}_H1_${s}.xml --input unconstrained_${t}_H1_${s}${m}${a}_CV.SBNspec.root --scalestring "lee" -v 0.0 -t unconstrained_${t}_H1_${s}${m}${a}_BKG
	    ./sbnfit_scale_spec -x const_${t}_H1_${s}.xml --input unconstrained_${t}_H1_${s}${m}${a}_detsys_CV.SBNspec.root --scalestring "lee" -v 0.0 -t unconstrained_${t}_H1_${s}${m}${a}_detsys_BKG
            
            echo "********* LABEL: ${s}${m}${a}"
            echo "********* PARAMS syst detsys gof mcstats zerobin: $syst $detsys $gof $mcstats $zerobin"
            #sensitivity test for stats only. Order: no additional error, with mc stats error only, with mc stats error + zero bin error, with zero bin error only
	    if [[ "$statsonly" == "true" ]]; then
	        ./sbnfit_lee_frequentist_study --xml const_${t}_H1_${s}.xml --signal constrained_${t}_H1_${s}_CV.SBNspec.root --background constrained_${t}_H1_${s}_BKG.SBNspec.root --poisson --tag constrained_${t}_H1_${s}_stats 
	        if [[ "$mcstats" == "true" && "$zerobin" == "false" ]]; then ./sbnfit_lee_frequentist_study --xml const_${t}_H1_${s}.xml --signal constrained_${t}_H1_${s}${m}_CV.SBNspec.root --background constrained_${t}_H1_${s}${m}_BKG.SBNspec.root --poisson --tag constrained_${t}_H1_${s}${m}_stats; fi 
	        if [[ "$mcstats" == "true" && "$zerobin" == "true" ]]; then ./sbnfit_lee_frequentist_study --xml const_${t}_H1_${s}.xml --signal constrained_${t}_H1_${s}${m}${a}_CV.SBNspec.root --background constrained_${t}_H1_${s}${m}${a}_BKG.SBNspec.root --poisson --tag constrained_${t}_H1_${s}${m}${a}_stats; fi 
	    fi
	    
	    #sensitivity test for flux+xsec only. Order: no additional error, with mc stats error only, with mc stats error + zero bin error, with zero bin error only
	    if [[ "$syst" == "true" && "$detsys" == "false" ]]; then
		./sbnfit_lee_frequentist_study --xml const_${t}_H1_${s}.xml --signal constrained_${t}_H1_${s}_CV.SBNspec.root --background constrained_${t}_H1_${s}_BKG.SBNspec.root --covariance constrained_${t}_H1_${s}.SBNcovar.root -e 1e-8 --tag constrained_${t}_H1_${s}_syst
		if [[ "$mcstats" == "true" && "$zerobin" == "false" ]]; then ./sbnfit_lee_frequentist_study --xml const_${t}_H1_${s}.xml --signal constrained_${t}_H1_${s}${m}_CV.SBNspec.root --background constrained_${t}_H1_${s}${m}_BKG.SBNspec.root --covariance constrained_${t}_H1_${s}${m}.SBNcovar.root -e 1e-8 --tag constrained_${t}_H1_${s}${m}_syst; fi
                if [[ "$gof" == true ]]; then
		if [[ "$mcstats" == "true" && "$zerobin" == "true"  ]]; then 
			./sbnfit_lee_frequentist_study --xml const_${t}_H1_${s}.xml --signal constrained_${t}_H1_${s}${m}${a}_CV.SBNspec.root --background constrained_${t}_H1_${s}${m}${a}_BKG.SBNspec.root --covariance constrained_${t}_H1_${s}${m}${a}.SBNcovar.root --fakedata constrained_${t}_H1_DATA_${s}_CV.SBNspec.root -e 1e-8 --tag constrained_${t}_H1_${s}${m}${a}_syst;
			./sbnfit_lee_frequentist_study --xml const_${t}_H1_${s}.xml --signal unconstrained_${t}_H1_${s}${m}${a}_CV.SBNspec.root --background unconstrained_${t}_H1_${s}${m}${a}_BKG.SBNspec.root --covariance unconstrained_${t}_H1_${s}${m}${a}.SBNcovar.root --fakedata constrained_${t}_H1_DATA_${s}_CV.SBNspec.root -e 1e-8 --tag unconstrained_${t}_H1_${s}${m}${a}_syst;
                fi
                fi
	    fi
	    #sensitivity test for flux+xsec+detsys. with mc stats error + zero bin error, with zero bin error only
	    if [[ "$detsys" == "true" ]]; then
                if [[ "$gof" == true ]]; then
	          if [[ "$mcstats" == "true" && "$zerobin" == "true" ]]; then ./sbnfit_lee_frequentist_study --xml const_${t}_H1_${s}.xml --signal constrained_${t}_H1_${s}${m}${a}_detsys_CV.SBNspec.root --background constrained_${t}_H1_${s}${m}${a}_detsys_BKG.SBNspec.root --covariance constrained_${t}_H1_${s}${m}${a}_detsys.SBNcovar.root --fakedata constrained_${t}_H1_DATA_${s}_CV.SBNspec.root -e 1e-8 --tag constrained_${t}_H1_${s}${m}${a}_syst_detsys; fi
	          #if [[ "$mcstats" == "true" && "$zerobin" == "true" ]]; then ./sbnfit_lee_frequentist_study --xml const_${t}_H1_${s}.xml --signal unconstrained_${t}_H1_${s}${m}${a}_CV.SBNspec.root --background unconstrained_${t}_H1_${s}${m}${a}_BKG.SBNspec.root --covariance unconstrained_${t}_H1_${s}${m}${a}_detsys.SBNcovar.root --fakedata constrained_${t}_H1_${s}_DATA_CV.SBNspec.root -e 1e-8 --tag unconstrained_${t}_H1_${s}${m}${a}_syst_detsys; fi
                fi
                if [[ "$gof" == false ]]; then
	          if [[ "$mcstats" == "true" && "$zerobin" == "true" ]]; then ./sbnfit_lee_frequentist_study --xml const_${t}_H1_${s}.xml --signal constrained_${t}_H1_${s}${m}${a}_detsys_CV.SBNspec.root --background constrained_${t}_H1_${s}${m}${a}_detsys_BKG.SBNspec.root --covariance constrained_${t}_H1_${s}${m}${a}_detsys.SBNcovar.root -e 1e-8 --tag constrained_${t}_H1_${s}${m}${a}_syst_detsys; fi
	          #if [[ "$mcstats" == "true" && "$zerobin" == "true" ]]; then ./sbnfit_lee_frequentist_study --xml const_${t}_H1_${s}.xml --signal unconstrained_${t}_H1_${s}${m}${a}_detsys_CV.SBNspec.root --background unconstrained_${t}_H1_${s}${m}${a}_detsys_BKG.SBNspec.root --covariance unconstrained_${t}_H1_${s}${m}${a}_detsys.SBNcovar.root -e 1e-8 --tag unconstrained_${t}_H1_${s}${m}${a}_syst_detsys; fi
                fi
	    fi
	fi #end run sensitivity test
	
    done
done
