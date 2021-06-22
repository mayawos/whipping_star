#!/bin/bash

usage="$(basename "$0") [-h] [-c -t -d -m -a -x -y -z -M -C -S] -- script to run the simple hypothesis test

where:
        -h  show this help text
        -c  input channel name of interest (options: '1eNp','1e0p','combined', default: running all channels)
    	-t  tag to use for the xml file (default: mc_collabOctFullMCunisim; multiple tags can be provided enclosed with single quotes)
        -d  add detsys error and perform the full systematics frequentist test by default
        -m  add mc stats error
        -a  add zero bin ext error
        -x  perform the stats only frequentist test
        -y  perform the flux+xsec systematics frequentist test
        -z  make the systematics tables (note:this includes creating the flux,genie,g4 only covariance table)
        -M  DO NOT make covariance matrix and spectrum
	-C  DO NOT run constraint procedure
	-S  DO NOT run sensitivity test"


#flags default
detsys=false; mcstats=false; zerobin=false; statsonly=false; syst=false; makesystable=false;
makecovariance=true; runconstraint=true; runsensitivity=true

while getopts ':dmaxyzMCSh:c:t' option; do
    case "$option" in
        h) echo "$usage"
            exit 1
            ;;
        c) channel=$OPTARG
            ;;
        t) tag=$OPTARG
           ;;
        d) detsys=true
            ;;
        m) mcstats=true
            ;;
        a) zerobin=true
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
elif [[ "`echo ${channel}`" == "1e0p" ]]
then
    sample="zp_numu_reco_e"
elif [[ "`echo ${channel}`" == "combined" ]]
then
    sample="np_zp_numu_reco_e"
fi

tag=""
if [[ "`echo ${tag}`" == "" ]];
then
    tag="mc_collabOctFullMCunisim"; 
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
	    ./sbnfit_make_covariance -x numu_reco_e_H1_DATA_${s}.xml -t numu_reco_e_H1_DATA_${s} #-p
	fi #end make covariance
	
	if [[ "$runconstraint" == "true" ]]; then
	    #go to the build/examples directory to run CreateConstrainedMatrix.cxx
	    cd ../examples/
	    make
            
            #Constraint procedures
            #Flags explanation:
	    #  -n, -z, -c : channels: 1eNp+numu, 1e0p+numu, 1eNp+1e0p+numu
	    #  -d : add PeLEE detsys error to the diagonals of full covariance matrix, then collapse, and constrain
	    #  -m : add mc stats error to the diagonal of the  full covariance matrix, then collapse, and constrain
	    #  -a : add zero ext bin error to the diagonal of the full covariance matrix, then collapse, and constrain
	    if [[ "`echo $t`" == "np_numu_reco_e" ]]; then
		./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_${s}.SBNspec.root -n  
		if [[ "$mcstats" == "true" ]]; then
                    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_${s}.SBNspec.root -n -m;
	            if [[ "$zerobin" == "true" ]]; then
			./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_${s}.SBNspec.root -n -m -a;
                    fi;
		fi;
		if [[ "$detsys" == "true"  ]]; then
                    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_${s}.SBNspec.root -n -d; 
	            if [[ "$mcstats" == "true" ]]; then
			./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_${s}.SBNspec.root -n -d -m;
	                if [[ "$zerobin" == "true" ]]; then
                            ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_${s}.SBNspec.root -n -d -m -a;
			fi;
                    fi;
		fi;
	    elif [[ "`echo $t`" == "zp_numu_reco_e" ]]; then
		./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_${s}.SBNspec.root -z 
		if [[ "$mcstats" == "true" ]]; then
                    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_${s}.SBNspec.root -z -m; 
	            if [[ "$zerobin" == "true" ]]; then
			./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_${s}.SBNspec.root -z -m -a;
                    fi;
		fi;
	        if [[ "$detsys" == "true"  ]]; then
                    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_${s}.SBNspec.root -z -d; 
	            if [[ "$mcstats" == "true" ]]; then
			./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_${s}.SBNspec.root -z -d -m; 
	                if [[ "$zerobin" == "true" ]]; then
                            ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_${s}.SBNspec.root -z -d -m -a;
			fi;
                    fi;
		fi;
	    elif [[ "`echo $t`" == "np_zp_numu_reco_e" ]]; then
		./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_${s}.SBNspec.root -c
		if [[ "$makesystable" == "true" ]]; then
		    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}_fluxonly.xml -t ${t}_H1_${s}_fluxonly -f ../bin/numu_reco_e_H1_DATA_${s}.SBNspec.root -c 
		    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}_genieonly.xml -t ${t}_H1_${s}_genieonly -f ../bin/numu_reco_e_H1_DATA_${s}.SBNspec.root -c 
		    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}_g4only.xml -t ${t}_H1_${s}_g4only -f ../bin/numu_reco_e_H1_DATA_${s}.SBNspec.root -c
		fi 
		if [[ "$mcstats" == "true" ]]; then
                    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_${s}.SBNspec.root -c -m;
	            if [[ "$zerobin" == "true" ]]; then
			./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_${s}.SBNspec.root -c -m -a;
                    fi;
		fi;
	       if [[ "$detsys" == "true"  ]]; then
                    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_${s}.SBNspec.root -c -d;
	            if [[ "$mcstats" == "true" ]]; then
			./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_${s}.SBNspec.root -c -d -m;
			if [[ "$zerobin" == "true" ]]; then
                            ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_${s}.SBNspec.root -c -d -m -a;
			    if [[ "$makesystable" == "true" ]]; then
				#create the systematics table
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
	    cp ../examples/unconstrained_${t}_H1_${s}.SBNspec.root unconstrained_${t}_H1_${s}_CV.SBNspec.root
	    cp ../examples/unconstrained_${t}_H1_${s}${m}.SBNspec.root unconstrained_${t}_H1_${s}${m}_CV.SBNspec.root
	    cp ../examples/unconstrained_${t}_H1_${s}${m}${a}.SBNspec.root unconstrained_${t}_H1_${s}${m}${a}_CV.SBNspec.root
	    cp ../examples/unconstrained_${t}_H1_${s}.SBNcovar.root .
	    cp ../examples/unconstrained_${t}_H1_${s}${m}.SBNcovar.root .
	    cp ../examples/unconstrained_${t}_H1_${s}${m}${a}.SBNcovar.root .
	    cp ../examples/unconstrained_${t}_H1_${s}_detsys.SBNcovar.root .
	    cp ../examples/unconstrained_${t}_H1_${s}${m}_detsys.SBNcovar.root .
	    cp ../examples/unconstrained_${t}_H1_${s}${m}${a}_detsys.SBNcovar.root .
	    cp ../examples/constrained_${t}_H1_${s}.SBNspec.root constrained_${t}_H1_${s}_CV.SBNspec.root
	    cp ../examples/constrained_${t}_H1_${s}${m}.SBNspec.root constrained_${t}_H1_${s}${m}_CV.SBNspec.root
	    cp ../examples/constrained_${t}_H1_${s}${m}${a}.SBNspec.root constrained_${t}_H1_${s}${m}${a}_CV.SBNspec.root
	    cp ../examples/constrained_${t}_H1_${s}.SBNcovar.root .
	    cp ../examples/constrained_${t}_H1_${s}${m}.SBNcovar.root .
	    cp ../examples/constrained_${t}_H1_${s}${m}${a}.SBNcovar.root .
	    cp ../examples/constrained_${t}_H1_${s}_detsys.SBNcovar.root .
	    cp ../examples/constrained_${t}_H1_${s}${m}_detsys.SBNcovar.root .
	    cp ../examples/constrained_${t}_H1_${s}${m}${a}_detsys.SBNcovar.root .
            
            #scale the lee channels to 0.0 for the different configuration
	    ./sbnfit_scale_spec -x const_${t}_H1_${s}.xml --input constrained_${t}_H1_${s}_CV.SBNspec.root --scalestring "lee" -v 0.0 -t constrained_${t}_H1_${s}_BKG
	    ./sbnfit_scale_spec -x const_${t}_H1_${s}.xml --input constrained_${t}_H1_${s}${m}_CV.SBNspec.root --scalestring "lee" -v 0.0 -t constrained_${t}_H1_${s}${m}_BKG
	    ./sbnfit_scale_spec -x const_${t}_H1_${s}.xml --input constrained_${t}_H1_${s}${m}${a}_CV.SBNspec.root --scalestring "lee" -v 0.0 -t constrained_${t}_H1_${s}${m}${a}_BKG
            
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
		if [[ "$mcstats" == "true" && "$zerobin" == "true"  ]]; then ./sbnfit_lee_frequentist_study --xml const_${t}_H1_${s}.xml --signal constrained_${t}_H1_${s}${m}${a}_CV.SBNspec.root --background constrained_${t}_H1_${s}${m}${a}_BKG.SBNspec.root --covariance constrained_${t}_H1_${s}${m}${a}.SBNcovar.root -e 1e-8 --tag constrained_${t}_H1_${s}${m}${a}_syst; fi
	    fi
	    #sensitivity test for flux+xsec+detsys. Order: no additional error, with mc stats error only, with mc stats error + zero bin error, with zero bin error only
	    if [[ "$detsys" == "true" ]]; then
		./sbnfit_lee_frequentist_study --xml const_${t}_H1_${s}.xml --signal constrained_${t}_H1_${s}_CV.SBNspec.root --background constrained_${t}_H1_${s}_BKG.SBNspec.root --covariance constrained_${t}_H1_${s}_detsys.SBNcovar.root -e 1e-8 --tag constrained_${t}_H1_${s}_syst_detsys; 
		if [[ "$mcstats" == "true" && "$zerobin" == "false" ]]; then ./sbnfit_lee_frequentist_study --xml const_${t}_H1_${s}.xml --signal constrained_${t}_H1_${s}${m}_CV.SBNspec.root --background constrained_${t}_H1_${s}${m}_BKG.SBNspec.root --covariance constrained_${t}_H1_${s}${m}_detsys.SBNcovar.root -e 1e-8 --tag constrained_${t}_H1_${s}${m}_syst_detsys; fi
		if [[ "$mcstats" == "true" && "$zerobin" == "true" ]]; then ./sbnfit_lee_frequentist_study --xml const_${t}_H1_${s}.xml --signal constrained_${t}_H1_${s}${m}${a}_CV.SBNspec.root --background constrained_${t}_H1_${s}${m}${a}_BKG.SBNspec.root --covariance constrained_${t}_H1_${s}${m}${a}_detsys.SBNcovar.root -e 1e-8 --tag constrained_${t}_H1_${s}${m}${a}_syst_detsys; fi
	    fi
	fi #end run sensitivity test
	
    done
done
