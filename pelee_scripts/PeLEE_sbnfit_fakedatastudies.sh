#!/bin/bash

usage="$(basename "$0") [-h] [-c -n -k -s -f] -- script to run the simple hypothesis test or FC calculation

where:
        -h  show this help text
        -c  input channel name of interest (options: '1eNp','1e0p','combined')
        -n  which fake data set? (if more than one set, need to be enclosed in single quotes)
        -k  keep the steps to make the SBNspec and SBNcovar
        -s  run the sbnfit_lee_frequentist test for the simple hypothesis test
        -f  run the sbnfit_uboone_scaling_fc for the Feldman-Cousins calculation"

feldman=false; keep=false; simple=false

while getopts ':hfks:c:n:' option; do
        case "$option" in
    h) echo "$usage"
            exit
            ;;
    c) channel=$OPTARG
            ;;
    n) number=$OPTARG
            ;;
    f) feldman=true
            ;;
    k) keep=true
            ;;
    s) simple=true
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

shift $((OPTIND - 1))

#channel names
if [[ "`echo ${channel}`" == "1eNp" ]]
then
    sample="Combined_1eNp_numu"
elif [[ "`echo ${channel}`" == "1e0p" ]] 
then
    sample="Combined_1e0p_numu"
elif [[ "`echo ${channel}`" == "combined" ]] 
then
    sample="Combined_1eNp_1e0p_numu"
fi

#variable names
var="Reconstructed Visible Energy [MeV]"
#if [[ "`echo ${tag}`" == "trueE" ]] 
#then
#    var="True Neutrino Energy [MeV]"
#elif [[ "`echo ${tag}`" == "leptonE" ]] 
#then
#    var="Reconstructed Lepton Energy [MeV]"
#fi

for n in $number; do

    mkdir -p /uboone/data/users/${USER}/logs

    if [[ "$keep" == "true" ]]
    then
        echo "Making the spectra and covariance matrices as inputs to the test statistics"
    fi
    if [[ "$simple" == "true" ]]
    then
        echo "Will run the simple LEE hypothesis test!"
    fi
    if [[ "$feldman" == "true" ]]
    then
        echo "Will run the Feldman-Cousins Signal Scaling!"
    fi

    make
    
    if [[ "$keep" == "true" ]]
    then
        ./sbnfit_make_covariance -x ../../xml/PELEE/Fakedata/${sample}_reco_energy_fakedata${n}_MC.xml -t ${sample}_reco_energy_fakedata${n}_MC
        ./sbnfit_make_covariance -x ../../xml/PELEE/Fakedata/${sample}_reco_energy_fakedata${n}_DATA.xml -t ${sample}_reco_energy_fakedata${n}_DATA
        ./sbnfit_make_covariance -x ../../xml/PELEE/Fakedata/${sample}_reco_energy_fakedata${n}_nolee_MC.xml -t ${sample}_reco_energy_fakedata${n}_nolee_MC
        ./sbnfit_make_spec -x ../../xml/PELEE/Fakedata/${sample}_reco_energy_fakedata${n}_MC.xml -t ${sample}_reco_energy_fakedata${n}_MC
        ./sbnfit_make_spec -x ../../xml/PELEE/Fakedata/${sample}_reco_energy_fakedata${n}_DATA.xml -t ${sample}_reco_energy_fakedata${n}_DATA
    
        ./sbnfit_scale_spec -x ../../xml/PELEE/Fakedata/${sample}_reco_energy_fakedata${n}_MC.xml --input ${sample}_reco_energy_fakedata${n}_MC.SBNspec.root --scalestring lee -v 0.0 -t ${sample}_reco_energy_fakedata${n}_MC_BKG
        ./sbnfit_scale_spec -x ../../xml/PELEE/Fakedata/${sample}_reco_energy_fakedata${n}_nolee_MC.xml --input ${sample}_reco_energy_fakedata${n}_nolee_MC.SBNspec.root --scalestring lee -v 0.0 -t ${sample}_reco_energy_fakedata${n}_nolee_MC_BKG
    fi

    if [[ "$simple" == "true" ]]
    then
        echo "will do simple hypothesis testing..."
        [ -e /uboone/data/users/${USER}/logs/${sample}_reco_energy_fakedata${n}_StatsOnly.log ] && rm /uboone/data/users/${USER}/logs/${sample}_reco_energy_fakedata${n}_StatsOnly.log
        ./sbnfit_lee_frequentist_study --xml ../../xml/PELEE/Fakedata/${sample}_reco_energy_fakedata${n}_MC.xml --signal ${sample}_reco_energy_fakedata${n}_MC.SBNspec.root --background ${sample}_reco_energy_fakedata${n}_MC_BKG.SBNspec.root --fakedata ${sample}_reco_energy_fakedata${n}_DATA.SBNspec.root --poisson --tag ${sample}_reco_energy_fakedata${n}_StatsOnly >& /uboone/data/users/${USER}/logs/${sample}_reco_energy_fakedata${n}_StatsOnly.log
        [ -e /uboone/data/users/${USER}/logs/${sample}_reco_energy_fakedata${n}_MC_H0matrix.log ] && rm /uboone/data/users/${USER}/logs/${sample}_reco_energy_fakedata${n}_MC_H0matrix.log
        ./sbnfit_lee_frequentist_study --xml ../../xml/PELEE/Fakedata/${sample}_reco_energy_fakedata${n}_MC.xml --signal ${sample}_reco_energy_fakedata${n}_MC.SBNspec.root --background ${sample}_reco_energy_fakedata${n}_MC_BKG.SBNspec.root --fakedata ${sample}_reco_energy_fakedata${n}_DATA.SBNspec.root --covariance ${sample}_reco_energy_fakedata${n}_nolee_MC.SBNcovar.root -e 1e-8 -f 20 --tag ${sample}_reco_energy_fakedata${n}_H0matrix >& /uboone/data/users/${USER}/logs/${sample}_reco_energy_fakedata${n}_MC_H0matrix.log
	    [ -e /uboone/data/users/${USER}/logs/${sample}_reco_energy_fakedata${n}_MC_H1matrix.log ] && rm /uboone/data/users/${USER}/logs/${sample}_reco_energy_fakedata${n}_MC_H1matrix.log
        ./sbnfit_lee_frequentist_study --xml ../../xml/PELEE/Fakedata/${sample}_reco_energy_fakedata${n}_MC.xml --signal ${sample}_reco_energy_fakedata${n}_MC.SBNspec.root --background ${sample}_reco_energy_fakedata${n}_MC_BKG.SBNspec.root --fakedata ${sample}_reco_energy_fakedata${n}_DATA.SBNspec.root --covariance ${sample}_reco_energy_fakedata${n}_MC.SBNcovar.root -e 1e-8 -f 20 --tag ${sample}_reco_energy_fakedata${n}_MC_H1matrix >& /uboone/data/users/${USER}/logs/${sample}_reco_energy_fakedata${n}_MC_H1matrix.log
    fi	
    if [[ "$feldman" == "true" ]]
    then
        [ -e /uboone/data/users/${USER}/logs/${sample}_reco_energy_fakedata${n}_MC_FC_StatsOnly.log ] && rm /uboone/data/users/${USER}/logs/${sample}_reco_energy_fakedata${n}_MC_FC_StatsOnly.log
        ./sbnfit_uboone_scaling_fc -x ../../xml/PELEE/Fakedata/${sample}_reco_energy_fakedata${n}_MC.xml -t ${sample}_reco_energy_fakedata${n}_MC -i lee -g "0 4 25" -m feldman --cnp --stat >& /uboone/data/users/${USER}/logs/${sample}_reco_energy_fakedata${n}_MC_FC_StatsOnly.log 
        [ -e /uboone/data/users/${USER}/logs/${sample}_reco_energy_fakedata${n}_MC_FC.log ] && rm /uboone/data/users/${USER}/logs/${sample}_reco_energy_fakedata${n}_MC_FC.log
        ./sbnfit_uboone_scaling_fc -x ../../xml/PELEE/Fakedata/${sample}_reco_energy_fakedata${n}_MC.xml -t ${sample}_reco_energy_fakedata${n}_MC -i lee -g "0 4 25" -m feldman --cnp --detsys >& /uboone/data/users/${USER}/logs/${sample}_reco_energy_fakedata${n}_MC_FC.log 
    fi
done
