#!/bin/bash

usage="$(basename "$0") [-h] [-n] -- script to build the xml files based on the original fakedata xml

where:
        -h  show this help text
        -n  which fake data set? (if more than one set, need to be enclosed in single quotes)"

while getopts ':h:n:' option; do
        case "$option" in
    h) echo "$usage"
            exit
            ;;
    n) number=$OPTARG
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

dir_in="/uboone/app/users/${USER}/SBNFitMarkBranch/whipping_star/build/bin"
dir_out="/uboone/app/users/${USER}/SBNFitToCommit/whipping_star/xml/PELEE/Fakedata"

for n in $number; do
## declare an array variable
    declare -a file_in=("nue_1e0p_numu_reco_e_H1_mc_fakedata${n}.xml" "nue_1e0p_numu_reco_e_H1_fakedata${n}.xml" "nue_1e0p_numu_reco_e_H1_mc_fakedata${n}_nolee.xml")
    declare -a file_out=("Combined_1eNp_1e0p_numu_reco_energy_fakedata${n}_MC.xml" "Combined_1eNp_1e0p_numu_reco_energy_fakedata${n}_DATA.xml" "Combined_1eNp_1e0p_numu_reco_energy_fakedata${n}_nolee_MC.xml")

    for ((i=0;i<${#file_in[@]};++i)); do
        printf "%s is an input to %s\n" "$dir_in/${file_in[i]}" "$dir_out/${file_out[i]}"
    
        [ -e $dir_out/${file_out[i]} ] && rm $dir_out/${file_out[i]}
        cat fakedata_header.xml >> $dir_out/${file_out[i]}
        if [[  ${file_in[i]} =~ "mc_fakedata" ]]
        then
            head -226 ${dir_in}/${file_in[i]} | tail -182 >> $dir_out/${file_out[i]}
            cat fakedata_tail.xml >> $dir_out/${file_out[i]}
        else
             head -60 ${dir_in}/${file_in[i]} | tail -33 >> $dir_out/${file_out[i]}
        fi
        sed -i 's/uBooNE_nue/uBooNE_1eNp/' $dir_out/${file_out[i]}
    done
done
