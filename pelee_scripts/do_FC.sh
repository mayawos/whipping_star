#!/bin/bash


GRID="0 4 100"

UNIV=5000

SAMPLE=("np" "zp" "np_zp")

DOLLR=$1 #( ' ' '--llr') #Option to use LLR

OUTPUTDIR='Unblinded_data'
mkdir $OUTPUTDIR

#make
i=0
for s in ${SAMPLE[@]}
	do
	
	FILES=constrained_${s}_numu_reco_e_H1_mc_collabOctFullMCunisim_with_mc_err_with_zerobin_err
	TAG=${FILES}
	#constrained_np_numu_reco_e_H1_DATA_mc_collabOctFullMCunisim_CV.SBNspec.root	
	#DATAFILE=${FILES}_DATA_CV.SBNspec.root
	DATAFILE=constrained_${s}_numu_reco_e_H1_DATA_mc_collabOctFullMCunisim_CV.SBNspec.root
	./sbnfit_uboone_scaling_fc -x const_${s}_numu_reco_e_H1_mc_collabOct.xml -t $TAG -i lee -g "${GRID}" -m feldman --detsys $DOLLR --cnp -n $UNIV
	./sbnfit_uboone_scaling_fc -x const_${s}_numu_reco_e_H1_mc_collabOct.xml -t $TAG -i lee -g "${GRID}" -m belt --detsys $DOLLR --cnp -n $UNIV  > Belts_${s}.txt
	 
	./sbnfit_uboone_scaling_fc -x const_${s}_numu_reco_e_H1_mc_collabOct.xml -t $TAG -i lee -g "${GRID}" -m beltdata --detsys $DOLLR --cnp -n $UNIV -d ${DATAFILE}  > Beltdata_fd${f}_${s}.txt
	
	./sbnfit_uboone_scaling_fc -x const_${s}_numu_reco_e_H1_mc_collabOct.xml -t $TAG -i lee -g "${GRID}" -m data --detsys $DOLLR --cnp -n $UNIV -d ${DATAFILE} > BF_fd${f}_${s}.txt
	
	cp SBNfeld_output_${TAG}.root $OUTPUTDIR
	mv DATACOMP_FC_confidence_belt_${TAG}.pdf $OUTPUTDIR
	mv FC_confidence_belt_${TAG}.pdf $OUTPUTDIR
	mv DeltaChi2Plot_DATA_Comp_${TAG}.pdf $OUTPUTDIR
	mv DeltaChi2Plot_${TAG}.pdf $OUTPUTDIR
	mv Belts_${s}.txt $OUTPUTDIR
	mv Beltdata_fd${f}_${s}.txt $OUTPUTDIR
	mv BF_fd${f}_${s}.txt $OUTPUTDIR
	
done
