for t in np_numu_reco_e zp_numu_reco_e np_zp_numu_reco_e
do
    for s in mc_collabOctFullMCunisim
    do
	make
	./sbnfit_make_covariance -x ${t}_H1_${s}.xml -t ${t}_H1_${s} #-p
	./sbnfit_make_covariance -x numu_reco_e_H1_DATA_mc_collabOct.xml -t numu_reco_e_H1_DATA_mc_collabOct #-p
	
	cd ../examples/
	
	make
	if [[ "`echo $t`" == "np_numu_reco_e" ]]; then
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_mc_collabOct.SBNspec.root -n 
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_mc_collabOct.SBNspec.root -n -d
	elif [[ "`echo $t`" == "zp_numu_reco_e" ]]; then
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_mc_collabOct.SBNspec.root -z 
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_mc_collabOct.SBNspec.root -z -d
	elif [[ "`echo $t`" == "np_zp_numu_reco_e" ]]; then
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_mc_collabOct.SBNspec.root -c 
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_mc_collabOct.SBNspec.root -c -d
	fi
	
	cd ../bin
	cp ../examples/collapsed_zero_ext_bin_${t}_H1_mc_${s}.SBNcovar.root .
	cp ../examples/unconstrained_${t}_H1_${s}.SBNspec.root unconstrained_${t}_H1_${s}_CV.SBNspec.root
	cp ../examples/constrained_${t}_H1_${s}.SBNspec.root constrained_${t}_H1_${s}_CV.SBNspec.root
	cp ../examples/constrained_${t}_H1_${s}.SBNcovar.root .
	cp ../examples/constrained_${t}_H1_${s}_detsys.SBNcovar.root .
	./sbnfit_scale_spec -x const_${t}_H1_mc_collabOct.xml --input constrained_${t}_H1_${s}_CV.SBNspec.root --scalestring "lee" -v 0.0 -t constrained_${t}_H1_${s}_BKG
	./sbnfit_lee_frequentist_study --xml const_${t}_H1_mc_collabOct.xml --signal constrained_${t}_H1_${s}_CV.SBNspec.root --background constrained_${t}_H1_${s}_BKG.SBNspec.root --poisson --tag constrained_${t}_H1_${s}_stats
	./sbnfit_lee_frequentist_study --xml const_${t}_H1_mc_collabOct.xml --signal constrained_${t}_H1_${s}_CV.SBNspec.root --background constrained_${t}_H1_${s}_BKG.SBNspec.root --covariance constrained_${t}_H1_${s}.SBNcovar.root -e 1e-8 --tag constrained_${t}_H1_${s}_syst
	./sbnfit_lee_frequentist_study --xml const_${t}_H1_mc_collabOct.xml --signal constrained_${t}_H1_${s}_CV.SBNspec.root --background constrained_${t}_H1_${s}_BKG.SBNspec.root --covariance constrained_${t}_H1_${s}_detsys.SBNcovar.root -e 1e-8 --tag constrained_${t}_H1_${s}_syst_detsys
    done
done
