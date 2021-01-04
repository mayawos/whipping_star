for t in np_numu_reco_e zp_numu_reco_e np_zp_numu_reco_e
do
    for n in 1 2 3 4 5
    do
	make
	./sbnfit_make_covariance -x ${t}_H1_mc_fakedata${n}.xml -t ${t}_H1_mc_fakedata${n} #-p
	./sbnfit_make_covariance -x ${t}_H1_fakedata${n}.xml -t ${t}_H1_fakedata${n} #-p
	
	cd ../examples/
	
	make
	if [[ "`echo $t`" == "np_numu_reco_e" ]]; then
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_mc_fakedata${n}.xml -t ${t}_H1_mc_fakedata${n} -f ../bin/${t}_H1_fakedata${n}.SBNspec.root -n 
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_mc_fakedata${n}.xml -t ${t}_H1_mc_fakedata${n} -f ../bin/${t}_H1_fakedata${n}.SBNspec.root -n -d
	elif [[ "`echo $t`" == "zp_numu_reco_e" ]]; then
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_mc_fakedata${n}.xml -t ${t}_H1_mc_fakedata${n} -f ../bin/${t}_H1_fakedata${n}.SBNspec.root -z 
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_mc_fakedata${n}.xml -t ${t}_H1_mc_fakedata${n} -f ../bin/${t}_H1_fakedata${n}.SBNspec.root -z -d
	elif [[ "`echo $t`" == "np_zp_numu_reco_e" ]]; then
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_mc_fakedata${n}.xml -t ${t}_H1_mc_fakedata${n} -f ../bin/${t}_H1_fakedata${n}.SBNspec.root -c 
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_mc_fakedata${n}.xml -t ${t}_H1_mc_fakedata${n} -f ../bin/${t}_H1_fakedata${n}.SBNspec.root -c -d
	fi
	
	cd ../bin

	cp ../examples/unconstrained_${t}_H1_mc_fakedata${n}_CV.SBNspec.root unconstrained_${t}_H1_mc_fakedata${n}_CV.SBNspec.root
	cp ../examples/constrained_${t}_H1_mc_fakedata${n}_CV.SBNspec.root constrained_${t}_H1_mc_fakedata${n}_CV.SBNspec.root
	cp ../examples/constrained_${t}_H1_mc_fakedata${n}.SBNcovar.root .
	cp ../examples/constrained_${t}_H1_mc_fakedata${n}_detsys.SBNcovar.root .

	./sbnfit_scale_spec -x const_${t}_H1_mc_fakedata${n}.xml --input constrained_${t}_H1_mc_fakedata${n}_CV.SBNspec.root --scalestring "lee" -v 0.0 -t constrained_${t}_H1_mc_fakedata${n}_BKG
	./sbnfit_lee_frequentist_study --xml const_${t}_H1_mc_fakedata${n}.xml --signal constrained_${t}_H1_mc_fakedata${n}_CV.SBNspec.root --background constrained_${t}_H1_mc_fakedata${n}_BKG.SBNspec.root --poisson --tag constrained_${t}_H1_mc_fakedata${n}_stats
	./sbnfit_lee_frequentist_study --xml const_${t}_H1_mc_fakedata${n}.xml --signal constrained_${t}_H1_mc_fakedata${n}_CV.SBNspec.root --background constrained_${t}_H1_mc_fakedata${n}_BKG.SBNspec.root --covariance constrained_${t}_H1_mc_fakedata${n}.SBNcovar.root -e 1e-8 --tag constrained_${t}_H1_mc_fakedata${n}_syst
	./sbnfit_lee_frequentist_study --xml const_${t}_H1_mc_fakedata${n}.xml --signal constrained_${t}_H1_mc_fakedata${n}_CV.SBNspec.root --background constrained_${t}_H1_mc_fakedata${n}_BKG.SBNspec.root --covariance constrained_${t}_H1_mc_fakedata${n}_detsys.SBNcovar.root -e 1e-8 --tag constrained_${t}_H1_mc_fakedata${n}_syst_detsys
    done
done
