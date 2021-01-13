for t in np_numu_reco_e #zp_numu_reco_e np_zp_numu_reco_e
do
    for s in mc_collabOctFullMCunisim #this tag identify the xml with the correct treatment of the unisim error
    do
	#make
	#./sbnfit_make_covariance -x ${t}_H1_${s}.xml -t ${t}_H1_${s} #-p
	#./sbnfit_make_covariance -x numu_reco_e_H1_DATA_mc_collabOct.xml -t numu_reco_e_H1_DATA_mc_collabOct #-p
	
	cd ../examples/
	make

        #Constraint procedures
        #Flags explanation:
	#  -n, -z, -c : channels: 1eNp+numu, 1e0p+numu, 1eNp+1e0p+numu
	#  -d : add PeLEE detsys error to the diagonals of full covariance matrix, then collapse, and constrain
	#  -m : add mc stats error to the diagonal of the  full covariance matrix, then collapse, and constrain
	#  -a : add zero ext bin error to the diagonal of the full covariance matrix, then collapse, and constrain
	if [[ "`echo $t`" == "np_numu_reco_e" ]]; then
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_mc_collabOct.SBNspec.root -n  
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_mc_collabOct.SBNspec.root -n -a 
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_mc_collabOct.SBNspec.root -n -m
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_mc_collabOct.SBNspec.root -n -m -a
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_mc_collabOct.SBNspec.root -n -d 
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_mc_collabOct.SBNspec.root -n -d -m
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_mc_collabOct.SBNspec.root -n -d -m -a
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_mc_collabOct.SBNspec.root -n -d -a
	elif [[ "`echo $t`" == "zp_numu_reco_e" ]]; then
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_mc_collabOct.SBNspec.root -z 
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_mc_collabOct.SBNspec.root -z -m 
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_mc_collabOct.SBNspec.root -z -m -a
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_mc_collabOct.SBNspec.root -z -d 
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_mc_collabOct.SBNspec.root -z -d -m 
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_mc_collabOct.SBNspec.root -z -d -m -a
	elif [[ "`echo $t`" == "np_zp_numu_reco_e" ]]; then
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_mc_collabOct.SBNspec.root -c 
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_mc_collabOct.SBNspec.root -c -m 
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_mc_collabOct.SBNspec.root -c -m -a
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_mc_collabOct.SBNspec.root -c -d
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_mc_collabOct.SBNspec.root -c -d -m
	    ./CreateConstrainedMatrix -x ../bin/${t}_H1_${s}.xml -t ${t}_H1_${s} -f ../bin/numu_reco_e_H1_DATA_mc_collabOct.SBNspec.root -c -d -m -a
       fi

        #tag used to identify different configuration for sensitivity plots and logfiles of sensitivity/pvals results
        m="_with_mc_err"
        a="_with_zerobin_err"
	
        #move to bin directory to write the sensitivity results
        cd ../bin

        #copy the constrained/unconstrained Covariance Matrices and Histograms/Spectras to ../bin for sensitivity test
	cp ../examples/collapsed_zero_ext_bin_${t}_H1_${s}${a}.SBNcovar.root .
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
	./sbnfit_scale_spec -x const_${t}_H1_mc_collabOct.xml --input constrained_${t}_H1_${s}_CV.SBNspec.root --scalestring "lee" -v 0.0 -t constrained_${t}_H1_${s}_BKG
	./sbnfit_scale_spec -x const_${t}_H1_mc_collabOct.xml --input constrained_${t}_H1_${s}${m}_CV.SBNspec.root --scalestring "lee" -v 0.0 -t constrained_${t}_H1_${s}${m}_BKG
	./sbnfit_scale_spec -x const_${t}_H1_mc_collabOct.xml --input constrained_${t}_H1_${s}${m}${a}_CV.SBNspec.root --scalestring "lee" -v 0.0 -t constrained_${t}_H1_${s}${m}${a}_BKG

        #sensitivity test for stats only. Order: no additional error, with mc stats error only, with mc stats error + zero bin error, with zero bin error only
	./sbnfit_lee_frequentist_study --xml const_${t}_H1_mc_collabOct.xml --signal constrained_${t}_H1_${s}_CV.SBNspec.root --background constrained_${t}_H1_${s}_BKG.SBNspec.root --poisson --tag constrained_${t}_H1_${s}_stats 
	./sbnfit_lee_frequentist_study --xml const_${t}_H1_mc_collabOct.xml --signal constrained_${t}_H1_${s}${m}_CV.SBNspec.root --background constrained_${t}_H1_${s}${m}_BKG.SBNspec.root --poisson --tag constrained_${t}_H1_${s}${m}_stats 
	./sbnfit_lee_frequentist_study --xml const_${t}_H1_mc_collabOct.xml --signal constrained_${t}_H1_${s}${a}_CV.SBNspec.root --background constrained_${t}_H1_${s}${a}_BKG.SBNspec.root --poisson --tag constrained_${t}_H1_${s}${a}_stats 
	./sbnfit_lee_frequentist_study --xml const_${t}_H1_mc_collabOct.xml --signal constrained_${t}_H1_${s}${m}_CV.SBNspec.root --background constrained_${t}_H1_${s}_BKG.SBNspec.root --poisson --tag constrained_${t}_H1_${s}${m}_stats  

        #sensitivity test for flux+xsec only. Order: no additional error, with mc stats error only, with mc stats error + zero bin error, with zero bin error only
	./sbnfit_lee_frequentist_study --xml const_${t}_H1_mc_collabOct.xml --signal constrained_${t}_H1_${s}_CV.SBNspec.root --background constrained_${t}_H1_${s}_BKG.SBNspec.root --covariance constrained_${t}_H1_${s}.SBNcovar.root -e 1e-8 --tag constrained_${t}_H1_${s}_syst
	./sbnfit_lee_frequentist_study --xml const_${t}_H1_mc_collabOct.xml --signal constrained_${t}_H1_${s}${m}.{m}_CV.SBNspec.root --background constrained_${t}_H1_${s}${m}.{m}_BKG.SBNspec.root --covariance constrained_${t}_H1_${s}${m}.SBNcovar.root -e 1e-8 --tag constrained_${t}_H1_${s}${m}.{m}_syst
	./sbnfit_lee_frequentist_study --xml const_${t}_H1_mc_collabOct.xml --signal constrained_${t}_H1_${s}${m}${a}.{m}${a}_CV.SBNspec.root --background constrained_${t}_H1_${s}${m}${a}.{m}${a}_BKG.SBNspec.root --covariance constrained_${t}_H1_${s}${m}${a}.SBNcovar.root -e 1e-8 --tag constrained_${t}_H1_${s}${m}${a}.{m}${a}_syst
	./sbnfit_lee_frequentist_study --xml const_${t}_H1_mc_collabOct.xml --signal constrained_${t}_H1_${s}_CV.SBNspec.root --background constrained_${t}_H1_${s}_BKG.SBNspec.root --covariance constrained_${t}_H1_${s}.SBNcovar.root -e 1e-7 --tag constrained_${t}_H1_${s}${m}_syst

        #sensitivity test for flux+xsec+detsys only. Order: no additional error, with mc stats error only, with mc stats error + zero bin error, with zero bin error only
	./sbnfit_lee_frequentist_study --xml const_${t}_H1_mc_collabOct.xml --signal constrained_${t}_H1_${s}_CV.SBNspec.root --background constrained_${t}_H1_${s}_BKG.SBNspec.root --covariance constrained_${t}_H1_${s}_detsys.SBNcovar.root -e 1e-8 --tag constrained_${t}_H1_${s}_syst_detsys
	./sbnfit_lee_frequentist_study --xml const_${t}_H1_mc_collabOct.xml --signal constrained_${t}_H1_${s}${m}_CV.SBNspec.root --background constrained_${t}_H1_${s}${m}_BKG.SBNspec.root --covariance constrained_${t}_H1_${s}${m}_detsys.SBNcovar.root -e 1e-8 --tag constrained_${t}_H1_${s}${m}_syst_detsys
	./sbnfit_lee_frequentist_study --xml const_${t}_H1_mc_collabOct.xml --signal constrained_${t}_H1_${s}${m}${a}_CV.SBNspec.root --background constrained_${t}_H1_${s}${m}${a}_BKG.SBNspec.root --covariance constrained_${t}_H1_${s}${m}${a}_detsys.SBNcovar.root -e 1e-8 --tag constrained_${t}_H1_${s}${m}${a}_syst_detsys
	./sbnfit_lee_frequentist_study --xml const_${t}_H1_mc_collabOct.xml --signal constrained_${t}_H1_${s}_CV.SBNspec.root --background constrained_${t}_H1_${s}_BKG.SBNspec.root --covariance constrained_${t}_H1_${s}_detsys.SBNcovar.root -e 1e-7 --tag constrained_${t}_H1_${s}${m}_syst_detsys 

    done
done
