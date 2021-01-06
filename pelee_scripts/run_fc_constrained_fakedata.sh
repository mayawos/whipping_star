for n in 1 2 3 4 5; 
do
    for s in np_zp_numu_reco_e np_numu_reco_e zp_numu_reco_e; 
    do
	./sbnfit_make_spec -x ${s}_H1_mc_fakedata${n}_DATA -t ${s}_H1_mc_fakedata${n}_DATA
	./sbnfit_uboone_scaling_fc -x const_${s}_H1_mc_fakedata${n}.xml -t constrained_${s}_H1_mc_fakedata${n} -i lee -g "0 4 80" -m feldman --detsys --cnp
	./sbnfit_uboone_scaling_fc -x const_${s}_H1_mc_fakedata${n}.xml -t constrained_${s}_H1_mc_fakedata${n} -i lee -g "0 4 80" -m belt --detsys --cnp
	./sbnfit_uboone_scaling_fc -x const_${s}_H1_mc_fakedata${n}.xml -t constrained_${s}_H1_mc_fakedata${n} -i lee -g "0 4 80" -m data -d ${s}_H1_mc_fakedata${n}_DATA.SBNspec.root --detsys --cnp
    done
done
