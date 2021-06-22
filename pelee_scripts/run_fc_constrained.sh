for s in np_zpnumu_reco_e np_numu_reco_e zp_numu_reco_e; 
do
  ./sbnfit_uboone_scaling_fc -x const_${s}_H1_mc_collabOct.xml -t constrained_${s}_H1_mc_collabOctFullMC -i lee -g "0 4 80" -m feldman --detsys --cnp
  ./sbnfit_uboone_scaling_fc -x const_${s}_H1_mc_collabOct.xml -t constrained_${s}_H1_mc_collabOctFullMC -i lee -g "0 4 80" -m belt --detsys --cnp
done
