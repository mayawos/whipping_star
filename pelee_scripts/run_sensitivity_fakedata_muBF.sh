for n in $1; do
for t in np_numu_reco_e_H1 zp_numu_reco_e_H1 np_zp_numu_reco_e_H1; do

 if [[ "$t" == "np_zp_numu_reco_e_H1" ]]; then 
       scale=0.95; #currently hardcoded to FD5 mu_BF
 elif [[ "$t" == "np_numu_reco_e_H1" ]]; then 
       scale=0.55;#0.5
 else 
       scale=3.4;#1.1
 fi
 #scale=1.0 #uncomment this to perform frequentist test wrt to H1 (1x LEE signal)

 echo "tag, scale = ${t}, ${scale}"
  cp ../examples/constrained_${t}_mc_fakedata${n}_CV.SBNspec.root .
  cp ../examples/constrained_${t}_mc_fakedata${n}_detsys.SBNcovar.root .
  ./sbnfit_scale_spec -x const_${t}_mc_fakedata.xml -i constrained_${t}_mc_fakedata${n}_CV.SBNspec.root -s lee -v 0.0  -t constrained_${t}_mc_fakedata${n}_BKG
  ./sbnfit_scale_spec -x const_${t}_mc_fakedata.xml -i constrained_${t}_mc_fakedata${n}_CV.SBNspec.root -s lee -v ${scale}  -t constrained_${t}_mc_fakedata${n}_SIG
  ./constrained_sbnfit_lee_frequentist_study --xml const_${t}_mc_fakedata.xml --signal constrained_${t}_mc_fakedata${n}_SIG.SBNspec.root --background constrained_${t}_mc_fakedata${n}_BKG.SBNspec.root --poisson -t constrained_${t}_mc_fakedata${n}_statsonly_BF --fakedata constrained_${t}_mc_fakedata${n}_DATA_CV.SBNspec.root
  ./constrained_sbnfit_lee_frequentist_study --xml const_${t}_mc_fakedata.xml --signal constrained_${t}_mc_fakedata${n}_SIG.SBNspec.root --background constrained_${t}_mc_fakedata${n}_BKG.SBNspec.root -c constrained_${t}_mc_fakedata${n}_no1st_bins.SBNcovar.root  -t constrained_${t}_mc_fakedata${n}_syst_BF  --fakedata constrained_${t}_mc_fakedata${n}_DATA_CV.SBNspec.root
  ./constrained_sbnfit_lee_frequentist_study --xml const_${t}_mc_fakedata.xml --signal constrained_${t}_mc_fakedata${n}_SIG.SBNspec.root --background constrained_${t}_mc_fakedata${n}_BKG.SBNspec.root -c constrained_${t}_mc_fakedata${n}_detsys.SBNcovar.root  -t constrained_${t}_mc_fakedata${n}_syst_detsys --fakedata ../examples/constrained_${t}_mc_fakedata${n}_DATA_CV.SBNspec.root
done
done
