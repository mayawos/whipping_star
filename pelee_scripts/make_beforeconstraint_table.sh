rm $1_$3_beforeconstraint_table.tex;
#nue all bins table
echo "\\begin{table}[H]" >> $1_$3_beforeconstraint_table.tex;
echo "\\centering" >> $1_$3_beforeconstraint_table.tex;
echo "\\begin{tabular}{| c | m{1.65cm} | m{1.65cm} | m{1.65cm} | m{1.65cm} | m{1.65cm} | m{1.65cm}| } " >> $1_$3_beforeconstraint_table.tex;
echo "\\hline" >> $1_$3_beforeconstraint_table.tex;
echo "Energy [GeV] & Flux Only & Genie Only  & G4 Only & Flux+ Genie+ G4 & Flux+ Genie+ G4+ Det.Syst. & Flux+ Genie+ G4+ Det.Syst.+ Sample Stats\\\\" >> $1_$3_beforeconstraint_table.tex;
echo "\hline" >> $1_$3_beforeconstraint_table.tex;
rows=`wc -l full_$1_systematicstable_beforeconstraint_$2_H1_$3.txt | awk '{print $1}'`
echo $rows
for (( l=1; l<=${rows}; l++ )) 
do
  col1=`head -${l} full_$1_systematicstable_beforeconstraint_$2_H1_$3_fluxonly.txt | tail -1 | awk '{print $2}'`; 
  col2=`head -${l} full_$1_systematicstable_beforeconstraint_$2_H1_$3_genieonly.txt | tail -1 | awk '{print $2}'`; 
  col3=`head -${l} full_$1_systematicstable_beforeconstraint_$2_H1_$3_g4only.txt | tail -1 | awk '{print $2}'`; 
  col4=`head -${l} full_$1_systematicstable_beforeconstraint_$2_H1_$3.txt | tail -1 | awk '{print $2}'`; 
  col5=`head -${l} full_$1_systematicstable_beforeconstraint_$2_H1_$3_with_detsys.txt | tail -1 | awk '{print $2}'`; 
  col6=`head -${l} full_$1_systematicstable_beforeconstraint_$2_H1_$3_with_mc_err_with_zerobin_err_with_detsys.txt | tail -1 | awk '{print $2}'`; 
  binstart=`ksh -c 'echo "$(( 0.15 + (('${l}' - 1.) * 0.1) ))"'`
  binend=`ksh -c 'echo "$(( 0.25 + (('${l}'- 1. ) * 0.1) ))"'`
  a=`ksh -c 'echo "$(('${col1}' * 100.))"'`; 
  b=`ksh -c 'echo "$(('${col2}' * 100.))"'`;
  c=`ksh -c 'echo "$(('${col3}' * 100.))"'`; 
  d=`ksh -c 'echo "$(('${col4}' * 100.))"'`;
  e=`ksh -c 'echo "$(('${col5}' * 100.))"'`; 
  f=`ksh -c 'echo "$(('${col6}' * 100.))"'`;
  echo "`printf "%.2f - %.2f" ${binstart} ${binend}` & `printf "%.2f" ${a}` & `printf "%.2f" ${b}` & `printf "%.2f" ${c}` & `printf "%.2f" ${d}` & `printf "%.2f" ${e}` & `printf "%.2f" ${f}`\\\\" >> $1_$3_beforeconstraint_table.tex; 
done
echo "\\hline" >> $1_$3_beforeconstraint_table.tex;
echo "\\end{tabular}" >> $1_$3_beforeconstraint_table.tex;
echo "\\label{tab:$1_errors}" >> $1_$3_beforeconstraint_table.tex;
echo "\\end{table}" >> $1_$3_beforeconstraint_table.tex;
echo "" >> $1_$3_beforeconstraint_table.tex;
echo "\\newpage" >> $1_$3_beforeconstraint_table.tex;
echo "" >> $1_$3_beforeconstraint_table.tex;
