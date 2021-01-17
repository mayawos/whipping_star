rm $1_$3_beforeafterconst_table.tex;
#nue all bins table
echo "\\begin{table}[H]" >> $1_$3_beforeafterconst_table.tex;
echo "\\centering" >> $1_$3_beforeafterconst_table.tex;
echo "\\begin{tabular}{| c | m{1.65cm} | m{1.65cm} | m{1.65cm} | m{1.65cm} | m{1.65cm} | m{1.65cm}|} " >> $1_$3_beforeafterconst_table.tex;
echo "\\cline{1-7}" >> $1_$3_beforeafterconst_table.tex;
echo "\\multirow{2}{*}{Energy [GeV]} &\multicolumn{3}{c|}{Flux+Genie+G4}&\multicolumn{3}{c|}{Flux+Genie+G4+Det.Systematics+Sample Stats}\\\\" >> $1_$3_beforeafterconst_table.tex;
echo "\\cline{2-7}" >> $1_$3_beforeafterconst_table.tex;
echo "{} &  before constraint & after constraint  & size of reduction & before constraint & after constraint & size of reduction \\\\" >> $1_$3_beforeafterconst_table.tex;
echo "\\hline" >> $1_$3_beforeafterconst_table.tex;
rows=`wc -l full_$1_systematicstable_beforeconstraint_$2_H1_$3$4$5.txt | awk '{print $1}'`
for (( l=1; l<=${rows}; l++ )) 
do 
  col1=`head -${l} full_$1_systematicstable_beforeconstraint_$2_H1_$3$4$5.txt | tail -1 | awk '{print $2}'`; 
  col2=`head -${l} full_$1_systematicstable_afterconstraint_$2_H1_$3$4$5.txt | tail -1 | awk '{print $2}'`; 
  col3=`head -${l} full_$1_systematicstable_beforeconstraint_$2_H1_$3$4$5_with_detsys.txt | tail -1 | awk '{print $2}'`; 
  col4=`head -${l} full_$1_systematicstable_afterconstraint_$2_H1_$3$4$5_with_detsys.txt | tail -1 | awk '{print $2}'`; 
  a=`ksh -c 'echo "$(('${col1}' * 100.))"'`; 
  b=`ksh -c 'echo "$(('${col2}' * 100.))"'`;
  c=`ksh -c 'echo "$(( ( 1. - ( '${col2}' / '${col1}') ) * 100. ))"'`;
  d=`ksh -c 'echo "$(('${col3}' * 100.))"'`;
  e=`ksh -c 'echo "$(('${col4}' * 100.))"'`;
  f=`ksh -c 'echo "$(( ( 1. - ( '${col4}' / '${col3}') ) * 100. ))"'`;
  binstart=`ksh -c 'echo "$(( 0.15 + (('${l}' - 1.) * 0.1) ))"'`
  binend=`ksh -c 'echo "$(( 0.25 + (('${l}'- 1. ) * 0.1) ))"'`
  echo "`printf "%.2f - %.2f" ${binstart} ${binend}` & `printf "%.2f" ${a}` & `printf "%.2f" ${b}` & `printf "%.2f" ${c}` & `printf "%.2f" ${d}` & `printf "%.2f" ${e}` & `printf "%.2f" ${f}` \\\\" >> $1_$3_beforeafterconst_table.tex; 
done
echo "\\hline" >> $1_$3_beforeafterconst_table.tex;
echo "\\end{tabular}" >> $1_$3_beforeafterconst_table.tex;
echo "\\label{tab:$1_numu_bdt}" >> $1_$3_beforeafterconst_table.tex;
echo "\\end{table}" >> $1_$3_beforeafterconst_table.tex;
echo "" >> $1_$3_beforeafterconst_table.tex;
echo "\\newpage" >> $1_$3_beforeafterconst_table.tex;
echo "" >> $1_$3_beforeafterconst_table.tex;
