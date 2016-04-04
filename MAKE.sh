rm R_acr_all_results.txt
if [ $# -eq 1 ]
  then
    if [[ $1 == "prepare" ]]
      then
        python C_patterns.py pcd multistate
        python C_patterns.py edit multistate
        python C_patterns.py simple multistate
        python C_patterns.py simple binary
    fi
fi

echo 'Analyzing ARBRE'
echo 'simple'
python C_analyze.py matrix=simple.binary acr tree=P_arbre.tre
echo 'simple multistate'
python C_analyze.py matrix=simple.multistate acr tree=P_arbre.tre
echo 'weighted multistate'
python C_analyze.py matrix=edit.multistate acr tree=P_arbre.tre
echo 'weighted directed multistate'
python C_analyze.py matrix=pcd.multistate acr tree=P_arbre.tre

echo 'Analyzing SHUXINGTU'
echo 'simple'
python C_analyze.py matrix=simple.binary acr tree=P_shuxingtu.tre
echo 'simple multistate'
python C_analyze.py matrix=simple.multistate acr tree=P_shuxingtu.tre
echo 'weighted multistate'
python C_analyze.py matrix=edit.multistate acr tree=P_shuxingtu.tre
echo 'weighted directed multistate'
python C_analyze.py matrix=pcd.multistate acr tree=P_shuxingtu.tre

echo 'Analyzing Southern Chinese'
echo 'simple'
python C_analyze.py matrix=simple.binary acr tree=P_southern_chinese.tre
echo 'simple multistate'
python C_analyze.py matrix=simple.multistate acr tree=P_southern_chinese.tre
echo 'weighted multistate'
python C_analyze.py matrix=edit.multistate acr tree=P_southern_chinese.tre
echo 'weighted directed multistate'
python C_analyze.py matrix=pcd.multistate acr tree=P_southern_chinese.tre

python C_analyze.py summarize
cp R_table.tex ../draft/includes/

if [ $# -eq 1 ]
  then
    if [[ $1 == "plot" ]]
      then python C_analyze.py matrix=pcd.multistate plot tree=P_arbre.tre
  fi
  else
    echo "Finished analysis."
fi
