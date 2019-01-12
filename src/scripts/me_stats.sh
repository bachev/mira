
echo -n "Reads used: " 
grep -c -v ^# step1_info_contigreadlist.txt
echo

echo "Results: disparate SNP sites"
echo "----------------------------"
echo -n "Contigs built: " 
cut -f 1 step1_info_contigreadlist.txt | grep '_c[0-9]\+$' | uniq | grep -c ^
echo -n "Singlets: " 
cut -f 1 step1_info_contigreadlist.txt | grep '_s[0-9]\+$' | uniq | grep -c ^
echo

echo "Results: unified SNP sites from above"
echo "-------------------------------------"
echo -n "Unified contigs with SNP sites: "
grep PAOS step3_info_consensustaglist.txt >/tmp/gna
grep PROS step3_info_consensustaglist.txt >>/tmp/gna
grep PIOS step3_info_consensustaglist.txt >>/tmp/gna
cut -f 1 /tmp/gna | uniq >/tmp/gna2
grep -c ^ /tmp/gna2

echo -n "Num. SNP positions: " 
grep -c ^ /tmp/gna

echo
echo


########## ...having # reads: grep -f gna2 step3_info_contigreadlist.txt  grep -c ^


