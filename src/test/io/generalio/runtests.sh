#!/bin/sh


rm tstgeneralio.bad ttt diff.txt 2>/dev/null
./tstgeneralio 1> ttt 2>&1
diff tstgeneralio.out ttt >diff.txt
if test -s diff.txt 
then
    mv ttt tstgeneralio.bad
    echo `pwd`
    echo "fasta tests differ! Original output is in tstgeneralio.out, actual in tstgeneralio.bad."
    echo "Differences are in file diff.txt"; 
else
    echo "generalio test ok."
    rm diff.txt ttt 2>/dev/null
fi