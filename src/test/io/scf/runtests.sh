#!/bin/sh


rm tstFASTA.bad ttt diff.txt 2>/dev/null
./tstFASTA 1> ttt 2>&1
diff tstFASTA.out ttt >diff.txt
if test -s diff.txt 
then
    mv ttt tstFASTA.bad
    echo `pwd`
    echo "fasta tests differ! Original output is in tstFASTA.out, actual in tstFASTA.bad."
    echo "Differences are in file diff.txt"; 
else
    echo "fasta test ok."
    rm diff.txt ttt 2>/dev/null
fi