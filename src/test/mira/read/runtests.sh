#!/bin/sh


rm tstRead.bad ttt diff.txt 2>/dev/null
./tstRead 1> ttt 2>&1
diff tstRead.out ttt >diff.txt
if test -s diff.txt 
then
    mv ttt tstRead.bad
    echo `pwd`
    echo "Read tests differ! Original output is in tstRead.out, actual in tstRead.bad."
    echo "Differences are in file diff.txt"; 
else
    echo "Read test ok."
    rm diff.txt ttt 2>/dev/null
fi