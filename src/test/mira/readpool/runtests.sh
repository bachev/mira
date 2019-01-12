#!/bin/sh


rm tstReadPool.bad ttt diff.txt 2>/dev/null
./tstReadPool 1> ttt 2>&1
diff tstReadPool.out ttt >diff.txt
if test -s diff.txt 
then
    mv ttt tstReadPool.bad
    echo `pwd`
    echo "ReadPool tests differ! Original output is in tstReadPool.out, actual in tstReadPool.bad."
    echo "Differences are in file diff.txt"; 
else
    echo "ReadPool test ok."
    rm diff.txt ttt 2>/dev/null
fi