#!/bin/sh


for name in *
do
    if test -d $name
    then
	if test -x $name/runtests.sh
	then
	    cd $name; ./runtests.sh; cd ..;
        fi
    fi
done