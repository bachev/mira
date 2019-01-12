SUBDIRS = template
EXTRA_DIST = README

clean:
	-mkdir saveMakefile
	mv Makefile* saveMakefile
	-rm -f *
	mv saveMakefile/* .
	rmdir saveMakefile

testconf: clean
	(cd template;  cp `ls -1 | grep -v Makefile| grep -v CVS` ..)
	ln -s ../src/progs/mira mira
