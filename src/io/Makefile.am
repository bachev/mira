AM_CPPFLAGS = -I$(top_srcdir)/src $(all_includes)

clean-local:
	rm -rf *flexer.cc *xxd.H

so.obo.xxd: so.obo
	@echo "Reducing so.obo for smaller so.obo.xxd"
	@grep -E '^\[|^id:|^name:' so.obo >so.obo.xxd

SUFFIXES= .xxd

%.xxd.H: %.xxd
	@echo "Making .H out of .xxd" $@
	@$(SHELL) -ec 'grep -v "^#" $< | xxd -i >$<.H'


libmirafio_a_XXDSRC = mira_ft_set.xxd gap4_ft_so_map.xxd so.obo.xxd so2xgap4map.xxd oldso2so.xxd
BUILT_SOURCES=$(libmirafio_a_XXDSRC:.xxd=.xxd.H)
nodist_libmirafio_a_SOURCES=$(libmirafio_a_XXDSRC:.xxd=.xxd.H)


noinst_LIBRARIES = libmirafio.a
libmirafio_a_SOURCES= \
	annotationmappings.C\
	exp.C\
	generalio.C\
	fasta.C\
	fastq-mira.C\
	phd.C\
	scf.C\
	ncbiinfoxml.C\
	exp_flexer.ll\
	gap4_ft_so_map.xxd\
	mira_ft_set.xxd\
	oldso2so.xxd\
	so2xgap4map.xxd\
	so.obo

noinst_HEADERS= \
	annotationmappings.H\
	exp.H\
	fasta.H\
	fastq-mira.H\
	fastq-lh.H\
	generalio.H\
	ncbiinfoxml.H\
	phd.H\
	scf.H\
	tokens.h
