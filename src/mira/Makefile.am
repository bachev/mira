AM_CPPFLAGS = -I$(top_srcdir)/src $(all_includes)

clean-local:
	rm -rf *flexer.cc *xxd.H

EXTRA_DIST = CHANGES.txt CHANGES_old.txt

SUFFIXES= .xxd

%.xxd.H: %.xxd
	@echo "Making .H out of .xxd" $@
	@$(SHELL) -ec 'grep -v "^#" $< | xxd -i >$<.H'


noinst_LIBRARIES = libmiraestass.a libmira.a
#noinst_LIBRARIES = libmirasupport.a libmira.a

libmira_a_XXDSRC = adaptorsforclip.454.xxd adaptorsforclip.solexa.xxd adaptorsforclip.iontor.xxd adaptorsregex.solexa.xxd adaptorsregex.iontor.xxd seqforfilter_phix174.solexa.xxd
BUILT_SOURCES=$(libmira_a_XXDSRC:.xxd=.xxd.H)
nodist_libmira_a_SOURCES=$(libmira_a_XXDSRC:.xxd=.xxd.H)

#dist_pkgdata_DATA = demoparam.prm
#pkgdatadir=$(prefix)/share/$(PACKAGE)

# put hashstats, assembly, read, contig, etc. in front: with make -j .... they'll get made first
#  and won't block the rest.

libmira_a_SOURCES= \
	$(libmira_a_XXDSRC) \
	hashstats_sdbg.C\
	hashstats.C\
	skim.C\
	assembly.C\
	contig.C\
	read.C\
	contig_analysis.C\
	ads.C\
	\
	adsfacts.C\
	align.C\
	assembly_info.C\
	assembly_buildcontigs.C\
	assembly_io.C\
	assembly_misc.C\
	assembly_output.C\
	assembly_pbcorrect.C\
	assembly_reduceskimhits.C\
	assembly_swalign.C\
	contig_consensus.C\
	contig_covanalysis.C\
	contig_edit.C\
	contig_featureinfo.C\
	contig_output.C\
	contig_pairconsistency.C\
	dataprocessing.C\
	dynamic.C\
	gbf_parse.C\
	gff_parse.C\
	gff_save.C\
	maf_parse.C\
	manifest.C\
	multitag.C\
	parameters_flexer.ll\
	parameters.C \
	pcrcontainer.C \
	ppathfinder.C\
	preventinitfiasco.C\
	readgrouplib.C\
	readlotstat.C\
	readpool.C\
	readpool_io.C\
	sam_collect.C\
	scaffolder.C\
	seqtohash.C\
	simplebloomfilter.C\
	skim_lowbph.C\
	warnings.C
noinst_HEADERS= adsfacts.H\
	ads.H\
	assembly_output.H\
	parameters.H\
	align.H\
	assembly_info.H\
	assembly.H\
	contig.H\
	dataprocessing.H\
	dynamic.H\
	enums.H\
	gbf_parse.H\
	gff_parse.H\
	gff_save.H\
	hashstats.H\
	hdeque.H\
	manifest.H\
	maf_parse.H\
	multitag.H\
	overlapedges.H\
	parameters_tokens.h\
	pcrcontainer.H \
	ppathfinder.H\
	read.H\
	readgrouplib.H\
	readlotstat.H\
	readpool.H\
	readpool_io.H\
	readseqtypes.H\
	sam_collect.H\
	scaffolder.H\
	seqtohash.H\
	simple_2Dsignalprocessing.H\
	simplebloomfilter.H\
	skim.H\
	stringcontainer.H\
	structs.H\
	timerestrict.H\
	types_basic.H\
	vhash.H\
	vluint.H\
	warnings.H

libmiraestass_a_SOURCES= estassembly.C
