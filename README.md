# MIRA - The Genome and Transcriptome Assembler and Mapper

For de-novo assemblies, MIRA contains integrated editors for all supported sequence technologies which iteratively remove many sequencing errors from the assembly project and improve the overall alignment quality.

MIRA can also be used for mapping assemblies and automatic tagging of difference site (SNPs, insertions or deletions) of mutant strains against a reference sequence.
For organisms where annotated files in GFF3 format are available (or for GenBank files without intron/exon structures), MIRA can generate tables which are ready to use for biologists as they show exactly which genes are hit and give a first estimate whether the function of the protein is perturbed by the change.

# Scope

A quick note at the start. Nowadays, the scope I recommend MIRA
for is:

- For genome de-novo assemblies or mapping projects, haploid organisms up
  to 20 to 40 megabases should be the limit.
- Do not use MIRA if you have PacBio or Oxford Nanopore reads. Then again,
  for polishing those assemblies with Illumina data, MIRA is really good.
- For mapping projects, do not use MIRA if you expect splicing like, e.g.,
  RNASeq against an eukaryotic genome. Genome to genome or RNASeq to
  transcript is fine.
- Lastly, Illumina projects with more than 40 to 60 million reads start to
  be so resource intensive that you might be better served with other
  assemblers or mapping programs. I know some people use MIRA for de-novo
  RNASeq with 300 million reads because they think it's worth it, but they
  wait more than a month then for the calculation to finish.

# Building

Before trying to build MIRA yourself, please consider using the [pre-built
binaries](https://github.com/bachev/mira/releases) which are often available
for Linux and Mac OSX. If your really must build yourself, consult the
[chapter handling
installation](http://mira-assembler.sourceforge.net/docs/DefinitiveGuideToMIRA.html#chap_installation)
in the "Definitive Guide to MIRA" for more information on pre-requisites,
available options and walkthrough for common systems. There's also a section
on building documentation in the same file.

# Need help?

Please consult the extensive online documentation (as
[HTML](http://mira-assembler.sourceforge.net/docs/DefinitiveGuideToMIRA.html)
or
[PDF](http://mira-assembler.sourceforge.net/docs/DefinitiveGuideToMIRA.pdf))
which covers more or less all aspects of MIRA. If questions persist, please
subscribe to [the MIRA talk list](https://www.freelists.org/list/mira_talk)
and mail general questions to the list via this address:
  mira_talk@freelists.org

To report bugs or ask for new features, please use
 [the GitHub issue system](https://github.com/bachev/mira/issues)

This ensures that requests do not get lost.
