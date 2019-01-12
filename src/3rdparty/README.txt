The scripts in this directory have been written by users of MIRA to tackle
specific tasks. I am thankful to the authors who made the scripts available
for public use and allowed me to package them for download on the MIRA site.

For thanks, encouragements and, well, bug reports, please contact the
respective authors directly.

  Bastien


INSTALL:
========
To install the scripts on your system, simply copy the files to some directory
which is in your $PATH.


SCRIPTS:
========


sff_extract
-----------
Written by Jose Blanca (University of Valencia) and Bastien Chevreux
(paired-end and clipping quirk support).

This tool needs a Python installation (version >= 2.4). It extracts sequences
and qualities from the 454 SFF files and writes them as FASTA and FASTA
quality files. Additionally, it generates traceinfo XML files with ancillary
information.

sff_extract is a replacement for the sff tools from the Roche OffInstruments
package. In short: now you don't need to bug Roche/454 anymore and sign NDAs
to have them send you one vital piece of software ... it performs everything
you need to get your sequences running with MIRA (or any other piece of
software).

The home of sff_extract is: http://bioinf.comav.upv.es/sff_extract/index.html
but I am thankful to Jose for giving permission to distribute the script in
the MIRA 3rd party package.


454pairedEnd2caf.pl
-------------------
Written by Jacqueline Weber-Lehmann, MWG Biotech AG

Transform fasta 454 paired read info to CAF format which MIRA can use for
paired end assembly: it splits up every read into a left and right part and
orders them so that MIRA can use them.

You might want to change the insert size in the code to fit your library size.

The script needs cross_match and the sequences must be already clipped. For a
version which uses BLAST, see below.


readpair2caf.pl
---------------
Written by Stephen Taylor, Oxford University

Transform fasta 454 paired read info to CAF format which MIRA can use for
paired end assembly: it splits up every read into a left and right part and
orders them so that MIRA can use them.

You might want to change the insert size in the code to fit your library size.

The script needs cross_match and the sequences must be already clipped. For a
version which uses cross_match, see above.


lucy2xml.pl
-----------
Written by Ross Whettten, NC State University Raleigh

A script to convert output from the Lucy quality & vector trimming program
into an XML file according to the NCBI Traceinfo RFC found at
http://www.ncbi.nlm.nih.gov/Traces/trace.cgi?cmd=show&f=rfc&m=main&s=rfc


bin_fasta_on_mid_primers.pl
---------------------------
Written by Gregory Harhay, United States Department of Agricultural,
Agricultural Research Service

A script to bin sequences coming from 454 sequencing according to their MID
primers.

The script needs BioPerl. A FASTA file containing example primers is also
provided (midi_screen.fasta).


qual2ball
---------
Written by Tony Travis, University of Aberdeen, Rowett Institute of Nutrition
and Health.

A script to create PHD balls for consed.



qual2ball
---------
Written by Tony Travis, University of Aberdeen, Rowett Institute of Nutrition
and Health.

A script to create PHD balls for consed.



caf2aceMiraConsed.pl
--------------------
Written by Lionel Guy, Uppsala University, Evolutionsbiologiskt centrum

It creates phdballs for solexa and 454 reads using consed, and then  
edits these (and possibly other phdball files) to take into account  
the edits done by mira, parsed from the caf file. It also edits the  
ace file produced by mira to make it a bit smaller, to set the correct  
time and dates, add the references to the phdball files, and also puts  
the tags produced by mira at the end of the ace file (so consed can  
read them).

