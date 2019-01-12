#!/bin/sh
# \
  exec tclsh "$0" ${1+"$@"}

#
# Written by Bastien Chevreux (BaCh)
#
# Copyright (C) 2007 and later by Bastien Chevreux
#
# All rights reserved.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the 
# Free Software Foundation, Inc., 
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
# 
#


# tool to convert the output of "sffinfo" into
#  - FASTA
#  - FASTA quality
#  - NCBI TRACEINFO XML

# the command line options are modeled after the ones found in the 
#  Roche Off-Instruments package "sffvolume" command

namespace eval sffi2mf {
    variable opts
    variable readinfo
}


proc sffi2mf::xmlencode {value} {
    return [string map {& &amp; ' &apos; < &lt; > &gt; \" &quot;} $value]
}

proc sffi2mf::xmlify {tout code value} {
    if {[string length $value]} {
	puts $tout "\t\t<$code>[xmlencode $value]</$code>"
    }
}

proc sffi2mf::writeTRACEINFOHeader {tout} {
    puts $tout "<?xml version=\"1.0\"?>"
    puts $tout "<trace_volume>"
}
proc sffi2mf::writeTRACEINFOTail {tout} {
    puts $tout "</trace_volume>"
}

proc sffi2mf::writeReadToFiles {fout qout tout} {
    variable opts
    variable readinfo

    #puts "Writing $readinfo(readname)"
    #puts $readinfo(readlen)
    #puts $readinfo(cql)
    #puts $readinfo(cqr)
    #puts $readinfo(cal)
    #puts $readinfo(car)
    #puts $readinfo(sequence)
    #puts $readinfo(qualities)

    puts $fout ">$readinfo(readname)"
    puts $fout "$readinfo(sequence)"

    puts $qout ">$readinfo(readname)"
    foreach qual $readinfo(qualities) {
	puts -nonewline $qout "$qual "
    }
    puts $qout ""

    puts $tout "\t<trace>"
    xmlify $tout trace_name $readinfo(readname)

    xmlify $tout center_name $opts(-center)
    xmlify $tout clone_id $opts(-clone)
    xmlify $tout seq_lib_id $opts(-seqlib)
    xmlify $tout source_type $opts(-source)
    xmlify $tout species_code $opts(-species)
    xmlify $tout strategy $opts(-strategy)

    xmlify $tout trace_type_code 454
    xmlify $tout program_id 454Basecaller

    if {$readinfo(cql) != 0 && $readinfo(cqr) != 0} {
	xmlify $tout clip_quality_left $readinfo(cql)
	xmlify $tout clip_quality_right $readinfo(cqr)
    }
    if {$readinfo(cal) != 0 && $readinfo(car) != 0} {
	xmlify $tout clip_vector_left $readinfo(cal)
	xmlify $tout clip_vector_right $readinfo(car)
    }
    puts $tout "\t</trace>"
}

proc sffi2mf::parseSFFInfoLine {line} {
    variable readinfo
    switch -regex $line {
	"^.*\# of Bases:" {
	    set readinfo(readlen) [lindex $line end]
	}
	"^.*Clip Qual Left:" {
	    set readinfo(cql) [lindex $line end]
	}
	"^.*Clip Qual Right:" {
	    set readinfo(cqr) [lindex $line end]
	}
	"^.*Clip Adap Left:" {
	    set readinfo(cal) [lindex $line end]
	}
	"^.*Clip Adap Right:" {
	    set readinfo(car) [lindex $line end]
	}
	"^Bases:" {
	    set readinfo(sequence) [lindex $line end]
	}
	"^Quality Scores:" {
	    set readinfo(qualities) [string map {"Quality Scores:" ""} $line]
	}
    }
}

proc sffi2mf::initNewRead {line} {
    variable readinfo
    array unset readinfo
    
    set readinfo(cql) 0
    set readinfo(cqr) 0
    set readinfo(cql) 0
    set readinfo(cqr) 0

    set readinfo(readname) [string map {> ""} [lindex $line 0]]
}

proc sffi2mf::processit {} {
    variable opts
    variable readinfo

    if {[string length $opts(-infile)] ==0} {
	set fin stdin
    } else {
	set fin [open $opts(-infile) r]
    }

    set fout [open $opts(-project)_in.454.fasta w]
    set qout [open $opts(-project)_in.454.fasta.qual w]
    set tout [open $opts(-project)_traceinfo_in.454.xml w]

    writeTRACEINFOHeader $tout

    while {[gets $fin line] != -1} {
	if {[string compare [string index $line 0] ">"] == 0} break
    }

    initNewRead $line

    while {[gets $fin line] != -1} {
	if {[string compare [string index $line 0] ">"] == 0} {
	    writeReadToFiles $fout $qout $tout
	    initNewRead $line
	} else {
	    parseSFFInfoLine $line
	}
    }
    writeReadToFiles $fout $qout $tout
    writeTRACEINFOTail $tout

    close $fin
    close $tout
    close $qout
    close $fout
}

proc usage {prgname emsg} {
    set prgname [lindex [split $prgname /] end]
    puts stderr "$prgname:
Converts the output of \"sffinfo\" into FASTA, FASTA quality and NCBI TRACEINFO
XML suited for assembly with the MIRA assembler.\n"

    puts stderr "Usage: $prgname ?options?

The command line options are modeled after the ones found in the Roche 
Off-Instruments package \"sffvolume\" command

Options:
\t-infile   filename      File to load sffinfo output from (stdin if ommited)
\t-project  string        Project base name

The following options (if set) must comply to NCBI TRACEINFO guidelines
\t-center   string        Sequencing center name
\t-clone    string        Clone name
\t-seqlib   string        SeqLib ID
\t-source   string        Sequencing center name
\t-species  string        Species name
\t-strategy string        Sequencing strategy

Examples:
$prgname -infile somefilefromsffinfo.txt
sffinfo EV10YMP01.sff | $prgname -project ecoli_k12 -species \"ESCHERICHIA COLI K12\"
"

    puts stderr $emsg

    exit
}

set sffi2mf::opts(-infile) ""
set sffi2mf::opts(-project) "mira"

set sffi2mf::opts(-center) ""
set sffi2mf::opts(-clone) ""
set sffi2mf::opts(-seqlib) ""
set sffi2mf::opts(-source) ""
set sffi2mf::opts(-species) ""
set sffi2mf::opts(-strategy) ""

set emsg ""
foreach {key val} $argv {
    if {![info exists sffi2mf::opts($key)]} {
	if {[string compare [string index $key 0] "-"] == 0} { 
	    append emsg "Bad option $key\n"
	}
	set val $key
	set key -infile
    }
    set sffi2mf::opts($key) $val
}
if {[string length $emsg]>0} {
    usage $argv0 $emsg
}

if {[string length $sffi2mf::opts(-project)] ==0} { usage $argv0 ; }

sffi2mf::processit
