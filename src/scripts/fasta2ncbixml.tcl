#!/bin/sh
# \
  exec tclsh "$0" ${1+"$@"}

#
# Written by Bastien Chevreux (BaCh)
#
# Copyright (C) 2004 and later by Bastien Chevreux
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

namespace eval fasta2ncbixml {
  variable opts

}


proc fasta2ncbixml::putspace {indent} {
    for {set i 0} {$i<$indent} {incr i} {
	puts -nonewline " "
    }
}

proc fasta2ncbixml::extractit {line} {
    variable opts

    set line [string range $line 1 end]

    set linefields [split $line]

#	    "name " {
#	    }
    foreach c $opts(-c) f $opts(-f) {
	switch -exact -- $f {
	    default {
		putspace 6
		puts "<$f>[lindex $linefields $c]</$f>"
	    }
	}
    }

}

proc fasta2ncbixml::processit {} {
    variable opts
    set fid [open $opts(-infile) r]

    puts "<?xml version=\"1.0\"?>"
    puts "  <trace_volume>"

    while {[gets $fid line] != -1} {
	if {[string compare [string index $line 0] ">"] == 0} {

	    putspace 4
	    puts "<trace>"
	    extractit $line
	    putspace 4
	    puts "</trace>"
	}
    }

    puts "  </trace_volume>"

    close $fid
}

proc usage {prgname} {
    puts stderr "$prgname: Extracts values from comments in fasta file and
writes them to NCBI xml format.\n"
    puts stderr "Usage: $prgname [-c ... -f ...]* ?options? infile
\t-c   int      Column number in fasta name line
\t-f   string   Field name of the above column for XML file
"
  exit
}

set fasta2ncbixml::opts(-c) 0
set fasta2ncbixml::opts(-f) trace_name
set fasta2ncbixml::opts(-infile) ""

foreach {key val} $argv {
    if {![info exists fasta2ncbixml::opts($key)]} {
	if {[string compare [string index $key 0] "-"] == 0} { 
	    puts stderr "Bad key $key\n"
	    usage $argv0
	}
	set val $key
	set key -infile
    }
    switch -exact -- $key {
	"-c" -
	"-f" {
	    lappend fasta2ncbixml::opts($key) $val
	}
	default {
	    set fasta2ncbixml::opts($key) $val
	}
    }
}

if {[string length $fasta2ncbixml::opts(-infile)] ==0} { usage $argv0 ; }

fasta2ncbixml::processit
