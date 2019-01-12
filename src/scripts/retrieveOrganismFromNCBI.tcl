#
# Written by Bastien Chevreux (BaCh)
#
# Copyright (C) 2005 and later by Bastien Chevreux
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

proc getNumReadsFromNCBI {species} {

    exec query_tracedb.pl "query count species_code='$species'" >/tmp/ttt

    set fid [open /tmp/ttt r]
    gets $fid numofreads
    close $fid
    
    return [expr {$numofreads}]
}

proc getReadListFromNCBI {species} {
    set maxgo 40000

    set numreads [getNumReadsFromNCBI $species]
    puts "$species has $numreads reads."
    
    exec echo -n "" >/tmp/ttt1
    set numpages [expr {$numreads/$maxgo}]
    incr numpages
    puts -nonewline "Must first retrieve $numpages information pages from NCBI:"

    for {set i 0} {$i < $numpages} {incr i} {
	puts -nonewline " $i ..."
	flush stdout
	exec query_tracedb.pl "query page_size $maxgo page_number $i text species_code='$species'" >>/tmp/ttt1
    }
    exec sort -g /tmp/ttt1 >/tmp/ttt2
    
    puts " done."
    set fid [open /tmp/ttt2 r]
    set data [read $fid]
    close $fid

    return $data
}

proc tilist2ranges {tilist} {
    set step 0
    lappend tilist endofgame
    foreach entry $tilist {
	switch $step {
	    0 {
		set from $entry
		set last $entry
		set step 1
	    }
	    1 { 
		if {$last+1 < $entry 
		    || $entry == "endofgame" 
		    || $entry-$from == 40000} {
		    set to $last
		    if {$to-$from == 0} {
			lappend data $to 1
		    } else {
			lappend data "$from-$to" [expr {$to-$from+1}]
		    }
		    set from $entry
		}
		set last $entry
	    }
	}
    }

    #puts "Chunks: $data"
    return $data
}


proc makeDownloadRanges {ranges} {
    set numindl 0
    set actchunk ""
    foreach {range numelem} $ranges {
	if {$numindl+$numelem > 40000} {
	    lappend data $actchunk $numindl
	    set actchunk ""
	} 
	if {[string length $actchunk] >0} {
	    append actchunk ",$range"
	    incr numindl $numelem
	} else {
	    set actchunk $range
	    set numindl $numelem
	}
    }
    lappend data $actchunk $numindl
    #puts "DLChunks: $data"
    return $data
}

proc downloadFromNCBI {species shortcode} {
    set tilist [getReadListFromNCBI $species]
    set ranges [tilist2ranges $tilist]
    set dlranges [makeDownloadRanges $ranges]

    set i 0
    foreach {dlrange numelem} $dlranges {
	puts "Now downloading $numelem reads from NCBI. TI numbers: $dlrange"
	set filename ${shortcode}_part$i.tar.gz
	exec query_tracedb.pl "retrieve_tar_gz fasta quality xml $dlrange" >$filename
	incr i
    }
}




#set species "BACILLUS ANTHRACIS STR. AMES"
#set shortcode bant_ames
#set species "BACILLUS ANTHRACIS STR. FRANCE"
#set shortcode bant_france
#set species "BACILLUS ANTHRACIS STR. KRUGER B"
#set shortcode bant_krugerb
#set species "BACILLUS ANTHRACIS STR. VOLLUM"
#set shortcode bant_vollum
set species "BACILLUS ANTHRACIS STR. WESTERN NORTH AMERICA"
set shortcode bant_wna


downloadFromNCBI $species $shortcode
