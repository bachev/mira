#!/bin/sh
# \
  exec tclsh "$0" ${1+"$@"}

#
# Written by Bastien Chevreux (BaCh)
#
# Copyright (C) 2002 and later by Bastien Chevreux
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


namespace eval phd2exp {
  variable opts

  variable printout
  variable sequence
  variable scyllaparamfile /tmp/phd2exp_scylladummyparam.[pid]

  set tags {
    ID TN EN LN LE LT QL QR SV SF SL SR SQ
  }
}

proc phd2exp::firstData {fn} {
  variable printout
  variable sequence

  catch {unset printout}

  set fid [open $fn r+]
  set data [read $fid]
  close $fid

  #########
  set printout(ID) ""
  regexp -line {^BEGIN_SEQUENCE (.*?)$} $data -> printout(ID)
  if {[string length $printout(ID)] == 0} {
    return -code error "Couldn't find BEGIN_SEQUENCE in phd file: $fn"
  }

  #########
  set sequence [extractSequence $data]
  set printout(LE) [expr {[string length $sequence] -1}]

  #########
}

proc phd2exp::extractSequence {data} {
  set begin [expr {[string first BEGIN_DNA $data] + [string length BEGIN_DNA]}]
  set end [expr {[string first END_DNA $data] -1}]

  set seqdata [string range $data $begin $end]

  set sequence ""
  foreach {char qual ppos} $seqdata {
    append sequence [string toupper $char]
  }
  return $sequence
}

proc phd2exp::clipSequence {fn reffile wordlen threshold} {
  variable scyllaparamfile
  variable sequence 

  set retvalues [list 0 [string length $sequence]]
  set maskfile /tmp/phd2exp_maskout.[pid]

  if {![file exists $scyllaparamfile]} {
    exec touch $scyllaparamfile
  }

  #puts "catch {exec scylla -Query=$fn -Reference=$reffile -Output=$maskfile -Param=$scyllaparamfile -align -MaskChar ~ -WordLen $wordlen -Threshold $threshold}"
  if {[catch {exec scylla -Query=$fn -Reference=$reffile -Output=$maskfile -Param=$scyllaparamfile -align -MaskChar ~ -WordLen $wordlen -Threshold $threshold} res]} {
    puts stderr "Scylla failed for $fn: $res"
    return $retvalues
  }

  if {[catch {
    set fid [open $maskfile r+]
    set data [read $fid]
    close $fid
  } res]} {
    puts stderrr "Failed reading masked file $maskfile?"
    return $retvalues
  }
  catch {file delete $maskfile}

  set maskseq [extractSequence $data]
  set seqlen [string length $maskseq]
  puts $maskseq
  set mask1ststart [string first ~ $maskseq]
  if {$mask1ststart >=0} {
    set mask1stend $mask1ststart
    while {$mask1stend < $seqlen \
	     && [string index $maskseq $mask1stend] == "~"} {incr mask1stend}
    #puts $mask1ststart
    #puts $mask1stend
    #puts [string range $maskseq 0 $mask1stend]
    set mask2ndend [string last ~ $maskseq]
    set hastwoblocks 0
    if {$mask2ndend > $mask1stend} {
      set hastwoblocks 1
      set mask2ndstart $mask2ndend
      while {$mask2ndstart > 0 \
	       && [string index $maskseq $mask2ndstart] == "~"} {incr mask2ndstart -1}
      incr mask2ndstart
    }
    set mid [expr {$seqlen/2}]
    if {$hastwoblocks} {
      set left 0

      # NONO! es kann auch ganz viel frei am Anfang sein :(
      #
      #if {$mask1ststart >50 && $mask1ststart > ($mask1stend-$mask1ststart)} {
      #	# hhhmmm ... flachfallen lassen, keine linke Grenze?
      #} else {
      #	set left $mask1stend
      #} 
      set left $mask1stend
      
      set right $seqlen

      # NONO! es kann auch ganz viel frei am Ende sein :(
      #
      #if {($seqlen-$mask2ndend) >50 && ($seqlen-$mask2ndend) > ($mask2ndend-$mask2ndstart)} {
      #	# hhhmmm ... flachfallen lassen, keine rechte Grenze?
      #} else {
      #	set right $mask2ndstart
      #}
      set right $mask2ndstart

      set retvalues [list $left $right]
    } else {
      # now, is this block on the left or the right?
      set isleft 1
      if {$mask1ststart>$mid} {
	set isleft 0
      } else {
	if {$mask1stend>$mid} {
	  # check on which side there are less bases and mask these
	  if {($seqlen - $mask1stend) < $mask1ststart} {
	    set isleft 0
	  }
	}
      }
      if {$isleft} {
	set retvalues [list $mask1stend [lindex $retvalues 1]]
      } else {
	set retvalues [list 0 $mask1ststart]
      }
    }
  }
  return $retvalues
}


proc phd2exp::setSVClips {fn} {
  variable printout
  variable opts

  if {[string lengt $opts(-seqvector)] == 0} {
    return
  }

  #########
  set clips [phd2exp::clipSequence \
	       $fn \
	       $opts(-seqvector) \
	       6 12]
  set printout(SF) $opts(-seqvector)
  set printout(SL) [lindex $clips 0]
  set printout(SR) [lindex $clips 1]
}

proc phd2exp::setQualClips {} {
  variable printout
  variable opts

  set scfn [file join $opts(-scfdir) $printout(LN)]
  if {[catch {set n [exec clipSCFbyQual $scfn 30 5 20 20 0]} res]} {
    puts stderr "Error loading $scfn"
    return
  }
  foreach {start end numbases avgqual} $n {}

  set printout(QL) $start
  set printout(QR) $end
}

proc phd2exp::extendData {} {
  variable printout
  variable sequence

  #########
  set printout(LN) $printout(ID).scf

  #########
  set printout(EN) $printout(ID)

  #########
  set sq \n
  set counter 0
  foreach char [split $sequence ""] {
    if {$counter == 0} {
      append sq "     "
    }
    append sq $char
    incr counter
    if {$counter == 60} {
      append sq \n
      set counter 0
      continue
    }
    if {[expr {$counter % 10}] == 0} {
      append sq " "
    }
  }
  if {$counter >0} {
    append sq \n
  }
  append sq "//"
  set printout(SQ) $sq

  #########
  set printout(TN) ""
  regexp {(.*?).[p|q][0-9][a-z]?[a-z]?.*} $printout(ID) -> printout(TN)

  #########
  set printout(LT) SCF

  #########
  #set printout(SV) p3t

}



proc phd2exp::writeEXP {} {
  variable printout
  variable tags
  variable opts

  set fon [file join $opts(-targetdir) $printout(ID).exp]
  if {[catch {
    set fod [open $fon w]
  } res]} {
    puts stderr "Could not open $fon for writing: $res"
    return
  }

  foreach tag $tags {
    if {[info exists printout($tag)]} {
      if {[string length $printout($tag)]} {
	puts $fod "$tag   $printout($tag)"
	unset printout($tag)
      }
    }
  }
  #puts [array get phd2exp::printout]
  close $fod
}

proc phd2exp::processFile {fn} {
  phd2exp::firstData $fn
  phd2exp::setSVClips $fn
  phd2exp::extendData
  phd2exp::setQualClips
  phd2exp::writeEXP
  #puts [array get phd2exp::printout]
}  


proc phd2exp::processList {fln} {

  if {[catch {
    set fid [open $fln r+]
    set data [read $fid]
    close $fid
  } res ]} {
    puts stderr "Error loading list of files $fln: $res"
    return
  }
  set lines [split $data \n]
  set numlines [llength $lines]
  set count 1
  foreach line $lines {
    set fn [string trim $line]
    if {[string length $fn] == 0} { 
      puts "Empty line"
      continue
    }
    puts "Processing $count/$numlines $fn"
    processFile $fn
    incr count
  }
}

proc usage {prgname} {
    puts stderr "Usage: $prgname ?options?
\t-phdfile   file      (Required) File containing names of phd files
\t                      to convert
\t-targetdir directory Directory to write resulting exp into. Will be
\t                      created if does not exist. Default: exp
\t-seqvector file      File containing sequencing vectors to screen 
\t                      against and use for clipping
\t-scfdir    directory Name of directory containing scf files.
"
  exit
}

set phd2exp::opts(-phdfile) ""
set phd2exp::opts(-targetdir) exp
set phd2exp::opts(-seqvector) ""
set phd2exp::opts(-scfdir) ""

foreach {key val} $argv {
  if {![info exists phd2exp::opts($key)]} {
    puts stderr "Bad key $key\n"
    usage $argv0
  }
  #puts "set $key"
  set phd2exp::opts($key) $val
}

if {[string length $phd2exp::opts(-phdfile)] == 0} {
  puts stderr "Missing -phdfile option!"
  usage $argv0
}
if {![file exists $phd2exp::opts(-targetdir)]} {
  if {[catch {file mkdir $phd2exp::opts(-targetdir)} res]} {
    puts stderr "Could not create target directory $phd2exp::opts(-targetdir):\n$res"
    exit
  }
}

set phd2exp::opts(-seqvector) [string trim $phd2exp::opts(-seqvector)]

phd2exp::processList $phd2exp::opts(-phdfile)

file delete $phd2exp::scyllaparamfile

