#!/bin/sh
# \
  exec tclsh "$0" ${1+"$@"}

#
# Written by Bastien Chevreux (BaCh)
#
# Copyright (C) 2009 and later by Bastien Chevreux
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

# Consed (up till at least version 19) has a bug: it can read
# consensus tags of an ACE only if they're at the end. This script
# reorganises an ACE so that consed can also read the consensus
# tags

# reads an ACE file and dumps it to STDOUT while
# placing all consensus tags to the end 

proc filterACE {filename} {

    set inct 0
    set fid [open $filename r]
    while {[gets $fid line] != -1} {
	set trimmedline [string trim $line]
	if {$inct} {
	    append ctags $line\n
	    if {$trimmedline == "\}"} {
		set inct 0
		append ctags \n
	    }
	} else {
	    if {$trimmedline == "CT\{"} {
		set inct 1
		append ctags $line\n
	    } else {
		puts $line
	    }
	}
    }
    close $fid

    puts $ctags
}

filterACE [lindex $argv 0]

