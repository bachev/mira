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

proc readExtractNames {filename} {
    upvar interestedin interestedin

    set fid [open $filename r]
    while {[gets $fid line] != -1} {
	set line [string map {" " ""} $line]
	set interestedin($line) 1
    }
    close $fid

}

proc doit {} {
    upvar interestedin interestedin

    set state 0
    while {[gets stdin line] != -1} {
	switch $state {
	    0 {
		set printout 0
		if {[regexp {<trace>} $line trace]} {
		    set state 1
		    set actstring $line
		    #puts "found trace"
		} else {
		    puts $line
		}
	    }
	    1 {
		if {[regexp {<trace_name>(.*)</trace_name>} $line . tracename]} {
		    set state 2
		    #puts "found $tracename"
		    if {[info exists interestedin($tracename)]} {
			#puts "found $tracename"
			set printout 1
		    }
		}
		append actstring "\n$line"
	    }
	    2 {
		append actstring "\n$line"
		if {[regexp {</trace>} $line trace]} {
		    set state 0
		    #puts "found /trace"
		    if {$printout>0} {
			puts $actstring
		    }
		}
	    }
	}
    }
}

readExtractNames [lindex $argv 0]
doit
