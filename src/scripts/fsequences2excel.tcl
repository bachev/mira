#!/bin/sh
# \
  exec tclsh "$0" ${1+"$@"}

#
# Written by Bastien Chevreux (BaCh)
#
# Copyright (C) 2016 and later by Bastien Chevreux
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

package require csv


proc readfile { filename } {
    variable entries

    set fin [open $filename r]

    set loopit 1
    while {$loopit} {
	set line [gets $fin]
	if {[string index $line 0] != "#"} {
	    set loopit 0
	} else {
	    puts "Header: $line"
	}
    }

    set loopit 1
    set entrynum 1
    while {![eof $fin]} {
	set linelist [::csv::split $line \t]
	if {[llength $linelist]} {
	    #puts "$entrynum $line"
	    set entries($entrynum) $linelist
	    incr entrynum
	} else {
	    puts "seen empty line???"
	}
	set line [gets $fin]
    }
    close $fin

    # last line
    set linelist [::csv::split $line \t]
    if {[llength $linelist]} {
	set entries($entrynum) $linelist
	incr entrynum
    }

    set entries(numentries) $entrynum
}

proc writeHeader { fout wsname } {
    variable entries

    puts $fout "<?xml version=\"1.0\"?>
<?mso-application progid=\"Excel.Sheet\"?>
<Workbook xmlns=\"urn:schemas-microsoft-com:office:spreadsheet\"
 xmlns:o=\"urn:schemas-microsoft-com:office:office\"
 xmlns:x=\"urn:schemas-microsoft-com:office:excel\"
 xmlns:ss=\"urn:schemas-microsoft-com:office:spreadsheet\"
 xmlns:html=\"http://www.w3.org/TR/REC-html40\">
 <DocumentProperties xmlns=\"urn:schemas-microsoft-com:office:office\">
  <LastAuthor>MIRA</LastAuthor>
  <Created>2008-08-16T00:00:00Z</Created>
  <LastSaved>2008-08-16T00:00:00Z</LastSaved>
  <Version>11.9999</Version>
 </DocumentProperties>
 <ExcelWorkbook xmlns=\"urn:schemas-microsoft-com:office:excel\">
  <WindowHeight>8580</WindowHeight>
  <WindowWidth>27660</WindowWidth>
  <WindowTopX>300</WindowTopX>
  <WindowTopY>60</WindowTopY>
  <ProtectStructure>False</ProtectStructure>
  <ProtectWindows>False</ProtectWindows>
 </ExcelWorkbook>
 <Styles>
  <Style ss:ID=\"Default\" ss:Name=\"Normal\">
   <Alignment ss:Vertical=\"Bottom\"/>
   <Borders/>
   <Font/>
   <Interior/>
   <NumberFormat/>
   <Protection/>
  </Style>
  <Style ss:ID=\"s21\">
   <Interior ss:Color=\"#CCFFCC\" ss:Pattern=\"Solid\"/>
  </Style>
  <Style ss:ID=\"s22\">
   <Alignment ss:Vertical=\"Bottom\" ss:WrapText=\"1\"/>
   <Interior ss:Color=\"#CCFFCC\" ss:Pattern=\"Solid\"/>
  </Style>
 </Styles>
 <Worksheet ss:Name=\"${wsname}_info_featuresequences\">
  <Table ss:ExpandedColumnCount=\"10\" ss:ExpandedRowCount=\"${entries(numentries)}\" x:FullColumns=\"1\"
   x:FullRows=\"1\">
   <Column ss:StyleID=\"s21\" ss:AutoFitWidth=\"0\"/>
   <Column ss:StyleID=\"s21\" ss:AutoFitWidth=\"0\" ss:Width=\"132.75\"/>
   <Column ss:StyleID=\"s21\" ss:AutoFitWidth=\"0\" ss:Span=\"2\"/>
   <Column ss:Index=\"6\" ss:AutoFitWidth=\"0\" ss:Width=\"141\" ss:Span=\"3\"/>
   <Row ss:AutoFitHeight=\"0\" ss:Height=\"30.75\" ss:StyleID=\"s22\">
    <Cell><Data ss:Type=\"String\"># Contig</Data></Cell>
    <Cell><Data ss:Type=\"String\">FType</Data></Cell>
    <Cell><Data ss:Type=\"String\">Locus</Data></Cell>
    <Cell><Data ss:Type=\"String\">Gene name</Data></Cell>
    <Cell><Data ss:Type=\"String\">Feature strain</Data></Cell>
    <Cell><Data ss:Type=\"String\">Mutant strain</Data></Cell>
    <Cell><Data ss:Type=\"String\">Protein in reference</Data></Cell>
    <Cell><Data ss:Type=\"String\">Protein in $wsname</Data></Cell>
    <Cell><Data ss:Type=\"String\">DNA in reference</Data></Cell>
    <Cell><Data ss:Type=\"String\">DNA in $wsname</Data></Cell>
   </Row>"
}

proc writeFooter { fout } {
    variable entries

    puts $fout "  </Table>
  <WorksheetOptions xmlns=\"urn:schemas-microsoft-com:office:excel\">
   <Selected/>
   <FreezePanes/>
   <FrozenNoSplit/>
   <SplitHorizontal>1</SplitHorizontal>
   <TopRowBottomPane>1</TopRowBottomPane>
   <SplitVertical>5</SplitVertical>
   <LeftColumnRightPane>5</LeftColumnRightPane>
   <ActivePane>0</ActivePane>
   <Panes>
    <Pane>
     <Number>3</Number>
    </Pane>
    <Pane>
     <Number>1</Number>
    </Pane>
    <Pane>
     <Number>2</Number>
    </Pane>
    <Pane>
     <Number>0</Number>
    </Pane>
   </Panes>
   <ProtectObjects>False</ProtectObjects>
   <ProtectScenarios>False</ProtectScenarios>
  </WorksheetOptions>
 </Worksheet>
</Workbook>"

}


proc writeCell { fout style type row index} {
    variable entries

    if {[string length $style]} {
	puts -nonewline $fout "    <Cell ss:StyleID=\"${style}\">"
    } else {
	puts -nonewline $fout "    <Cell>"
    }
    puts $fout "<Data ss:Type=\"$type\">[lindex $entries($row) $index]</Data></Cell>"
}

proc writeBody { fout } {
    variable entries

    #for {set row 1} {$row < 3} {incr row}
    for {set row 1} {$row < $entries(numentries)} {incr row} {
	set interesting [lindex $entries($row) 4]
	puts $fout "   <Row>"

	set index 0
	writeCell $fout "" String $row $index
	incr index
	writeCell $fout "" String $row $index
	incr index
	writeCell $fout "" String $row $index
	incr index
	writeCell $fout "" String $row $index
	incr index
	writeCell $fout "" String $row $index
	incr index
	writeCell $fout "" String $row $index
	incr index
	writeCell $fout "" String $row $index
	incr index
	writeCell $fout "" String $row $index
	incr index
	writeCell $fout "" String $row $index
	incr index
	writeCell $fout "" String $row $index
	incr index

	puts $fout "   </Row>"
    }
}


proc writefile { filename wsname } {
    set fout [open $filename w]

    writeHeader $fout $wsname
    writeBody $fout
    writeFooter $fout
}


set project [lindex $argv 0]

puts "Working on $project"

readfile ${project}_info_featuresequences.txt
writefile ${project}_info_featuresequences.xml ${project}
