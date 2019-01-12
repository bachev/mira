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

    set worksheetname ${wsname}_info_featuresummary
    if {[string length $worksheetname]>16} {
	set worksheetname [string range $worksheetname 0 15]
    }

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
  <WindowHeight>8385</WindowHeight>
  <WindowWidth>27420</WindowWidth>
  <WindowTopX>720</WindowTopX>
  <WindowTopY>225</WindowTopY>
  <ProtectStructure>False</ProtectStructure>
  <ProtectWindows>False</ProtectWindows>
 </ExcelWorkbook>


 <Styles>
  <Style ss:ID=\"Default\" ss:Name=\"Normal\">
   <Alignment ss:Vertical=\"Bottom\"/>
   <Borders/>
   <Font ss:FontName=\"Calibri\" x:Family=\"Swiss\" ss:Size=\"11\" ss:Color=\"#000000\"/>
   <Interior/>
   <NumberFormat/>
   <Protection/>
  </Style>
  <Style ss:ID=\"s62\">
   <Alignment ss:Horizontal=\"Right\" ss:Vertical=\"Bottom\"/>
  </Style>
  <Style ss:ID=\"s68\">
   <Alignment ss:Horizontal=\"Center\" ss:Vertical=\"Bottom\"/>
  </Style>
  <Style ss:ID=\"s77\">
   <Alignment ss:Horizontal=\"Center\" ss:Vertical=\"Bottom\" ss:WrapText=\"1\"/>
   <Font ss:FontName=\"Calibri\" x:Family=\"Swiss\" ss:Color=\"#000000\"/>
   <Interior ss:Color=\"#D9E1F2\" ss:Pattern=\"Solid\"/>
   <NumberFormat ss:Format=\"@\"/>
  </Style>
  <Style ss:ID=\"s78\">
   <Interior ss:Color=\"#D9E1F2\" ss:Pattern=\"Solid\"/>
  </Style>
  <Style ss:ID=\"s89\">
   <Alignment ss:Horizontal=\"Left\" ss:Vertical=\"Bottom\"/>
   <Font ss:FontName=\"Calibri\" x:Family=\"Swiss\" ss:Color=\"#000000\"/>
  </Style>
 </Styles>
 <Worksheet ss:Name=\"${worksheetname}\">
  <Names>
   <NamedRange ss:Name=\"_FilterDatabase\"
    ss:RefersTo=\"=${worksheetname}!R1C1:R${entries(numentries)}C26\" ss:Hidden=\"1\"/>
  </Names>
  <Table ss:ExpandedColumnCount=\"28\" ss:ExpandedRowCount=\"${entries(numentries)}\" x:FullColumns=\"1\"
   x:FullRows=\"1\" ss:DefaultRowHeight=\"15\">
   <Column ss:StyleID=\"s78\" ss:AutoFitWidth=\"0\" ss:Width=\"155.25\"/>
   <Column ss:StyleID=\"s78\" ss:AutoFitWidth=\"0\" ss:Span=\"1\"/>
   <Column ss:Index=\"4\" ss:StyleID=\"s62\" ss:AutoFitWidth=\"0\"/>
   <Column ss:Index=\"6\" ss:StyleID=\"s68\" ss:AutoFitWidth=\"0\"/>
   <Column ss:AutoFitWidth=\"0\" ss:Width=\"43.5\"/>
   <Column ss:StyleID=\"s89\" ss:AutoFitWidth=\"0\" ss:Width=\"46.5\"/>
   <Column ss:StyleID=\"s68\" ss:AutoFitWidth=\"0\" ss:Width=\"46.5\" ss:Span=\"5\"/>
   <Column ss:Index=\"15\" ss:AutoFitWidth=\"0\" ss:Width=\"41.25\" ss:Span=\"8\"/>
   <Column ss:Index=\"24\" ss:StyleID=\"s78\"/>
   <Column ss:StyleID=\"s78\"/>
   <Column ss:StyleID=\"s78\" ss:AutoFitWidth=\"0\" ss:Width=\"177.75\"/>
   <Column ss:StyleID=\"s78\" ss:AutoFitWidth=\"0\" ss:Width=\"135\"/>
   <Column ss:StyleID=\"s78\" ss:AutoFitWidth=\"0\" ss:Width=\"622.5\"/>
   <Row ss:AutoFitHeight=\"0\" ss:Height=\"38.25\" ss:StyleID=\"s77\">
    <Cell><Data ss:Type=\"String\">Locus</Data><NamedCell ss:Name=\"_FilterDatabase\"/></Cell>
    <Cell><Data ss:Type=\"String\">Gene name</Data><NamedCell
      ss:Name=\"_FilterDatabase\"/></Cell>
    <Cell><Data ss:Type=\"String\">Feature type</Data><NamedCell
      ss:Name=\"_FilterDatabase\"/></Cell>
    <Cell><Data ss:Type=\"String\">Feature checksum</Data><NamedCell
      ss:Name=\"_FilterDatabase\"/></Cell>
    <Cell><Data ss:Type=\"String\">Genome map pos</Data><NamedCell
      ss:Name=\"_FilterDatabase\"/></Cell>
    <Cell><Data ss:Type=\"String\">Interesting?</Data><NamedCell
      ss:Name=\"_FilterDatabase\"/></Cell>
    <Cell><Data ss:Type=\"String\">Your own filter</Data><NamedCell
      ss:Name=\"_FilterDatabase\"/></Cell>
    <Cell><Data ss:Type=\"String\">Coverage status</Data><NamedCell
      ss:Name=\"_FilterDatabase\"/></Cell>
    <Cell><Data ss:Type=\"String\">First codon is start</Data><NamedCell
      ss:Name=\"_FilterDatabase\"/></Cell>
    <Cell><Data ss:Type=\"String\">Changed start codon</Data><NamedCell
      ss:Name=\"_FilterDatabase\"/></Cell>
    <Cell><Data ss:Type=\"String\">Destroyed start codon</Data><NamedCell
      ss:Name=\"_FilterDatabase\"/></Cell>
    <Cell><Data ss:Type=\"String\">Changed stop codon</Data><NamedCell
      ss:Name=\"_FilterDatabase\"/></Cell>
    <Cell><Data ss:Type=\"String\">Destroyed stop codon</Data><NamedCell
      ss:Name=\"_FilterDatabase\"/></Cell>
    <Cell><Data ss:Type=\"String\">Premature stop codon</Data><NamedCell
      ss:Name=\"_FilterDatabase\"/></Cell>
    <Cell><Data ss:Type=\"String\">Inter-genic mutation</Data><NamedCell
      ss:Name=\"_FilterDatabase\"/></Cell>
    <Cell><Data ss:Type=\"String\">Insertion in locus</Data><NamedCell
      ss:Name=\"_FilterDatabase\"/></Cell>
    <Cell><Data ss:Type=\"String\">Deletion in locus</Data><NamedCell
      ss:Name=\"_FilterDatabase\"/></Cell>
    <Cell><Data ss:Type=\"String\">Silent in locus</Data><NamedCell
      ss:Name=\"_FilterDatabase\"/></Cell>
    <Cell><Data ss:Type=\"String\">AA change in locus</Data><NamedCell
      ss:Name=\"_FilterDatabase\"/></Cell>
    <Cell><Data ss:Type=\"String\">Insertion untrans-lated</Data><NamedCell
      ss:Name=\"_FilterDatabase\"/></Cell>
    <Cell><Data ss:Type=\"String\">Deletion untrans-lated</Data><NamedCell
      ss:Name=\"_FilterDatabase\"/></Cell>
    <Cell><Data ss:Type=\"String\">Silent untrans-lated</Data><NamedCell
      ss:Name=\"_FilterDatabase\"/></Cell>
    <Cell><Data ss:Type=\"String\">AA change untranslated</Data><NamedCell
      ss:Name=\"_FilterDatabase\"/></Cell>
    <Cell><Data ss:Type=\"String\">GO process</Data><NamedCell ss:Name=\"_FilterDatabase\"/></Cell>
    <Cell><Data ss:Type=\"String\">GO function</Data><NamedCell ss:Name=\"_FilterDatabase\"/></Cell>
    <Cell><Data ss:Type=\"String\">Product</Data><NamedCell ss:Name=\"_FilterDatabase\"/></Cell>
    <Cell><Data ss:Type=\"String\">Function</Data><NamedCell
      ss:Name=\"_FilterDatabase\"/></Cell>
    <Cell><Data ss:Type=\"String\">Note</Data><NamedCell ss:Name=\"_FilterDatabase\"/></Cell>
   </Row>"
}


proc writeFooter { fout } {
    variable entries

    puts $fout "  </Table>
  <WorksheetOptions xmlns=\"urn:schemas-microsoft-com:office:excel\">
   <Selected/>
   <FreezePanes/>
   <FrozenNoSplit/>
   <FilterOn/>
   <SplitHorizontal>1</SplitHorizontal>
   <TopRowBottomPane>1</TopRowBottomPane>
   <SplitVertical>7</SplitVertical>
   <LeftColumnRightPane>7</LeftColumnRightPane>
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
     <ActiveRow>0</ActiveRow>
     <ActiveCol>0</ActiveCol>
    </Pane>
   </Panes>
   <ProtectObjects>False</ProtectObjects>
   <ProtectScenarios>False</ProtectScenarios>
  </WorksheetOptions>
  <AutoFilter x:Range=\"R1C1:R${entries(numentries)}C26\"
   xmlns=\"urn:schemas-microsoft-com:office:excel\">
   <AutoFilterColumn x:Index=\"3\" x:Type=\"Custom\">
    <AutoFilterOr>
     <AutoFilterCondition x:Operator=\"Equals\" x:Value=\"CDS\"/>
     <AutoFilterCondition x:Operator=\"Equals\" x:Value=\"Figr\"/>
    </AutoFilterOr>
   </AutoFilterColumn>
   <AutoFilterColumn x:Index=\"6\" x:Type=\"Custom\">
    <AutoFilterOr>
     <AutoFilterCondition x:Operator=\"Equals\" x:Value=\"perhaps\"/>
     <AutoFilterCondition x:Operator=\"Equals\" x:Value=\"yes\"/>
    </AutoFilterOr>
   </AutoFilterColumn>
  </AutoFilter>
  <ConditionalFormatting xmlns=\"urn:schemas-microsoft-com:office:excel\">
   <Range>R2C6:R${entries(numentries)}C6</Range>
   <Condition>
    <Qualifier>Equal</Qualifier>
    <Value1>&quot;yes&quot;</Value1>
    <Format Style='color:red;font-weight:700;text-line-through:none'/>
   </Condition>
   <Condition>
    <Qualifier>Equal</Qualifier>
    <Value1>&quot;perhaps&quot;</Value1>
    <Format Style='color:#ED7D31;font-weight:700'/>
   </Condition>
  </ConditionalFormatting>
  <ConditionalFormatting xmlns=\"urn:schemas-microsoft-com:office:excel\">
   <Range>R2C8:R${entries(numentries)}C8</Range>
   <Condition>
    <Qualifier>Equal</Qualifier>
    <Value1>&quot;ok&quot;</Value1>
    <Format Style='color:#548235;font-weight:700;text-line-through:none'/>
   </Condition>
   <Condition>
    <Qualifier>NotEqual</Qualifier>
    <Value1>&quot;ok&quot;</Value1>
    <Format Style='color:red;font-weight:700;text-line-through:none'/>
   </Condition>
  </ConditionalFormatting>
  <ConditionalFormatting xmlns=\"urn:schemas-microsoft-com:office:excel\">
   <Range>R2C9:R${entries(numentries)}C14</Range>
   <Condition>
    <Qualifier>Equal</Qualifier>
    <Value1>&quot;n/a&quot;</Value1>
    <Format Style='color:white'/>
   </Condition>
  </ConditionalFormatting>
  <ConditionalFormatting xmlns=\"urn:schemas-microsoft-com:office:excel\">
   <Range>R2C15:R${entries(numentries)}C23</Range>
   <Condition>
    <Qualifier>Equal</Qualifier>
    <Value1>0</Value1>
    <Format Style='color:white'/>
   </Condition>
   <Condition>
    <Qualifier>NotEqual</Qualifier>
    <Value1>0</Value1>
    <Format Style='color:red;font-weight:700'/>
   </Condition>
  </ConditionalFormatting>
 </Worksheet>
</Workbook>
"
}


proc writeCell { fout style type row index} {
    variable entries

    if {[string length $style]} {
	puts -nonewline $fout "    <Cell ss:StyleID=\"${style}\">"
    } else {
	puts -nonewline $fout "    <Cell>"
    }
    puts $fout "<Data ss:Type=\"$type\">[lindex $entries($row) $index]</Data><NamedCell ss:Name=\"_FilterDatabase\"/></Cell>"
}

proc writeBody { fout } {
    variable entries

    #for {set row 1} {$row < 3} {incr row}
    for {set row 1} {$row < $entries(numentries)} {incr row} {
	set interesting [lindex $entries($row) 4]
	if {$interesting != "yes"} {
	    puts $fout "   <Row ss:Hidden=\"1\">"
	} else {
	    puts $fout "   <Row>"
	}
	set index 0
	writeCell $fout "" String $row $index
	incr index
	writeCell $fout "" String $row $index
	incr index
	writeCell $fout "" String $row $index
	incr index
	writeCell $fout "" String $row $index
	incr index
	writeCell $fout "" Number $row $index
	incr index
	writeCell $fout "" String $row $index
	incr index
	writeCell $fout "" String $row $index
	incr index

	# coverage status
	writeCell $fout "" String $row $index
	incr index

	# first codon is start
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

	# intergenic mutation<
	writeCell $fout "" Number $row $index
	incr index

	# xxx in locus
	writeCell $fout "" Number $row $index
	incr index
	writeCell $fout "" Number $row $index
	incr index
	writeCell $fout "" Number $row $index
	incr index
	writeCell $fout "" Number $row $index
	incr index

	# xxx untranslated
	writeCell $fout "" Number $row $index
	incr index
	writeCell $fout "" Number $row $index
	incr index
	writeCell $fout "" Number $row $index
	incr index
	writeCell $fout "" Number $row $index
	incr index

	# GO product & function
	writeCell $fout "" String $row $index
	incr index
	writeCell $fout "" String $row $index
	incr index

	# product, function, note
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

regsub _info_featuresummary.txt $project {} project

puts "Working on $project"

readfile ${project}_info_featuresummary.txt
writefile ${project}_info_featuresummary.xml ${project}
