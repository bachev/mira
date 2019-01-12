proc protectExpMatch {up1 up2} {
    upvar $up1 exp
    upvar $up2 mat
    set exp [string map {
	\\ \\\\ 
	\{ \\\{ 
	\} \\\} 
	\[ \\\[ 
	\] \\\] 
	( \\( 
	) \\) 
	| \\| 
	* \\* 
	+ \\+ 
	. \. 
	? \\?} $exp]
    set mat [string map {& \\\&} $mat]
}

proc rawdocbookify {up1} {
    upvar $up1 what
    set what [string map {
	& &amp;
	<= &le;
	< &lt;
	>= &ge;
	> &gt;
	' &apos;
	\" &quot;
    } $what]
}

proc closeSect {seclevel tolevel} {
    while {$seclevel>=$tolevel} {
	puts "</sect$seclevel>"
	incr seclevel -1
    }
}

proc openSect {up1 newsec sectitle} {
    upvar $up1 seclevel
    if {$seclevel>0} {
	closeSect $seclevel $newsec
    }
    set id [string tolower [string trim $sectitle]]
    # get internal id clean of funny things
    set id [string map {
	" " _ 
	" / " _
	/ _
	= "" 
	& "" 
	. "" 
	, "" 
	( "" 
	) "" 
	{ "" 
	} "" 
	[ "" 
	] "" 
	< "" 
	> "" 
	- ""
	\" ""
    } $id]
    puts "<sect$newsec id=\"sect${newsec}_$id\">\n<title>\n$sectitle\n</title>\n<para>"
    set seclevel $newsec
}

proc closeInpara {} {
    upvar inpara inpara
    if {$inpara} {
	puts "</para>"
	set inpara 0
    }
}

proc doit {filename} {
    set fin [open $filename r]

    set inenums 0
    set closeenums 0
    set initems 0
    set closeitems 0
    set indescs 0
    set closedescs 0
    set newlinepara 0
    set inpara 0
    #set lastsect ""
    set seclevel 0
    while {[gets $fin line]>=0} {
	set strim [string trim $line]
	if {[string length $strim]==0} {
	    set newlinepara 1
	    continue
	}
	if {[string index $line 0] == "%"} {
	    puts "<!--$line-->"
	    continue
	}
	#set line [string map {\\_ _ $>$ > $>=$ >= $<$ < $<=$ <= $|$ ORDEFINE1 | ORDEFINE2} $line]
	set line [string map {\\_ _ $>$ > $>=$ >= $<$ < $<=$ <= $|$ | $=$ =} $line]
	#puts "NEWLINE: $line"

	if {[regexp {\\begin\{verbatim\}} $line]} {
	    closeInpara
	    puts "<screen>"
	    while {[gets $fin line]>=0} {
		if {[regexp {\\end\{verbatim\}} $line]} {
		    puts "</screen>\n<para>"
		    set inpara 1
		    break
		} else {
		    rawdocbookify line
		    puts $line
		}
	    }
	    continue
	}

	if {[regexp {\\section\{([^\}]+)\}} $line wholeexp sectitle]} {
	    closeInpara
	    openSect seclevel 1 $sectitle
	    set newlinepara 0
	    set inpara 1
	    continue
	}
	if {[regexp {\\subsection\{([^\}]+)\}} $line wholeexp sectitle]} {
	    closeInpara
	    openSect seclevel 2 $sectitle
	    set newlinepara 0
	    set inpara 1
	    continue
	}
	if {[regexp {\\subsubsection\{([^\}]+)\}} $line wholeexp sectitle]} {
	    closeInpara
	    openSect seclevel 3 $sectitle
	    set newlinepara 0
	    set inpara 1
	    continue
	}
	while {[regexp {\\URL\{([^\}]+)\}} $line wholeexp urlname]} {
	    protectExpMatch wholeexp urlname
	    set urlname [string map {\\\\& \\&amp;} $urlname]
	    regsub -all $wholeexp $line "<ulink url=\"$urlname\"/>" line
	}
	while {[regexp {\\Email\{([^\}]+)\}} $line wholeexp name]} {
	    protectExpMatch wholeexp name
	    regsub -all $wholeexp $line "<email>$name</email>" line
	}
	while {[regexp {\\File\{([^\}]+)\}} $line wholeexp filename]} {
	    protectExpMatch wholeexp filename
	    regsub -all $wholeexp $line "<filename>$filename</filename>" line
	}
	while {[regexp {\\emph\{([^\}]+)\}} $line wholeexp emphwhat]} {
	    protectExpMatch wholeexp emphwhat
	    regsub -all $wholeexp $line "<emphasis>$emphwhat</emphasis>" line
	}
	while {[regexp {\\textbf\{([^\}]+)\}} $line wholeexp emphwhat]} {
	    protectExpMatch wholeexp emphwhat
	    regsub -all $wholeexp $line "<emphasis role=\"bold\">$emphwhat</emphasis>" line
	}
	while {[regexp {\\underline\{([^\}]+)\}} $line wholeexp emphwhat]} {
	    protectExpMatch wholeexp emphwhat
	    regsub -all $wholeexp $line "<emphasis role=\"underline\">$emphwhat</emphasis>" line
	}
	while {[regexp {\\texttt\{([^\}]+)\}} $line wholeexp ttwhat]} {
	    protectExpMatch wholeexp ttwhat
	    regsub -all $wholeexp $line "<literal>$ttwhat</literal>" line
	}
	while {[regexp {\\Prog\{([^\}]+?)\}} $line wholeexp progwhat]} {
	    protectExpMatch wholeexp progwhat
	    regsub -all $wholeexp $line "$progwhat" line
	}
	while {[regexp {\\Cmd\{([^\}]+?)\}\{([^\}]+?)\}} $line wholeexp cmdwhat cmdsec]} {
	    protectExpMatch wholeexp cmdwhat
	    regsub -all $wholeexp $line "<command>$cmdwhat</command>" line
	}
	if {[regexp {\\begin\{enumerate\}} $line]} {
	    set closeenums 0
	    set inenums 1
	    closeInpara
	    set line "<orderedlist>"
	}
	if {[regexp {\\end\{enumerate\}} $line]} {
	    if {$closeenums} {
		set line "</para>\n</listitem>\n"
	    }
	    append line "</orderedlist>\n<para>"
	    set closeenums 0
	    set inenums 0
	    set inpara 1
	}
	if {[regexp {\\begin\{itemize\}} $line]} {
	    set closeitems 0
	    set initems 1
	    closeInpara
	    set line "<itemizedlist>"
	}
	if {[regexp {\\end\{itemize\}} $line]} {
	    if {$closeitems} {
		set line "</para>\n</listitem>\n"
	    }
	    append line "</itemizedlist>\n<para>"
	    set closeitems 0
	    set initems 0
	    set inpara 1
	}
	if {[regexp {\\begin\{description\}} $line]} {
	    if {$indescs} {
		error "Eeek5 desc: $line"
	    }
	    if {$inenums} {
		error "Eeek5 enum: $line"
	    }
	    if {$initems} {
		error "Eeek5 item: $line"
	    }
	    set closedescs 0
	    set indescs 1
	    closeInpara
	    set line "<variablelist>"
	}
	if {[regexp {\\end\{description\}} $line]} {
	    if {!$closedescs} {
		error "Eeeek6"
	    }
	    set line "    </para>\n    </listitem>\n  </varlistentry>\n</variablelist>\n<para>"
	    set indescs 0
	    set closedescs 0
	    set inpara 1
	}
	set initemline 0
	if {[regexp -all {\\item\[([^\]]+?)\]} $line wholeexp match]} {
	    protectExpMatch wholeexp match
	    if {$closedescs} {
		regsub -all $wholeexp $line "</para>\n</listitem>\n  </varlistentry>\n  <varlistentry>\n    <term>\n$match\n</term>\n    <listitem>\n<para>" line
	    } else {
		regsub -all $wholeexp $line "  <varlistentry>\n    <term>\n$match\n</term>\n    <listitem>\n<para>" line
	    }
	    set closedescs 1
	    set initemline 1
	}
	if {[regexp -all {\\item} $line]} {
	    if {$indescs && ($inenums || $initems)} {
		error "Eeek2 line: $line"
	    } elseif {$closeenums || $closeitems} {
		regsub -all {\\item} $line "</para>\n</listitem>\n<listitem>\n<para>\n" line
		if {$closeitems && $closeenums} {
		    error "Eeek"
		}
	    } else {
		if {$closeitems && $closeenums} {
		    error "Eeek3"
		} elseif {$initems} {
		    set closeitems 1
		} elseif {$inenums} {
		    set closeenums 1
		} else {
		    error "Eeeek4"
		}
		regsub -all {\\item} $line "<listitem>\n<para>\n" line
	    }
	}

	while {[regexp {\\OptArg\{([^\}]*?)\}\{([^\}]*?)\}} $line wholeexp m1 m2]} {
	    protectExpMatch wholeexp m1
	    if {$initemline} {
		regsub -all $wholeexp $line "<arg>$m1<replaceable>$m2</replaceable></arg>" line
	    } else {
		regsub -all $wholeexp $line "<arg>$m1$m2</arg>" line
	    }
	}
	# only in synopsis anyway
	if {[regexp {\\oOptArg\{([^\}]*?)\}\{([^\}]*?)\}} $line wholeexp m1 m2]} {
	    protectExpMatch wholeexp m1
	    regsub -all $wholeexp $line "<arg>$m1<replaceable>$m2</replaceable></arg><sbr/>" line
	}
	if {[regexp {\\oOpt\{([^\}]+)\}} $line wholeexp m1]} {
	    protectExpMatch wholeexp m1
	    regsub -all $wholeexp $line "<arg>$m1</arg>" line
	}

	if {[regexp {\\\\} $line]} {
	    if {$inpara} {
		set line [string map {\\\\ GNAGNU1\n</para>\n<para>} $line]
	    }
	}

	# special rules: mirdocs sometimes used <filename> and similar as userinput
	set line [string map {
	    <projectname> &lt;projectname&gt;
	    <extension> &lt;extension&gt;
	    <seqtype> &lt;seqtype&gt;
	    <filetype> &lt;filetype&gt;
	    <directoryname> &lt;directoryname&gt;
	    <strainname> &lt;strainname&gt;
	    "<empty string>" "&lt;empty string&gt;"
	    "<single character>" "&lt;single character&gt;"
	    <ssaha2options> &lt;ssaha2options&gt;
	    integer>0 "integer &gt; 0"
	    <=100 "&le;100"
	    0<=integer<=0 "0 &le; integer &le; 0"
	    \\\\ GNAGNU2\n</para>\n<para>
	    \\copyright &copy;
	    \\SP " "
	    <caf,fasta,gbf> &lt;caf,fasta,gbf&lt;
	    <nameOriginalStrain> &lt;nameOriginalStrain&lt;
	  } $line]
	if {[regexp {<filename>} $line]} {
	    if {![regexp {</filename>} $line]} {
		set line [string map {<filename> &lt;filename&gt;} $line]
	    }
	}
	if {![regexp {<name>[^<]+?</name>} $line]} {
	    set line [string map {<name> &lt;name&gt;} $line]
	}
	if {![regexp {<type>[^<]+?</type>} $line]} {
	    set line [string map {<type> &lt;type&gt;} $line]
	}

	if {$newlinepara} {
	    if {$inpara} {
		puts "</para>\n<para>"
	    }
	    set newlinepara 0
	}

	set line [string map {
	    "\\%" "%"
	    " \\& " " &amp; "
	    " <= " " &le; "
	    " < " " &lt; "
	    " >= " " &ge; "
	    " > " " &gt; "} $line]

	puts $line

    }

    closeInpara
    closeSect $seclevel 1

    close $fin
}


puts "<?xml version=\"1.0\" ?>
<!DOCTYPE book PUBLIC \"-//OASIS//DTD DocBook XML V4.5//EN\" \"http://www.docbook.org/xml/4.5/docbookx.dtd\">
<chapter>
<title>MIRA3 reference</title>
"

doit ../../solexadev_main.tex


puts "</chapter>
"