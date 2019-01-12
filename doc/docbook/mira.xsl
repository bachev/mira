<?xml version='1.0'?> 
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">

<xsl:import href="http://docbook.sourceforge.net/release/xsl/current/html/docbook.xsl"/>

<xsl:param name="html.stylesheet" select="'doccss/miradocstyle.css'"/>
<xsl:param name="admon.graphics" select="1"/>

<!-- no "Chapter X" written, but section numbering-->
<xsl:param name="chapter.autolabel" select="0"/>
<xsl:param name="section.autolabel" select="1"/>


<!-- toc level -->
<xsl:param name="toc.max.depth">4</xsl:param>
<xsl:param name="toc.section.depth">4</xsl:param>

<!-- that switches of chapter tocs, only book toc remains -->
<xsl:param name="generate.toc">
appendix  toc,title
article/appendix  nop
article   toc,title
book      toc,title,figure,table,example,equation
chapter   toc,title
part      toc,title
preface   toc,title
qandadiv  toc
qandaset  toc
reference toc,title
sect1     toc
sect2     toc
sect3     toc
sect4     toc
sect5     toc
section   toc
set       toc,title
</xsl:param>


<!-- book: only book toc
<xsl:template match="preface|chapter|appendix|article" mode="toc">
  <xsl:param name="toc-context" select="."/>

  <xsl:choose>
    <xsl:when test="local-name($toc-context) = 'book'">
      <xsl:call-template name="subtoc">
        <xsl:with-param name="toc-context" select="$toc-context"/>
        <xsl:with-param name="nodes" select="foo"/>
      </xsl:call-template>
    </xsl:when>
    <xsl:otherwise>
      <xsl:call-template name="subtoc">
        <xsl:with-param name="toc-context" select="$toc-context"/>
        <xsl:with-param name="nodes"
              select="section|sect1|glossary|bibliography|index
                     |bridgehead[$bridgehead.in.toc != 0]"/>
      </xsl:call-template>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>
-->






</xsl:stylesheet>

<!--
mira parameter type

<xsl:template match="mpt">
  <xsl:call-template name="inline.boldseq"/>
  <xsl:call-template name="inline.italicseq" >
    <xsl:with-param name="content"> = on|yes|1, off|no|0</xsl:with-param>
  </xsl:call-template> 
</xsl:template>

-->
