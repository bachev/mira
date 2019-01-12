Documentation is provided

1. As full manual in PDF and HTML format in the "doc" directory
   - For the HTML documentation: xsltproc needs to be functional
   - For the PDF documentation: xsltproc & dblatex need to be functional
2. As barebone man pages in the "man" directory

The man pages are mostly meant as simple help for the small utilities
(miraconvert and mirabait). They are located one directory up to make it
easier to just compile the MIRA binaries and install them alongside with just
the man pages as the PDF and HTML documentation has pretty heavy dependencies
(dblatex etc.).
