EPStoPDF = epstopdf

ALL: Applications.ps Applications.pdf

eps_source = $(wildcard *.eps)

pdf_source = $(eps_source:.eps=.pdf)

Applications.dvi: Applications.tex Applications.bib $(eps_source)
	latex Applications.tex < /dev/null
	bibtex Applications < /dev/null
	latex Applications.tex < /dev/null
	latex Applications.tex < /dev/null
	latex Applications.tex < /dev/null

Applications.pdf: Applications.tex Applications.bib $(pdf_source)
	pdflatex Applications.tex < /dev/null
	bibtex Applications < /dev/null
	pdflatex Applications.tex < /dev/null
	pdflatex Applications.tex < /dev/null
	pdflatex Applications.tex < /dev/null

pdf:	Applications.pdf

Applications.ps: Applications.dvi
	dvips -t letter -o $@ $<

%.pdf: %.eps
	$(EPStoPDF) $<

clean:
	$(RM) $(pdf_source) Applications.dvi 
	$(RM) Applications.blg Applications.log Applications.pdf
	$(RM) Applications.aux Applications.ps Applications.bbl
	$(RM) Applications.out *~

.PHONY: clean
