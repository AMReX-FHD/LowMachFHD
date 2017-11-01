EPStoPDF = epstopdf

ALL: Electroosmosis.ps Electroosmosis.pdf

eps_source = $(wildcard *.eps)

pdf_source = $(eps_source:.eps=.pdf)

Electroosmosis.dvi: Electroosmosis.tex Electroosmosis.bib $(eps_source)
	latex Electroosmosis.tex < /dev/null
	bibtex Electroosmosis < /dev/null
	latex Electroosmosis.tex < /dev/null
	latex Electroosmosis.tex < /dev/null
	latex Electroosmosis.tex < /dev/null

Electroosmosis.pdf: Electroosmosis.tex Electroosmosis.bib $(pdf_source)
	pdflatex Electroosmosis.tex < /dev/null
	bibtex Electroosmosis < /dev/null
	pdflatex Electroosmosis.tex < /dev/null
	pdflatex Electroosmosis.tex < /dev/null
	pdflatex Electroosmosis.tex < /dev/null

pdf:	Electroosmosis.pdf

Electroosmosis.ps: Electroosmosis.dvi
	dvips -t letter -o $@ $<

%.pdf: %.eps
	$(EPStoPDF) $<

clean:
	$(RM) $(pdf_source) Electroosmosis.dvi 
	$(RM) Electroosmosis.blg Electroosmosis.log Electroosmosis.pdf
	$(RM) Electroosmosis.aux Electroosmosis.ps Electroosmosis.bbl
	$(RM) Electroosmosis.out *~

.PHONY: clean
