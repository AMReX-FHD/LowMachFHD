EPStoPDF = epstopdf

ALL: Stokes_Preconditioners.ps Stokes_Preconditioners.pdf

eps_source = $(wildcard *.eps)

pdf_source = $(eps_source:.eps=.pdf)

Stokes_Preconditioners.dvi: Stokes_Preconditioners.tex $(eps_source)
	latex Stokes_Preconditioners.tex < /dev/null
	latex Stokes_Preconditioners.tex < /dev/null
	latex Stokes_Preconditioners.tex < /dev/null
	latex Stokes_Preconditioners.tex < /dev/null

Stokes_Preconditioners.pdf: Stokes_Preconditioners.tex $(pdf_source)
	pdflatex Stokes_Preconditioners.tex < /dev/null
	pdflatex Stokes_Preconditioners.tex < /dev/null
	pdflatex Stokes_Preconditioners.tex < /dev/null
	pdflatex Stokes_Preconditioners.tex < /dev/null

pdf:	Stokes_Preconditioners.pdf

Stokes_Preconditioners.ps: Stokes_Preconditioners.dvi
	dvips -t letter -o $@ $<

%.pdf: %.eps
	$(EPStoPDF) $<

clean:
	$(RM) $(pdf_source) Stokes_Preconditioners.dvi 
	$(RM) Stokes_Preconditioners.blg Stokes_Preconditioners.log Stokes_Preconditioners.pdf
	$(RM) Stokes_Preconditioners.aux Stokes_Preconditioners.ps Stokes_Preconditioners.bbl
	$(RM) Stokes_Preconditioners.out *~

.PHONY: clean
