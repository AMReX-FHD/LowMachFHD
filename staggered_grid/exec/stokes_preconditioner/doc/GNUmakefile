EPStoPDF = epstopdf

ALL: PreconditionerNotes.ps PreconditionerNotes.pdf

eps_source = $(wildcard *.eps)

pdf_source = $(eps_source:.eps=.pdf)

PreconditionerNotes.dvi: PreconditionerNotes.tex $(eps_source)
	latex PreconditionerNotes.tex < /dev/null
	latex PreconditionerNotes.tex < /dev/null
	latex PreconditionerNotes.tex < /dev/null
	latex PreconditionerNotes.tex < /dev/null

PreconditionerNotes.pdf: PreconditionerNotes.tex $(pdf_source)
	pdflatex PreconditionerNotes.tex < /dev/null
	pdflatex PreconditionerNotes.tex < /dev/null
	pdflatex PreconditionerNotes.tex < /dev/null
	pdflatex PreconditionerNotes.tex < /dev/null

pdf:	PreconditionerNotes.pdf

PreconditionerNotes.ps: PreconditionerNotes.dvi
	dvips -t letter -o $@ $<

%.pdf: %.eps
	$(EPStoPDF) $<

clean:
	$(RM) $(pdf_source) PreconditionerNotes.dvi 
	$(RM) PreconditionerNotes.blg PreconditionerNotes.log PreconditionerNotes.pdf
	$(RM) PreconditionerNotes.aux PreconditionerNotes.ps PreconditionerNotes.bbl
	$(RM) PreconditionerNotes.out *~

.PHONY: clean
