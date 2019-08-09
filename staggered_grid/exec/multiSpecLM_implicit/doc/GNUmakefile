EPStoPDF = epstopdf

ALL: DesignDocument.ps DesignDocument.pdf

eps_source = $(wildcard *.eps)

pdf_source = $(eps_source:.eps=.pdf)

DesignDocument.dvi: DesignDocument.tex $(eps_source)
	latex DesignDocument.tex < /dev/null
	latex DesignDocument.tex < /dev/null
	latex DesignDocument.tex < /dev/null
	latex DesignDocument.tex < /dev/null

DesignDocument.pdf: DesignDocument.tex $(pdf_source)
	pdflatex DesignDocument.tex < /dev/null
	pdflatex DesignDocument.tex < /dev/null
	pdflatex DesignDocument.tex < /dev/null
	pdflatex DesignDocument.tex < /dev/null

pdf:	DesignDocument.pdf

DesignDocument.ps: DesignDocument.dvi
	dvips -t letter -o $@ $<

%.pdf: %.eps
	$(EPStoPDF) $<

clean:
	$(RM) $(pdf_source) DesignDocument.dvi 
	$(RM) DesignDocument.blg DesignDocument.log DesignDocument.pdf
	$(RM) DesignDocument.aux DesignDocument.ps DesignDocument.bbl
	$(RM) DesignDocument.out *~

.PHONY: clean
