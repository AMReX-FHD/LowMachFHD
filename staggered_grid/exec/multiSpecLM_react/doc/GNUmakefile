ALL: scheme_LM.pdf scheme_Bousq.pdf

scheme_LM.pdf: scheme_LM.tex
	pdflatex scheme_LM.tex
	pdflatex scheme_LM.tex

scheme_Bousq.pdf: scheme_Bousq.tex
	pdflatex scheme_Bousq.tex
	pdflatex scheme_Bousq.tex

clean:
	$(RM) scheme_*.aux scheme_*.log scheme*Notes.bib scheme_*.out scheme_*.pdf *~
