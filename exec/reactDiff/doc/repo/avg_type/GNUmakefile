EPStoPDF = epstopdf

ALL: repo1.pdf repo2.pdf

repo1.pdf: repo1.tex
	pdflatex repo1.tex < /dev/null
	pdflatex repo1.tex < /dev/null

repo2.pdf: repo2.tex
	pdflatex repo2.tex < /dev/null
	pdflatex repo2.tex < /dev/null

clean:
	$(RM) repo*.blg repo*.log repo*.pdf repo*.aux

.PHONY: clean
