EPStoPDF = epstopdf

ALL: repo1.pdf repo2.pdf repo3.pdf repo4.pdf repo5.pdf

repo1.pdf: repo1.tex
	pdflatex repo1.tex < /dev/null
	pdflatex repo1.tex < /dev/null

repo2.pdf: repo2.tex
	pdflatex repo2.tex < /dev/null
	pdflatex repo2.tex < /dev/null

repo3.pdf: repo3.tex
	pdflatex repo3.tex < /dev/null
	pdflatex repo3.tex < /dev/null

repo4.pdf: repo4.tex
	pdflatex repo4.tex < /dev/null
	pdflatex repo4.tex < /dev/null

repo5.pdf: repo5.tex
	pdflatex repo5.tex < /dev/null
	pdflatex repo5.tex < /dev/null

clean:
	$(RM) repo*.blg repo*.log repo*.pdf repo*.aux fig*/*-eps-converted-to.pdf

.PHONY: clean
