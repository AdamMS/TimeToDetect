
all:
	pdflatex Draft.tex
	bibtex Draft
	bibtex Draft
	pdflatex Supplement.tex
	pdflatex Supplement.tex
	pdflatex Draft.tex
	pdflatex Draft.tex