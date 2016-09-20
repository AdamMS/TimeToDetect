
all: code.zip title-page.pdf Draft.pdf Supplement.pdf
	pdflatex Draft.tex
	pdflatex Draft.tex
	
Draft.pdf: Draft.tex Supplement.tex masterbib.bib
	pdflatex Draft.tex
	bibtex Draft
	pdflatex Draft.tex
	pdflatex Draft.tex
	
Supplement.pdf: Supplement.tex Draft.tex 
	pdflatex Supplement.tex
	pdflatex Supplement.tex
	
code.zip: 
	zip -r code.zip supplemental/
	
title-page.pdf:
	pdflatex title-page.tex
	