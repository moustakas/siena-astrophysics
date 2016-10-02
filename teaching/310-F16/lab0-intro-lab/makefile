ALL=	example-lab-report.pdf

LATEX =	pdflatex
BIBTEX =bibtex

%.pdf: %.tex
	$(LATEX) $*
	$(LATEX) $*
	$(BIBTEX) $*
	$(LATEX) $*
	$(LATEX) $*
	rm -f *.log *.dvi *.aux *.bbl *.blg *.ps

all: $(ALL)

clean:
	rm -f *.log *.dvi *.aux *.bbl *.blg *.ps

dummy:

