ALL = \
	physics_bs_schedule.pdf

LATEX = xelatex

%.pdf: %.tex
	$(LATEX) $*; $(LATEX) $*
	rm -f $*.aux $*.log $*.out $*.tex~

all: $(ALL)

clean:
	rm -f .log *.dvi *.aux *.bbl *.blg

fullclean:
	rm -f *.pdf *.log *.dvi *.aux *.bbl *.blg

dummy:
