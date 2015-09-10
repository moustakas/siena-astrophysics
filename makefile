ALL = \
	dayone_handout.pdf \
	physics_bs_schedule.pdf \
	template_schedule_conceptual_astro.pdf
#	template_schedule_conceptual.pdf \

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
