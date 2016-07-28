TEX = $(wildcard *.tex)
PDFTEX = $(wildcard *.pdftex)
PSTEX = $(wildcard *.pstex)
PSTEX_T = $(wildcard *.pstex_t)
EPS = $(wildcard *.eps)
DVI = $(patsubst %.tex,%.dvi, $(TEX))
LOG = $(patsubst %.tex,%.log, $(TEX))
AUX = $(patsubst %.tex,%.aux, $(TEX))
TOC = $(patsubst %.tex,%.toc, $(TEX))
OUT = $(patsubst %.tex,%.out, $(TEX))
PS  = $(patsubst %.dvi,%.ps, $(DVI))
PDF = $(patsubst %.ps,%.pdf, $(PS))

%.dvi: %.tex $(EPS) $(PSTEX) $(PSTEX_T)
	( \
	\latex $<; \
	while \grep -q "Rerun to get cross-references right." $(<:.tex=.log); \
	do \
		\latex $<; \
	done \
	)

%.ps: %.dvi ; dvips $< -o $@

%.pdf: %.dvi ; dvi2pdf_Perso $<

all: $(DVI) $(PS) $(PDF)

allpdf: $(DVI) $(PS) $(PDF)

# reason for *.aux is that if a FOO.tex file is inexistent but in
# \input{FOO.tex} then a FOO.aux is created and is not removed if we put
# $(AUX) in the target...

clean:
	rm -f $(PS) $(DVI) $(PDF) $(LOG) *.aux $(TOC) $(OUT)
