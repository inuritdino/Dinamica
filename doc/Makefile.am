main = dinamica
figures = *.eps
TEX = latex

$(main).dvi: $(figures) $(main).tex
	while ($(TEX) $(main) ; \
		grep -q "Rerun to get cross" $(main).log ) do true; \
		done

$(main).pdf: $(main).ps
	ps2pdf $(main).ps $(main).pdf

$(main).ps: $(main).dvi
	dvips $(main)

dist_doc_DATA = $(main).pdf $(figures) $(main).tex $(main).bib
