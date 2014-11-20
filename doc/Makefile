tex = $(wildcard *.tex)
pdf = $(patsubst %.tex, %.pdf, $(tex))
html = $(patsubst %.tex, %.html, $(tex))

all: $(pdf)

$(pdf) : %.pdf: %.tex
	pdflatex $<
	pdflatex $<

$(html) : %.html: %.tex
	tth -i $<

# -i: means to use italic font for equations

IPO=cz1@129.109.88.204
#IPO=cz1@he.host.utmb.edu

Bossman: $(pdf)
	rsync -avzrL *.pdf *.tex $(IPO):/Bossman/cz1/notes/singer1985mu

# upload to github
github:
	git push origin master

push: github

clean:
	rm -f *~ *.out *.dvi *.aux *.log *.idv *.lg .*.un~
	rstrip.py

.PHONY: clean all

