prefix = hrp1983
epsfigs = $(prefix)fig1.eps \
	  $(prefix)fig2.eps \
	  $(prefix)fig3.eps \
	  $(prefix)fig4.eps \
	  $(prefix)fig5.eps \
	  $(prefix)fig6.eps \
	  $(prefix)fig7.eps
pngfigs = $(patsubst %.eps,%.png, $(epsfigs))

default:
	./run.sh
	$(MAKE) png

png: $(pngfigs)

$(pngfigs) : %.png : %.eps Makefile
	convert -density 300 -resize 640x640 $< $@

clean:
	rm -rf *~ .*~ .*.un~ *.dat

deepclean: clean
	rm -rf *.eps *.png
