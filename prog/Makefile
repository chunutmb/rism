# debugging
CC = gcc -O0 -Wall -Wunused -g
LM = -lm -lfftw3
# optimized compiling
CX = icc -O3 -Wall -Wremarks
LMX = -lfftw3

bins = rism0
bins_d = $(patsubst %,%_d, $(bins))

all: $(bins)

rism0_deps = util.h models.h Makefile

rism0 : rism0.c $(rism0_deps)
	$(CX) -o $@ $< $(LMX)

rism0_d : rism0.c $(rism0_deps)
	$(CX) -o $@ $< $(LMX)
clean:
	rm -f $(bins) $(bins_d) *~ a.out *.dat MTSEED
	rm -rf r[0-9]*hs

