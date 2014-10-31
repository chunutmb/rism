# debugging
CC = gcc -O0 -Wall -Wunused -g
LM = -lm -lfftw3
# optimized compiling
# Note: -O3 makes it slower
CX = icc -O2 -Wall -Wremarks
LMX = -lfftw3
# alternative
# CX = gcc -Ofast -Wall -Wextra
# LMX = -lfftw3 -lm

bins = rism0
bins_d = $(patsubst %,%_d, $(bins))

all: $(bins)

rism0_deps = util.h models.h lmv.h mdiis.h uv.h Makefile

rism0 : rism0.c $(rism0_deps)
	$(CX) -o $@ $< $(LMX)

rism0_d : rism0.c $(rism0_deps)
	$(CX) -o $@ $< $(LMX)

clean:
	rm -f $(bins) $(bins_d) *~ a.out *.dat gmon.out .*.un~ MTSEED
	rm -rf r[0-9]*hs
	rstrip.py -Rlv

