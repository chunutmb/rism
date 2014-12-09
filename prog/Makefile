# debugging
CC = gcc -O0 -Wall -Wunused -Wextra -Wno-missing-field-initializers -g
LM = -lm

# optimized compiling
CX = gcc -O3 -Wall -Wextra -Wno-missing-field-initializers
LMX = -lm

# alternative using intel C compiler
# Note: -O3 makes it slower
#CX = icc -O2 -Wall -Wremarks
#LMX =

bins = rism0
bins_d = $(patsubst %, %_d, $(bins))
bins_nofftw = $(patsubst %, %_nofftw, $(bins))
deps = $(wildcard *.h) Makefile

all: $(bins)

nofftw: $(bins_nofftw)

$(bins): % : %.c $(deps)
	$(CX) -o $@ $< -lfftw3 $(LMX)

$(bins_d): %_d : %.c $(deps)
	$(CC) -o $@ $< -lfftw3 $(LM)

$(bins_nofftw): %_nofftw : %.c $(deps)
	$(CX) -DNOFFTW -o $@ $< $(LMX)

clean:
	rm -rf $(bins) $(bins_d) $(bins_nofftw) *~ a.out *.dat* gmon.out .*.un~ */*~ */.*.un~ .#* MTSEED
	rm -rf r[0-9]*hs
	rstrip.py -Rlv

