# debugging
CC = gcc -O0 -Wall -Wunused -g
LM = -lm -lfftw3
# optimized compiling
CX = icc -O3 -Wall -Wremarks
LMX = -lfftw3

deps = util.h Makefile
bins = rism0
bins_d = $(patsubst %,%_d, $(bins))

all: $(bins)

$(bins) : %: %.c $(deps)
	$(CX) -o $@ $< $(LMX)

$(bins_d) : %_d: %.c $(deps)
	$(CC) -o $@ $< $(LM)

clean:
	rm -f $(bins) $(bins_d) *~ a.out *.dat MTSEED
	rm -rf r[0-9]*hs

