CC      ?= cc
CFLAGS  ?= -O3 -std=c11 -Wall -Wextra -pedantic

# Do not add -ffast-math: the prover's interval arithmetic assumes
# ordinary IEEE-754 semantics and nextafter padding.

euclid_prover: euclid_prover.c
	$(CC) $(CFLAGS) -o $@ $< -lm

clean:
	rm -f euclid_prover

.PHONY: clean
