FC = gfortran
FFLAGS = -O2 -Wall -Wextra -std=f2008
BIN = shocktube
SRC = src/shock_tube_solver.f90

all: $(BIN)

$(BIN): $(SRC)
	$(FC) $(FFLAGS) $(SRC) -o $(BIN)

run: $(BIN)
	mkdir -p results
	./$(BIN) godunov 201 0.25 results/godunov_n201.csv
	./$(BIN) muscl   201 0.25 results/muscl_n201.csv
	./$(BIN) hybrid  201 0.25 results/hybrid_n201.csv
	./$(BIN) godunov 1001 0.25 results/godunov_n1001.csv
	./$(BIN) muscl   1001 0.25 results/muscl_n1001.csv
	./$(BIN) hybrid  1001 0.25 results/hybrid_n1001.csv

plots:
	python3 scripts/plot_comparison.py

clean:
	rm -f $(BIN) *.o *.mod

clean-results:
	rm -f results/*.csv figures/*.png
