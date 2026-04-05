# Makefile for 2D Heat Equation Solver
# Builds parallel MPI program using OpenMPI compiler wrapper
#
# On macOS: ensure TMPDIR is set to /tmp to avoid shared memory errors:
#   export TMPDIR=/tmp
# or add to your .zshrc or .bashrc

# Default target: build the Heat executable
all: Heat

# Link object files to create the executable
Heat: HeatUtils.o Heat.o Interfaces.o
	mpic++ -Wall Heat.o HeatUtils.o Interfaces.o -o Heat

# Compile Heat.cpp
Heat.o: Heat.cpp
	mpic++ -Wall -g -c Heat.cpp

# Compile HeatUtils.cpp
HeatUtils.o: HeatUtils.cpp
	mpic++ -Wall -g -c HeatUtils.cpp

# Compile Interfaces.cpp
Interfaces.o: Interfaces.cpp
	mpic++ -Wall -g -c Interfaces.cpp

# Remove object files
clean:
	rm -f *.o

# Clean all generated files (object files, output files, executable)
mrproper: clean
	rm -f *.txt *~ Heat
