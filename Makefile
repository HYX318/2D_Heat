# On Osx export TMPDIR=/tmp to avoid
# A system call failed during shared memory initialization that should
# not have.  It is likely that your MPI job will now either abort or
# experience performance degradation.

all : Heat


Heat : HeatUtils.o Heat.o Interfaces.o
	mpic++ -Wall Heat.o HeatUtils.o Interfaces.o -o Heat

Heat.o : Heat.cpp
	mpic++ -Wall -g -c Heat.cpp

HeatUtils.o : HeatUtils.cpp
	mpic++ -Wall -g -c HeatUtils.cpp

Interfaces.o : Interfaces.cpp
	mpic++ -Wall -g -c Interfaces.cpp

clean :
	rm -f *.o 

mrproper: clean
	rm -f *.txt *~ Heat
	
	

