SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

# ---------------------------------------------------------------------
# Compiler selection
# ---------------------------------------------------------------------

CCC = g++ -O3

# ---------------------------------------------------------------------
# Compiler options
# ---------------------------------------------------------------------

CCOPT = -m64 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD -DLINUX -DCPLEX_MILP -DMILP -fopenmp

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------
CPLEXDIR      = ./CPLEX/CPLEXDir/cplex
CONCERTDIR    = ./CPLEX/CPLEXDir/concert
CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)
CCLNFLAGS = -lconcert -lilocplex -lcplex -lm -lpthread
CFLAGS    = -I./ -I$(CPLEXDIR)/include -I$(CONCERTDIR)/include -w -std=c++11 -fpermissive $(CCOPT)

obj = MPSchedulingFR.o MultiProcessorSystem.o SchedulabilityAnalysis.o RandomSystemsTest.o StatisticSet.o Util.o BitMask.o 

BitMask.o: BitMask.cpp BitMask.h
	$(CCC) $(CFLAGS) $(CCLNDIRS) -c  BitMask.cpp

MPSchedulingFR.o: MPSchedulingFR.cpp MPSchedulingFR.h
	$(CCC) $(CFLAGS) $(CCLNDIRS) -c  MPSchedulingFR.cpp

MultiProcessorSystem.o: MultiProcessorSystem.cpp
	$(CCC) $(CFLAGS) $(CCLNDIRS) -c  MultiProcessorSystem.cpp

RandomSystemsTest.o: RandomSystemsTest.cpp RandomSystemsTest.h
	$(CCC) $(CFLAGS) $(CCLNDIRS) -c  RandomSystemsTest.cpp

SchedulabilityAnalysis.o: SchedulabilityAnalysis.cpp SchedulabilityAnalysis.h
	$(CCC) $(CFLAGS) $(CCLNDIRS) -c  SchedulabilityAnalysis.cpp

StatisticSet.o: StatisticSet.cpp StatisticSet.h
	$(CCC) $(CFLAGS) $(CCLNDIRS) -c  StatisticSet.cpp

Util.o: Util.cpp Util.h
	$(CCC) $(CFLAGS) $(CCLNDIRS) -c  Util.cpp

main.o: main.cpp
	$(CCC) $(CFLAGS) $(CCLNDIRS) -c main.cpp

#unit_test.o: unit_test.cpp
#	$(CCC) $(CFLAGS) $(CCLNDIRS) -c unit_test.cpp

allObj: $(obj)

exec: $(obj) main.o
	$(CCC) $(CFLAGS) $(CCLNDIRS) -o GlobalFixedPriorityScheduling $(obj) main.o $(CCLNFLAGS)

clean :
	rm $(obj) main.o 