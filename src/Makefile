CXX=g++
CXXFLAGS= -O3 -g -Wall -std=c++11 -march=native -fopenmp
CC=gcc
avx:=$(shell if [ `grep -c "avx2" /proc/cpuinfo` -ne 0 ] ; then echo "-D__AVX"; fi)
CCFLAGS= -O3 -march=native -D__SIM_SSE3 -D_GNU_SOURCE -D_FMA -lm $(avx)
#ifeq($(grep -c "avx2" /proc/cpuinfo), 0)
#CCFLAGS= -O3 -march=native -D__SIM_SSE3 -D_GNU_SOURCE -D_FMA -lm
#else
#CCFLAGS= -O3 -march=native -D__AVX -D__SIM_SSE3 -D_GNU_SOURCE -D_FMA -lm
#endif

#var:=$(shell grep -c "avx2" /proc/cpuinfo )
#ifeq($$(grep -c "avx2" /proc/cpuinfo), 0)
#CCFLAGS = $(CCFLAGS) -D__AVX
#endif

DEPS = word.hpp sequence.hpp match.hpp pattern.hpp pseudoalignment.hpp pseudoaligner.hpp \
	util.hpp rhasbari/patternset.h rhasbari/variance.h options.hpp paindex.hpp stats.hpp randommatchfinder.hpp

OBJS = multi-SpaM.o word.o options.o sequence.o pseudoalignment.o \
component.o  stats.o randommatchfinder.o\
raxml/axml_extract.o raxml/mem_alloc.o raxml/models.o raxml/evaluateGenericSpecial.o raxml/searchAlgo.o \
raxml/newviewGenericSpecial.o raxml/makenewzGenericSpecial.o raxml/evaluatePartialGenericSpecial.o \
raxml/optimizeModel.o raxml/avxLikelihood.o rasbhari/patternset.o rasbhari/variance.o

all: multi-SpaM

# raxml
%.o : %.c
	$(CC) $(CCFLAGS) -c -g -o $@ $<

%.o : %.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

multi-SpaM: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

clean:  
	find . -name "*.o" -exec rm {} ';'
