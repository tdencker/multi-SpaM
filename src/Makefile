CXX=g++
CXXFLAGS= -O3 -g -Wall -std=c++11 -march=native -fopenmp
CC=gcc
avx:=$(shell if [ `grep -c "avx2" /proc/cpuinfo` -ne 0 ] ; then echo "-D__AVX"; fi)
avx_obj:=$(shell if [ `grep -c "avx2" /proc/cpuinfo` -ne 0 ] ; then echo "raxml/avxLikelihood.o"; fi)
CCFLAGS= -O3 -march=native -D__SIM_SSE3 -D_GNU_SOURCE -D_FMA -lm $(avx)

DEPS = word.hpp sequence.hpp match.hpp pattern.hpp pseudoalignment.hpp pseudoaligner.hpp \
	util.hpp rhasbari/patternset.h rhasbari/variance.h options.hpp paindex.hpp stats.hpp randommatchfinder.hpp

OBJS = multi-SpaM.o word.o options.o sequence.o quartetblock.o \
component.o  stats.o randommatchfinder.o pattern.o\
raxml/axml_extract.o raxml/mem_alloc.o raxml/models.o raxml/evaluateGenericSpecial.o raxml/searchAlgo.o \
raxml/newviewGenericSpecial.o raxml/makenewzGenericSpecial.o raxml/evaluatePartialGenericSpecial.o \
raxml/optimizeModel.o $(avx_obj) rasbhari/patternset.o rasbhari/variance.o

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
