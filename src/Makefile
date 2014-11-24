CFLAGS = -o2 $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs-without-cblas) -lopenblas -lgfortran -lrt

executables = pca fastpca evec2weights project_weights

.PHONY: all debug clean

all: $(executables)

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CFLAGS += -g
endif

debug: CFLAGS += -g
debug: all

clean:
	rm -rf *.o $(executables)

pca: kjg_geno.o kjg_genoIO.o kjg_util.o kjg_GRM.o kjg_gsl.o

fastpca: kjg_geno.o kjg_genoIO.o kjg_bedIO.o kjg_gsl.o kjg_util.o kjg_fpca.o kjg_geno_gsl.o glue.o

evec2weights: kjg_geno.o kjg_gsl.o kjg_genoIO.o kjg_geno_gsl.o glue.o

evec2weights-bed: kjg_geno.o kjg_gsl.o kjg_bedIO.o kjg_geno_gsl.o

project_weights: kjg_gsl.o kjg_geno.o kjg_genoIO.o kjg_fpca.o

sim_star: kjg_geno.o kjg_geno_rand.o
