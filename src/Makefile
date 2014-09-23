CFLAGS = $(shell gsl-config --cflags) -I/home/Kevin/local/include
LDLIBS = $(shell gsl-config --libs-without-cblas) -L/home/Kevin/local/lib -lopenblas -lgfortran -lrt

executables = pca fastpca evec2weights project_weights

.PHONY: all debug clean

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CFLAGS += -g
endif

all: $(executables)

debug: CFLAGS += -g
debug: all

clean:
	rm -rf *.o $(executables)

pca: kjg_geno.o kjg_genoIO.o kjg_util.o kjg_GRM.o kjg_gsl.o

fastpca: kjg_geno.o kjg_genoIO.o kjg_gsl.o kjg_util.o kjg_fpca.o kjg_geno_gsl.o

fastpca-bed: kjg_geno.o kjg_bedIO.o kjg_gsl.o kjg_util.o kjg_fpca.o kjg_geno_gsl.o

evec2weights: kjg_geno.o kjg_gsl.o  kjg_genoIO.o kjg_fpca.o

project_weights: kjg_gsl.o kjg_geno.o kjg_genoIO.o kjg_fpca.o
