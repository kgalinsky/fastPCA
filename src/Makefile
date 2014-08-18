CFLAGS=$(shell gsl-config --cflags)
LDLIBS=$(shell gsl-config --libs-without-cblas) -lopenblas -lgfortran -lrt

executables = pca fastpca evec2weights project_weights

.PHONY: all debug clean

all: $(executables)

debug: CFLAGS+=-g
debug: $(executables)

clean:
	rm -rf *.o $(executables)

pca: pca.o kjg_2bit.o kjg_geno.o kjg_genoIO.o kjg_util.o kjg_GRM.o kjg_gsl.o

fastpca: fastpca.o kjg_2bit.o kjg_geno.o kjg_genoIO.o kjg_gsl.o kjg_util.o kjg_fpca.o

evec2weights: evec2weights.o kjg_gsl.o kjg_2bit.o kjg_geno.o kjg_genoIO.o kjg_fpca.o

project_weights: project_weights.o kjg_gsl.o kjg_2bit.o kjg_geno.o kjg_genoIO.o kjg_fpca.o
