CFLAGS=$(shell gsl-config --cflags)
LDLIBS=$(shell gsl-config --libs-without-cblas) -lopenblas -lgfortran -lrt

executables = pca fastpca fastpca-iterative evec2weights

all: $(executables)

debug: CFLAGS+=-g
debug: $(executables)

pca: pca.o kjg_geno.o kjg_genoIO.o kjg_util.o kjg_GRM.o

fastpca: fastpca.o kjg_geno.o kjg_genoIO.o kjg_gsl.o kjg_util.o kjg_fpca.o

fastpca-iterative: fastpca-iterative.o kjg_geno.o kjg_genoIO.o kjg_gsl.o kjg_util.o kjg_fpca.o

evec2weights: evec2weights.o kjg_gsl.o kjg_geno.o kjg_genoIO.o kjg_fpca.o

clean:
	rm -rf *.o $(executables)
