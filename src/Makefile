CFLAGS=$(shell gsl-config --cflags)
LDLIBS=$(shell gsl-config --libs)

all: pca fastpca fastpca-iterative

pca: pca.o kjg_geno.o kjg_genoIO.o kjg_gsl.o kjg_util.o

fastpca: fastpca.o kjg_geno.o kjg_genoIO.o kjg_gsl.o kjg_util.o

fastpca-iterative: fastpca-iterative.o kjg_geno.o kjg_genoIO.o kjg_gsl.o kjg_util.o

clean:
	rm -rf *.o *.exe
