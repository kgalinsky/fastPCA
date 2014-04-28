CFLAGS=$(shell gsl-config --cflags)
LDLIBS=$(shell gsl-config --libs-without-cblas) -L/home/Kevin/local/lib -lopenblas

all: pca fastpca fastpca-iterative

pca: pca.o kjg_geno.o kjg_genoIO.o kjg_gsl.o kjg_util.o

fastpca: fastpca.o kjg_geno.o kjg_genoIO.o kjg_gsl.o kjg_util.o kjg_fpca.o

fastpca-iterative: fastpca-iterative.o kjg_geno.o kjg_genoIO.o kjg_gsl.o kjg_util.o kjg_fpca.o

clean:
	rm -rf *.o *.exe

dgemm_gsl_vs_blas: foo.o kjg_gsl.o
