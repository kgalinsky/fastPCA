CFLAGS=$(shell gsl-config --cflags)
LDLIBS=$(shell gsl-config --libs-without-cblas) -lopenblas -lgfortran

executables = pca fastpca fastpca-iterative

all: $(executables)

debug: CFLAGS+=-g
debug: $(executables)

pca: pca.o kjg_geno.o kjg_genoIO.o kjg_util.o kjg_GRM.o

fastpca: fastpca.o kjg_geno.o kjg_genoIO.o kjg_gsl.o kjg_util.o kjg_fpca.o

fastpca-iterative: fastpca-iterative.o kjg_geno.o kjg_genoIO.o kjg_gsl.o kjg_util.o kjg_fpca.o

clean:
	rm -rf *.o *.exe $(executables)

dgemm_gsl_vs_blas: dgemm_gsl_vs_blas.o kjg_gsl.o

submatrix: submatrix.o kjg_gsl.o
