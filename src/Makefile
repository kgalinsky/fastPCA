CFLAGS=$(shell gsl-config --cflags)
LDLIBS=$(shell gsl-config --libs-without-cblas) -L/home/Kevin/local/lib -lopenblas -lgfortran

all: pca fastpca fastpca-iterative

pca: pca.o kjg_geno.o kjg_genoIO.o kjg_util.o ~/local/lib/libopenblas.a

fastpca: fastpca.o kjg_geno.o kjg_genoIO.o kjg_gsl.o kjg_util.o kjg_fpca.o

fastpca-iterative: fastpca-iterative.o kjg_geno.o kjg_genoIO.o kjg_gsl.o kjg_util.o kjg_fpca.o

clean:
	rm -rf *.o *.exe

dgemm_gsl_vs_blas: foo.o kjg_gsl.o
