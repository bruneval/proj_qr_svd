#SCALAPACK= -m64  -L${MKLROOT}/lib -lmkl_scalapack_lp64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lmkl_blacs_openmpi_lp64 -lgomp -lpthread -lm -ldl
#SCALAPACK= -m64  -L/ccc/products/mkl-23.2.0/system/default/23.2.0//lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lpthread -lm
#SCALAPACK= ${SCALAPACK_LDFLAGS}
#SCALAPACK= -L/ccc/products/mkl-24.0.0/system/default/24.0.0//lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_sequential -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lpthread -lm \
	   -Wl,-rpath,/ccc/products/mkl-24.0.0/system/default/24.0.0//lib/intel64/
#SCALAPACK= -L/ccc/products/mkl-24.0.0/system/default/24.0.0//lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lpthread -lm \
	   -Wl,-rpath,/ccc/products/mkl-24.0.0/system/default/24.0.0//lib/intel64/
#SCA= -L/ccc/products/mkl-20.0.0/intel--20.0.0__openmpi--4.0.1/default/20.0.0/mkl/lib/intel64/ -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lpthread -lm

SCA=-L/home/spack/spack-test/opt/spack/linux-rocky8-skylake_avx512/gcc-13.1.0/netlib-scalapack-2.2.0-zrsvgssreu7x7llbji23iwmj3drbfa5m/lib/ -lscalapack  \
        -Wl,-rpath,/home/spack/spack-test/opt/spack/linux-rocky8-skylake_avx512/gcc-13.1.0/netlib-scalapack-2.2.0-zrsvgssreu7x7llbji23iwmj3drbfa5m/lib/ \
	-L/home/spack/spack-test/opt/spack/linux-rocky8-skylake_avx512/gcc-13.1.0/openblas-0.3.25-uoc3cz2ar6fhirx3udbxdoek6dwnilhh/lib/ -lopenblas \
	-Wl,-rpath,/home/spack/spack-test/opt/spack/linux-rocky8-skylake_avx512/gcc-13.1.0/openblas-0.3.25-uoc3cz2ar6fhirx3udbxdoek6dwnilhh/lib/

SCA= ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_openmpi_lp64.a -Wl,--end-group -lgomp -lpthread -lm -ldl

LAPACK= -m64  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl

LAPACK = ${SCA}


FC=mpifort -fcheck=all -fopenmp -cpp

all: proj_qr_svd_scalapack proj_qr_svd_lapack

proj_qr_svd_scalapack: low_level.o input_param.o proj_qr_svd_scalapack.o
	$(FC) -o proj_qr_svd_scalapack proj_qr_svd_scalapack.o low_level.o input_param.o $(SCA)

proj_qr_svd_scalapack.o: proj_qr_svd_scalapack.f90
	$(FC) -c proj_qr_svd_scalapack.f90

proj_qr_svd_lapack: low_level.o input_param.o proj_qr_svd_lapack.o
	$(FC) -o proj_qr_svd_lapack proj_qr_svd_lapack.o low_level.o input_param.o  $(LAPACK)

proj_qr_svd_lapack.o: low_level.f90 input_param.f90 proj_qr_svd_lapack.f90
	$(FC) -c proj_qr_svd_lapack.f90 

input_param.o: low_level.f90 input_param.f90
	$(FC) -c input_param.f90

low_level.o: low_level.f90
	$(FC) -c low_level.f90


clean:
	$(RM)  proj_qr_svd_scalapack proj_qr_svd_lapack *.o *.mod
