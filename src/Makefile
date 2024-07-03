#SCALAPACK= -m64  -L${MKLROOT}/lib -lmkl_scalapack_lp64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lmkl_blacs_openmpi_lp64 -lgomp -lpthread -lm -ldl
#SCALAPACK= -m64  -L/ccc/products/mkl-23.2.0/system/default/23.2.0//lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lpthread -lm
#SCALAPACK= ${SCALAPACK_LDFLAGS}
#SCALAPACK= -L/ccc/products/mkl-24.0.0/system/default/24.0.0//lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_sequential -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lpthread -lm \
	   -Wl,-rpath,/ccc/products/mkl-24.0.0/system/default/24.0.0//lib/intel64/
#SCALAPACK= -L/ccc/products/mkl-24.0.0/system/default/24.0.0//lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lpthread -lm \
	   -Wl,-rpath,/ccc/products/mkl-24.0.0/system/default/24.0.0//lib/intel64/
#MKLROOT=/ccc/products/mkl-24.0.0/system/default/24.0.0/
#MKLROOT=/ccc/products/mkl-20.0.0/intel--20.0.0__openmpi--4.0.1/default/20.0.0/mkl/lib/intel64/
#MKLROOT=/ccc/products/mkl-20.0.0/gcc--11.1.0__openmpi--4.0.1/default/20.0.0/mkl/
#SCA= -L/ccc/products/mkl-20.0.0/intel--20.0.0__openmpi--4.0.1/default/20.0.0/mkl/lib/intel64/ -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lpthread -lm

SCA=-L/home/spack/spack-test/opt/spack/linux-rocky8-skylake_avx512/gcc-13.1.0/netlib-scalapack-2.2.0-zrsvgssreu7x7llbji23iwmj3drbfa5m/lib/ -lscalapack  \
        -Wl,-rpath,/home/spack/spack-test/opt/spack/linux-rocky8-skylake_avx512/gcc-13.1.0/netlib-scalapack-2.2.0-zrsvgssreu7x7llbji23iwmj3drbfa5m/lib/ \
	-L/home/spack/spack-test/opt/spack/linux-rocky8-skylake_avx512/gcc-13.1.0/openblas-0.3.25-uoc3cz2ar6fhirx3udbxdoek6dwnilhh/lib/ -lopenblas \
	-Wl,-rpath,/home/spack/spack-test/opt/spack/linux-rocky8-skylake_avx512/gcc-13.1.0/openblas-0.3.25-uoc3cz2ar6fhirx3udbxdoek6dwnilhh/lib/

SCA= ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_openmpi_lp64.a -Wl,--end-group -lgomp -lpthread -lm -ldl



FC=mpifort -fcheck=all -fopenmp -cpp

all: proj_qr_svd_scalapack

proj_qr_svd_scalapack: proj_qr_svd_scalapack.f90
	$(FC) -o proj_qr_svd_scalapack proj_qr_svd_scalapack.f90 $(SCA)

proj_qr_svd_fortran: proj_qr_svd_fortran.f90
	$(FC) -o proj_qr_svd_fortran proj_qr_svd_fortran.f90 $(LAPACK)

clean:
	$(RM)  proj_qr_svd_scalapack proj_qr_svd_fortran 
