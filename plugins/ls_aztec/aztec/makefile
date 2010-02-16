
# compilers/archive programs and compilation flags
CC = gcc-4.3
FC = gfortran-4.3
AR = ar -ruvs

CFLAGS = -O3
FFLAGS = -O3
CFLAGS += -fPIC -Dappend_ -DAZ_SERIAL
FFLAGS += -fPIC

# Aztec/Netlib (BLAS, Y12M and Lapack) sources and objects/library
AZ_SRC = az_bilu.c              az_cg.c               az_cgs.c           \
         az_cgstab.c            az_check.c            az_comm.c          \
         az_converge.c          az_dd_overlap.c       az_dgemv2.c        \
         az_dgemv3.c            az_domain_decomp.c    az_fortran_wrap.c  \
         az_scaling.c           az_flop_cnt.c         az_gmres.c         \
         az_gmresr.c            az_ilu_util.c         az_ilut.c          \
         az_interface.c         az_lu_y12.c           az_matrix_util.c   \
         az_matvec_mult.c       az_old_matvec_mult.c  az_pad_utils.c     \
         az_poly.c              az_precond.c          az_qmrcgs.c        \
         az_rilu.c              az_solve.c            az_sort.c          \
         az_subdomain_solver.c  az_tools.c            az_util.c          \
         az_icc.c               az_fix_pt.c                              \
         md_wrap_scalar_c.c     md_timer_generic.c    # az_reorder.f

NETLIB = blas/blas1.f           blas/blas2.f          blas/blas3.f        \
         blas/xerbla.f          blas/lsame.f                              \
         y12m/y12m.f            y12m/y12mae.f         y12m/y12maf.f       \
         y12m/y12mbe.f          y12m/y12mbf.f         y12m/y12mce.f       \
         y12m/y12mcf.f          y12m/y12mde.f         y12m/y12mdf.f       \
         y12m/y12mfe.f          y12m/y12mge.f         y12m/y12mhe.f       \
         y12m/y12cck.f                                                    \
         lapack/la_dgetrf.f     lapack/la_dgetri.f    lapack/la_ilaenv.f  \
         lapack/la_dgetf2.f     lapack/la_dlaswp.f    lapack/la_dtrtri.f  \
         lapack/la_dtrti2.f     lapack/la_dgetrs.f    lapack/la_dgeco.f   \
         lapack/la_dgedi.f      lapack/la_dgefa.f     lapack/la_dpotf2.f  \
         lapack/la_dpotrf.f     lapack/la_dlamch.f    lapack/la_dlaic1.f  \
         lapack/la_ieeeck.f     lapack/la_iparmq.f

OBJ = $(AZ_SRC:.c=.o) az_reorder.o $(NETLIB:.f=.o)
LIB = libaztec.a

# targets
all : $(LIB)

$(LIB) : $(OBJ) makefile
	@echo "Building $@"
	@rm -f $@
	$(AR) $@ $(OBJ)

clean :
	@rm -f $(OBJ) $(LIB)

