AC_INIT([RASPA], [2.0])
AC_CONFIG_SRCDIR([src/main.c])
AM_INIT_AUTOMAKE

AC_CHECK_LIB(blas,ddot_,[LIBS="-lblas ${LIBS}"],,${LIBS})
AC_CHECK_LIB(lapack,dsyevd_,[LIBS="-llapack ${LIBS}"] HAVE_LAPACK="yes",,${LIBS})

if test "$HAVE_LAPACK" = yes; then
  AC_DEFINE(HAVE_LAPACK, 0, [Define if you have LAPACK library.])
fi

AC_PROG_CC
AC_PROG_CXX
AC_PROG_F77
AC_PROG_RANLIB

dnl Checks for header files.
AC_STDC_HEADERS
AC_HAVE_HEADERS(string.h)


AC_OUTPUT(Makefile \
src/Makefile \
forcefield/Makefile \
forcefield/GenericMOFs/Makefile \
forcefield/GarciaPerez2006/Makefile \
forcefield/Dubbeldam2007FlexibleIRMOF-1/Makefile \
forcefield/Dubbeldam2007FlexibleIRMOF-10/Makefile \
forcefield/Dubbeldam2007FlexibleIRMOF-16/Makefile \
forcefield/CoreShellCatlow/Makefile \
forcefield/CoreShellSchroderSauer/Makefile \
forcefield/CoreShellSastreGale2005/Makefile \
forcefield/RigidIonKramer1991/Makefile \
forcefield/Pascual2004/Makefile \
forcefield/AmirjalayerSchmid/Makefile \
forcefield/Nicholas/Makefile \
forcefield/Dubbeldam2012MIL-100/Makefile \
forcefield/Dubbeldam2012MIL-101/Makefile \
forcefield/McKendrick2012FlexibleMIL-53/Makefile \
forcefield/GarciaPerez2012/Makefile \
forcefield/Dzubak2012-IRMOF-1/Makefile \
forcefield/Dzubak2012-MgMOF-74/Makefile \
forcefield/Dzubak2012-MgMOF-74-extended/Makefile \
forcefield/Dzubak2012-ZnMOF-74/Makefile \
forcefield/MgMOF-74-Yazaydin/Makefile \
forcefield/CastilloVlugtCalero2009/Makefile \
forcefield/TraPPE/Makefile \
molecules/Makefile \
molecules/TraPPE/Makefile \
molecules/TraPPE-UA/Makefile \
molecules/CastilloVlugtCalero2009/Makefile \
framework/Makefile \
framework/Dubbeldam2007FlexibleIRMOF-1/Makefile \
framework/Dubbeldam2007FlexibleIRMOF-10/Makefile \
framework/Dubbeldam2007FlexibleIRMOF-16/Makefile \
framework/CoreShellCatlow/Makefile \
framework/CoreShellSchroderSauer/Makefile \
framework/CoreShellSastreGale2005/Makefile \
framework/NicholasModified/Makefile \
framework/AmirjalayerSchmid/Makefile \
framework/McKendrick2012FlexibleMIL-53/Makefile \
structures/Makefile \
structures/zeolites/Makefile \
structures/mofs/Makefile \
structures/minerals/Makefile \
structures/ceramics/Makefile \
structures/zeolites/cif/Makefile \
structures/zeolites/block/Makefile \
structures/zeolites/ions/Makefile \
structures/mofs/cif/Makefile  \
structures/mofs/block/Makefile \
structures/mofs/ions/Makefile \
structures/minerals/cif/Makefile  \
structures/ceramics/cif/Makefile  \
) 
