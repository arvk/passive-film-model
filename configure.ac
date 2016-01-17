AC_PREREQ(2.59)

dnl INITIALIZE autoconf and automake
AC_INIT([Passive Film Model],[1.1.1],[aravindk@mit.edu],[pfm])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_YACC

dnl IDENTIFY FORTRAN COMPILERS
AC_PROG_FC([ftn mpif90 mpf90 mpxlf90 mpxlf95 mpxlf_r], [Fortran 90])


dnl TEST FOR SPPARKS - READ FROM STANDARD INPUT. THEN TEST FOR CUSTOM APPS
AC_ARG_WITH(spparks, [AS_HELP_STRING([--with-spparks=DIR], [Directory where SPPARKS was installed.])],[SPPARKS_LIB_DIR="$withval"],[])
AC_SUBST(SPPARKS_LIB_DIR)
AC_CHECK_FILE(["$withval/libspparks.so"],[AC_MSG_NOTICE([Found SPPARKS library])],[AC_MSG_ERROR([Didn't find SPPARKS at $SPPARKS_LIB_DIR])])
AC_CHECK_FILES([$withval/app_metfilm.cpp $withval/app_filmenv.cpp],[AC_MSG_NOTICE([Found custom SPPARKS apps])],[AC_MSG_ERROR([Didn't find custom SPPARKS apps at $SPPARKS_LIB_DIR])])



dnl TEST FOR PETSC - FIRST CHECK ENVIRONMENTAL VARIABLES, OR ELSE READ FROM STANDARD INPUT
if test x"$PETSC_DIR" != x""; then
  AC_MSG_NOTICE([PETSC directory is pre-defined. Using $PETSC_DIR])
else
  AC_ARG_WITH(petscdir, [AS_HELP_STRING([--with-petscdir=DIR], [PETSC base directory.])],[PETSC_DIR="$withval"],[AC_MSG_ERROR([PETSC directory is undefined])])
fi

if test x"$PETSC_ARCH" != x""; then
  AC_MSG_NOTICE([PETSC architecture is pre-defined. Using $PETSC_ARCH])
else
  AC_ARG_WITH(petscarch, [AS_HELP_STRING([--with-petscarch=DIR], [PETSC architecture])],[PETSC_ARCH="$withval"],[AC_MSG_ERROR([PETSC architecture is undefined])])
fi

AC_CHECK_FILE(["$PETSC_DIR/$PETSC_ARCH/lib/libpetsc.so"],[AC_MSG_NOTICE([Found PETSC library])],[AC_MSG_ERROR([Didn't find PETSC at $PETSC_DIR])])


dnl WRITE OUTPUT
AC_OUTPUT(Makefile)
