AC_PREREQ(2.59)

AC_INIT([Passive Film Model],[1.1.1],[aravindk@mit.edu],[pfm])

AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_PROG_FC([ftn mpif90 mpf90 mpxlf90 mpxlf95 mpxlf_r], [Fortran 90])

AC_ARG_WITH(spparks, [AS_HELP_STRING([--with-spparks=DIR], [Directory where SPPARKS was installed.])],[SPPARKS_LIB_DIR="$withval"],[])
AC_SUBST(SPPARKS_LIB_DIR)

AC_OUTPUT(Makefile)
