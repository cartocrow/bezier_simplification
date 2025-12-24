/* spiro-config.h.  Generated from spiro-config.h.in by configure.  */
/* spiro-config.h.in.  Generated from configure.ac by autoheader.  */

#ifndef _SPIRO_CONFIG_H
#define _SPIRO_CONFIG_H 1

/* Check if input values are realistic and not infinite in value. */
/* #undef CHECK_INPUT_FINITENESS */

/* Use 'gettimeofday()==true, else use older sys/timeb.h */
#define DO_TIME_DAY 1

/* Define to 1 if you have the <dlfcn.h> header file. */
/* #undef HAVE_DLFCN_H */

/* Define to 1 if you have the 'finite' function. */
/* #undef HAVE_FINITE */

/* Define to 1 if you have the 'hypot' function. */
#define HAVE_HYPOT 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the 'isfinite' function. */
#define HAVE_ISFINITE 1

/* Have pthreads.h. Do multi-user check in call-test. */
#define HAVE_PTHREADS 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdio.h> header file. */
#define HAVE_STDIO_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Libspiro version major value */
#define LS_VERSION_MJ 1

/* Libspiro version minor value */
#define LS_VERSION_MN 5

/* Libspiro Major.Minor value. Report this back */
#define LS_VERSION_STR "1.5"

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* Name of package */
#define PACKAGE "libspiro"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "fontforge-devel@lists.sourceforge.net"

/* Define to the full name of this package. */
#define PACKAGE_NAME "spiro"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "spiro 20240903"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "libspiro"

/* Define to the home page for this package. */
#define PACKAGE_URL "https://github.com/fontforge/libspiro"

/* Define to the version of this package. */
#define PACKAGE_VERSION "20240903"

/* Define to 1 if all of the C89 standard headers exist (not just the ones
   required in a freestanding environment). This macro is provided for
   backward compatibility; new code need not use it. */
#define STDC_HEADERS 1

/* Do 'S_TESTS' multi-thread 'call-testm' checks when you run 'make check'. */
#define S_TESTP 100

/* Verbose library printf() output enabled for debugging. */
/* #undef VERBOSE */

/* Version number of package */
#define VERSION "20240903"

/* Define IS_FINITE(x) to isfinite(x) or finite(x) */
#if HAVE_ISFINITE
#define IS_FINITE(x) isfinite(x)
#else
#if HAVE_FINITE
#define IS_FINITE(x) finite(x)
#endif
#endif

#endif
