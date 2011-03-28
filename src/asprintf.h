


#ifndef ASPRINTF_INC
#define ASPRINTF_INC

#ifdef __cplusplus
extern "C" {
#endif

#ifndef NEED_ASPRINTF
#include <stdio.h>
#else

#include <stdarg.h>

/* defined in gnulib/asprintf.c */
int asprintf (char **resultp, const char *format, ...);

/* defined in gnulib/vasprintf.c */
int vasprintf(char **strp, const char *fmt, va_list ap);

#endif


#ifdef __cplusplus
}
#endif

#endif



