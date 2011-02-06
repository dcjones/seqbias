


#ifndef ASPRINTF_INC
#define ASPRINTF_INC

#ifdef __cplusplus
extern "C" {
#endif

#ifndef NEED_ASPRINTF
#include <stdio.h>
#else
int asprintf (char **resultp, const char *format, ...);
/* defined in gnulib/asprintf.c */
#endif


#ifdef __cplusplus
}
#endif

#endif



