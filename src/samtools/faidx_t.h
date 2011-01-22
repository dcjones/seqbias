/* 
 * I moved some definitions from faidx.c to here so I can manipulated the
 * faidx_t structure.
 *
 * (Daniel Jones / 2011.01.05)
 */

#ifndef __AC_FAIDX_T_H
#define __AC_FAIDX_T_H

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "faidx.h"
#include "khash.h"



typedef struct {
	uint64_t len:32, line_len:16, line_blen:16;
	uint64_t offset;
} faidx1_t;
KHASH_MAP_INIT_STR(s, faidx1_t)

#ifndef _NO_RAZF
#include "razf.h"
#else
#ifdef _WIN32
#define ftello(fp) ftell(fp)
#define fseeko(fp, offset, whence) fseek(fp, offset, whence)
#else
extern off_t ftello(FILE *stream);
extern int fseeko(FILE *stream, off_t offset, int whence);
#endif
#define RAZF FILE
#define razf_read(fp, buf, size) fread(buf, 1, size, fp)
#define razf_open(fn, mode) fopen(fn, mode)
#define razf_close(fp) fclose(fp)
#define razf_seek(fp, offset, whence) fseeko(fp, offset, whence)
#define razf_tell(fp) ftello(fp)
#endif
#ifdef _USE_KNETFILE
#include "knetfile.h"
#endif

struct __faidx_t {
	RAZF *rz;
	int n, m;
	char **name;
	khash_t(s) *hash;
};



#endif


