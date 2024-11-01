#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdbool.h>
#include <limits.h>
#include <unistd.h>
#include <inttypes.h>
#include <sys/ioctl.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/tbx.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>

#include "ketopt.h"
#include "cgranges.h"
#include "thpool.h"
#include "dplot.h"
#include "version.h"

extern const char *__progname;

#define VERSION "0.1.0"
#define CHUNK 0xFFFF
#define MAX_IS 0x400

#define ARR "\e[90m\xE2\x97\x82\e[0m"
#define BUL "\e[90m\xE2\x97\x8F\e[0m"
#define ERR "\e[31m\xE2\x9C\x97\e[0m"
#define PP fprintf(stderr, "%s\t%d\t<%s>\n", __FILE__, __LINE__, __func__);

#define error(msg, ...) do {                                 \
    fputs("\e[31m\xE2\x9C\x97\e[0m ", stderr);               \
    fprintf(stderr, msg, ##__VA_ARGS__);                     \
    if (errno) fprintf(stderr, ": %s", strerror(errno));     \
    fflush(stderr);                                          \
    exit(EXIT_FAILURE);                                      \
} while (0)

// argument struct
typedef struct
{
	char *in, *out, *plot, *sub, *ann, *dep, *ctg;
	int mis;
} arg_t;

static ko_longopt_t long_options[] = {
	{ "in",                        ko_required_argument, 'i' },
	{ "out",                       ko_required_argument, 'o' },
	{ "mis",                       ko_required_argument, 'm' },
	{ "plot",                      ko_required_argument, 'p' },
	{ "sub",                       ko_required_argument, 's' },
	{ "ctg",                       ko_required_argument, 'c' },
	{ "help",                      ko_no_argument, 'h' },
	{ "version",                   ko_no_argument, 'v' },
	{ NULL, 0, 0 }
};

void prs_arg(int argc, char **argv, arg_t *arg);
void ld_os(bam_hdr_t *hdr, int ci, kh_t *os, uint64_t *gl);
void ld_gr(samFile *fp, bam_hdr_t *h, int ci, int mis, cgranges_t *gr);
void ld_dp(const cgranges_t *cr, bam_hdr_t *hdr, int ci, dp_t **dp, uint32_t *md,
		uint64_t *nd);
void out_bed(const cgranges_t *gr, bam_hdr_t *h, int ci, const char *out);
void out_bgzf(const cgranges_t *gr, bam_hdr_t *h, int ci, const char *out);
void dump_dp(bam_hdr_t *hdr, dp_t *dp, uint64_t nd, const char *out);
void prep_an(const dp_t *dp, uint64_t nd, uint64_t gl, char *an);
void draw_canvas(cairo_surface_t *sf, cairo_t *cr, bam_hdr_t *hdr, int ci,
		const kh_t *os, const char *tt, const char *st, const char *an, uint32_t md,
		uint64_t gl);
void draw_ped1(cairo_t *cr, kh_t *os, uint32_t md, uint64_t gl, dp_t *dp);
void draw_axis(cairo_t *cr, uint32_t md, uint64_t gl);
int is_gzip(const char *fn);
bool ends_with(const char *str, const char *sfx);
int strlen_wo_esc(const char *str);
void horiz_em(const int _n);
void horiz(const int _n);
static void usage();
