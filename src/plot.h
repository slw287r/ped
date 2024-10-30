#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include <inttypes.h>
#include <cairo/cairo.h>
#include <cairo/cairo-svg.h>
#include "khashl.h"

/*
 * Standard Cell height is 15, which is equal to 20 pixels, and that's about 1/4 inch.
 * And standard Cell width is is 8.43 which is equal to 64 pixels and that's about 3/4 inch.
 */
// depth plots
#define MARGIN 35*3
#define DIM_X 605*2.5
#define DIM_Y 165*2.5
#define WIDTH (DIM_X + MARGIN)
#define HEIGHT (DIM_Y + MARGIN)

// map of unsigned long longs
KHASHL_MAP_INIT(KH_LOCAL, kh_t, kh, uint64_t, uint64_t, kh_hash_uint64, kh_eq_generic)

typedef struct
{
	int tid, pos, len, dep;
} dp_t;

void kh_ins(kh_t *h, uint64_t n, uint64_t v);
uint64_t kh_xval(kh_t *h, const uint64_t n);
void ld_pd(const char *fn, kh_t *pd);

//static void draw_rrect(cairo_t *cr);
void draw_arrow(cairo_t *cr, double start_x, double start_y, double end_x, double end_y);

/**
 * @brief Draw hist of depth
 *
 * @param cr  cairo sufface
 * @param pd  hash to pos-depth
 * @param gs  genome size
 * @param md  max depth
 */
void draw_hist(cairo_t *cr, kh_t *pd, const uint64_t gs, const int md);
