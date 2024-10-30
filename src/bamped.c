#include "bamped.h"

int main(int argc, char *argv[])
{
	if (argc == 1)
		usage();
	setenv("FONTCONFIG_PATH", "/etc/fonts", 1);
	arg_t *arg = calloc(1, sizeof(arg_t));
	arg->mis = MAX_IS;
	prs_arg(argc, argv, arg);
	cgranges_t *gr = cr_init();
	samFile *fp = sam_open(arg->in, "r");
	if (!fp)
		error("Error: failed to read input bam [%s]\n", arg->in);
	bam_hdr_t *hdr = sam_hdr_read(fp);
	// multiple ctg offsets
	int i;
	khint_t k;
	uint64_t gl = 0; // genome length in total
	kh_t *os = kh_init();
	ld_os(hdr, os, &gl);
	
	/*
	cairo_save(cr);
	cairo_scale(cr, DIM_X, DIM_Y);
	draw_hist(cr, pd, op->len, md);
	cairo_restore(cr);
	*/
	/* dbg os
	for (i = 0; i < hdr->n_targets; ++i)
		printf("%s\t%d\t%"PRIu64"\n", hdr->target_name[i], hdr->target_len[i], kh_xval(os, i));
	printf("%"PRIu64"\n", gl);
	*/
	// PE regions
	ld_gr(fp, hdr, arg->mis, gr);
	// TODO save dp to array and get max dep
	if (arg->out && ends_with(arg->out, ".gz"))
		out_bgzf(gr, hdr, arg->out);
	else
		out_bed(cg, hdr, arg->out);
	// dep plot
	char *tt = NULL, *st = NULL, *an = NULL;
	cairo_t *cr = NULL;
	cairo_surface_t *sf = NULL;
	if (arg->plot)
		draw_canvas(sf, cr, hdr, tt, st, an, md, gl);
	cr_destroy(gr);
	// TODO iterate dp array and draw to cr
	bam_hdr_destroy(hdr);
	hts_close(fp);
	if (arg->plot)
	{
		cairo_surface_write_to_png(cairo_get_target(cr), arg->plot);
		cairo_surface_destroy(sf);
		cairo_destroy(cr);
	}
	free(arg);
	return 0;
}

void ld_os(bam_hdr_t *hdr, kh_t *os, uint64_t *gl)
{
	int i, shift = 0;
	for (i = 0; i < hdr->n_targets; ++i)
	{
		kh_ins(os, i, shift);
		shift += hdr->target_len[i];
	}
	*gl = shift;
}

void ld_cr(samFile *fp, bam_hdr_t *h, int mis, cgranges_t *cr)
{
	bam1_t *b = bam_init1();
	while (sam_read1(fp, h, b) >= 0)
	{
		bam1_core_t *c = &b->core;
		if (c->flag & (BAM_FQCFAIL | BAM_FUNMAP | BAM_FREVERSE | BAM_FSECONDARY | BAM_FSUPPLEMENTARY))
			continue;
		if ((c->flag & BAM_FPAIRED) && (c->tid == c->mtid) && c->isize <= mis)
			cr_add(cr, h->target_name[c->tid], c->pos, c->pos + c->isize, 0);
	}
	bam_destroy1(b);
	cr_sort(cr);
	cr_index(cr);
}

void draw_canvas(cairo_surface_t *sf, cairo_t *cr, bam_hdr_t *hdr, const char *tt,
		const char *st, const char *an, int md, int gl)
{
	sf = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, WIDTH, HEIGHT);
	cv = cairo_create(sf);
	cairo_set_antialias(cr, CAIRO_ANTIALIAS_BEST);
	draw_rrect(cr);
	cairo_set_source_rgb (cr, 0, 0, 0);
	cairo_translate(cr, MARGIN / 1.25, MARGIN / 2.0);
	// axis labels
	double x, y;
	cairo_text_extents_t ext;
	//cairo_set_source_rgb(cr, 0.25, 0.25, 0.25);
	if (tt) // title
	{
		cairo_set_font_size(cr, 24.0);
		cairo_select_font_face(cr, "serif", CAIRO_FONT_SLANT_ITALIC, CAIRO_FONT_WEIGHT_NORMAL);
		cairo_text_extents(cr, tt, &ext);
		x = DIM_X / 2.0 - (ext.width / 2.0 + ext.x_bearing);
		y = ext.height / 2 + ext.y_bearing;
		cairo_move_to(cr, x, y);
		cairo_show_text(cr, tt);
	}
	if (st) // sub-title
	{
		cairo_set_font_size(cr, 24.0);
		cairo_select_font_face(cr, "serif", CAIRO_FONT_SLANT_ITALIC, CAIRO_FONT_WEIGHT_NORMAL);
		cairo_text_extents(cr, st, &ext);
		x = DIM_X / 2.0 - (ext.width / 2.0 + ext.x_bearing);
		y = ext.height + ext.y_bearing;
		cairo_move_to(cr, x, y);
		cairo_show_text(cr, st);
	}
	if (an) // zlab
	{
		char zlab[NAME_MAX];
		sprintf(zlab, "%s", an);
		cairo_select_font_face(cr, "Open Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
		cairo_text_extents(cr, zlab, &ext);
		x = DIM_X - (ext.width + ext.x_bearing);
		y = ext.height / 2 + ext.y_bearing;
		cairo_move_to(cr, x, y);
		cairo_show_text(cr, zlab);
	}
	// xlab
	char xlab[] = "Genome coordinates";
	char ylab[] = "Depth (PE)";
	cairo_select_font_face(cr, "Open Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
	cairo_text_extents(cr, xlab, &ext);
	x = DIM_X / 2.0 - (ext.width / 2.0 + ext.x_bearing);
	y = DIM_Y + MARGIN / 4.0 - (ext.height / 2 + ext.y_bearing);
	cairo_move_to(cr, x, y);
	cairo_show_text(cr, xlab);
	cairo_save(cr);
	// ylab
	cairo_translate(cr, MARGIN / 1.25, HEIGHT / 2.0); // translate origin to the center
	cairo_rotate(cr, 3 * M_PI / 2.0);
	cairo_text_extents(cr, ylab, &ext);
	cairo_move_to(cr, MARGIN / 2.5, -MARGIN * 1.25);
	cairo_show_text(cr, ylab);
	cairo_restore(cr);
	// axis
	draw_arrow(cr, 0, DIM_Y, DIM_X, DIM_Y); // xaxis
	double w1 = 1.0, w2 = 1.0;
	cairo_device_to_user_distance(cr, &w1, &w2);
	cairo_set_line_width(cr, fmin(w1, w2));
	cairo_set_line_cap(cr, CAIRO_LINE_CAP_SQUARE);
	cairo_set_source_rgb(cr, 0, 0, 0);
	cairo_move_to(cr, 0, 0);
	cairo_line_to(cr, 0, DIM_Y); // yaxis
	draw_yticks(cr, md);
	char buf[sizeof(uint64_t) * 8 + 1];
	if (gl >= 1e9)
		sprintf(buf, "%.2fG", gl * 1.0e-9);
	else if (gl >= 1e6)
		sprintf(buf, "%.2fM", gl * 1.0e-6);
	else if (gl >= 1e3)
		sprintf(buf, "%.2fK", gl * 1.0e-3);
	else
		sprintf(buf, "%"PRIu64, gl);
	cairo_select_font_face(cr, "Open Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
	cairo_text_extents(cr, buf, &ext);
	x = DIM_X - ext.width - ext.x_bearing;
	y = DIM_Y + MARGIN / 4.0 - (ext.height / 2 + ext.y_bearing);
	cairo_move_to(cr, x, y);
	cairo_show_text(cr, buf);
}

void draw_ped(cairo_t *cr, kh_t *os, int md, dp_t *dp)
{
	;
}

void out_bed(const cgranges_t *cr, bam_hdr_t *h, const char *out)
{
	int64_t i, j, d = 0, n, p = 0, *cr_b = 0, max_b = 0;
	FILE *fp = out ? fopen(out, "w") : stdout;
	for (i = 0; i < h->n_targets; ++i)
	{
		for (j = 0; j < h->target_len[i]; ++j)
		{
			if ((n = cr_overlap(cr, h->target_name[i], j, j, &cr_b, &max_b)))
			{
				if (n != d)
				{
					if (p)
					{
						fprintf(fp, "%"PRId64"\t%"PRId64"\n", p + 1, d);
						d = p = 0;
					}
					fprintf(fp, "%s\t%"PRId64"\t", h->target_name[i], j);
				}
				d = n;
				p = j;
			}
			else
			{
				if (d)
					fprintf(fp, "%"PRId64"\t%"PRId64"\n", p + 1, d);
				d = p = 0;
			}
		}
		if (d)
			fprintf(fp, "%"PRId64"\t%"PRId64"\n", p + 1, d);
		d = p = 0;
	}
	fclose(fp);
	free(cr_b);
}

void out_bgzf(const cgranges_t *cr, bam_hdr_t *h, const char *out)
{
	int64_t i, j, d = 0, n, p = 0, *cr_b = 0, max_b = 0;
	kstring_t ks = {0, 0, NULL};
	for (i = 0; i < h->n_targets; ++i)
	{
		for (j = 0; j < h->target_len[i]; ++j)
		{
			if ((n = cr_overlap(cr, h->target_name[i], j, j, &cr_b, &max_b)))
			{
				if (n != d)
				{
					if (p)
					{
						ksprintf(&ks, "%"PRId64"\t%"PRId64"\n", p + 1, d);
						d = p = 0;
					}
					ksprintf(&ks, "%s\t%"PRId64"\t", h->target_name[i], j);
				}
				d = n;
				p = j;
			}
			else
			{
				if (d)
					ksprintf(&ks, "%"PRId64"\t%"PRId64"\n", p + 1, d);
				d = p = 0;
			}
		}
		if (d)
			ksprintf(&ks, "%"PRId64"\t%"PRId64"\n", p + 1, d);
		d = p = 0;
	}
	BGZF *fp = bgzf_open(out, "w");
	if (ks.l != bgzf_write(fp, ks.s, ks.l))
		error("Error writing to file [%s]\n", out);
	bgzf_close(fp);
    tbx_conf_t conf = tbx_conf_bed;
	if(tbx_index_build3(out, NULL, 14, 8, &conf))
		error("Error building index for [%s]\n", out);
	ks_release(&ks);
}

int is_gzip(const char *fn)
{
	char buf[2];
	FILE *fp;
	int gzip = 0;
	if ((fp = fopen(fn, "rb")) == NULL)
		error("[ERROR] Unable to open file: %s\n", fn);
	if (fread(buf, 1, 2, fp) == 2)
		if (((int)buf[0] == 0x1f) && ((int)(buf[1]&0xFF) == 0x8b))
			gzip = 1;
	fclose(fp);
	return gzip;
}

bool ends_with(const char *str, const char *sfx)
{
	int ret = 0;
	int str_len = strlen(str);
	int sfx_len = strlen(sfx);
	if ((str_len >= sfx_len) && (0 == strcasecmp(str + (str_len-sfx_len), sfx)))
		ret = 1;
	return ret;
}

void prs_arg(int argc, char **argv, arg_t *arg)
{
	int c = 0;
	ketopt_t opt = KETOPT_INIT;
	const char *opt_str = "i:o:m:p:s:hv";
	while ((c = ketopt(&opt, argc, argv, 1, opt_str, long_options)) >= 0)
	{
		switch (c)
		{
			case 'i': arg->in = opt.arg; break;
			case 'o': arg->out = opt.arg; break;
			case 'm': arg->mis = atoi(opt.arg); break;
			case 'p': arg->plot = opt.arg; break;
			case 's': arg->sub = opt.arg; break;
			case 'h': usage(); break;
			case 'v':
				if (strlen(BRANCH_COMMIT))
					printf("%s (%s)\n", VERSION, BRANCH_COMMIT);
				else
					puts(VERSION);
				exit(EXIT_SUCCESS);
			case '?':
				printf("Invalid option: [%c]\n", opt.opt); exit(EXIT_SUCCESS);
			default:
				printf("Invalid option: [%s]", opt.arg); exit(EXIT_SUCCESS);
		}
	}
	if (!arg->in || access(arg->in, R_OK))
		error("Error: input bam is unspecified or inaccessible!\n");
	if (!(ends_with(arg->in, ".bam") && is_gzip(arg->in)))
		error("Oops! only bam input is supported. Invalid bam file [%s]\n", arg->in);
	if (arg->mis <= 0)
		error("Error: invalid max insert size value specified [%d]\n", arg->mis);
	char bai[PATH_MAX];
	snprintf(bai, PATH_MAX, "%s.bai", arg->in);
	if (access(bai, R_OK))
		error("Error: bam's index file (.bai) is required, please use samtools sort and index to create it.\n");
	/*
	 char *svg
	if (arg->plot && !ends_with(arg->plot, ".svg"))
	{
		asprintf(&svg, "%s.svg", arg->plot);
		arg->plot = svg;
	}
	*/
}

int strlen_wo_esc(const char *str)
{
	int length = 0;
	while (*str != '\0')
	{
		// Check if this is the start of an ANSI escape sequence (ESC [ ... m)
		if (*str == '\033' && *(str + 1) == '[')
		{
			// Move past \033[
			str += 2;
			// Skip until we reach 'm', which ends the ANSI sequence
			while (*str != '\0' && *str != 'm')
				str++;
			// Move past 'm' if we found it
			if (*str == 'm')
				str++;
			continue;
		}
		// Count visible characters
		length++;
		str++;
	}
	return length;
}

void horiz(const int _n)
{
	struct winsize w;
	ioctl(0, TIOCGWINSZ, &w);
	int i, n = (w.ws_col >= _n) ? _n : w.ws_col;
	for (i = 0; i < n; ++i) fputs("\e[90m\xe2\x94\x80\e[0m", stdout);
	fputc('\n', stdout);
}

void horiz_em(const int _n)
{
	struct winsize w;
	ioctl(0, TIOCGWINSZ, &w);
	int i, n = (w.ws_col >= _n) ? _n : w.ws_col;
	for (i = 0; i < n; ++i) fputs("\e[1;90m\xe2\x94\x80\e[0;0m", stdout);
	fputc('\n', stdout);
}


static void usage()
{
	int w = 58;
	horiz_em(w);
	char title[] = "\e[1mCalculate the Paired-End depth of microbial genome\e[0m";
	int title_len = strlen_wo_esc(title);
	printf("%*.*s\n", (int)((w - title_len) / 2 + strlen(title)), (int)strlen(title), title);
	horiz(w);
	printf("%s \e[1mUsage\e[0m: \e[1;31m%s\e[0;0m \e[1;90m[options]\e[0;0m --in <bam> --out <bed>\n", BUL, __progname);
	putchar('\n');
	puts(BUL " \e[1mOptions\e[0m:");
	puts("  -i, --in  \e[3mFILE\e[0m   Input BAM file with bai index");
	puts("  -o, --out \e[3mSTR\e[0m    Output PE depth bed\e[90m[.gz]\e[0m file \e[90m[stdout]\e[0m");
	printf("  -m, --mis \e[3mINT\e[0m    Maximum insert size allowed \e[90m[%d]\e[0m\n", MAX_IS);
	puts("  -p, --plot \e[3mFILE\e[0m  Depth plot svg file \e[90m[none]\e[0m (TODO)");
	puts("  -s, --sub \e[3mFILE\e[0m   Sub-title of depth plot \e[90m[none]\e[0m (TODO)");
	putchar('\n');
	puts("  -h               Show help message");
	puts("  -v               Display program version");
	putchar('\n');
	puts(BUL " \e[1mContact\e[0m: \e[4mmeta@geneplus.cn\e[0m for support and bug report");
	horiz_em(w);
	exit(1);
}
