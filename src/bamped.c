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
	int i, ci = -1;
	// check ctg
	if (arg->ctg)
	{
		for (i = 0; i < hdr->n_targets; ++i)
		{
			if (!strcmp(hdr->target_name[i], arg->ctg))
			{
				ci = i;
				break;
			}
		}
		if (ci == -1)
		{
			fprintf(stderr, "%s Error: specified contig [%s] not found in bam\n", ERR, arg->ctg);
			fprintf(stderr, "Please use samtools view -H %s to view contigs contained\n", arg->in);
			exit(EXIT_FAILURE);
		}
	}
	// multiple ctg offsets
	uint64_t gl = 0, nd = 0; // genome length in total
	kh_t *os = kh_init();
	ld_os(hdr, ci, os, &gl);
	/* dbg os
	for (i = 0; i < hdr->n_targets; ++i)
		printf("%s\t%d\t%"PRIu64"\n", hdr->target_name[i], hdr->target_len[i], kh_xval(os, i));
	printf("%"PRIu64"\n", gl);
	*/
	// PE regions
	ld_gr(fp, hdr, ci, arg->mis, gr);
	cairo_t *cr = NULL;
	cairo_surface_t *sf = NULL;
	if (!arg->plot)
	{
		if (arg->out && ends_with(arg->out, ".gz"))
			out_bgzf(gr, hdr, ci, arg->out);
		else
			out_bed(gr, hdr, ci, arg->out);
	}
	else
	{
		uint64_t i;
		uint32_t md = 0;
		dp_t *dp = NULL;
		ld_dp(gr, hdr, ci, &dp, &md, &nd);
		/* debug dp
		for (i = 0; i < nd; ++i)
			printf("%s\t%d\t%d\t%d\n", hdr->target_name[dp[i].tid], dp[i].pos, dp[i].pos + dp[i].len, dp[i].dep);
		*/
		// dep plot
		if (ends_with(arg->plot, ".svg"))
			sf = cairo_svg_surface_create(arg->plot, WIDTH * 1.02, HEIGHT);
		else if (ends_with(arg->plot, ".png"))
			sf = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, WIDTH, HEIGHT);
		else
			error("Error: unsupported plot format: [%s]\n", arg->plot);
		cr = cairo_create(sf);
		char tt[PATH_MAX] = {'\0'}, an[NAME_MAX] = {'\0'};
		// prepare title
		char *p = strrchr(arg->in, '/');
		snprintf(tt, PATH_MAX, "%s", p ? p + 1 : arg->in);
		*strstr(tt, ".bam") = '\0';
		if (ci != -1)
			snprintf(tt, PATH_MAX, "%s: %s", tt, arg->ctg);
		// prepare anno
		prep_an(dp, nd, gl, an);
		draw_canvas(sf, cr, hdr, ci, os, tt, arg->sub, an, md, gl);
		cairo_save(cr);
		cairo_scale(cr, DIM_X, DIM_Y);
		// iterate dp array and draw to cr
		for (i = 0; i < nd; ++i)
			draw_ped1(cr, os, md, gl, dp + i);
		cairo_restore(cr);
		draw_axis(cr, md, gl);
		dump_dp(hdr, dp, nd, arg->out);
		free(dp);
	}
	cr_destroy(gr);
	bam_hdr_destroy(hdr);
	hts_close(fp);
	if (arg->plot)
	{
		if (ends_with(arg->plot, ".png"))
			cairo_surface_write_to_png(cairo_get_target(cr), arg->plot);
		cairo_surface_destroy(sf);
		cairo_destroy(cr);
	}
	kh_destroy(os);
	free(arg);
	return 0;
}

void ld_os(bam_hdr_t *hdr, int ci, kh_t *os, uint64_t *gl)
{
	int i;
	uint64_t shift = 0;
	if (ci != -1)
	{
		kh_ins(os, ci, 0);
		*gl = hdr->target_len[ci];
	}
	else
	{
		for (i = 0; i < hdr->n_targets; ++i)
		{
			kh_ins(os, i, shift);
			shift += hdr->target_len[i];
		}
		*gl = shift;
	}
}

void ld_gr(samFile *fp, bam_hdr_t *hdr, int ci, int mis, cgranges_t *cr)
{
	bam1_t *b = bam_init1();
	while (sam_read1(fp, hdr, b) >= 0)
	{
		bam1_core_t *c = &b->core;
		if (c->flag & (BAM_FQCFAIL | BAM_FUNMAP | BAM_FREVERSE | BAM_FSECONDARY | BAM_FSUPPLEMENTARY))
			continue;
		if (ci != -1 && c->tid != ci)
			continue;
		if ((c->flag & BAM_FPAIRED) && (c->tid == c->mtid) && c->isize <= mis)
			cr_add(cr, hdr->target_name[c->tid], c->pos, c->pos + c->isize, 0);
	}
	bam_destroy1(b);
	cr_sort(cr);
	cr_index(cr);
}

void prep_an(const dp_t *dp, uint64_t nd, uint64_t gl, char *an)
{
	uint64_t i;
	double cl = 0.0f;
	for (i = 0; i < nd; ++i)
		cl += dp[i].len;
	cl /= gl / 100;
	if (cl == 0)
		snprintf(an, NAME_MAX, "Genome coverage: N/A");
	else if (cl < 0.001)
		snprintf(an, NAME_MAX, "Genome coverage: <0.01‰");
	else if (cl < 0.01)
		snprintf(an, NAME_MAX, "Genome coverage: <0.01%%");
	else
		snprintf(an, NAME_MAX, "Genome coverage: %.*f%%", fmin(100, cl) == 100 ? 0 : 2, fmin(100, cl));
}

void draw_canvas(cairo_surface_t *sf, cairo_t *cr, bam_hdr_t *hdr, int ci,
		const kh_t *os, const char *tt, const char *st, const char *an, uint32_t md,
		uint64_t gl)
{
	cairo_set_antialias(cr, CAIRO_ANTIALIAS_BEST);
	draw_rrect(cr);
	cairo_set_source_rgb(cr, 0, 0, 0);
	cairo_translate(cr, MARGIN / 1.25, MARGIN / 2.0);
	// axis labels
	double x, y;
	cairo_text_extents_t ext;
	if (tt && strlen(tt)) // title
	{
		cairo_set_font_size(cr, 24.0);
		cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
		cairo_text_extents(cr, tt, &ext);
		x = DIM_X / 2.0 - (ext.width / 2.0 + ext.x_bearing);
		y = (st && strlen(st)) ? -ext.height / 2 + ext.y_bearing * 1.25 : ext.height + ext.y_bearing * 2.5;
		cairo_move_to(cr, x, y);
		cairo_show_text(cr, tt);
	}
	if (st && strlen(st)) // sub-title
	{
		cairo_set_font_size(cr, 20.0);
		cairo_select_font_face(cr, "serif", CAIRO_FONT_SLANT_ITALIC, CAIRO_FONT_WEIGHT_BOLD);
		cairo_text_extents(cr, st, &ext);
		x = DIM_X / 2.0 - (ext.width / 2.0 + ext.x_bearing);
		y = ext.height + ext.y_bearing * 2;
		cairo_move_to(cr, x, y);
		cairo_show_text(cr, st);
	}
	if (an && strlen(an)) // zlab
	{
		char zlab[NAME_MAX];
		sprintf(zlab, "%s", an);
		cairo_set_font_size(cr, 18.0);
		cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
		cairo_text_extents(cr, zlab, &ext);
		x = DIM_X - (ext.width + ext.x_bearing);
		y = ext.height / 2 + ext.y_bearing;
		cairo_move_to(cr, x, y);
		cairo_show_text(cr, zlab);
	}
	// xlab
	char xlab[] = "Genome coordinates";
	char ylab[] = "Depth (PE)";
	cairo_set_font_size(cr, 18.0);
	cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
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
	// contig shades
	int i;
	cairo_set_source_rgb(cr, 0.9, 0.9, 0.9);
	if (ci == -1)
	{
		for (i = 1; i < hdr->n_targets; i += 2)
		{
			x = kh_xval(os, i);
			cairo_rectangle(cr, (double)DIM_X * x / gl, 0, (double)DIM_X * hdr->target_len[i] / gl, DIM_Y);
			cairo_fill(cr);
		}
	}
}

void draw_axis(cairo_t *cr, uint32_t md, uint64_t gl)
{
	double x, y;
	cairo_text_extents_t ext;
	draw_arrow(cr, 0, DIM_Y, DIM_X, DIM_Y); // xaxis
	double w1 = 1.0, w2 = 1.0;
	cairo_set_font_size(cr, 16.0);
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
	cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
	cairo_text_extents(cr, buf, &ext);
	x = DIM_X - ext.width - ext.x_bearing;
	y = DIM_Y + MARGIN / 4.0 - (ext.height / 2 + ext.y_bearing);
	cairo_move_to(cr, x, y);
	cairo_show_text(cr, buf);
}

void draw_ped1(cairo_t *cr, kh_t *os, uint32_t md, uint64_t gl, dp_t *dp)
{
	double w1 = 1.0, w2 = 1.0, x, y, w, h;
	cairo_device_to_user_distance(cr, &w1, &w2);
	double lw = fmin(w1, w2) / 2;
	cairo_set_line_width(cr, lw);
	cairo_set_line_cap(cr, CAIRO_LINE_CAP_BUTT);
	cairo_set_source_rgb(cr, 87 / 255.0, 122 / 255.0, 166 / 255.0);
	cairo_set_antialias(cr, CAIRO_ANTIALIAS_NONE);
	double ymx = ceil(log10(md)) + 1;
	cairo_set_source_rgb(cr, 87 / 255.0, 122 / 255.0, 166 / 255.0);
	if (dp->len <= 5) // use hist instead of rectangle to make it visible
	{
		cairo_set_line_width(cr, lw);
		x = (double)(dp->pos + kh_xval(os, dp->tid) + (double)dp->len / 2) / gl;
		y = 1 - (log10(dp->dep) + 1) / ymx;
		cairo_move_to(cr, x, 1);
		cairo_line_to(cr, x, y);
		cairo_stroke(cr);
	}
	else
	{
		x = (double)(dp->pos + kh_xval(os, dp->tid)) / gl;
		y = 1 - (log10(dp->dep) + 1) / ymx;
		w = (double)dp->len / gl;
		h = (log10(dp->dep) + 1) / ymx;
		cairo_rectangle(cr, x, y, w, h);
		cairo_fill(cr);
	}
}

void ld_dp(const cgranges_t *cr, bam_hdr_t *hdr, int ci, dp_t **dp, uint32_t *md,
		uint64_t *nd)
{
	int64_t i, j, d = 0, m = CHUNK, n, p = 0, *cr_b = 0, max_b = 0;
	*dp = malloc(m * sizeof(dp_t));
	for (i = 0; i < hdr->n_targets; ++i)
	{
		if (ci != -1 && i != ci)
			continue;
		for (j = 0; j < hdr->target_len[i]; ++j)
		{
			if ((n = cr_overlap(cr, hdr->target_name[i], j, j, &cr_b, &max_b)))
			{
				if (n != d)
				{
					if (p)
					{
						(*dp)[*nd].len = j - (*dp)[*nd].pos + 1;
						d = p = 0;
						if (++*nd == m)
						{
							m <<= 1;
							*dp = realloc(*dp, m * sizeof(dp_t));
						}
					}
					(*dp)[*nd].tid = i;
					(*dp)[*nd].pos = j;
					(*dp)[*nd].dep = n;
					*md = fmax(*md, n);
				}
				d = n;
				p = j;
			}
			else
			{
				if (d)
				{
					(*dp)[*nd].len = j - (*dp)[*nd].pos + 1;
					if (++*nd == m)
					{
						m <<= 1;
						*dp = realloc(*dp, m * sizeof(dp_t));
					}
				}
				d = p = 0;
			}
		}
		if (d)
		{
			(*dp)[*nd].len = j - (*dp)[*nd].pos + 1;
			if (++*nd == m)
			{
				m <<= 1;
				*dp = realloc(*dp, m * sizeof(dp_t));
			}
		}
		d = p = 0;
	}
	*dp = realloc(*dp, *nd * sizeof(dp_t));
	free(cr_b);
}

void out_bed(const cgranges_t *cr, bam_hdr_t *hdr, int ci, const char *out)
{
	int64_t i, j, d = 0, n, p = 0, *cr_b = 0, max_b = 0;
	FILE *fp = out ? fopen(out, "w") : stdout;
	for (i = 0; i < hdr->n_targets; ++i)
	{
		if (ci != -1 && i != ci)
			continue;
		for (j = 0; j < hdr->target_len[i]; ++j)
		{
			if ((n = cr_overlap(cr, hdr->target_name[i], j, j, &cr_b, &max_b)))
			{
				if (n != d)
				{
					if (p)
					{
						fprintf(fp, "%"PRId64"\t%"PRId64"\n", p + 1, d);
						d = p = 0;
					}
					fprintf(fp, "%s\t%"PRId64"\t", hdr->target_name[i], j);
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

void out_bgzf(const cgranges_t *cr, bam_hdr_t *hdr, int ci, const char *out)
{
	int64_t i, j, d = 0, n, p = 0, *cr_b = 0, max_b = 0;
	kstring_t ks = {0, 0, NULL};
	BGZF *fp = bgzf_open(out, "w");
	for (i = 0; i < hdr->n_targets; ++i)
	{
		if (ci != -1 && i != ci)
			continue;
		for (j = 0; j < hdr->target_len[i]; ++j)
		{
			if ((n = cr_overlap(cr, hdr->target_name[i], j, j, &cr_b, &max_b)))
			{
				if (n != d)
				{
					if (p)
					{
						ksprintf(&ks, "%"PRId64"\t%"PRId64"\n", p + 1, d);
						if (ks.l >= BGZF_BLOCK_SIZE)
						{
							if (ks.l != bgzf_write(fp, ks.s, ks.l))
								error("Error writing to file [%s]\n", out);
							ks.l = 0;
						}
						d = p = 0;
					}
					ksprintf(&ks, "%s\t%"PRId64"\t", hdr->target_name[i], j);
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
	if (ks.l && ks.l != bgzf_write(fp, ks.s, ks.l))
		error("Error writing to file [%s]\n", out);
	ks.l = 0;
	bgzf_close(fp);
    tbx_conf_t conf = tbx_conf_bed;
	if(tbx_index_build3(out, NULL, 14, 8, &conf))
		error("Error building index for [%s]\n", out);
	ks_release(&ks);
}

void dump_dp(bam_hdr_t *hdr, dp_t *dp, uint64_t nd, const char *out)
{
	uint64_t i;
	kstring_t ks = {0, 0, NULL};
	if (out && ends_with(out, ".gz"))
	{
		BGZF *fp = bgzf_open(out, "w");
		for (i = 0; i < nd; ++i)
		{
			ksprintf(&ks, "%s\t%"PRIu64"\t%"PRIu64"\t%d\n", hdr->target_name[dp[i].tid],
					dp[i].pos, dp[i].pos + dp[i].len, dp[i].dep);
			if (ks.l >= BGZF_BLOCK_SIZE)
			{
				if (ks.l != bgzf_write(fp, ks.s, ks.l))
					error("Error writing to file [%s]\n", out);
				ks.l = 0;
			}
		}
		if (ks.l && ks.l != bgzf_write(fp, ks.s, ks.l))
			error("Error writing to file [%s]\n", out);
		ks.l = 0;
		bgzf_close(fp);
		tbx_conf_t conf = tbx_conf_bed;
		if(tbx_index_build3(out, NULL, 14, 8, &conf))
			error("Error building index for [%s]\n", out);
	}
	else
	{
		FILE *fp = out ? fopen(out, "w") : stdout;
		for (i = 0; i < nd; ++i)
			fprintf(fp, "%s\t%"PRIu64"\t%"PRIu64"\t%d\n", hdr->target_name[dp[i].tid],
					dp[i].pos, dp[i].pos + dp[i].len, dp[i].dep);
		fclose(fp);
	}
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
	const char *opt_str = "i:o:m:p:s:c:hv";
	while ((c = ketopt(&opt, argc, argv, 1, opt_str, long_options)) >= 0)
	{
		switch (c)
		{
			case 'i': arg->in = opt.arg; break;
			case 'o': arg->out = opt.arg; break;
			case 'm': arg->mis = atoi(opt.arg); break;
			case 'p': arg->plot = opt.arg; break;
			case 's': arg->sub = opt.arg; break;
			case 'c': arg->ctg = opt.arg; break;
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

void pe_dep_sch()
{
	puts("\e[90m                            ┌─┐\n"
		"                    ┌───────┘ └───┐\n"
		"              ──────┴─────────────┴──────\e[0m\n"
		"                    \e[31m——→\e[0m    \e[34m←——\e[0m\n"
		"                            \e[31m——\e[0m\e[34m←\e[0m\e[31m→\e[0m\e[34m——\e[0m");
}

static void usage()
{
	int w = 58;
	horiz_em(w);
	char title[] = "\e[1mCalculate the Paired-End depth of microbial genome\e[0m";
	int title_len = strlen_wo_esc(title);
	printf("%*.*s\n", (int)((w - title_len) / 2 + strlen(title)), (int)strlen(title), title);
	pe_dep_sch();
	horiz(w);
	printf("%s \e[1mUsage\e[0m: \e[1;31m%s\e[0;0m \e[1;90m[options]\e[0;0m --in <bam> --out <bed>\n", BUL, __progname);
	putchar('\n');
	puts(BUL " \e[1mOptions\e[0m:");
	puts("  -i, --in  \e[3mFILE\e[0m   Input BAM file with bai index");
	puts("  -o, --out \e[3mSTR\e[0m    Output PE depth bed\e[90m[.gz]\e[0m file \e[90m[stdout]\e[0m");
	printf("  -m, --mis \e[3mINT\e[0m    Maximum insert size allowed \e[90m[%d]\e[0m\n", MAX_IS);
	puts("  -p, --plot \e[3mFILE\e[0m  Depth plot png file \e[90m[none]\e[0m");
	puts("  -s, --sub \e[3mFILE\e[0m   Sub-title of depth plot \e[90m[none]\e[0m");
	puts("  -c, --ctg \e[3mSTR\e[0m    Restrict analysis to this contig \e[90m[none]\e[0m");
	putchar('\n');
	puts("  -h               Show help message");
	puts("  -v               Display program version");
	putchar('\n');
	puts(BUL " \e[1mContact\e[0m: \e[4mmeta@geneplus.cn\e[0m for support and bug report");
	horiz_em(w);
	exit(1);
}
