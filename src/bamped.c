#include "bamped.h"

int main(int argc, char *argv[])
{
	if (argc == 1)
		usage();
	arg_t *arg = calloc(1, sizeof(arg_t));
	arg->mis = MAX_IS;
	prs_arg(argc, argv, arg);
	samFile *fp = sam_open(arg->in, "r");
	if (!fp)
		error("Error: failed to read input bam [%s]\n", arg->in);
	bam_hdr_t *h = sam_hdr_read(fp);
	bam1_t *b = bam_init1();
	cgranges_t *cr = cr_init();
	while (sam_read1(fp, h, b) >= 0)
	{
		bam1_core_t *c = &b->core;
		if (c->flag & (BAM_FQCFAIL | BAM_FUNMAP | BAM_FREVERSE | BAM_FSECONDARY | BAM_FSUPPLEMENTARY))
			continue;
		if ((c->flag & BAM_FPAIRED) && (c->tid == c->mtid) && c->isize <= arg->mis)
			cr_add(cr, h->target_name[c->tid], c->pos, c->pos + c->isize, 0);
	}
	cr_sort(cr);
	cr_index(cr);
	if (arg->out && ends_with(arg->out, ".gz"))
		out_bgzf(cr, h, arg->out);
	else
		out_bed(cr, h, arg->out);
	cr_destroy(cr);
	bam_hdr_destroy(h);
	bam_destroy1(b);
	hts_close(fp);
	free(arg);
	return 0;
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
	char *bai;
	asprintf(&bai, "%s.bai", arg->in);
	if (access(bai, R_OK))
	{
		free(bai);
		error("Error: bam's index file (.bai) is required, please use samtools sort and index to create it.\n");
	}
	free(bai);
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
	puts("  -p, --plot \e[3mFILE\e[0m  Depth plot svg file \e[90m[none] (TODO)\e[0m");
	puts("  -s, --sub \e[3mFILE\e[0m   Sub-title of depth plot \e[90m[none] (TODO)\e[0m");
	putchar('\n');
	puts("  -h               Show help message");
	puts("  -v               Display program version");
	putchar('\n');
	puts(BUL " \e[1mContact\e[0m: \e[4mmeta@geneplus.cn\e[0m for support and bug report");
	horiz_em(w);
	exit(1);
}
