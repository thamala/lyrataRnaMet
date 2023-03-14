/*
 Copyright (C) 2023 Tuomas Hamala

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 For any other inquiries send an email to tuomas.hamala@gmail.com

 ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

 Program for compiling methylation data to create meta-plots. Output is location scaled in relation to a gene/TE and average methylation propotion.
 Locations are shown as follows: -1 to 0 upstream, 0 to 1 gene/TE body, 1 to 2 downstream.
 Methylation proportions are assumed to be in combined BEDGRAPH format, created with bedtools unionbedg (example: bedtools unionbedg -header -filler . -names ind0 ind1 ind2 -i met0.bg met1.bg met2.bg > out.bg).

 Compiling: gcc bg2meta_plot.c -o bg2meta_plot -lm

 Usage:
 -bg [file] Methylation propotions in BEDGRAPH format. Needs to be sorted based on chrom and start position.
 -bed [file] Bed file listing regions to use (required fields: chrom, start, end, name, score, strand). Needs to be sorted based on chrom and start position.
 -inds [file] File listing individuals to include. Optional.
 -bp [int] Distance around regions to include. Default 1000.
 -min [int] Minimum number of individuals required to consider a site. Default 1.

 Example:
 ./bg2meta_plot -bg test.bg -bed genes.bed -inds inds.txt -bp 1000 -min 2 > out.txt
*/

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#define merror "\nERROR: System out of memory\n\n"

typedef struct {
    double start, end;
    char str, chr[50], id[100];
} bed_s;

void openFiles(int argc, char *argv[]);
bed_s *readBed(FILE *bed_file, int *n);
char **readInds(FILE *inds_file, int *n);
void readBg(FILE *bg_file, bed_s *beds, char **inds, int min, int bed_n, int ind_n, double bp);
int isNumeric(const char *s);
void lineTerminator(char *line);

int main(int argc, char *argv[]) {
    int second = 0, minute = 0, hour = 0;
    time_t timer = 0;

    timer = time(NULL);
    openFiles(argc, argv);
    second = time(NULL) - timer;
    minute = second / 60;
    hour = second / 3600;

    fprintf(stderr, "\nDone!");
    if(hour > 0)
        fprintf(stderr, "\nElapsed time: %i h, %i min & %i sec\n\n", hour, minute - hour * 60, second - minute * 60);
    else if(minute > 0)
        fprintf(stderr, "\nElapset time: %i min & %i sec\n\n", minute, second - minute * 60);
    else if(second > 5)
        fprintf(stderr, "\nElapsed time: %i sec\n\n", second);
    else
        fprintf(stderr, "\n\n");

    return 0;
}

void openFiles(int argc, char *argv[]) {
    int i, bed_n = 0, ind_n = 0, min = 1;
    double bp = 1000;
    char **inds;
    bed_s *beds = NULL;
    FILE *bg_file = NULL, *bed_file = NULL, *ind_file = NULL;

    fprintf(stderr, "\nParameters:\n");

    for(i = 1; i < argc; i++) {
        if(strcmp(argv[i], "-bg") == 0) {
            if((bg_file = fopen(argv[++i], "r")) == NULL) {
                fprintf(stderr, "\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-bg %s\n", argv[i]);
        }

        else if(strcmp(argv[i], "-bed") == 0) {
            if((bed_file = fopen(argv[++i], "r")) == NULL) {
                fprintf(stderr, "\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-bed %s\n", argv[i]);
        }

        else if(strcmp(argv[i], "-inds") == 0) {
            if((ind_file = fopen(argv[++i], "r")) == NULL) {
                fprintf(stderr, "\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-inds %s\n", argv[i]);
        }

        else if(strcmp(argv[i], "-bp") == 0) {
            if(isNumeric(argv[++i]))
                bp = atof(argv[i]);
            fprintf(stderr, "\t-bp %s\n", argv[i]);
        }

        else if(strcmp(argv[i], "-min") == 0) {
            if(isNumeric(argv[++i]))
                min = atoi(argv[i]);
            fprintf(stderr, "\t-min %s\n", argv[i]);
        }

        else {
            fprintf(stderr, "\nERROR: Unknown argument '%s'\n\n", argv[i]);
            exit(EXIT_FAILURE);
        }
    }

    fprintf(stderr, "\n");

    if(bg_file == NULL || bed_file == NULL) {
        fprintf(stderr, "\nERROR: -bg [file] and -bed [file] are required!\n");
        exit(EXIT_FAILURE);
    }

    beds = readBed(bed_file, &bed_n);
    if(ind_file != NULL)
        inds = readInds(ind_file, &ind_n);

    readBg(bg_file, beds, inds, min, bed_n, ind_n, bp);
}

bed_s *readBed(FILE *bed_file, int *n) {
    int i, char_i = 0, line_i = 0, maxchar = 0;
    char c, *line = NULL;
    bed_s *list = NULL;

    while((c = fgetc(bed_file)) != EOF) {
        char_i++;
        if(c == '\n') {
            line_i++;
            if(char_i > maxchar)
                maxchar = char_i;
            char_i = 0;
        }
    }

    rewind(bed_file);

    if((line = malloc((maxchar + 1) * sizeof(char))) == NULL) {
        fprintf(stderr, merror);
        exit(EXIT_FAILURE);
    }
    if((list = malloc(line_i * sizeof(bed_s))) == NULL) {
        fprintf(stderr, merror);
        exit(EXIT_FAILURE);
    }

    while(fgets(line, maxchar + 1, bed_file) != NULL) {
        if(line[0] == '\n')
            continue;
        lineTerminator(line);
        strncpy(list[*n].chr, strtok(line, "\t"), 49);
        list[*n].start = atof(strtok(NULL, "\t")) + 1;
        list[*n].end = atof(strtok(NULL, "\t"));
        strncpy(list[*n].id, strtok(NULL, "\t"), 99);
        strtok(NULL, "\t");
        list[*n].str = strtok(NULL, "\t")[0];
        *n = *n + 1;
    }

    free(line);
    fclose(bed_file);

    return list;
}

char **readInds(FILE *ind_file, int *n) {
    int i = 0, char_i = 0, maxchar = 0, line_i = 0;
    char c, *line, **list;

    while((c = fgetc(ind_file)) != EOF) {
        char_i++;
        if(c == '\n') {
            line_i++;
            if(char_i > maxchar)
                maxchar = char_i;
            char_i = 0;
        }
    }

    rewind(ind_file);

    if((line = malloc((maxchar + 1) * sizeof(char))) == NULL) {
        fprintf(stderr, merror);
        exit(EXIT_FAILURE);
    }
    if((list = malloc((line_i + 1) * sizeof(char *))) == NULL) {
        fprintf(stderr, merror);
        exit(EXIT_FAILURE);
    }
    for(i = 0; i < line_i + 1; i++) {
        if((list[i] = malloc((maxchar + 1) * sizeof(char))) == NULL) {
            fprintf(stderr, merror);
            exit(EXIT_FAILURE);
        }
    }
    i = 0;
    while(fgets(line, maxchar + 1, ind_file) != NULL) {
        if(line[0] == '\n')
            continue;
        lineTerminator(line);
        strcpy(list[i], line);
        i++;
        *n = *n + 1;
    }
    list[i][0] = '\0';

    free(line);
    fclose(ind_file);

    return list;
}

void readBg(FILE *bg_file, bed_s *beds, char **inds, int min, int bed_n, int ind_n, double bp) {
    int i, j, k = 0, l = 0, n = 0, ok = 0, pos = 0, bed_i = 0, *mlist;
    double dist = 0, met = 0, met_i = 0, bmet = 0, bmet_i = 0;
    char chr[50], *line = NULL, *temp = NULL;
    size_t len = 0;
    ssize_t read;

    while((read = getline(&line, &len, bg_file)) != -1) {
        if(line[0] == '\n')
            continue;
        lineTerminator(line);
        temp = strtok(line, "\t");
        i = 1;
        n = 0;
        ok = 0;
        met = 0;
        met_i = 0;
        if(strcmp(temp, "chrom") == 0) {
            if(ind_n == 0)
                continue;
            if((mlist = malloc(len * sizeof(int))) == NULL) {
                fprintf(stderr, merror);
                exit(EXIT_FAILURE);
            }
            while(temp != NULL) {
                if(i > 3) {
                    for(j = 0; j < ind_n; j++) {
                        if(strcmp(temp, inds[j]) == 0)
                            break;
                    }
                    if(j == ind_n)
                        mlist[n] = 0;
                    else
                        mlist[n] = 1;
                    n++;
                }
                temp = strtok(NULL, "\t");
                i++;
            }
            continue;
        }
        while(temp != NULL) {
            if(i == 1)
                strncpy(chr, temp, 49);
            else if(i == 3) {
                pos = atoi(temp);
                while(bed_i < bed_n) {
                    if(strcmp(chr, beds[bed_i].chr) == 0) {
                        if(pos <= beds[bed_i].end + bp && pos >= beds[bed_i].start - bp) {
                            ok = 1;
                            break;
                        } else if(pos < beds[bed_i].start - bp)
                            break;
                    } else if(strcmp(chr, beds[bed_i].chr) < 0)
                        break;
                    bed_i++;
                }
                if(ok == 0)
                    break;
            } else if(i > 3) {
                if(temp[0] != '.') {
                    if(ind_n > 0) {
                        if(mlist[n] == 1) {
                            met += atof(temp) / 100;
                            met_i++;
                        }
                    } else {
                        met += atof(temp) / 100;
                        met_i++;
                    }
                }
                n++;
            }
            temp = strtok(NULL, "\t");
            i++;
        }
        if(temp == NULL && met_i >= min) {
            for(i = l; i < bed_n; i++) {
                if(strcmp(chr, beds[i].chr) == 0) {
                    if(pos <= beds[i].end + bp && pos >= beds[i].start - bp) {
                        if(pos < beds[i].start && beds[i].str == '+')
                            dist = (pos - beds[i].start) / bp;
                        else if(pos < beds[i].start && beds[i].str == '-')
                            dist = 1 + (beds[i].start - pos) / bp;
                        else if(pos > beds[i].end && beds[i].str == '+')
                            dist = 1 + (pos - beds[i].end) / bp;
                        else if(pos > beds[i].end && beds[i].str == '-')
                            dist = (beds[i].end - pos) / bp;
                        else {
                            dist = (beds[i].end - pos) / (beds[i].end - beds[i].start + 1);
                            bmet += met;
                            bmet_i += met_i;
                        }
                        printf("%f\t%f\t%s\n", dist, met / met_i, beds[i].id);
                    } else if(pos < beds[i].start - bp) {
                        l = i;
                        for(j = 1; j <= i; j++) {
                            if(pos <= beds[i - j].end + bp && pos >= beds[i - j].start - bp)
                                l = i - j;
                            else if(pos > beds[i - j].end + bp) {
                                if(i - j - 1 >= 0) {
                                    if(pos > beds[i - j - 1].end + bp)
                                        break;
                                } else
                                    break;
                            }
                        }
                        break;
                    }
                } else if(strcmp(chr, beds[i].chr) < 0) {
                    l = i;
                    for(j = 1; j <= i; j++) {
                        if(strcmp(chr, beds[i - j].chr) == 0) {
                            if(pos <= beds[i - j].end + bp && pos >= beds[i - j].start - bp)
                                l = i - j;
                            else if(pos > beds[i - j].end + bp) {
                                if(i - j - 1 >= 0) {
                                    if(pos > beds[i - j - 1].end + bp)
                                        break;
                                } else
                                    break;
                            }
                        } else if(strcmp(chr, beds[i - j].chr) > 0)
                            break;
                    }
                    break;
                }
            }
        }
    }
    fprintf(stderr, "Average body methylation = %.2f\n", bmet / bmet_i);
}

int isNumeric(const char *s) {
    char *p;

    if(s == NULL || *s == '\0' || isspace(*s))
        return 0;
    strtod(s, &p);
    return *p == '\0';
}

void lineTerminator(char *line) {
    int i;

    for(i = 0; line[i] != 0; i++) {
        if(line[i] == '\n' || line[i] == '\r')
            line[i] = '\0';
    }
}
