/*
 Copyright (C) 2021 Tuomas Hamala

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

 Program for estimating pairwise nucleotide diversity (pi) using genotype probabilities.
 The probability file is expected to include both variant and invariant sites.

 Compiling: gcc probs2pi.c -o probs2pi -lm

 Usage:
 -beagle [file] Posterior genotype probabilities in Beagle format (generated e.g., with Angsd or PCAngsd).
 -genes [file] Tab delimited file listing genes (format chr, start, end, strand [+ or -], id). Optional.
 -bp [int] Distance around genes to calculate pi for up- and downstream areas. Optional.
 -min [int] Minimum number of individuals required to consider a site. Default 2.

 Example:
 ./probs2pi -beagle postprobs.beagle -genes genes.txt -bp 1000 -min 6 > test.txt
*/

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#define merror "\nERROR: System out of memory\n\n"
#define NA 0.333333

typedef struct {
    int L;
    double tP;
} Theta_s;

typedef struct {
    int start, end;
    char str, chr[101], id[101];
    Theta_s up, cds, down;
} Gene_s;

void openFiles(int argc, char *argv[]);
Gene_s *readGenes(FILE *gene_file, int *n);
void readBeagle(FILE *beagle_file, Gene_s *genes, int bp, int gene_n, int min);
void printOut(Gene_s gene, int bp, int i);
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
    int i, j, gene_n = 0, min = 2, bp = 0;
    Gene_s *genes = NULL;
    FILE *beagle_file = NULL, *gene_file = NULL;

    fprintf(stderr, "\nParameters:\n");

    for(i = 1; i < argc; i++) {
        if(strcmp(argv[i], "-beagle") == 0) {
            if((beagle_file = fopen(argv[++i], "r")) == NULL) {
                fprintf(stderr, "\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-beagle %s\n", argv[i]);
        }

        else if(strcmp(argv[i], "-genes") == 0) {
            if((gene_file = fopen(argv[++i], "r")) == NULL) {
                fprintf(stderr, "\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-genes %s\n", argv[i]);
        }

        else if(strcmp(argv[i], "-bp") == 0) {
            if(isNumeric(argv[++i]))
                bp = atoi(argv[i]);
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

    if(beagle_file == NULL) {
        fprintf(stderr, "\nERROR: -beagle [file] is required!\n");
        exit(EXIT_FAILURE);
    }

    if(gene_file != NULL)
        genes = readGenes(gene_file, &gene_n);

    readBeagle(beagle_file, genes, bp, gene_n, min);
}

Gene_s *readGenes(FILE *gene_file, int *n) {
    int i, char_i = 0, line_i = 0, maxchar = 0;
    char c, *line = NULL;
    Gene_s *list = NULL;

    while((c = fgetc(gene_file)) != EOF) {
        char_i++;
        if(c == '\n') {
            line_i++;
            if(char_i > maxchar)
                maxchar = char_i;
            char_i = 0;
        }
    }

    rewind(gene_file);

    if((line = malloc((maxchar + 1) * sizeof(char))) == NULL) {
        fprintf(stderr, merror);
        exit(EXIT_FAILURE);
    }
    if((list = malloc(line_i * sizeof(Gene_s))) == NULL) {
        fprintf(stderr, merror);
        exit(EXIT_FAILURE);
    }
    memset(list, 0, sizeof(Gene_s) * line_i);

    while(fgets(line, maxchar + 1, gene_file) != NULL) {
        if(line[0] == '\n')
            continue;
        lineTerminator(line);
        strncpy(list[*n].chr, strtok(line, "\t"), 100);
        list[*n].start = atoi(strtok(NULL, "\t"));
        list[*n].end = atoi(strtok(NULL, "\t"));
        list[*n].str = strtok(NULL, "\t")[0];
        strncpy(list[*n].id, strtok(NULL, "\t"), 100);
        *n = *n + 1;
    }

    free(line);
    fclose(gene_file);

    return list;
}

void readBeagle(FILE *beagle_file, Gene_s *genes, int bp, int gene_n, int min) {
    int i, j = 0, k = 0, pos = 0, ok = 0, gene_i = 0, kept_i = 0, site_i = 0;
    double prob[3] = {0}, p = 0, n = 0;
    char chr[101], *line = NULL, *temp = NULL, *end = NULL;
    size_t len = 0;
    ssize_t read;

    while((read = getline(&line, &len, beagle_file)) != -1) {
        if(line[0] == '\n')
            continue;
        lineTerminator(line);
        temp = strtok_r(line, "\t", &end);
        i = 1;
        j = 0;
        p = 0;
        n = 0;
        if(strcmp(temp, "marker") == 0)
            continue;
        site_i++;
        strncpy(chr, strtok(temp, "_"), 100);
        pos = atoi(strtok(NULL, "_"));
        if(gene_n > 0) {
            ok = 0;
            while(gene_i < gene_n) {
                if(strcmp(chr, genes[gene_i].chr) == 0) {
                    if(pos <= genes[gene_i].end + bp && pos >= genes[gene_i].start - bp) {
                        ok = 1;
                        break;
                    } else if(pos < genes[gene_i].start - bp)
                        break;
                } else if(strcmp(chr, genes[gene_i].chr) < 0)
                    break;
                gene_i++;
            }
            if(ok == 0)
                continue;
        }
        while(temp != NULL) {
            if(i > 3) {
                prob[j] = atof(temp);
                j++;
                if(j == 3) {
                    if(prob[0] != NA | prob[1] != NA | prob[2] != NA) {
                        p += prob[1] + 2 * prob[2];
                        n += 2;
                    }
                    j = 0;
                }
            }
            temp = strtok_r(NULL, "\t", &end);
            i++;
        }
        if(temp == NULL && n / 2 >= min) {
            p /= n;
            kept_i++;
            if(gene_n == 0)
                printf("%s\t%i\t%f\n", chr, pos, 2 * p * (1 - p));
            else {
                for(i = k; i < gene_n; i++) {
                    if(strcmp(chr, genes[i].chr) == 0) {
                        if(pos <= genes[i].end + bp && pos >= genes[i].start - bp) {
                            if(pos < genes[i].start && genes[i].str == '+') {
                                genes[i].up.tP += 2 * p * (1 - p);
                                genes[i].up.L++;
                            } else if(pos < genes[i].start && genes[i].str == '-') {
                                genes[i].down.tP += 2 * p * (1 - p);
                                genes[i].down.L++;
                            } else if(pos > genes[i].end && genes[i].str == '+') {
                                genes[i].down.tP += 2 * p * (1 - p);
                                genes[i].down.L++;
                            } else if(pos > genes[i].end && genes[i].str == '-') {
                                genes[i].up.tP += 2 * p * (1 - p);
                                genes[i].up.L++;
                            } else {
                                genes[i].cds.tP += 2 * p * (1 - p);
                                genes[i].cds.L++;
                            }
                        } else if(pos < genes[i].start - bp) {
                            k = i;
                            for(j = 1; j <= i; j++) {
                                if(pos <= genes[i - j].end + bp && pos >= genes[i - j].start - bp)
                                    k = i - j;
                                else if(pos > genes[i - j].end + bp) {
                                    if(i - j - 1 >= 0) {
                                        if(pos > genes[i - j - 1].end + bp)
                                            break;
                                    } else
                                        break;
                                }
                            }
                            break;
                        }
                    } else if(strcmp(chr, genes[i].chr) < 0) {
                        k = i;
                        for(j = 1; j <= i; j++) {
                            if(strcmp(chr, genes[i - j].chr) == 0) {
                                if(pos <= genes[i - j].end + bp && pos >= genes[i - j].start - bp)
                                    k = i - j;
                                else if(pos > genes[i - j].end + bp) {
                                    if(i - j - 1 >= 0) {
                                        if(pos > genes[i - j - 1].end + bp)
                                            break;
                                    } else
                                        break;
                                }
                            } else if(strcmp(chr, genes[i - j].chr) > 0)
                                break;
                        }
                        break;
                    }
                }
            }
        }
    }

    if(gene_n > 0) {
        if(isatty(1))
            fprintf(stderr, "\n");
        for(i = 0; i < gene_n; i++)
            printOut(genes[i], bp, i);
    }

    if(isatty(1))
        fprintf(stderr, "\n");
    fprintf(stderr, "Kept %i out of %i sites\n", kept_i, site_i);

    free(line);
    if(gene_n > 0)
        free(genes);
    fclose(beagle_file);
}

void printOut(Gene_s gene, int bp, int i) {
    if(bp == 0) {
        if(i == 0)
            printf("id\tcoding_tP\tcoding_n\n");
        printf("%s\t%f\t%i\n", gene.id, gene.cds.tP, gene.cds.L);
    } else {
        if(i == 0)
            printf("id\tup_tP\tup_n\tcoding_tP\tcoding_n\tdown_tP\tdown_n\n");
        printf("%s\t%f\t%i\t%f\t%i\t%f\t%i\n", gene.id, gene.up.tP, gene.up.L, gene.cds.tP, gene.cds.L, gene.down.tP, gene.down.L);
    }
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