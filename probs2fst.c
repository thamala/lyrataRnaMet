/*
 Copyright (C) 2022 Tuomas Hamala

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

 Program for estimating Weir & Cockerham's Fst across arbitrary number of populations using genotype probabilities.

 Compiling: gcc probs2fst.c -o probs2fst -lm

 Usage:
 -beagle [file] Posterior genotype probabilities in Beagle format (generated e.g. with Angsd or PCAngsd).
 -pop [file] File listing individuals from a single population. Can be used >= 2 times.
 -genes [file] Tab delimited file listing genes (format chr, start, end, strand [+ or -], id). Optional.
 -bp [int] Distance around genes to calculate Fst for up- and downstream areas. Optional.
 -min [int] Minimum number of individuals per population required to consider a site. Default 1.
 -maf [double] Minimum minor allele frequency required to consider a site. Default 0.

 Example:
 ./probs2fst -beagle postprobs.beagle -pop list1.txt -pop list2.txt -pop list3.txt -genes genes.txt -bp 1000 -min 6 -maf 0.05 > test.txt
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
    int n;
    double hw, hb;
} Var_s;

typedef struct {
    int start, end;
    char str, chr[101], id[101];
    Var_s up, cds, down;
} Gene_s;

void openFiles(int argc, char *argv[]);
char **readPop(FILE *pop_file, int *n);
Gene_s *readGenes(FILE *gene_file, int *n);
void readBeagle(FILE *beagle_file, char ***pops, Gene_s *genes, int bp, int pop_n, int gene_n, int ind_n, int min, double maf);
Var_s estVars(double **dosaga, int **plist, int pop_n, int ind_n, int plist_n, int min, double maf);
double estFst(Var_s vars);
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
    int i, j, gene_n = 0, min = 1, pop_n = 0, ind_n = 0, bp = 0;
    double maf = 0;
    char ***pops = NULL;
    Gene_s *genes = NULL;
    FILE *beagle_file = NULL, *pop_file = NULL, *gene_file = NULL;

    fprintf(stderr, "\nParameters:\n");

    if((pops = malloc(argc * sizeof(char **))) == NULL) {
        fprintf(stderr, merror);
        exit(EXIT_FAILURE);
    }

    for(i = 1; i < argc; i++) {
        if(strcmp(argv[i], "-beagle") == 0) {
            if((beagle_file = fopen(argv[++i], "r")) == NULL) {
                fprintf(stderr, "\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-beagle %s\n", argv[i]);
        }

        else if(strcmp(argv[i], "-pop") == 0) {
            if((pop_file = fopen(argv[++i], "r")) == NULL) {
                fprintf(stderr, "\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-pop %s\n", argv[i]);

            pops[pop_n] = readPop(pop_file, &ind_n);
            pop_n++;
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

        else if(strcmp(argv[i], "-maf") == 0) {
            if(isNumeric(argv[++i]))
                maf = atof(argv[i]);
            fprintf(stderr, "\t-maf %s\n", argv[i]);
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

    if(pop_n < 2) {
        fprintf(stderr, "\nERROR: at least two population files (-pop [file) are required!\n");
        exit(EXIT_FAILURE);
    }

    if(gene_file != NULL)
        genes = readGenes(gene_file, &gene_n);

    readBeagle(beagle_file, pops, genes, bp, pop_n, gene_n, ind_n, min, maf);
}

char **readPop(FILE *pop_file, int *n) {
    int i = 0, char_i = 0, maxchar = 0, line_i = 0;
    char c, *line, **list;

    while((c = fgetc(pop_file)) != EOF) {
        char_i++;
        if(c == '\n') {
            line_i++;
            if(char_i > maxchar)
                maxchar = char_i;
            char_i = 0;
        }
    }

    rewind(pop_file);

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
    while(fgets(line, maxchar + 1, pop_file) != NULL) {
        if(line[0] == '\n')
            continue;
        lineTerminator(line);
        strcpy(list[i], line);
        i++;
        *n = *n + 1;
    }
    list[i][0] = '\0';

    free(line);
    fclose(pop_file);

    return list;
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

void readBeagle(FILE *beagle_file, char ***pops, Gene_s *genes, int bp, int pop_n, int gene_n, int ind_n, int min, double maf) {
    int i, j = 0, k = 0, l = 0, m = 0, n = 0, pos = 0, ok = 0, p_i = 0, gene_i = 0, kept_i = 0, site_i = 0, **plist = NULL;
    double prob[3] = {0}, **dosage;
    char chr[101], *line = NULL, *temp = NULL, *end = NULL;
    Var_s vars = {0};
    size_t len = 0;
    ssize_t read;

    if((plist = malloc(ind_n * sizeof(int *))) == NULL) {
        fprintf(stderr, merror);
        exit(EXIT_FAILURE);
    }
    for(i = 0; i < ind_n; i++) {
        if((plist[i] = malloc(2 * sizeof(int))) == NULL) {
            fprintf(stderr, merror);
            exit(EXIT_FAILURE);
        }
    }

    while((read = getline(&line, &len, beagle_file)) != -1) {
        if(line[0] == '\n')
            continue;
        lineTerminator(line);
        temp = strtok_r(line, "\t", &end);
        i = 1;
        j = 0;
        k = 0;
        if(strcmp(temp, "marker") == 0) {
            while(temp != NULL) {
                if(i > 3) {
                    if(j == 0) {
                        for(k = 0; k < pop_n; k++) {
                            for(l = 0; pops[k][l][0] != '\0'; l++) {
                                if(strcmp(temp, pops[k][l]) == 0) {
                                    plist[p_i][0] = n;
                                    plist[p_i][1] = k;
                                    p_i++;
                                }
                            }
                        }
                        n++;
                    }
                    j++;
                    if(j == 3)
                        j = 0;
                }
                temp = strtok_r(NULL, "\t", &end);
                i++;
            }
            if(p_i == 0) {
                fprintf(stderr, "ERROR: Individuals in pop files were not found in the Beagle file!\n\n");
                exit(EXIT_FAILURE);
            }
            if(p_i < ind_n)
                fprintf(stderr, "Warning: Pop files contain individuals that are not in the Beagle file\n");
            fprintf(stderr, "Kept %i individuals from %i populations\n", p_i, pop_n);
            if((dosage = malloc(n * sizeof(double *))) == NULL) {
                fprintf(stderr, merror);
                exit(EXIT_FAILURE);
            }
            for(j = 0; j < n; j++) {
                if((dosage[j] = malloc(2 * sizeof(double))) == NULL) {
                    fprintf(stderr, merror);
                    exit(EXIT_FAILURE);
                }
            }
            if(gene_n == 0) {
                if(isatty(1))
                    fprintf(stderr, "\n");
                printf("chr\tbp\tfst\n");
            }
            continue;
        }
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
                        dosage[k][0] = prob[1] + 2 * prob[2];
                        dosage[k][1] = prob[1];
                    } else
                        dosage[k][0] = 9;
                    j = 0;
                    k++;
                }
            }
            temp = strtok_r(NULL, "\t", &end);
            i++;
        }
        if(temp == NULL) {
            vars = estVars(dosage, plist, pop_n, n, p_i, min, maf);
            if(isnan(vars.hw) == 1)
                continue;
            kept_i++;
            if(gene_n == 0)
                printf("%s\t%i\t%f\n", chr, pos, estFst(vars));
            else {
                for(i = m; i < gene_n; i++) {
                    if(strcmp(chr, genes[i].chr) == 0) {
                        if(pos <= genes[i].end + bp && pos >= genes[i].start - bp) {
                            if(pos < genes[i].start && genes[i].str == '+') {
                                genes[i].up.hw += vars.hw;
                                genes[i].up.hb += vars.hb;
                                genes[i].up.n++;
                            } else if(pos < genes[i].start && genes[i].str == '-') {
                                genes[i].down.hw += vars.hw;
                                genes[i].down.hb += vars.hb;
                                genes[i].down.n++;
                            } else if(pos > genes[i].end && genes[i].str == '+') {
                                genes[i].down.hw += vars.hw;
                                genes[i].down.hb += vars.hb;
                                genes[i].down.n++;
                            } else if(pos > genes[i].end && genes[i].str == '-') {
                                genes[i].up.hw += vars.hw;
                                genes[i].up.hb += vars.hb;
                                genes[i].up.n++;
                            } else {
                                genes[i].cds.hw += vars.hw;
                                genes[i].cds.hb += vars.hb;
                                genes[i].cds.n++;
                            }
                        } else if(pos < genes[i].start - bp) {
                            m = i;
                            for(j = 1; j <= i; j++) {
                                if(pos <= genes[i - j].end + bp && pos >= genes[i - j].start - bp)
                                    m = i - j;
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
                        m = i;
                        for(j = 1; j <= i; j++) {
                            if(strcmp(chr, genes[i - j].chr) == 0) {
                                if(pos <= genes[i - j].end + bp && pos >= genes[i - j].start - bp)
                                    m = i - j;
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
        if(bp == 0) {
            printf("id\tcoding_fst\tcoding_n\n");
            for(i = 0; i < gene_n; i++)
                printf("%s\t%f\t%i\n", genes[i].id, estFst(genes[i].cds), genes[i].cds.n);
        } else {
            printf("id\tup_fst\tup_n\tcoding_fst\tcoding_n\tdown_fst\tdown_n\n");
            for(i = 0; i < gene_n; i++)
                printf("%s\t%f\t%i\t%f\t%i\t%f\t%i\n", genes[i].id, estFst(genes[i].up), genes[i].up.n, estFst(genes[i].cds), genes[i].cds.n, estFst(genes[i].down), genes[i].down.n);
        }
    }

    if(isatty(1))
        fprintf(stderr, "\n");
    fprintf(stderr, "Kept %i out of %i sites\n", kept_i, site_i);

    free(line);
    free(plist);
    free(pops);
    if(gene_n > 0)
        free(genes);
    fclose(beagle_file);
}

Var_s estVars(double **dosage, int **plist, int pop_n, int ind_n, int plist_n, int min, double maf) {
    int i, j, ok = 1;
    double a = 0, b = 0, c = 0, pbar = 0, nbar = 0, hbar = 0, n_sum = 0, n_sum2 = 0, nc = 0, r = 0, s2 = 0, *p = NULL, *n = NULL;
    Var_s vars;

    if((p = malloc(pop_n * sizeof(double))) == NULL) {
        fprintf(stderr, merror);
        exit(EXIT_FAILURE);
    }
    if((n = malloc(pop_n * sizeof(double))) == NULL) {
        fprintf(stderr, merror);
        exit(EXIT_FAILURE);
    }
    memset(p, 0, sizeof(double) * pop_n);
    memset(n, 0, sizeof(double) * pop_n);

    for(i = 0; i < ind_n; i++) {
        if(dosage[i][0] == 9)
            continue;
        for(j = 0; j < plist_n; j++) {
            if(plist[j][0] == i) {
                p[plist[j][1]] += dosage[i][0];
                n[plist[j][1]]++;
                pbar += dosage[i][0];
                hbar += dosage[i][1];
                n_sum++;
            }
        }
    }
    for(i = 0; i < pop_n; i++) {
        p[i] /= n[i] * 2;
        n_sum2 += (n[i] * n[i]);
        if(n[i] < min)
            ok = 0;
    }
    r = (double)pop_n;
    nbar = n_sum / r;
    pbar /= n_sum * 2;
    hbar /= n_sum;
    if(ok == 0 || pbar < maf || 1 - pbar < maf) {
        free(p);
        free(n);
        vars.hw = 0.0 / 0.0;
        return vars;
    }
    for(i = 0; i < pop_n; i++)
        s2 += n[i] * (p[i] - pbar) * (p[i] - pbar);
    s2 /= (r - 1.0) * nbar;
    nc = n_sum - (n_sum2 / n_sum) / (r - 1.0);
    a = (s2 - (pbar * (1.0 - pbar) - (((r - 1.0) * s2) / r) - (hbar / 4.0)) / (nbar - 1.0)) * nbar / nc;
    b = (pbar * (1.0 - pbar) - (s2 * (r - 1.0) / r) - hbar * (((2.0 * nbar) - 1.0) / (4.0 * nbar))) * nbar / (nbar - 1.0);
    c = hbar / 2.0;
    if(isnan(a) || isnan(b) || isnan(c))
        vars.hw = 0.0 / 0.0;
    else {
        vars.hw = a;
        vars.hb = a + b + c;
    }
    free(p);
    free(n);
    return vars;
}

double estFst(Var_s vars) {
    return vars.hw / vars.hb;
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
