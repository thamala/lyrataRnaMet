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

 Program for generating input files for est-sfs (Keightley & Jackson, 2018), using genotype probabilities and whole-genome alignments conducted with MUMmer.
 Requires coordinates files produced by 'show-coords' program from MUMmer (settings -T -H) and substitution files produced by 'show-snps' program from MUMmer (settings -C -H -T).
 Polymorphism data are assumed to be genotype probabilities in Beagle format, generated e.g. with ANGSD (-doMaf 2 -doMajorMinor 4 -doPost 1 -beagleProb 1).

 Only sites that have outgroup information in two of the three species are used.
 If the Beagle file contains missing data, missing alleles are imputed by drawing them from a Bernoulli distribution. 

 Compiling: gcc make_est-sfs.c -o make_est-sfs -lm

 Usage:
 -coord1 [file] coordinates file from outgroup 1 (closet outgroup)
 -coord2 [file] coordinates file from outgroup 2 (mid outgroup)
 -coord3 [file] coordinates file from outgroup 3 (distant outgroup)
 -div1 [file] substitution file from outgroup 1 (closet outgroup)
 -div2 [file] substitution file from outgroup 2 (mid outgroup)
 -div3 [file] substitution file from outgroup 3 (distant outgroup)
 -beagle [file] genotype probabilities in Beagle format
 -region [file] tab-delimited file defining regions to include (chr, start, end)
 -sites [file] tab-delimited file defining sites to include (chr, pos)

 Example:
 ./make_est-sfs \
    -coord1 lyrata-thaliana.coord \
    -coord2 lyrata-capsella.coord \
    -coord3 lyrata-arabis.coord \
    -div1 lyrata-thaliana.snps \
    -div2 lyrata-capsella.snps \
    -div3 lyrata-arabis.snps \
    -beagle <(zcat J1.beagle.gprobs.gz) \
    -region DEG_field.txt \
    -sites 0fold.sites > J1_DEG_field_est-sfs.txt
*/

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#define merror "\nERROR: System out of memory\n"
#define NA 0.333333

typedef struct {
    int chr, start, stop;
    char id[50];
} Region_s;

typedef struct {
    int chr, pos;
    char ref, alt;
} Site_s;

void openFiles(int argc, char *argv[]);
Region_s *readRegions(FILE *region_file, int *n);
Site_s *readSites(FILE *site_file, int *n);
Region_s *readCoord(FILE *coord_file, int *n);
Site_s *readDiv(FILE *div_file, int *n);
void readBeagle(FILE *beagle_file, Region_s *coord1, Region_s *coord2, Region_s *coord3, Region_s *regions, Site_s *div1, Site_s *div2, Site_s *div3, Site_s *sites, int co_n1, int co_n2, int co_n3, int rg_n, int div_n1, int div_n2, int div_n3, int site_n);
char defOut(Region_s *coord, Site_s *div, char ref, int chr, int pos, int co_n, int div_n, int *co_i, int *div_i);
void printOut(char nuc, int end);
void lineTerminator(char *line);

int main(int argc, char *argv[]) {
    int second = 0, minute = 0, hour = 0;
    time_t timer = 0;

    timer = time(NULL);

    openFiles(argc, argv);

    second = time(NULL) - timer;

    minute = second / 60;

    hour = second / 3600;

    if(isatty(1))
        fprintf(stderr, "\n");

    if(hour > 0)
        fprintf(stderr, "Run finished in %i h, %i min & %i sec\n\n", hour, minute - hour * 60, second - minute * 60);
    else if(minute > 0)
        fprintf(stderr, "Run finished in %i min & %i sec\n\n", minute, second - minute * 60);
    else if(second > 5)
        fprintf(stderr, "Run finished in %i sec\n\n", second);
    else
        fprintf(stderr, "\n");

    return 0;
}

void openFiles(int argc, char *argv[]) {
    int i, co_n1 = 0, co_n2 = 0, co_n3 = 0, rg_n = 0, div_n1 = 0, div_n2 = 0, div_n3 = 0, site_n = 0;
    FILE *coord_file1 = NULL, *coord_file2 = NULL, *coord_file3 = NULL, *div_file1 = NULL, *div_file2 = NULL, *div_file3 = NULL, *beagle_file = NULL, *region_file = NULL, *site_file = NULL;
    Region_s *coord1, *coord2, *coord3, *regions;
    Site_s *div1, *div2, *div3, *sites;

    fprintf(stderr, "\nParameters:\n");

    for(i = 1; i < argc; i++) {
        if(strcmp(argv[i], "-coord1") == 0) {
            if((coord_file1 = fopen(argv[++i], "r")) == NULL) {
                fprintf(stderr, "\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-coord1 %s\n", argv[i]);
        }

        else if(strcmp(argv[i], "-coord2") == 0) {
            if((coord_file2 = fopen(argv[++i], "r")) == NULL) {
                fprintf(stderr, "\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-coord2 %s\n", argv[i]);
        }

        else if(strcmp(argv[i], "-coord3") == 0) {
            if((coord_file3 = fopen(argv[++i], "r")) == NULL) {
                fprintf(stderr, "\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-coord3 %s\n", argv[i]);
        }

        else if(strcmp(argv[i], "-div1") == 0) {
            if((div_file1 = fopen(argv[++i], "r")) == NULL) {
                fprintf(stderr, "\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-div1 %s\n", argv[i]);
        }

        else if(strcmp(argv[i], "-div2") == 0) {
            if((div_file2 = fopen(argv[++i], "r")) == NULL) {
                fprintf(stderr, "\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-div2 %s\n", argv[i]);
        }

        else if(strcmp(argv[i], "-div3") == 0) {
            if((div_file3 = fopen(argv[++i], "r")) == NULL) {
                fprintf(stderr, "\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-div3 %s\n", argv[i]);
        }

        else if(strcmp(argv[i], "-beagle") == 0) {
            if((beagle_file = fopen(argv[++i], "r")) == NULL) {
                fprintf(stderr, "\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-beagle %s\n", argv[i]);
        }

        else if(strcmp(argv[i], "-region") == 0) {
            if((region_file = fopen(argv[++i], "r")) == NULL) {
                fprintf(stderr, "\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-region %s\n", argv[i]);
        }

        else if(strcmp(argv[i], "-sites") == 0) {
            if((site_file = fopen(argv[++i], "r")) == NULL) {
                fprintf(stderr, "\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-sites %s\n", argv[i]);
        }

        else {
            fprintf(stderr, "\nERROR: Unknown argument '%s'\n\n", argv[i]);
            exit(EXIT_FAILURE);
        }
    }

    fprintf(stderr, "\n");

    if(div_file1 == NULL || div_file2 == NULL || div_file3 == NULL || coord_file1 == NULL || coord_file2 == NULL || coord_file3 == NULL || beagle_file == NULL) {
        fprintf(stderr, "ERROR: The following parameters are required: -coord1 [file] -coord2 [file] -coord3 [file] -div1 [file] -div2 [file] -div3 [file] -beagle [file]\n\n");
        exit(EXIT_FAILURE);
    }

    if(region_file != NULL)
        regions = readRegions(region_file, &rg_n);
    if(site_file != NULL)
        sites = readSites(site_file, &site_n);

    coord1 = readCoord(coord_file1, &co_n1);
    coord2 = readCoord(coord_file2, &co_n2);
    coord3 = readCoord(coord_file3, &co_n3);
    div1 = readDiv(div_file1, &div_n1);
    div2 = readDiv(div_file2, &div_n2);
    div3 = readDiv(div_file3, &div_n3);
    readBeagle(beagle_file, coord1, coord2, coord3, regions, div1, div2, div3, sites, co_n1, co_n2, co_n3, rg_n, div_n1, div_n2, div_n3, site_n);
}

Region_s *readRegions(FILE *region_file, int *n) {
    int i, char_i = 0, line_i = 0, maxchar = 0;
    char c, *line, *temp;
    Region_s *list;

    while((c = fgetc(region_file)) != EOF) {
        char_i++;
        if(c == '\n') {
            line_i++;
            if(char_i > maxchar)
                maxchar = char_i;
            char_i = 0;
        }
    }

    rewind(region_file);

    if((list = malloc((line_i) * sizeof(Region_s))) == NULL) {
        fprintf(stderr, merror);
        exit(EXIT_FAILURE);
    }

    if((line = malloc((maxchar + 1) * sizeof(char))) == NULL) {
        fprintf(stderr, merror);
        exit(EXIT_FAILURE);
    }

    while(fgets(line, maxchar + 1, region_file) != NULL) {
        if(isdigit(line[0])) {
            lineTerminator(line);
            list[*n].chr = atoi(strtok(line, "\t"));
            list[*n].start = atoi(strtok(NULL, "\t"));
            list[*n].stop = atoi(strtok(NULL, "\t"));
            *n = *n + 1;
        }
    }

    free(line);
    fclose(region_file);

    return list;
}

Site_s *readSites(FILE *site_file, int *n) {
    int i, j, char_i = 0, line_i = 0, maxchar = 0;
    char c, *line, *temp;
    Site_s *list;

    while((c = fgetc(site_file)) != EOF) {
        char_i++;
        if(c == '\n') {
            line_i++;
            if(char_i > maxchar)
                maxchar = char_i;
            char_i = 0;
        }
    }

    rewind(site_file);

    if((list = malloc(line_i * sizeof(Site_s))) == NULL) {
        fprintf(stderr, merror);
        exit(EXIT_FAILURE);
    }

    if((line = malloc((maxchar + 1) * sizeof(char))) == NULL) {
        fprintf(stderr, merror);
        exit(EXIT_FAILURE);
    }

    while(fgets(line, maxchar + 1, site_file) != NULL) {
        if(isdigit(line[0])) {
            lineTerminator(line);
            list[*n].chr = atoi(strtok(line, "\t"));
            list[*n].pos = atoi(strtok(NULL, "\t"));
            *n = *n + 1;
        }
    }

    free(line);
    fclose(site_file);

    return list;
}

Region_s *readCoord(FILE *coord_file, int *n) {
    int i, char_i = 0, line_i = 0, maxchar = 0;
    char c, *line, *temp;
    Region_s *list;

    while((c = fgetc(coord_file)) != EOF) {
        char_i++;
        if(c == '\n') {
            line_i++;
            if(char_i > maxchar)
                maxchar = char_i;
            char_i = 0;
        }
    }

    rewind(coord_file);

    if((list = malloc((line_i) * sizeof(Region_s))) == NULL) {
        fprintf(stderr, merror);
        exit(EXIT_FAILURE);
    }

    if((line = malloc((maxchar + 1) * sizeof(char))) == NULL) {
        fprintf(stderr, merror);
        exit(EXIT_FAILURE);
    }

    while(fgets(line, maxchar + 1, coord_file) != NULL) {
        lineTerminator(line);
        list[*n].start = atoi(strtok(line, "\t"));
        list[*n].stop = atoi(strtok(NULL, "\t"));
        for(i = 0; i < 6; i++)
            temp = strtok(NULL, "\t");
        if(isdigit(temp[0])) {
            list[*n].chr = temp[0] - '0';
            *n = *n + 1;
        }
    }

    free(line);
    fclose(coord_file);

    return list;
}

Site_s *readDiv(FILE *div_file, int *n) {
    int i, j, char_i = 0, line_i = 0, maxchar = 0;
    char c, *line, *temp;
    Site_s *list;

    while((c = fgetc(div_file)) != EOF) {
        char_i++;
        if(c == '\n') {
            line_i++;
            if(char_i > maxchar)
                maxchar = char_i;
            char_i = 0;
        }
    }

    rewind(div_file);

    if((list = malloc(line_i * sizeof(Site_s))) == NULL) {
        fprintf(stderr, merror);
        exit(EXIT_FAILURE);
    }

    if((line = malloc((maxchar + 1) * sizeof(char))) == NULL) {
        fprintf(stderr, merror);
        exit(EXIT_FAILURE);
    }

    while(fgets(line, maxchar + 1, div_file) != NULL) {
        lineTerminator(line);
        list[*n].pos = atoi(strtok(line, "\t"));
        temp = strtok(NULL, "\t");
        list[*n].ref = temp[0];
        temp = strtok(NULL, "\t");
        list[*n].alt = temp[0];
        for(i = 0; i < 6; i++)
            temp = strtok(NULL, "\t");
        if(isdigit(temp[0])) {
            list[*n].chr = temp[0] - '0';
            *n = *n + 1;
        }
    }

    free(line);
    fclose(div_file);

    return list;
}

void readBeagle(FILE *beagle_file, Region_s *coord1, Region_s *coord2, Region_s *coord3, Region_s *regions, Site_s *div1, Site_s *div2, Site_s *div3, Site_s *sites, int co_n1, int co_n2, int co_n3, int rg_n, int div_n1, int div_n2, int div_n3, int site_n) {
    int i, j = 0, chr = 0, pos = 0, co_i1 = 0, co_i2 = 0, co_i3 = 0, rg_i = 0, div_i1 = 0, div_i2 = 0, div_i3 = 0, site_i = 0, ok = 0;
    double ref_i = 0, alt_i = 0, mis_i = 0, p = 0, ind[3] = {0};
    char ref = 'N', alt = 'N', out1 = 'N', out2 = 'N', out3 = 'N', *line = NULL, *loc = NULL, *temp = NULL, *end = NULL;
    size_t len = 0;
    ssize_t read;
    FILE *out_file;

    srand(time(NULL));

    if((out_file = fopen("info.txt", "w")) == NULL) {
        fprintf(stderr, "\nERROR: Cannot create file 'info.txt' \n\n");
        exit(EXIT_FAILURE);
    }

    while((read = getline(&line, &len, beagle_file)) != -1) {
        if(line[0] == '\n' || line[0] == 'm')
            continue;
        lineTerminator(line);
        temp = strtok_r(line, "\t", &end);
        chr = atoi(strtok(temp, "_"));
        pos = atoi(strtok(NULL, "_"));
        if(rg_n > 0) {
            while(rg_i < rg_n) {
                if(chr == regions[rg_i].chr) {
                    if(pos <= regions[rg_i].stop && pos >= regions[rg_i].start) {
                        ok = 1;
                        break;
                    } else if(pos < regions[rg_i].start)
                        break;
                } else if(chr < regions[rg_i].chr)
                    break;
                rg_i++;
            }
            if(ok == 0)
                continue;
        }
        if(site_n > 0) {
            ok = 0;
            while(site_i < site_n) {
                if(chr == sites[site_i].chr) {
                    if(pos == sites[site_i].pos) {
                        ok = 1;
                        break;
                    } else if(pos < sites[site_i].pos)
                        break;
                } else if(chr < sites[site_i].chr)
                    break;
                site_i++;
            }
            if(ok == 0)
                continue;
        }
        i = 1;
        ref_i = 0;
        alt_i = 0;
        mis_i = 0;
        ok = 0;
        out1 = 'N';
        out2 = 'N';
        out3 = 'N';
        while(temp != NULL) {
            if(i == 2) {
                if(temp[0] == '0')
                    ref = 'A';
                else if(temp[0] == '1')
                    ref = 'C';
                else if(temp[0] == '2')
                    ref = 'G';
                else
                    ref = 'T';
            } else if(i == 3) {
                if(temp[0] == '0')
                    alt = 'A';
                else if(temp[0] == '1')
                    alt = 'C';
                else if(temp[0] == '2')
                    alt = 'G';
                else
                    alt = 'T';
                ok = 0;
                out1 = defOut(coord1, div1, ref, chr, pos, co_n1, div_n1, &co_i1, &div_i1);
                out2 = defOut(coord2, div2, ref, chr, pos, co_n2, div_n2, &co_i2, &div_i2);
                out3 = defOut(coord3, div3, ref, chr, pos, co_n3, div_n3, &co_i3, &div_i3);
                if((out1 == 'N' && out2 == 'N') || (out1 == 'N' && out3 == 'N') || (out2 == 'N' && out3 == 'N'))
                    break;
            } else if(i > 3) {
                ind[j] = atof(temp);
                j++;
                if(j == 3) {
                    if(ind[0] != NA | ind[1] != NA | ind[2] != NA) {
                        ref_i += 2 * ind[0] + ind[1];
                        alt_i += ind[1] + 2 * ind[2];
                    } else
                        mis_i += 2;
                    j = 0;
                }
            }
            temp = strtok_r(NULL, "\t", &end);
            i++;
        }
        if(temp != NULL)
            continue;
        if(mis_i > 0) {
            p = round(alt_i) / round(ref_i + alt_i);
            if(p == 0)
                ref_i += mis_i;
            else if(p == 1)
                alt_i += mis_i;
            else {
                for(i = 0; i < mis_i; i++) {
                    if((double)rand() / RAND_MAX < p)
                        alt_i++;
                    else
                        ref_i++;
                }
            }
        }
        fprintf(out_file, "%i\t%i\n", chr, pos);
        ref_i = round(ref_i);
        alt_i = round(alt_i);
        if(ref == 'A')
            printf("%.0f,", ref_i);
        else if(alt == 'A')
            printf("%.0f,", alt_i);
        else
            printf("0,");
        if(ref == 'C')
            printf("%.0f,", ref_i);
        else if(alt == 'C')
            printf("%.0f,", alt_i);
        else
            printf("0,");
        if(ref == 'G')
            printf("%.0f,", ref_i);
        else if(alt == 'G')
            printf("%.0f,", alt_i);
        else
            printf("0,");
        if(ref == 'T')
            printf("%.0f\t", ref_i);
        else if(alt == 'T')
            printf("%.0f\t", alt_i);
        else
            printf("0\t");
        printOut(out1, 0);
        printOut(out2, 0);
        printOut(out3, 1);
    }
    free(line);
    free(coord1);
    free(coord2);
    free(coord3);
    free(div1);
    free(div2);
    free(div3);
    if(rg_n > 0)
        free(regions);
    if(site_n > 0)
        free(sites);
    fclose(beagle_file);
}

char defOut(Region_s *coord, Site_s *div, char ref, int chr, int pos, int co_n, int div_n, int *co_i, int *div_i) {
    int ok = 0;

    while(*co_i < co_n) {
        if(chr == coord[*co_i].chr) {
            if(pos <= coord[*co_i].stop && pos >= coord[*co_i].start) {
                ok = 1;
                while(*div_i < div_n) {
                    if(chr == div[*div_i].chr) {
                        if(pos == div[*div_i].pos) {
                            ok = 2;
                            break;
                        } else if(pos < div[*div_i].pos)
                            break;
                    } else if(chr < div[*div_i].chr)
                        break;
                    *div_i = *div_i + 1;
                }
                if(ok == 1)
                    return ref;
                else if(ok == 2) {
                    if(div[*div_i].alt == '.')
                        return 'N';
                    else
                        return div[*div_i].alt;
                }
                break;
            } else if(pos < coord[*co_i].start)
                break;
        } else if(chr < coord[*co_i].chr)
            break;
        *co_i = *co_i + 1;
    }

    return 'N';
}

void printOut(char nuc, int end) {
    if(end == 0) {
        if(nuc == 'A')
            printf("1,0,0,0 ");
        else if(nuc == 'C')
            printf("0,1,0,0 ");
        else if(nuc == 'G')
            printf("0,0,1,0 ");
        else if(nuc == 'T')
            printf("0,0,0,1 ");
        else
            printf("0,0,0,0 ");
    } else {
        if(nuc == 'A')
            printf("1,0,0,0\n");
        else if(nuc == 'C')
            printf("0,1,0,0\n");
        else if(nuc == 'G')
            printf("0,0,1,0\n");
        else if(nuc == 'T')
            printf("0,0,0,1\n");
        else
            printf("0,0,0,0\n");
    }
}

void lineTerminator(char *line) {
    int i;
    for(i = 0; line[i] != 0; i++) {
        if(line[i] == '\n' || line[i] == '\r')
            line[i] = '\0';
    }
}
