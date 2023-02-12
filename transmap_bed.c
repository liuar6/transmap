/* The MIT License (MIT)

   Copyright (c) 2023 Anrui Liu <liuar6@gmail.com>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   “Software”), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "transmap_bed.h"

static char ** strsplit(char * line, char ** results, int length,char c){
    char *start=line;
    char *end=NULL;
    int i=0;
    while ((end=strchr(start, c))!=NULL && i<length){
        end[0]='\0';
        results[i]=start;
        start=end+1;
        i=i+1;
    }
    if (i<length && start[0]!='\0') {
        results[i]=start;
        i=i+1;
    }
    for (;i<length;++i) results[i]=NULL;
    return results;
}

void record_free(bed_t* record){
    if (record){
        if (record->chrom) free(record->chrom);
        if (record->name) free(record->name);
        free(record);
    }
}

void bed_free(bed_dict_t *bed){
    if (!bed) return;
    if (bed->record) {
        for (int i = 0; i < bed->size; ++i){
            bed_t *record = bed->record[i];
            record_free(record);
        }
        free(bed->record);
    }
    free(bed);
}

bed_dict_t *bed_parse(const char* fname){
    char buffer[1024];
    char *items[7];
    int new_tid = 0;
    FILE *f = fopen(fname, "r");
    if (!f) return NULL;
    bed_dict_t *bed;
    if (!(bed = malloc(sizeof(*bed)))) return NULL;
    bed->record = malloc(sizeof(*(bed->record)));
    if (bed->record == NULL) goto clean_up;
    bed->size = 0;
    bed->capacity = 1;
    bed->idx = bioidx_init();
    if (!bed->idx) goto clean_up;
    while (fgets(buffer, 1024, f)){
        strsplit(buffer, items, 7, '\t');
        bed_t *record = calloc(1, sizeof(bed_t));
        if (!record) goto clean_up;
        record->new_tid = new_tid++;
        record->chrom = strdup(items[0]);
        if (!record->chrom)  goto clean_up;
        record->start = strtol(items[1], NULL, 0);
        record->end = strtol(items[2], NULL, 0);
        record->name = strdup(items[3]);
        if (!record->name) goto clean_up;
        record->strand = items[5][0];
        if (bed->size == bed->capacity) {
            void *new_record = realloc(bed->record, ((bed->capacity)<<=1) * sizeof(*(bed->record)));
            if (!new_record) goto clean_up;
            bed->record = new_record;
        }
        bed->record[bed->size++] = record;
    }
    fclose(f);
    return bed;
    clean_up:
    fclose(f);
    bed_free(bed);
    return NULL;

}

#define TRANSMAP_UNMAPPED_ORIGINAL 5
#define TRANSMAP_UNMAPPED_NO_OVERLAP 4
#define TRANSMAP_UNMAPPED_PARTIAL 3
#define TRANSMAP_UNMAPPED_NO_MATCH 2
#define TRANSMAP_MAPPED 0




