#include "stdlib.h"
#include "stdio.h"
#include "htslib/sam.h"
#include "htslib/khash.h"
#include "vector.h"
#include "bioidx/bioidx.h"
#include "transmap_gtf.h"

int exon_comp(const void *a, const void *b){
    return (*(exon_t **)a)->start - (*(exon_t **)b)->start;
}



void transcript_free(transcript_t *tr){
    exon_t *exon;
    int j;
    if (tr->exons) {
        for (j = 0; j < tr->exons->size; ++j) {
            if ((exon = tr->exons->data[j])) free(exon);
        }
        vec_destroy(exon, tr->exons);
    }
    if (tr->chrom) free(tr->chrom);
    if (tr->name) free(tr->name);
    free(tr);
}

void gtf_free(gtf_dict_t *gtf){
    transcript_t *tr;
    khiter_t i;
    int j;
    if (gtf->idx) bioidx_destroy(gtf->idx);
    if (gtf->record) {
        for (i = 0; i < kh_end(gtf->record); ++i) {
            if (kh_exist(gtf->record, i) && (tr = kh_val(gtf->record, i)) != NULL)
                transcript_free(tr);
        }
        kh_destroy(transcript, gtf->record);
    }
    free(gtf);
}
gtf_dict_t *gtf_parse(const char* fname, const char* used_feature, const char* used_attribute){
    int used_attribute_len = strlen(used_attribute);
    char buffer[2048];
    char *items[10];
    int new_tid = 0;
    FILE *f = fopen(fname, "r");
    if (!f) return NULL;
    gtf_dict_t *gtf;
    if (!(gtf = malloc(sizeof(*gtf)))) return NULL;
    gtf->record = NULL;
    gtf->idx = NULL;
    gtf->record = kh_init(transcript);
    if (!gtf->record) goto clean_up;
    gtf->idx = bioidx_init();
    if (!gtf->idx) goto clean_up;
    int ret;
    char *attr_begin, *attr_end;
    int line_count = 0;
    while (fgets(buffer, 2048, f)){
        line_count++;
        if (buffer[0] == '#') continue;
        strsplit(buffer, items, 10, '\t');
        if (strcmp(items[2], used_feature) != 0) continue;
        attr_begin = items[8];
        while (attr_begin && (strncmp(attr_begin, used_attribute, used_attribute_len) != 0 || attr_begin[used_attribute_len] != ' ')){
            attr_begin = strpbrk(attr_begin, ";\"");
            if (attr_begin && *attr_begin == '\"') {
                attr_begin = strchr(attr_begin + 1, '\"');
                if (attr_begin) attr_begin = strchr(attr_begin + 1, ';');
            }
            if (attr_begin){
                attr_begin++;
                while (*attr_begin == ' ') attr_begin++;
            }

        }
        if (!attr_begin) {
            fprintf(stderr, "[gtf parse] attribute \"%s\" not found for line %d.\n", used_attribute, line_count);
            continue;
        }
        attr_begin += strlen(used_attribute);
        while (*attr_begin == ' ') attr_begin++;
        if (*attr_begin == '\"') attr_end = strchr(++attr_begin, '\"');
        else attr_end = strchr(attr_begin, ';');
        if (!attr_end) {
            fprintf(stderr, "[gtf parse] the attribute field seems to be incomplete for line %d.\n", line_count);
            continue;
        }
        *attr_end = '\0';
        khiter_t i = kh_get(transcript, gtf->record, attr_begin);
        transcript_t *tr;
        char *new_str;
        if (i == kh_end(gtf->record)){
            new_str = strdup(attr_begin);
            if (!new_str) goto clean_up;
            i = kh_put(transcript, gtf->record, new_str, &ret);
            if (ret < 0) goto clean_up;
            tr = calloc(1, sizeof(*tr));
            if (!tr) {kh_val(gtf->record, i) = NULL; goto clean_up;}
            tr->new_tid = new_tid++;
            tr->chrom = strdup(items[0]);
            tr->name = new_str;
            tr->strand = items[6][0];
            if (!tr->chrom) goto clean_up;
            tr->exons = vec_init(exon);
            if (!tr->exons) goto clean_up;
            kh_val(gtf->record, i) = tr;
        } else tr = kh_val(gtf->record, i);
        exon_t *exon = calloc(1, sizeof(*exon));
        if (!exon) goto clean_up;
        exon->new_tid = tr->new_tid;
        exon->chrom = tr->chrom;
        exon->start = strtol(items[3], NULL, 0) - 1;
        exon->end = strtol(items[4], NULL, 0);
        exon->strand = items[6][0];
        exon->tr = tr;
        if (vec_add(exon, tr->exons, exon) != 0) goto clean_up;
    }
    fclose(f);
    khiter_t i;
    int j;
    for (i = 0; i < kh_end(gtf->record); ++i){
        if (!kh_exist(gtf->record, i)) continue;
        transcript_t *tr = kh_val(gtf->record, i);
        exon_t **exons = tr->exons->data;
        size_t exon_size = tr->exons->size;
        exon_t *exon, *prev_exon;
        qsort(exons, exon_size, sizeof(*exons), exon_comp);
        tr->len = 0;
        tr->start = exons[0]->start;
        tr->end = exons[exon_size - 1]->end;
        for (j = 0; j < exon_size; ++j) {
            exon = exons[j];
            exon->idx = j;
            exon->tstart = tr->len;
            tr->len += exon->end - exon->start;
            if (j != 0){
                prev_exon = exons[j - 1];
                if (prev_exon->end > exon->start){
                    fprintf(stderr, "[gtf parse] exons overlap for transcript %s.\n", tr->name);
                    break;
                }
                if (strcmp(prev_exon->chrom, exon->chrom) != 0) {
                    fprintf(stderr, "[gtf parse] different chromosomes for exons from transcript %s.\n",tr->name);
                    break;
                }
                if (prev_exon->strand != exon->strand){
                    fprintf(stderr, "[gtf parse] inconsistent strand for exons from transcript %s.\n", tr->name);
                    break;
                }
            }
        }
        if (j != exon_size){
            kh_del(transcript, gtf->record, i);
            transcript_free(tr);
        }

    }
    return gtf;
    clean_up:
    fclose(f);
    gtf_free(gtf);
    return NULL;
}




