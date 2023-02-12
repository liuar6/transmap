#include <stdlib.h>
#include <string.h>
#include "htslib/sam.h"
#include "htslib/khash.h"
#include "vector.h"
#include "bioidx/bioidx.h"

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

typedef struct exon_t{
    //int32_t tid;
    int32_t new_tid;
    char *chrom;
    hts_pos_t start;
    hts_pos_t tstart;
    hts_pos_t end;
    //char *name;
    char strand;
    int idx;
    struct transcript_t *tr;
} exon_t;
VEC_INIT(exon, exon_t *);
typedef struct transcript_t{
    char* chrom;
    char *name;
    char strand;
    hts_pos_t start;
    hts_pos_t end;
    int32_t len;
    int32_t new_tid;
    int32_t tid;
    vec_t(exon) *exons;
} transcript_t;
KHASH_MAP_INIT_STR(transcript, transcript_t *);

typedef struct gtf_dict_t{
    khash_t (transcript) *record;
    bioidx_t *idx;
} gtf_dict_t;

gtf_dict_t *gtf_parse(const char* fname, const char *used_feature, const char *used_attribute);
void gtf_free(gtf_dict_t *);

static int exon_search_comp(const void *a, const void *b){
    if ((*(exon_t **)a)->tr->new_tid !=  (*(exon_t **)b)->tr->new_tid) return ((*(exon_t **)a)->tr->new_tid - (*(exon_t **)b)->tr->new_tid);
    else return (*(exon_t **)a)->idx - (*(exon_t **)b)->idx;
}

static inline void gtf_search(bioidx_t *idx, int32_t key, int32_t start, int32_t end, vec_t(exon) *hits){
    bioidx_itr_t itr;
    exon_t *hit;
    size_t init_index;
    int i, j;
    init_index = hits->size;
    if (bioidx_search(idx, &itr, key, start, end) != 0) return;
    while ((hit = bioidx_itr_next(&itr))) vec_add(exon, hits, hit);
    if (hits->size > init_index) {
        qsort(hits->data + init_index, hits->size - init_index, sizeof(*(hits->data)), exon_search_comp);
        for (i = init_index, j = init_index + 1; j < hits->size; ++j){
            if (hits->data[j]->new_tid != hits->data[i]->new_tid) hits->data[++i] = hits->data[j];
        }
        hits->size = i + 1;
    }
}

static inline void gtf_search_one(bioidx_t *idx, int32_t key, int32_t start, int32_t end, vec_t(exon) *hits){
    vec_clear(exon, hits);
    gtf_search(idx, key, start, end, hits);
}

static inline void gtf_search_any(bioidx_t *idx, int32_t key, int32_t start1, int32_t end1, int32_t start2, int32_t end2, vec_t(exon) *hits){
    vec_clear(exon, hits);
    gtf_search(idx, key, start1, end1, hits);
    gtf_search(idx, key, start2, end2, hits);
    if (hits->size > 0) {
        qsort(hits->data, hits->size, sizeof(*(hits->data)), exon_search_comp);
        int i, j;
        for (i = 0, j = 1; j < hits->size; ++j){
            if (hits->data[j]->new_tid != hits->data[i]->new_tid) hits->data[++i] = hits->data[j];
        }
        hits->size = i + 1;
    }
};

static inline void gtf_search_both(bioidx_t *idx, int32_t key, int32_t start1, int32_t end1, int32_t start2, int32_t end2, vec_t(exon) *hits){
    vec_clear(exon, hits);
    gtf_search(idx, key, start1, end1, hits);
    gtf_search(idx, key, start2, end2, hits);
    qsort(hits->data, hits->size, sizeof(*(hits->data)), exon_search_comp);
    int i, j;

    for (i = 0, j = 0; j + 1 < hits->size; ++j) {
        if (hits->data[j]->new_tid == hits->data[j + 1]->new_tid) hits->data[i++] = hits->data[j];
    }
    hits->size = i;
}