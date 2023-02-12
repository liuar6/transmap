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
#include "htslib/sam.h"
#include "transmap_bam.h"

bam_vector_t *bam_vector_init(){
    bam_vector_t *bv;
    bv = calloc(1, sizeof(*bv));
    return bv;
}

bam1_t *bam_vector_next(bam_vector_t *bv){
    if (bv->size == bv->capacity) {
        size_t new_capacity = bv->capacity<<1u;
        if (new_capacity < 8) new_capacity = 8;
        bam1_t **new_data = realloc(bv->data, new_capacity * sizeof(*bv->data));
        if (!new_data) return NULL;
        for (int i = bv->capacity; i < new_capacity; ++i) {
            if ((new_data[i] = bam_init1()) == NULL) {
                for (int j = i - 1; j >= bv->capacity; --j) bam_destroy1(new_data[i]);
                return NULL;
            }
        }
        bv->capacity = new_capacity;
        bv->data = new_data;
    }
    return bv->data[bv->size];

}

void bam_vector_destroy(bam_vector_t *bv) {
    for (int i = 0; i < bv->capacity; ++i) bam_destroy1(bv->data[i]);
    free(bv);
}

int sam_parser_comp(const void *_b1, const void *_b2){
    bam1_t *b1 = *(bam1_t **)_b1;
    bam1_t *b2 = *(bam1_t **)_b2;
    void *aux = NULL;
    int ni1 = INT_MAX, ni2 = INT_MAX;
    if ((aux = bam_aux_get(b1, "HI"))) ni1 = bam_aux2i(aux);
    if ((aux = bam_aux_get(b2, "HI"))) ni2 = bam_aux2i(aux);
    if (ni1 != ni2) return ni1 - ni2;
    if ((b1->core.flag & BAM_FSUPPLEMENTARY) != (b2->core.flag & BAM_FSUPPLEMENTARY)) return (b1->core.flag & BAM_FSUPPLEMENTARY) - (b2->core.flag & BAM_FSUPPLEMENTARY);
    return (b2->core.flag & BAM_FREAD1) - (b1->core.flag & BAM_FREAD1);
}

sam_parser_t *sam_parser_open(const char* fn){
    sam_parser_t *p = malloc(sizeof(sam_parser_t));
    p->fn = strdup(fn);
    p->fp = sam_open(fn, "r");
    if (!p->fp) return NULL;
    p->hdr = sam_hdr_read(p->fp);
    if (!p->hdr){free(p->fn); sam_close(p->fp) ;return NULL;}
    bam1_t *b = bam_init1();
    if (!b) {free(p->fn); sam_hdr_destroy(p->hdr);sam_close(p->fp) ;return NULL;}
    int ret = sam_read1(p->fp, p->hdr, b);
    if (ret < 0) {free(p->fn); bam_destroy1(b); sam_hdr_destroy(p->hdr);sam_close(p->fp) ;return NULL;}
    p->b = b;
    return p;
}

int sam_parser_close(sam_parser_t *p) {
    free(p->fn);
    bam_hdr_destroy(p->hdr);
    sam_close(p->fp);
    if (p->b) bam_destroy1(p->b);
    free(p);
    return 0;
}

int sam_parser_next(sam_parser_t *p, bam_vector_t *bv){
    if (!p->b) return 0;
    bam1_t *b, *b1;
    int ret;
    if (!(b1 = bam_vector_next(bv))) return -1;
    bv->data[bv->size] = p->b;
    size_t init_index = bv->size++;
    while ((b = bam_vector_next(bv)) && (ret = sam_read1(p->fp, p->hdr, b)) >= 0){
        if (strcmp(bam_get_qname(b), bam_get_qname(p->b)) == 0) {
            if (b->core.flag & BAM_FSUPPLEMENTARY) continue;
            bv->size++;
        } else {
            p->b = b;
            bv->data[bv->size] = b1;
            break;
        }
    }
    if (!b) return -1;
    if (ret < -1) return -1;
    if (ret == -1) {
        p->b = NULL;
        bam_destroy1(b1);
    }
    if (bv->size - init_index > 1) qsort(bv->data + init_index, bv->size - init_index, sizeof(bam1_t *), sam_parser_comp);
    return bv->size - init_index;
}

static int sam_realloc_bam_data(bam1_t *b, size_t desired) /* from htslib-1.15 */
{
    uint32_t new_m_data;
    uint8_t *new_data;
    new_m_data = desired;
    kroundup32(new_m_data);
    if (new_m_data < desired) {
        errno = ENOMEM; // Not strictly true but we can't store the size
        return -1;
    }
    if ((bam_get_mempolicy(b) & BAM_USER_OWNS_DATA) == 0) {
        new_data = realloc(b->data, new_m_data);
    } else {
        if ((new_data = malloc(new_m_data)) != NULL) {
            if (b->l_data > 0)
                memcpy(new_data, b->data,
                       b->l_data < b->m_data ? b->l_data : b->m_data);
            bam_set_mempolicy(b, bam_get_mempolicy(b) & (~BAM_USER_OWNS_DATA));
        }
    }
    if (!new_data) return -1;
    b->data = new_data;
    b->m_data = new_m_data;
    return 0;
}
int bam_set_cigar(bam1_t *b, uint32_t *new_cigars, uint32_t new_n_cigar){
    int32_t byte_shift = ((int32_t)new_n_cigar - (int32_t)b->core.n_cigar)<<2u;
    if (byte_shift != 0) {
        size_t desired = b->l_data + byte_shift;
        if (desired > b->m_data && (sam_realloc_bam_data(b, desired) < 0)) return -1;
        uint8_t *s = (uint8_t *)bam_get_seq(b);
        memmove(s + byte_shift, s,  b->l_data - (s - b->data));
        b->l_data += byte_shift;
        b->core.n_cigar = new_n_cigar;
    }
    memmove(bam_get_cigar(b), new_cigars, new_n_cigar<<2u);
    return 0;
}