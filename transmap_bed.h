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

#include "bioidx/bioidx.h"
#include "vector.h"
#include "htslib/sam.h"

#ifndef __TRANSCRIPT_BED_H
#define __TRANSCRIPT_BED_H
typedef struct bed_t{
    int32_t tid;
    int32_t new_tid;
    char *chrom;
    hts_pos_t start;
    hts_pos_t end;
    char *name;
    char strand;
} bed_t;

typedef struct bed_dict_t{
    bed_t ** record;
    bioidx_t *idx;
    int64_t size;
    int64_t capacity;
} bed_dict_t;

VEC_INIT(bed, bed_t *)

bed_dict_t *bed_parse(const char* fname);
void bed_free(bed_dict_t *bed);
void bed_search1(bed_dict_t *bed, bam1_t *b, vec_t(bed) *hits);
void bed_search2(bed_dict_t *bed, bam1_t *r1, bam1_t *r2, vec_t(bed) *hits, int mode);
#endif /* __TRANSCRIPT_BED_H */