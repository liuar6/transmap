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

typedef struct bam_vector_s{
    bam1_t **data;
    size_t size;
    size_t capacity;
} bam_vector_t;

typedef struct sam_parser_s{
    char *fn;
    samFile *fp;
    sam_hdr_t *hdr;
    bam1_t *b;
} sam_parser_t;

bam_vector_t *bam_vector_init();
bam1_t *bam_vector_next(bam_vector_t *bv);
void bam_vector_destroy(bam_vector_t *bv);

sam_parser_t *sam_parser_open(const char* fn);
int sam_parser_close(sam_parser_t *p);
int sam_parser_next(sam_parser_t *p, bam_vector_t *bv);