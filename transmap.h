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

#include "vector.h"
#include "transmap_bam.h"
#include "transmap_bed.h"

#define TRANSMAP_VERSION "1.0.0"

#define OPTION_ALLOW_PARTIAL 1u
#define OPTION_ALLOW_END_DEL 2u
#define OPTION_FIX_MD 4u
#define OPTION_FIX_NM 8u
#define OPTION_KEEP_UNMAPPED 16u
#define OPTION_FIX_NH 32u
#define OPTION_REQUIRE_BOTH_MATE 64u

struct transmap_option {
    const char* sam_file;
    const char* bed_file;
    const char *out_file;
    int index_cutoff;
    int show_help;
    int show_version;
    uint64_t others;
};

void transmap_option(struct transmap_option *options, int argc, char *argv[]);
void transmap_usage(const char* msg);
void transmap_version();

#define TRANSMAP_UNMAPPED_ORIGINAL 6
#define TRANSMAP_ONE_MATE_UNMAPPED_ORIGINAL 5
#define TRANSMAP_UNMAPPED_NO_OVERLAP 4
#define TRANSMAP_UNMAPPED_PARTIAL 3
#define TRANSMAP_UNMAPPED_NO_MATCH 2
#define TRANSMAP_MAPPED 0

struct transmap_statistic {
    int n_align_processed;
    int align_statistics[10];
    int n_read_processed;
    int read_statistics[10];
};

int hdrmap(sam_hdr_t *hdr1, sam_hdr_t *hdr, bed_dict_t *bed);
int transmap_single(bam1_t **bam, int count, bed_dict_t *bed_dict, bam_vector_t *r1v, bam_vector_t *r2v, vec_t(bed) *bed_hit, uint8_t **buffer, size_t *buffer_size, struct transmap_statistic *statistics, struct transmap_option *options);
int transmap_paired(bam1_t **bam, int count, bed_dict_t *bed_dict, bam_vector_t *r1v, bam_vector_t *r2v, vec_t(bed) *bed_hit, uint8_t **buffer, size_t *buffer_size, struct transmap_statistic *statistics, struct transmap_option *options);
int transmap(bam1_t *b, bam1_t *b1, bed_t *bed, uint32_t options, uint8_t **buffer, size_t *buffer_size);





#define is_unmap(b) ((b)->core.flag & (uint16_t)BAM_FUNMAP)
#define is_paired(b) ((b)->core.flag & (uint16_t)BAM_FPAIRED)
#define is_read1(b) ((b)->core.flag & (uint16_t)BAM_FREAD1)
#define is_read2(b) ((b)->core.flag & (uint16_t)BAM_FREAD2)
#define is_same_HI(b1, b2) (bam_aux2i(bam_aux_get(b1, "HI")) == bam_aux2i(bam_aux_get(b2, "HI")))

static void set_mate_unmapped(bam1_t *b){
    b->core.flag |= BAM_FMUNMAP;
    b->core.flag &= ~(BAM_FPROPER_PAIR | BAM_FMREVERSE);
    b->core.mtid = -1;
    b->core.mpos = -1;
    b->core.isize = 0;
}

#define max(a, b) (((a) > (b))?(a):(b))
#define min(a, b) (((a) < (b))?(a):(b))

static uint8_t **need_buffer(size_t needed_size, uint8_t **buffer, size_t *buffer_size){
    if (*buffer == NULL || needed_size > *buffer_size){
        void *new_buffer = realloc(*buffer, needed_size);
        if (!new_buffer) return NULL;
        *buffer = new_buffer;
        *buffer_size = needed_size;
    }
    return buffer;
}



