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
#include <getopt.h>
#include <stdio.h>
#include <stdint.h>
#include "htslib/sam.h"
#include "bioidx/bioidx.h"
#include "transmap.h"

int main(int argc, char *argv[]) {
    struct transmap_option options;
    struct transmap_statistic statistics;
    memset(&statistics, 0, sizeof(struct transmap_statistic));
    transmap_option(&options, argc, argv);
    if (options.show_help || options.show_version) return 0;

    int ret;
    sam_parser_t *sam = NULL;
    samFile *out = NULL;
    sam_hdr_t *hdr = NULL, *new_hdr = NULL;
    bam1_t **record;
    bam_vector_t *r1v = NULL, *r2v = NULL;
    bam_vector_t *bv = NULL;
    vec_t(bed) *bed_hit = NULL;

    if (!(bv = bam_vector_init())) {ret = 1; goto clean_up;}
    if (!(r1v = bam_vector_init())) {ret = 1; goto clean_up;}
    if (!(r2v = bam_vector_init())) {ret = 1; goto clean_up;}
    if (!(bed_hit = vec_init(bed))) {ret = 1; goto clean_up;}

    bed_dict_t *bed = NULL;
    uint8_t *buffer = NULL;
    size_t buffer_size = 0;
    int ret_val, count;



    new_hdr = sam_hdr_init();
    if (!new_hdr){
        ret = 1;
        goto clean_up;
    }
    if ((sam = sam_parser_open(options.sam_file)) == NULL){
        transmap_usage("[transmap] Error: can not open the input bam file.");
        ret = 1;
        goto clean_up;
    }
    char out_mode[3];
    strncpy(out_mode, "w\0\0", 3);
    if (strcmp(options.out_file + strlen(options.out_file) - 4, ".bam") == 0) out_mode[1] = 'b';
    if ((out = sam_open(options.out_file, out_mode)) == NULL){
        transmap_usage("[transmap] Error: can not open the output bam file.");
        ret = 1;
        goto clean_up;
    }
    if (!(bed = bed_parse(options.bed_file))){
        transmap_usage("[transmap] Error: can not open the bed file.");
        ret = 1;
        goto clean_up;
    }
    hdr = sam->hdr;
    hdrmap(new_hdr, hdr, bed);
    char *s = NULL;
    if (!(s = stringify_argv(argc, argv)) || sam_hdr_add_pg(new_hdr, "transmap", "VN", "0.1", "CL", s, NULL) != 0){
        ret = 1;
        free(s);
        goto clean_up;
    }
    free(s);
    if (sam_hdr_write(out, new_hdr) != 0){
        ret = 1;
        goto clean_up;
    };

    while ((count = sam_parser_next(sam, bv)) > 0){
        record = bv->data + bv->size - count;
        if (is_paired(record[0])){
            ret_val = transmap_paired(record, count, bed, r1v, r2v, bed_hit, &buffer, &buffer_size, &statistics, &options);
        } else ret_val = transmap_single(record, count, bed, r1v, r2v, bed_hit, &buffer, &buffer_size, &statistics, &options);
        if (ret_val != 0) {ret = 1; goto clean_up;}
        if (r1v->size >= 1000) {
            for (int i = 0; i < r1v->size; ++i){
                if (r1v->data[i]->core.tid != -1) if (sam_write1(out, new_hdr, r1v->data[i]) < 0) {ret = 1; goto clean_up;}
                if (r2v->data[i]->core.tid != -1) if (sam_write1(out, new_hdr, r2v->data[i]) < 0) {ret = 1; goto clean_up;}
            }
            r1v->size = 0;
            r2v->size = 0;
        }
        if (bv->size >= 1000) bv->size = 0;
    }
    for (int i = 0; i < r1v->size; ++i){
        if (r1v->data[i]->core.tid != -1) if (sam_write1(out, new_hdr, r1v->data[i]) < 0) {ret = 1; goto clean_up;}
        if (r2v->data[i]->core.tid != -1) if (sam_write1(out, new_hdr, r2v->data[i]) < 0) {ret = 1; goto clean_up;}
    }
    fprintf(stderr, "[Read statistics]\n");
    fprintf(stderr, "Total:\t\t\t\t%d\n", statistics.n_read_processed);
    fprintf(stderr, "Mapped:\t\t\t\t%d\n", statistics.read_statistics[TRANSMAP_MAPPED]);
    fprintf(stderr, "Unmapped original:\t\t%d\n", statistics.read_statistics[TRANSMAP_UNMAPPED_ORIGINAL]);
    if (options.others & OPTION_REQUIRE_BOTH_MATE)     fprintf(stderr, "Unmapped mate:\t\t%d\n", statistics.read_statistics[TRANSMAP_ONE_MATE_UNMAPPED_ORIGINAL]);
    fprintf(stderr, "Unmapped no overlap:\t\t%d\n", statistics.read_statistics[TRANSMAP_UNMAPPED_NO_OVERLAP]);
    if (options.others & OPTION_ALLOW_PARTIAL) fprintf(stderr, "Unmapped no match:\t\t%d\n", statistics.read_statistics[TRANSMAP_UNMAPPED_NO_MATCH]);
    else fprintf(stderr, "Unmapped partial:\t\t%d\n", statistics.read_statistics[TRANSMAP_UNMAPPED_PARTIAL]);
    fprintf(stderr, "\n[Alignment statistics]\n");
    fprintf(stderr, "Total:\t\t\t\t%d\n", statistics.n_align_processed);
    fprintf(stderr, "Mapped:\t\t\t\t%d\n", statistics.align_statistics[TRANSMAP_MAPPED]);
    //fprintf(stderr, "Unmapped original:\t\t%d\n", statistics.align_statistics[TRANSMAP_UNMAPPED_ORIGINAL]);
    if (options.others & OPTION_REQUIRE_BOTH_MATE)     fprintf(stderr, "Unmapped mate:\t\t%d\n", statistics.align_statistics[TRANSMAP_ONE_MATE_UNMAPPED_ORIGINAL]);
    fprintf(stderr, "Unmapped no overlap:\t\t%d\n", statistics.align_statistics[TRANSMAP_UNMAPPED_NO_OVERLAP]);
    if (options.others & OPTION_ALLOW_PARTIAL) fprintf(stderr, "Unmapped no match:\t\t%d\n", statistics.align_statistics[TRANSMAP_UNMAPPED_NO_MATCH]);
    else fprintf(stderr, "Unmapped partial:\t\t%d\n", statistics.align_statistics[TRANSMAP_UNMAPPED_PARTIAL]);
    ret = 0;
    clean_up:
    if (bed_hit) vec_destroy(bed, bed_hit);
    if (bv) bam_vector_destroy(bv);
    if (r1v) bam_vector_destroy(r1v);
    if (r2v) bam_vector_destroy(r2v);
    if (new_hdr) sam_hdr_destroy(new_hdr);
    if (bed) bed_free(bed);
    if (buffer) free(buffer);
    if (sam) sam_parser_close(sam);
    if (out) sam_close(out);
    return ret;
}

void transmap_version(){
    fprintf(stderr, "transmap-%s\n\n", TRANSMAP_VERSION);
}

void transmap_usage(const char* msg){
    const char *usage_info = "\
transmap: convert genomic alignments to transcriptome.\n\
Usage:  transmap [options] --fi <alignment file> --fo <output file> --bed <bed file>\n\
[options]\n\
-i/--fi             : input bam file sorted (or grouped) by read name. [required]\n\
-o/--fo             : output bam file. [required]\n\
-b/--bed            : bed file providing the coordinates of target regions. [required]\n\
-h/--help           : show help informations.\n\
--partial           : also process the bam records not fully (partially) contained by the target regions.\n\
--no-trim           : do not trim the marginal I or N for bam records when --partial is specified.\n\
--both-mate         : require both mate of paired-end alignments to be mapped for reporting.\n\
--fix-NH            : fix the NH and HI tag.\n\
--fix-MD            : fix the MD tag if exists.\n\
--fix-NM            : fix the NM tag when --fix-MD is specified. \n\n";
    if (msg==NULL || msg[0] == '\0') fprintf(stderr, "%s", usage_info);
    else fprintf(stderr, "%s\n\n%s", msg, usage_info);
    exit(1);
}

void transmap_option(struct transmap_option *options, int argc, char *argv[]){
    char c;

    options->sam_file = NULL;
    options->bed_file = NULL;
    options->out_file = "-";
    options->show_help = 0;
    options->show_version = 0;
    options->index_cutoff = 20;
    options->others = 0;
    if (argc == 1) transmap_usage("");
    const char *short_options = "hvo:i:b:PTMN";
    const struct option long_options[] =
            {
                    { "help" , no_argument , NULL, 'h' },
                    { "fi" , required_argument , NULL, 'i' },
                    { "fo" , required_argument, NULL, 'o' },
                    { "bed" , required_argument, NULL, 'b' },
                    { "partial" , no_argument, NULL, 'P' },
                    { "no-trim" , no_argument, NULL, 'T' },
                    { "fix-NH" , no_argument, NULL, 'N' },
                    { "fix-MD" , no_argument, NULL, 'D' },
                    { "fix-NM" , no_argument, NULL, 'M' },
                    { "index" , required_argument, NULL, 'B' },
                    {NULL, 0, NULL, 0} ,
            };

    while ((c = getopt_long(argc, argv, short_options, long_options, NULL)) >= 0)
    {
        switch (c)
        {
            case 'h':
                options->show_help = 1;
                transmap_usage(NULL);
                break;
            case 'v':
                options->show_version = 1;
                transmap_version();
                return;
                break;
            case 'o':
                options->out_file = optarg;
                break;
            case 'i':
                options->sam_file = optarg;
                break;
            case 'b':
                options->bed_file = optarg;
                break;
            case 'P':
                options->others |= OPTION_ALLOW_PARTIAL;
                break;
            case 'T':
                options->others |= OPTION_ALLOW_END_DEL;
                break;
            case 'N':
                options->others |= OPTION_FIX_NH;
                break;
            case 'D':
                options->others |= OPTION_FIX_MD;
                break;
            case 'M':
                options->others |= OPTION_FIX_NM;
                break;
            case 'B':
                options->index_cutoff = strtol(optarg, NULL, 10);
                break;
            default:
                transmap_usage("[transmap] Error:unrecognized parameter");
        }
    }
    if (argc != optind) transmap_usage("[transmap] Error:unrecognized parameter");
};

int fix_NH(bam1_t **b, int size){
    int i;
    for (i = 0; i < size; ++i){
        if (b[i]->core.tid != -1){
            if (bam_aux_update_int(b[i], "NH", size) != 0) return -1;
            if (bam_aux_update_int(b[i], "HI", i + 1) != 0) return -1;
        }
    }
    return 0;
}

int transmap_single(bam1_t **bam, int count, bed_dict_t *bed_dict, bam_vector_t *r1v, bam_vector_t *r2v, vec_t(bed) *bed_hit, uint8_t **buffer, size_t *buffer_size, struct transmap_statistic *statistics, struct transmap_option *options) {
    bam1_t *r1, *t1, *t2;
    bed_t **bed_list;
    int bed_size;
    int read_status, align_status;
    int init_index = r1v->size;
    int i = 0, j = 0;
    int ret;
    read_status = TRANSMAP_UNMAPPED_ORIGINAL;
    statistics->n_read_processed++;
    while (i < count) {
        r1 = bam[i++];
        if (is_unmap(r1)) continue;
        statistics->n_align_processed++;
        align_status = TRANSMAP_UNMAPPED_NO_OVERLAP;
        if (bed_dict->size > options->index_cutoff) {
            bed_search1(bed_dict, r1, bed_hit);
            bed_list = bed_hit->data;
            bed_size = bed_hit->size;
        } else {
            bed_list = bed_dict->record;
            bed_size = bed_dict->size;
        }
        for (j = 0; j < bed_size; ++j) {
            if (!(t1 = bam_vector_next(r1v))) return -1;
            if (!(t2 = bam_vector_next(r2v))) return -1;
            ret = transmap(r1, t1, bed_list[j], options->others, buffer, buffer_size);
            if (ret < 0) return -1;
            align_status = min(align_status, ret);
            if (ret != TRANSMAP_MAPPED) continue;
            t2->core.tid = -1;
            r1v->size++;
            r2v->size++;
        }
        statistics->align_statistics[align_status]++;
        read_status = min(read_status, align_status);
    }
    if (options->others & OPTION_FIX_NH) if (fix_NH(r1v->data + init_index, r1v->size - init_index) != 0) return -1;
    statistics->read_statistics[read_status]++;
    return 0;
}


int transmap_paired(bam1_t **bam, int count, bed_dict_t *bed_dict, bam_vector_t *r1v, bam_vector_t * r2v, vec_t(bed) *bed_hit, uint8_t **buffer, size_t *buffer_size, struct transmap_statistic *statistics, struct transmap_option *options){
    bam1_t *r1, *r2, *t1, *t2;
    bed_t **bed_list;
    int bed_size;
    int read_status, align_status;
    int i = 0, j = 0;
    int init_index = r1v->size;
    int ret, ret1, ret2;
    read_status = TRANSMAP_UNMAPPED_ORIGINAL;
    statistics->n_read_processed++;
    while (i < count) {
        r1 = NULL;
        r2 = NULL;
        if (is_unmap(bam[i])) {i++; continue;}
        if (is_read1(bam[i])) r1 = bam[i++];
        if (i < count && !is_unmap(bam[i]) && is_read2(bam[i]) && (!r1 || is_same_HI(r1, bam[i])))
            r2 = bam[i++];
        if (r1 == NULL && r2 == NULL) continue; /* currently impoosibile case */
        /* finishing this stage indicates an alignment is extracted */
        statistics->n_align_processed++;
        if ((options->others & OPTION_REQUIRE_BOTH_MATE) && (r1 == NULL || r2 == NULL)){
            statistics->align_statistics[TRANSMAP_ONE_MATE_UNMAPPED_ORIGINAL]++;
            read_status = min(read_status, TRANSMAP_ONE_MATE_UNMAPPED_ORIGINAL);
            continue;
        }
        /* finishing previous stage indicates the alignment is at least "no overlap" status */
        align_status = TRANSMAP_UNMAPPED_NO_OVERLAP;
        /* bed_search filter the beds that produce status no better than "no overlap" */
        if (bed_dict->size > options->index_cutoff) {
            bed_search2(bed_dict, r1, r2, bed_hit, (options->others & OPTION_REQUIRE_BOTH_MATE)?2:1);
            bed_list = bed_hit->data;
            bed_size = bed_hit->size;
        } else {
            bed_list = bed_dict->record;
            bed_size = bed_dict->size;
        }
        for (j = 0; j < bed_size; ++j) {
            if (!(t1 = bam_vector_next(r1v))) return -1;
            if (!(t2 = bam_vector_next(r2v))) return -1;
            if (r1) ret1 = transmap(r1, t1, bed_list[j], options->others, buffer, buffer_size);
            else ret1 = TRANSMAP_UNMAPPED_ORIGINAL;
            if (r2) ret2 = transmap(r2, t2, bed_list[j], options->others, buffer, buffer_size);
            else ret2 = TRANSMAP_UNMAPPED_ORIGINAL;
            if (ret1 < 0 || ret2 < 0) return -1;
            if (options->others & OPTION_REQUIRE_BOTH_MATE) ret = max(ret1, ret2);
            else ret = min(ret1, ret2);
            align_status = min(align_status, ret);
            if (ret != TRANSMAP_MAPPED) continue;
            r1v->size++;
            r2v->size++;
            if (ret1 > 0){set_mate_unmapped(t2); t1->core.tid = -1;}
            if (ret2 > 0){set_mate_unmapped(t1); t2->core.tid = -1;}
        }
        statistics->align_statistics[align_status]++;
        read_status = min(read_status, align_status);
    }
    /* fix NH and HI tag */
    if (options->others & OPTION_FIX_NH){
        if (fix_NH(r1v->data + init_index, r1v->size - init_index) != 0) return -1;
        if (fix_NH(r2v->data + init_index, r1v->size - init_index) != 0) return -1;
    }
    statistics->read_statistics[read_status]++;
    /* if (options->others & OPTION_KEEP_UNMAPPED) TODO: remember to reverse the sequence and quality for certain cases */
    return 0;
}

int hdrmap(sam_hdr_t *hdr1, sam_hdr_t *hdr, bed_dict_t *bed){
        int i, j;
        hdr1->n_targets = bed->size;
        hdr1->target_name = malloc(sizeof(char *) * bed->size);
        memset(hdr1->target_name, '\0',sizeof(char *) * bed->size);
        hdr1->target_len = malloc(sizeof(uint32_t) * bed->size);
        if (!hdr1->target_name || !hdr1->target_len) return -1;
        for (i = 0; i < bed->size; ++i){
            bed_t *record = bed->record[i];
            record->new_tid = i;
            record->tid = sam_hdr_name2tid(hdr, record->chrom);
            hdr1->target_name[i] = strdup(record->name);
            if (!hdr1->target_name[i]) return -1;
            hdr1->target_len[i] = record->end - record->start;
            bioidx_insert(bed->idx, record->tid, record->start, record->end, record);
        }
        i = 0;
        const char *hdr_lines = sam_hdr_str(hdr);
        int hdr_size = sam_hdr_length(hdr);
        while (i < hdr_size) {
            j = strchr(hdr_lines + i, '\n') - hdr_lines +1;
            if (strncmp(hdr_lines + i, "@SQ", 3) != 0) {
                if (sam_hdr_add_lines(hdr1, hdr_lines + i, j - i) != 0) return -1;
            }
            i = j;
        }
        return 0;
    }

uint8_t comp_base[] = {15, 8, 4, 15, 2, 15, 15, 15, 1, 15, 15, 15, 15, 15, 15, 15, 15};
void bam_rev_seq(bam1_t *b){
    int32_t l;
    uint8_t *s;
    uint8_t v1, v2;
    l = b->core.l_qseq;
    s = bam_get_seq(b);
    for (int i = 0, j = l - 1; i <= j; ++i, --j){
        v1 = comp_base[bam_seqi(s, i)];
        v2 = comp_base[bam_seqi(s, j)];
        bam_set_seqi(s, i, v2);
        bam_set_seqi(s, j, v1);
    }
}

void bam_rev_qual(bam1_t *b){
    int32_t l;
    uint8_t *q;
    uint8_t v;
    l = b->core.l_qseq;
    q = bam_get_qual(b);
    for (int i = 0, j = l - 1; i < j; ++i, --j){
        v = q[i];
        q[i] = q[j];
        q[j] = v;
    }
}
void bam_rev_cigar(bam1_t *b){
    int32_t l;
    uint32_t *cigar;
    uint32_t v;
    l = b->core.n_cigar;
    cigar = bam_get_cigar(b);
    for (int i = 0, j = l - 1; i < j; ++i, --j){
        v = cigar[i];
        cigar[i] = cigar[j];
        cigar[j] = v;
    }
}

static inline char rev_base(char b){
    switch(b){
        case 'A':
            return 'T';
        case 'T':
            return 'A';
        case 'C':
            return 'G';
        case 'G':
            return 'C';
        case 'a':
            return 't';
        case 't':
            return 'a';
        case 'c':
            return 'g';
        case 'g':
            return 'c';
        default:
            return 'N';
    }
}
static inline void bam_rev_aux_md(uint8_t *md, uint8_t *buffer, size_t md_len){
    off_t index, start, i;
    buffer[md_len] = '\0';
    index = md_len;
    i = 0;
    while (i < md_len){
        if (md[i] == '^') {
            start = i++;
            while (i < md_len && (md[i] < '0' || md[i] > '9') && md[i] != '^') ++i;
            for (int j = start + 1; j < i; ++j) buffer[--index] = rev_base(md[j]);
            buffer[--index] = '^';
        } else if (md[i] >= '0' && md[i] <= '9'){
            start = i++;
            while (i < md_len && md[i] >= '0' && md[i] <= '9') ++i;
            strncpy(buffer + (index -= i - start), md + start, i - start);
        } else {
            buffer[--index] = rev_base(md[i]);
            ++i;
        }
    }
    strcpy(md, buffer);
}

static inline int allow_clip(uint32_t mode, uint32_t cigar_type){
    return (mode & cigar_type) == mode;
}

static inline int sub_md(uint8_t* new_md, uint8_t* md, uint32_t clip, uint32_t len){
    uint8_t *s, *s1, *s2;
    uint32_t count, first_count, last_count, sub_count;
    int del_prefix = 0;

    s1 = md;
    count = 0;
    while(count < clip){
        if (s1[0] >= '1' && s1[0] <= '9'){
            count += strtol((char *)s1, (char **)&s1, 10);
            del_prefix = 0;
        } else if (s1[0] == '0'){
            s1++;
            del_prefix = 0;
        } else if (s1[0] == '^'){
            s1++;
            del_prefix = 1;
        } else {
            s1++;
            count++;
        }
    }
    first_count = count - clip;
    if (s1[0] >= '0' && s1[0] <= '9') del_prefix = 0;
    while (s1[0] == '0') ++s1;

    s2 = s1;
    count = first_count;
    last_count = 0;
    while (count < len){
        if (s2[0] >= '0' && s2[0] <= '9'){
            sub_count = strtol((char *)s2, (char **)&s, 10);
            if (count + sub_count >= len){
                last_count = len - count;
            } else s2 = s;
            count += sub_count;
        } else if (s2[0] == '^') {
            s2++;
        } else {
            s2++;
            count++;
        }
    }
    new_md[0] = '\0';
    if (first_count > len) {
        sprintf(new_md, "%d", len);
    } else {
        s = new_md;
        if (first_count > 0) s += sprintf(s, "%d", first_count);
        else if (s1[0] < '0' || s1[0] > '9') s += sprintf(s, "%d", 0);
        if (del_prefix) s += sprintf(s, "%c", '^');
        strncpy(s, s1, s2 - s1);
        s += s2 - s1;
        s[0] = '\0';
        if (last_count > 0) s += sprintf(s, "%d", last_count);
        else if (s2[-1] < '0' || s2[-1] > '9') s += sprintf(s, "%d", 0);
    }
    return 0;
}

int sub_alignment(bam1_t *b, bam1_t *b1, hts_pos_t start, hts_pos_t end, uint32_t options, uint8_t **buffer, size_t *buffer_size){
    uint32_t n_cigar = b->core.n_cigar;
    uint32_t *cigar = bam_get_cigar(b);
    uint32_t cigar_op = 0, cigar_type = 0, cigar_len = 0;
    uint32_t l_cigar_index = 0, l_cigar_new_len, r_cigar_index = 0, r_cigar_new_len, l_clip = 0, r_clip = 0;
    uint32_t l_md_clip = 0, r_md_clip = 0;
    uint32_t new_n_cigar = 0;
    uint32_t *new_cigar;
    hts_pos_t pos, new_pos;
    hts_pos_t end_pos = bam_endpos(b);
    int clip_mode = (options & OPTION_ALLOW_END_DEL)?2:3;
    int i;

    pos = b->core.pos;
    for (i = 0; i < n_cigar; ++i){
        cigar_op = bam_cigar_op(cigar[i]);
        cigar_type = bam_cigar_type(cigar_op);
        cigar_len = bam_cigar_oplen(cigar[i]);
        if (cigar_type & 1u) l_clip += cigar_len;
        if (cigar_type & 2u) {
            pos += cigar_len;
            if (cigar_op != BAM_CREF_SKIP) l_md_clip += cigar_len;
        }
        if (allow_clip(clip_mode, cigar_type) && pos > start) {
            if (pos - cigar_len >= end) return 1;
            l_cigar_index = i;
            l_cigar_new_len = min(cigar_len, pos - start);
            if (cigar_type & 1u) l_clip -= l_cigar_new_len;
            if (cigar_op != BAM_CREF_SKIP) l_md_clip -= l_cigar_new_len;
            new_pos = pos - l_cigar_new_len;
            break;
        }
    }
    if (i == n_cigar) return 1;

    pos = end_pos;
    for (i = n_cigar - 1; i >= 0; --i){
        cigar_op = bam_cigar_op(cigar[i]);
        cigar_type = bam_cigar_type(cigar_op);
        cigar_len = bam_cigar_oplen(cigar[i]);
        if (cigar_type & 1u) r_clip += cigar_len;
        if (cigar_type & 2u) {
            pos -= cigar_len;
            if (cigar_op != BAM_CREF_SKIP) r_md_clip += cigar_len;
        }
        if (allow_clip(clip_mode, cigar_type) && pos < end){
            if (pos + cigar_len <= start) return 1;
            r_cigar_index = i;
            r_cigar_new_len = min(end - pos, cigar_len);
            if (cigar_type & 1u) r_clip -= r_cigar_new_len;
            if (cigar_op != BAM_CREF_SKIP) r_md_clip -= r_cigar_new_len;
            break;
        }
    }
    if (i < 0) return 1;

    /* find temperary buffer for new CIGAR */
    if (!(need_buffer((b->core.n_cigar + 2) * 4, buffer, buffer_size))) return -1;
    new_cigar = (uint32_t *)*buffer;

    /* generate the new CIGAR array */
    new_n_cigar = 0;
    if (bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP) new_cigar[new_n_cigar++] = cigar[0];
    if (l_clip > 0) new_cigar[new_n_cigar++] = (l_clip << BAM_CIGAR_SHIFT) | BAM_CSOFT_CLIP;
    if (l_cigar_index == r_cigar_index)
        new_cigar[new_n_cigar++] = ((r_cigar_new_len + l_cigar_new_len - bam_cigar_oplen(cigar[l_cigar_index])) << BAM_CIGAR_SHIFT) | bam_cigar_op(cigar[l_cigar_index]);
    else {
        new_cigar[new_n_cigar++] =(l_cigar_new_len << BAM_CIGAR_SHIFT) | bam_cigar_op(cigar[l_cigar_index]);
        for (i = l_cigar_index + 1; i < r_cigar_index; ++i) new_cigar[new_n_cigar++] = cigar[i];
        new_cigar[new_n_cigar++] =(r_cigar_new_len << BAM_CIGAR_SHIFT) | bam_cigar_op(cigar[r_cigar_index]);
    }
    if (r_clip > 0) new_cigar[new_n_cigar++] = (r_clip << BAM_CIGAR_SHIFT) | BAM_CSOFT_CLIP;
    if (bam_cigar_op(cigar[n_cigar - 1]) == BAM_CHARD_CLIP) new_cigar[new_n_cigar++] = cigar[n_cigar - 1];

    /* generate the new bam record */
    if (bam_set1(b1, strlen(bam_get_qname(b)), bam_get_qname(b), b->core.flag, b->core.tid, new_pos, b->core.qual,
                 new_n_cigar, new_cigar, b->core.mtid, b->core.mpos, b->core.isize, b->core.l_qseq, bam_get_qual(b),
                 bam_get_qual(b), bam_get_l_aux(b)) < 0) return -1;
    memcpy(bam_get_seq(b1), bam_get_seq(b), (b->core.l_qseq + 1) / 2); /* copy the query sequence */
    memcpy(bam_get_aux(b1), bam_get_aux(b), bam_get_l_aux(b)); /* copy the originary aux data */
    b1->l_data += bam_get_l_aux(b);

    /* regenerate the MD and NM tag */
    if (options & OPTION_FIX_MD){
        uint8_t *md;
        uint8_t *new_md;
        hts_pos_t md_base = 0;
        int64_t new_nm = 0;
        if ((md = bam_aux_get(b, "MD")) != NULL){
            size_t needed_len = strlen(md);
            if (!(need_buffer(needed_len, buffer, buffer_size))) return -1;
            new_md = *buffer;
            new_md[0] = '\0';
            for (i = 0; i < n_cigar; ++i){
                if ((bam_cigar_type(bam_cigar_op(cigar[i])) & 2u) && bam_cigar_op(cigar[i])!= BAM_CREF_SKIP)
                    md_base += bam_cigar_oplen(cigar[i]);
            }
            sub_md(new_md, md + 1, l_md_clip, md_base - l_md_clip - r_md_clip);
            bam_aux_update_str(b1, "MD", strlen(new_md), new_md);
            if (options & OPTION_FIX_NM) {
                for (i = 0; i < new_n_cigar; ++i){
                    if (bam_cigar_op(bam_get_cigar(b1)[i]) == BAM_CINS) new_nm+= bam_cigar_oplen(bam_get_cigar(b1)[i]);
                }
                uint8_t *s = new_md;
                while(*s != '\0'){
                    if (!(s[0] >= '0' && s[0] <= '9') && s[0] != '^') new_nm++;
                    s++;
                }
                bam_aux_update_int(b1, "NM", new_nm);
            }
        }
    }
    return 0;
}

int transmap(bam1_t *b, bam1_t *b1, bed_t *bed, uint32_t options, uint8_t **buffer, size_t *buffer_size){
    hts_pos_t pos = b->core.pos;
    hts_pos_t end_pos = bam_endpos(b);
    if (b->core.tid != bed->tid) return TRANSMAP_UNMAPPED_NO_OVERLAP;
    if (end_pos <= bed->start || pos >= bed->end) return TRANSMAP_UNMAPPED_NO_OVERLAP;
    if (pos < bed->start || end_pos > bed->end){;
        if (!(options & OPTION_ALLOW_PARTIAL)) return TRANSMAP_UNMAPPED_PARTIAL;
        if (sub_alignment(b, b1, bed->start, bed->end, options, buffer, buffer_size) != 0) return TRANSMAP_UNMAPPED_NO_MATCH;
    } else if (!bam_copy1(b1, b)) return -1;

    b1->core.tid = bed->new_tid;
    if (bed->strand == '-') {
        b1->core.pos = bed->end - bam_endpos(b1);
        b1->core.flag^=16u;
        bam_rev_cigar(b1);
        bam_rev_seq(b1);
        bam_rev_qual(b1);
        uint8_t *md;
        if ((options & OPTION_FIX_MD) && (md = bam_aux_get(b1, "MD")) != NULL) {
            size_t md_len = strlen(++md);
            if (need_buffer(md_len + 1, buffer, buffer_size) == NULL) return -1;
            bam_rev_aux_md(md, *buffer, md_len);
        }
    } else {
        b1->core.pos = b1->core.pos - bed->start;
    }
    return TRANSMAP_MAPPED;
};
