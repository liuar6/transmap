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
    sam_hdr_t *new_hdr = NULL;
    bam1_t **record;
    bam_vector_t *r1v = NULL, *r2v = NULL;
    bam_vector_t *bv = NULL;
    vec_t(bed) *bed_hit = NULL;
    vec_t(exon) *exon_hit = NULL;

    if (!(bv = bam_vector_init())) {ret = 1; goto clean_up;}
    if (!(r1v = bam_vector_init())) {ret = 1; goto clean_up;}
    if (!(r2v = bam_vector_init())) {ret = 1; goto clean_up;}
    if (!(bed_hit = vec_init(bed))) {ret = 1; goto clean_up;}
    if (!(exon_hit = vec_init(exon))) {ret = 1; goto clean_up;}

    bed_dict_t *bed = NULL;
    gtf_dict_t *gtf = NULL;
    uint8_t *buffer = NULL;
    size_t buffer_size = 0;
    int ret_val, count;

    new_hdr = sam_hdr_init();
    if (!new_hdr){
        ret = 1;
        goto clean_up;
    }
    if ((sam = sam_parser_open(options.sam_file)) == NULL){
        fprintf(stderr, "[transmap] Error: can not open the input bam file.\n");
        ret = 1;
        goto clean_up;
    }
    char out_mode[3];
    strncpy(out_mode, "w\0\0", 3);
    if (strcmp(options.out_file + strlen(options.out_file) - 4, ".bam") == 0) out_mode[1] = 'b';
    if ((out = sam_open(options.out_file, out_mode)) == NULL){
        fprintf(stderr, "[transmap] Error: can not open the output bam file.");
        ret = 1;
        goto clean_up;
    }

    if (options.others & OPTION_GTF_MODE){
        if (!(gtf = gtf_parse(options.in_file, options.gtf_feature, options.gtf_attribute))){
            fprintf(stderr, "[transmap] Error: can not open the gtf file.");
            ret = 1;
            goto clean_up;
        }
        if (!(new_hdr = hdrmap_gtf(sam->hdr, gtf))){
            fprintf(stderr, "[transmap] Error: can not generate the new bam header.");
            ret = 1;
            goto clean_up;
        }
        if (gtf->record->size > options.index_cutoff) options.others |= OPTION_USE_INDEX;
    } else {
        if (!(bed = bed_parse(options.in_file))){
            fprintf(stderr, "[transmap] Error: can not open the bed file.");
            ret = 1;
            goto clean_up;
        }
        if (!(new_hdr = hdrmap_bed(sam->hdr, bed))){
            fprintf(stderr, "[transmap] Error: can not generate the new bam header.");
            ret = 1;
            goto clean_up;
        }
        if (bed->size > options.index_cutoff) options.others |= OPTION_USE_INDEX;
    }

    char *s = NULL;
    if (!(s = stringify_argv(argc, argv)) || sam_hdr_add_pg(new_hdr, "transmap", "VN", "0.1", "CL", s, NULL) != 0){
        fprintf(stderr, "[transmap] Error: can not generate the new bam header.");
        ret = 1;
        if (s) free(s);

        goto clean_up;
    }
    free(s);
    if (sam_hdr_write(out, new_hdr) != 0){
        fprintf(stderr, "[transmap] Error: can not write the bam header.");
        ret = 1;
        goto clean_up;
    };

    while ((count = sam_parser_next(sam, bv)) > 0){
        record = bv->data + bv->size - count;
        if (is_paired(record[0])){
            ret_val = transmap_paired(record, count, gtf, r1v, r2v, exon_hit, &buffer, &buffer_size, &statistics, &options);
        } else ret_val = transmap_single(record, count, gtf, r1v, r2v, exon_hit, &buffer, &buffer_size, &statistics, &options);
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
    fprintf(stderr, "Total:                      %d\n", statistics.n_read_processed);
    fprintf(stderr, "Mapped unique:              %d\n", statistics.read_statistics[TRANSMAP_MAPPED]);
    fprintf(stderr, "Mapped multiple:            %d\n", statistics.read_statistics[TRANSMAP_MULTI_MAPPED]);
    fprintf(stderr, "Unmapped unaligned:         %d\n", statistics.read_statistics[TRANSMAP_UNALIGNED]);
    if (options.others & OPTION_REQUIRE_BOTH_MATE){
        fprintf(stderr, "Unmapped mate unaligned:    %d\n", statistics.read_statistics[TRANSMAP_MATE_UNALIGNED]);
        fprintf(stderr, "Unmapped mate missing:      %d\n", statistics.read_statistics[TRANSMAP_MATE_MISSING]);
        fprintf(stderr, "Unmapped improper pair:     %d\n", statistics.read_statistics[TRANSMAP_PAIR_IMPROPER]);
    }
    fprintf(stderr, "Unmapped no overlap:        %d\n", statistics.read_statistics[TRANSMAP_UNMAPPED_NO_OVERLAP]);
    if (!(options.others & OPTION_ALLOW_PARTIAL)) fprintf(stderr, "Unmapped partial:           %d\n", statistics.read_statistics[TRANSMAP_UNMAPPED_PARTIAL]);
    if (options.others & OPTION_GTF_MODE) fprintf(stderr, "Unmapped exon imcompatible: %d\n", statistics.read_statistics[TRANSMAP_EXON_IMCOMPATIBLE]);
    if ((options.others & OPTION_ALLOW_PARTIAL && !(options.others & OPTION_GTF_MODE)) || options.others & OPTION_IRREGULAR) fprintf(stderr, "Unmapped no match:          %d\n", statistics.read_statistics[TRANSMAP_UNMAPPED_NO_MATCH]);

    fprintf(stderr, "\n[Alignment statistics]\n");
    fprintf(stderr, "Total:                      %d\n", statistics.n_align_processed);
    fprintf(stderr, "Mapped unique:              %d\n", statistics.align_statistics[TRANSMAP_MAPPED]);
    fprintf(stderr, "Mapped multiple:            %d\n", statistics.align_statistics[TRANSMAP_MULTI_MAPPED]);
    if (options.others & OPTION_REQUIRE_BOTH_MATE){
        fprintf(stderr, "Unmapped mate unaligned:    %d\n", statistics.align_statistics[TRANSMAP_MATE_UNALIGNED]);
        fprintf(stderr, "Unmapped mate missing:      %d\n", statistics.align_statistics[TRANSMAP_MATE_MISSING]);
        fprintf(stderr, "Unmapped improper pair:     %d\n", statistics.align_statistics[TRANSMAP_PAIR_IMPROPER]);
    }
    fprintf(stderr, "Unmapped no overlap:        %d\n", statistics.align_statistics[TRANSMAP_UNMAPPED_NO_OVERLAP]);
    if (!(options.others & OPTION_ALLOW_PARTIAL)) fprintf(stderr, "Unmapped partial:           %d\n", statistics.align_statistics[TRANSMAP_UNMAPPED_PARTIAL]);
    if (options.others & OPTION_GTF_MODE) fprintf(stderr, "Unmapped exon imcompatible: %d\n", statistics.align_statistics[TRANSMAP_EXON_IMCOMPATIBLE]);
    if ((options.others & OPTION_ALLOW_PARTIAL && !(options.others & OPTION_GTF_MODE)) || options.others & OPTION_IRREGULAR) fprintf(stderr, "Unmapped no match:          %d\n", statistics.align_statistics[TRANSMAP_UNMAPPED_NO_MATCH]);

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
    exit(0);
}

void transmap_usage(const char* msg){
    const char *usage_info = "\
transmap: convert genomic alignments to transcriptome.\n\
Usage:  transmap [options] --fi <alignment file> --fo <output file> --bed <bed file>\n\
[options]\n\
-i/--fi             : input bam file sorted (or grouped) by read name. [required]\n\
-o/--fo             : output bam file. [required]\n\
-b/--bed            : bed file providing the coordinates of target regions. [required]\n\
-g/--gtf            : gtf file providing the coordinates of target transcripts. [required]\n\
--gtf-feature       : gtf feature used to define regions of transcripts. default: exon.\n\
--gtf-attribute     : gtf attribute used as the name of transcripts. default: transcript_id.\n\
--partial           : also process the bam records not fully (partially) contained by the target regions.\n\
--no-trim           : do not trim the marginal I or N for bam records when --partial is specified.\n\
--both-mate         : require both mate of paired-end alignments to be mapped for reporting.\n\
--fix-NH            : fix the NH and HI tag.\n\
--fix-MD            : fix the MD tag if exists.\n\
--fix-NM            : fix the NM tag when --fix-MD is specified. \n\
-h                  : show help informations.\n\n";
    if (msg==NULL || msg[0] == '\0') fprintf(stderr, "%s", usage_info);
    else fprintf(stderr, "%s\n\n%s", msg, usage_info);
    exit(1);
}

void transmap_option(struct transmap_option *options, int argc, char *argv[]){
    char c;
    options->sam_file = NULL;
    options->in_file = NULL;
    options->out_file = "-";
    options->gtf_feature = "exon";
    options->gtf_attribute = "transcript_id";
    options->show_help = 0;
    options->show_version = 0;
    options->index_cutoff = 0;
    options->others = 0;
    if (argc == 1) transmap_usage("");
    const char *short_options = "hvo:i:b:g:F:A:OPTNDMIB:";
    const struct option long_options[] =
            {
                    { "help" , no_argument , NULL, 'h' },
                    { "version" , no_argument , NULL, 'v' },
                    { "fi" , required_argument , NULL, 'i' },
                    { "fo" , required_argument, NULL, 'o' },
                    { "bed" , required_argument, NULL, 'b' },
                    { "gtf" , required_argument, NULL, 'g' },
                    { "gtf-feature" , required_argument, NULL, 'F' },
                    { "gtf-attribute" , required_argument, NULL, 'A' },
                    { "both-mate" , no_argument, NULL, 'O' },
                    { "partial" , no_argument, NULL, 'P' },
                    { "no-trim" , no_argument, NULL, 'T' },
                    { "fix-NH" , no_argument, NULL, 'N' },
                    { "fix-MD" , no_argument, NULL, 'D' },
                    { "fix-NM" , no_argument, NULL, 'M' },
                    { "irregular" , no_argument, NULL, 'I' },
                    { "index-cutoff" , required_argument, NULL, 'B' },
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
                break;
            case 'o':
                options->out_file = optarg;
                break;
            case 'i':
                options->sam_file = optarg;
                break;
            case 'b':
                options->in_file = optarg;
                options->others |= OPTION_BED_MODE;
                break;
            case 'g':
                options->in_file = optarg;
                options->others |= OPTION_GTF_MODE;
                break;
            case 'F':
                options->gtf_feature = optarg;
                break;
            case 'A':
                options->gtf_attribute = optarg;
                break;
            case 'O':
                options->others |= OPTION_REQUIRE_BOTH_MATE;
                break;
            case 'P':
                options->others |= OPTION_ALLOW_PARTIAL;
                break;
            case 'T':
                options->others |= OPTION_NO_POLISH;
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
            case 'I':
                options->others |= OPTION_IRREGULAR;
                break;
            case 'B':
                options->index_cutoff = strtol(optarg, NULL, 10);
                break;
            default:
                transmap_usage("[transmap] Error:unrecognized parameter");
        }
    }
    if (argc != optind) transmap_usage("[transmap] Error:unrecognized parameter");
    if ((options->others & OPTION_BED_MODE) && (options->others & OPTION_GTF_MODE))
        transmap_usage("[transmap] Error: you can only provide one of --bed or --gtf.");
    if (options->in_file == NULL) transmap_usage("[transmap] Error: you should specify either --bed or --gtf.");
    if (options->sam_file == NULL) transmap_usage("[transmap] Error: you should provide the input bam file via --bam.");
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


int transmap_single(bam1_t **bam, int count, void *dict, bam_vector_t *r1v, bam_vector_t *r2v, void *candidate, uint8_t **buffer, size_t *buffer_size, struct transmap_statistic *statistics, struct transmap_option *options) {
    bam1_t *r1, *t1, *t2;
    uint64_t others = options->others;
    int read_status, align_status;
    int init_index = r1v->size;
    int i = 0, j = 0, align_n_mapped = 0, read_n_mapped = 0;
    int cand_size;
    int ret;
    read_status = TRANSMAP_UNALIGNED;
    statistics->n_read_processed++;
    while (i < count) {
        r1 = bam[i++];
        if (is_unmap(r1)) continue;
        statistics->n_align_processed++;
        align_status = TRANSMAP_UNMAPPED_NO_OVERLAP;
        if (others & OPTION_GTF_MODE) {
            if (others & OPTION_USE_INDEX) gtf_search_one(((gtf_dict_t *)dict)->idx, r1->core.tid, r1->core.pos, bam_endpos(r1), (vec_t(exon) *)candidate);
            cand_size = ((vec_t(exon) *)candidate)->size;
        } else {
            if (others & OPTION_USE_INDEX) bed_search_one(((bed_dict_t *)dict)->idx, r1->core.tid, r1->core.pos, bam_endpos(r1), (vec_t(bed) *)candidate);
            cand_size = ((vec_t(bed) *)candidate)->size;
        }
        align_n_mapped = 0;
        for (j = 0; j < cand_size; ++j) {
            if (!(t1 = bam_vector_next(r1v))) return -1;
            if (!(t2 = bam_vector_next(r2v))) return -1;
            if (others & OPTION_GTF_MODE) ret = transmap_gtf(r1, t1, ((vec_t(exon) *)candidate)->data[j], others, buffer, buffer_size);
            else  ret = transmap_bed(r1, t1, ((vec_t(bed) *)candidate)->data[j], others, buffer, buffer_size);
            if (ret < 0) return -1;
            align_status = min(align_status, ret);
            if (ret != TRANSMAP_MAPPED) continue;
            align_n_mapped++;
            t2->core.tid = -1;
            r1v->size++;
            r2v->size++;
        }
        if (align_n_mapped > 1) statistics->align_statistics[TRANSMAP_MULTI_MAPPED]++;
        else statistics->align_statistics[align_status]++;
        read_n_mapped += align_n_mapped;
        read_status = min(read_status, align_status);
    }
    if (others & OPTION_FIX_NH) if (fix_NH(r1v->data + init_index, r1v->size - init_index) != 0) return -1;
    if (read_n_mapped > 1) statistics->read_statistics[TRANSMAP_MULTI_MAPPED]++;
    else statistics->read_statistics[read_status]++;
    return 0;
}


int transmap_paired(bam1_t **bam, int count, void *dict, bam_vector_t *r1v, bam_vector_t * r2v, void *candidate, uint8_t **buffer, size_t *buffer_size, struct transmap_statistic *statistics, struct transmap_option *options){
    bam1_t *r1, *r2, *t1, *t2;
    uint64_t others = options->others;
    int read_status, align_status;
    int i = 0, j = 0, align_n_mapped = 0, read_n_mapped = 0;
    int init_index = r1v->size;
    int ret, ret1, ret2;
    int cand_size;
    read_status = TRANSMAP_UNALIGNED;
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
        if (others & OPTION_REQUIRE_BOTH_MATE) {
            if (r1 == NULL) {
                if (r2->core.flag & (uint16_t) BAM_FMUNMAP) align_status = TRANSMAP_MATE_UNALIGNED;
                else align_status = TRANSMAP_MATE_MISSING;
                statistics->align_statistics[align_status]++;
                read_status = min(read_status, align_status);
                continue;
            } else if (r2 == NULL) {
                if (r1->core.flag & (uint16_t) BAM_FMUNMAP) align_status = TRANSMAP_MATE_UNALIGNED;
                else align_status = TRANSMAP_MATE_MISSING;
                statistics->align_statistics[align_status]++;
                read_status = min(read_status, align_status);
                continue;
            } else if (!(r1->core.flag & (uint16_t)BAM_FPROPER_PAIR) || r1->core.tid != r2->core.tid) {
                align_status = TRANSMAP_PAIR_IMPROPER;
                read_status = min(read_status, align_status);
                continue;
            }
        }
        /* finishing previous stage indicates the alignment is at least "no overlap" status */
        align_status = TRANSMAP_UNMAPPED_NO_OVERLAP;
        if (others & OPTION_GTF_MODE) {
            if (others & OPTION_USE_INDEX){
                if (others & OPTION_REQUIRE_BOTH_MATE)
                    gtf_search_both(((gtf_dict_t *)dict)->idx, r1->core.tid, r1->core.pos, bam_endpos(r1), r2->core.pos, bam_endpos(r2), (vec_t(exon) *)candidate);
                else {
                    if (!r1) gtf_search_one(((gtf_dict_t *)dict)->idx, r2->core.tid, r2->core.pos, bam_endpos(r2), (vec_t(exon) *)candidate);
                    else if (!r2) gtf_search_one(((gtf_dict_t *)dict)->idx, r1->core.tid, r1->core.pos, bam_endpos(r1), (vec_t(exon) *)candidate);
                    else gtf_search_any(((gtf_dict_t *)dict)->idx, r1->core.tid, r1->core.pos, bam_endpos(r1), r2->core.pos, bam_endpos(r2), (vec_t(exon) *)candidate);
                }

            }
            cand_size = ((vec_t(exon) *)candidate)->size;
        } else {
            if (others & OPTION_USE_INDEX){
                if (others & OPTION_REQUIRE_BOTH_MATE)
                    bed_search_both(((bed_dict_t *)dict)->idx, r1->core.tid, r1->core.pos, bam_endpos(r1), r2->core.pos, bam_endpos(r2), (vec_t(bed) *)candidate);
                else {
                    if (!r1) bed_search_one(((bed_dict_t *)dict)->idx, r2->core.tid, r2->core.pos, bam_endpos(r2), (vec_t(bed) *)candidate);
                    else if (!r2) bed_search_one(((bed_dict_t *)dict)->idx, r1->core.tid, r1->core.pos, bam_endpos(r1), (vec_t(bed) *)candidate);
                    else bed_search_any(((bed_dict_t *)dict)->idx, r1->core.tid, r1->core.pos, bam_endpos(r1), r2->core.pos, bam_endpos(r2), (vec_t(bed) *)candidate);
                }
            }
            cand_size = ((vec_t(bed) *)candidate)->size;
        }
        align_n_mapped = 0;
        for (j = 0; j < cand_size; ++j) {
            if (!(t1 = bam_vector_next(r1v))) return -1;
            if (!(t2 = bam_vector_next(r2v))) return -1;
            ret1 = TRANSMAP_UNALIGNED;
            ret2 = TRANSMAP_UNALIGNED;
            if (others & OPTION_GTF_MODE) {
                exon_t *hit =  ((vec_t(exon) *)candidate)->data[j];
                if (r1) ret1 = transmap_gtf(r1, t1, hit, others, buffer, buffer_size);
                if (r2) ret2 = transmap_gtf(r2, t2, hit, others, buffer, buffer_size);
            } else {
                bed_t *hit =  ((vec_t(bed) *)candidate)->data[j];
                if (r1) ret1 = transmap_bed(r1, t1, hit, others, buffer, buffer_size);
                if (r2) ret2 = transmap_bed(r2, t2, hit, others, buffer, buffer_size);
            }
            if (ret1 < 0 || ret2 < 0) return -1;
            if (others & OPTION_REQUIRE_BOTH_MATE) ret = max(ret1, ret2);
            else ret = min(ret1, ret2);
            align_status = min(align_status, ret);
            if (ret != TRANSMAP_MAPPED) continue;
            align_n_mapped++;
            r1v->size++;
            r2v->size++;
            if (ret1 == TRANSMAP_MAPPED && ret2 == TRANSMAP_MAPPED) {fix_mate(t1, t2);}
            else {
                if (ret1 != TRANSMAP_MAPPED){set_mate_unmapped(t2); t1->core.tid = -1;}
                if (ret2 != TRANSMAP_MAPPED){set_mate_unmapped(t1); t2->core.tid = -1;}
            }
        }
        if (align_n_mapped > 1) statistics->align_statistics[TRANSMAP_MULTI_MAPPED]++;
        else statistics->align_statistics[align_status]++;
        read_n_mapped += align_n_mapped;
        read_status = min(read_status, align_status);
    }
    /* fix NH and HI tag */
    if (others & OPTION_FIX_NH){
        if (fix_NH(r1v->data + init_index, read_n_mapped) != 0) return -1;
        if (fix_NH(r2v->data + init_index, read_n_mapped) != 0) return -1;
    }
    if (read_n_mapped > 1) statistics->read_statistics[TRANSMAP_MULTI_MAPPED]++;
    else statistics->read_statistics[read_status]++;
    /* if (options->others & OPTION_KEEP_UNMAPPED) TODO: remember to reverse the sequence and quality for certain cases */
    return 0;
}

sam_hdr_t *hdrmap_bed(sam_hdr_t *hdr, bed_dict_t *bed){
    int i, j;
    sam_hdr_t *new_hdr;
    if (!(new_hdr = sam_hdr_init())) return NULL;
    new_hdr->n_targets = bed->size;
    new_hdr->target_name = calloc(bed->size, sizeof(char *));
    if (!new_hdr->target_name) goto clean_up;
    new_hdr->target_len = calloc(bed->size, sizeof(uint32_t));
    if (!new_hdr->target_len) goto clean_up;
    for (i = 0; i < bed->size; ++i){
        bed_t *record = bed->record[i];
        record->tid = sam_hdr_name2tid(hdr, record->chrom);
        if (!(new_hdr->target_name[i] = strdup(record->name))) goto clean_up;
        if (record->tid < 0) continue;
        new_hdr->target_len[i] = record->end - record->start;
        bioidx_insert(bed->idx, record->tid, record->start, record->end, record);
    }
    i = 0;
    const char *hdr_lines = sam_hdr_str(hdr);
    int hdr_size = sam_hdr_length(hdr);
    while (i < hdr_size) {
        j = strchr(hdr_lines + i, '\n') - hdr_lines +1;
        if (strncmp(hdr_lines + i, "@SQ", 3) != 0) {
            if (sam_hdr_add_lines(new_hdr, hdr_lines + i, j - i) != 0) goto clean_up;
        }
        i = j;
    }
    return new_hdr;

    clean_up:
    if (new_hdr->target_name){
        for (i = 0; i < new_hdr->n_targets; ++i)
            if (new_hdr->target_name[i]) free(new_hdr->target_name[i]);
        free(new_hdr->target_name);
        new_hdr->target_name = NULL;
    }
    if (new_hdr->target_len) free(new_hdr->target_len);
    free(new_hdr);
    return NULL;
}

sam_hdr_t *hdrmap_gtf(sam_hdr_t *hdr, gtf_dict_t *gtf){
    khash_t (transcript) *record = gtf->record;
    vec_t(exon) *exons;
    transcript_t *tr;
    khiter_t k;
    int  i, j;
    sam_hdr_t *new_hdr;
    if (!(new_hdr = sam_hdr_init())) return NULL;
    new_hdr->n_targets = kh_size(record);
    new_hdr->target_name = calloc(new_hdr->n_targets, sizeof(char *));
    if (!new_hdr->target_name) goto clean_up;
    new_hdr->target_len = malloc(sizeof(uint32_t) * kh_size(record));
    if (!new_hdr->target_len) goto clean_up;
    i = 0;
    for (k = 0; k < kh_end(record); ++k){
        if (!kh_exist(record, k)) continue;
        tr = kh_val(record, k);
        tr->tid = sam_hdr_name2tid(hdr, tr->exons->data[0]->chrom);
        if (!(new_hdr->target_name[tr->new_tid] = strdup(tr->name))) goto clean_up;
        new_hdr->target_len[i++] = tr->len;
        if (tr->tid < 0) continue;
        exons = tr->exons;
        for (j = 0; j < exons->size; ++j){
            if (bioidx_insert(gtf->idx, tr->tid, exons->data[j]->start, exons->data[j]->end, exons->data[j]) != 0)
                goto clean_up;
        }
    }
    i = 0;
    const char *hdr_lines = sam_hdr_str(hdr);
    int hdr_size = sam_hdr_length(hdr);
    while (i < hdr_size) {
        j = strchr(hdr_lines + i, '\n') - hdr_lines +1;
        if (strncmp(hdr_lines + i, "@SQ", 3) != 0) {
            if (sam_hdr_add_lines(new_hdr, hdr_lines + i, j - i) != 0) goto clean_up;
        }
        i = j;
    }
    return new_hdr;

    clean_up:
    if (new_hdr->target_name){
        for (i = 0; i < new_hdr->n_targets; ++i)
            if (new_hdr->target_name[i]) free(new_hdr->target_name[i]);
        free(new_hdr->target_name);
        new_hdr->target_name = NULL;
    }
    if (new_hdr->target_len) free(new_hdr->target_len);
    free(new_hdr);
    return NULL;
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
void trim_cigar(hts_pos_t start, hts_pos_t end, hts_pos_t *_pos, hts_pos_t *_end_pos, uint32_t *cigar, uint32_t n_cigar, uint32_t *new_cigar, uint32_t *_new_n_cigar, uint32_t md_clip[4], uint32_t options){
    uint32_t cigar_op = 0, cigar_type = 0, cigar_len = 0;
    uint32_t l_cigar_index = 0, l_cigar_new_len, r_cigar_index = 0, r_cigar_new_len, l_clip = 0, r_clip = 0;
    hts_pos_t pos, end_pos;
    uint8_t trim_mode = (options & OPTION_NO_POLISH)?2:3;
    int i;

    pos = *_pos;
    for (i = 0; i < n_cigar && pos < end; ++i){
        cigar_op = bam_cigar_op(cigar[i]);
        cigar_len = bam_cigar_oplen(cigar[i]);
        cigar_type = bam_cigar_type(cigar_op);
        if (((trim_mode & cigar_type) == trim_mode) && (end_pos = pos + cigar_len) > start) {
            *_pos = max(pos, start);
            l_cigar_index = i;
            l_cigar_new_len = end_pos - *_pos;
            if (cigar_type & 1u) l_clip += cigar_len - l_cigar_new_len;

            break;
        }
        if (cigar_type & 1u) l_clip += cigar_len;
        if (cigar_type & 2u) pos += cigar_len;
    }
    if (i == n_cigar || pos >= end) {*_new_n_cigar = 0; return;}

    end_pos = *_end_pos;
    for (i = n_cigar - 1; i >= 0 ; --i){
        cigar_op = bam_cigar_op(cigar[i]);
        cigar_type = bam_cigar_type(cigar_op);
        cigar_len = bam_cigar_oplen(cigar[i]);
        if (((trim_mode & cigar_type) == trim_mode) && (pos = end_pos - cigar_len)  < end){
            *_end_pos = min(end_pos, end);
            r_cigar_index = i;
            r_cigar_new_len = *_end_pos - pos;
            if (cigar_type & 1u) r_clip += cigar_len - r_cigar_new_len;
            break;
        }
        if (cigar_type & 1u) r_clip += cigar_len;
        if (cigar_type & 2u) end_pos -= cigar_len;
    }

    /* generate the new CIGAR array */
    uint32_t new_n_cigar = 0;
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
    *_new_n_cigar = new_n_cigar;

    /* prepare information for fixing MD and NM field */
    if (options & OPTION_FIX_MD) {
        memset(md_clip, 0, sizeof(uint32_t) << 2u);
        for (i = 0; i < n_cigar; ++i){
            cigar_op = bam_cigar_op(cigar[i]);
            cigar_type = bam_cigar_type(cigar_op);
            cigar_len = bam_cigar_oplen(cigar[i]);
            if (cigar_op != BAM_CREF_SKIP && (cigar_type & 2u)) {
                md_clip[2] += cigar_len;
                if (i <= l_cigar_index) md_clip[0] += cigar_len;
                if (i >= r_cigar_index) md_clip[1] += cigar_len;
            }
            if ((options & OPTION_FIX_NM) && cigar_op == BAM_CINS && i > l_cigar_index && i < r_cigar_index) md_clip[3] += cigar_len;
        }
        cigar_op = bam_cigar_op(cigar[l_cigar_index]);
        cigar_type = bam_cigar_type(cigar_op);
        if (cigar_op != BAM_CREF_SKIP && (cigar_type & 2u)) md_clip[0] -= l_cigar_new_len;
        cigar_op = bam_cigar_op(cigar[r_cigar_index]);
        cigar_type = bam_cigar_type(cigar_op);
        if (cigar_op != BAM_CREF_SKIP && (cigar_type & 2u)) md_clip[1] -= r_cigar_new_len;
    }
}

void stitch_cigar(hts_pos_t *_pos, uint32_t *cigar, uint32_t n_cigar, uint32_t *_new_n_cigar, int *need_stitch_md){
    uint32_t new_n_cigar = 0; /* refer to current index, not count */
    int i;
    if (bam_cigar_op(cigar[0]) == BAM_CREF_SKIP) {/* TODO: This is an irregular case */
        for (i = 0; i < n_cigar && bam_cigar_op(cigar[i]) == BAM_CREF_SKIP; ++i)
            *_pos += bam_cigar_oplen(cigar[i]);
        if (i == n_cigar) { *_new_n_cigar = 0; return;}
        cigar[0] = cigar[i++];
        for ( ;i < n_cigar; ++i) {
            if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) continue;
            if (bam_cigar_op(cigar[i]) == bam_cigar_op(cigar[new_n_cigar])){
                cigar[new_n_cigar] += bam_cigar_oplen(cigar[i]) << BAM_CIGAR_SHIFT;
                if (bam_cigar_op(cigar[i]) == BAM_CDEL) *need_stitch_md = 1; /* TODO: This is an irregular case */
            } else cigar[++new_n_cigar] = cigar[i];
        }
        new_n_cigar++;
    } else {
        for (i = 1 ;i < n_cigar; ++i) {
            if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) continue;
            if (bam_cigar_op(cigar[i]) == bam_cigar_op(cigar[new_n_cigar])){
                cigar[new_n_cigar] += bam_cigar_oplen(cigar[i]) << BAM_CIGAR_SHIFT;
                if (bam_cigar_op(cigar[i]) == BAM_CDEL) *need_stitch_md = 1; /* TODO: This is an irregular case */
            } else cigar[++new_n_cigar] = cigar[i];
        }
        new_n_cigar++;
    }
    *_new_n_cigar = new_n_cigar;
}

void stitch_MD(uint8_t *md){
    uint8_t *s1, *s2;
    int del_prefix = 0;
    for (s1 = s2 = md; s1[0]; ++s1, ++s2) {
        if (s1[0] == '^') del_prefix = 1;
        if (s1[0] >= '0' && s1[0] <= '9') {
            if (del_prefix && s1[0] == '0' && s1[1] == '^') s1 += 2;
            else del_prefix = 0;
        }
        s2[0] = s1[0];
    }
    s2[0] = '\0';
}

int fix_MD(bam1_t *b, uint8_t **buffer, size_t *buffer_size, uint32_t md_clip[4], int stitch_md, int fix_nm){
    if (md_clip[0] == 0 && md_clip[1] == 0) return 0;
    uint8_t *md, *new_md;
    int64_t new_nm = 0;
    if ((md = bam_aux_get(b, "MD")) != NULL){
        size_t needed_len = strlen((char *)md) + 1;
        if (!(need_buffer(needed_len, buffer, buffer_size))) return -1;
        new_md = *buffer;
        new_md[0] = '\0';
        sub_md(new_md, md + 1, md_clip[0], md_clip[2] - md_clip[0] - md_clip[1]);
        if (stitch_md) stitch_MD(new_md); /* TODO: This is an irregular case */
        bam_aux_update_str(b, "MD", strlen(new_md), new_md);
        if (fix_nm) {
            new_nm = md_clip[3];
            uint8_t *s = new_md;
            while(*s != '\0'){
                if ((s[0] < '0' || s[0] > '9') && s[0] != '^') new_nm++;
                s++;
            }
            bam_aux_update_int(b, "NM", new_nm);
        }
    }
    return 0;
}

int check_exon_compatible(int32_t pos, int32_t end_pos, const uint32_t *cigars, int32_t n_cigar, exon_t *exon){
    int32_t block_start, block_end = pos;
    transcript_t *tr = exon->tr;
    exon_t **exons = tr->exons->data;
    int exon_count = tr->exons->size;
    int exon_index = - 1;
    int i = 0;
    int pass = 1;
    do{
        block_start = block_end;
        if (bam_cigar_op(cigars[i]) == BAM_CREF_SKIP) block_start += bam_cigar_oplen(cigars[i++]);
        for (block_end = block_start ; i < n_cigar && bam_cigar_op(cigars[i]) != BAM_CREF_SKIP; ++i)
            if (bam_cigar_type(bam_cigar_op(cigars[i])) & 2) block_end += bam_cigar_oplen(cigars[i]);
        if (block_start == block_end){ /* TODO: This is an irregular case */
            if (block_start == pos && (bam_cigar_op(cigars[i - 1]) == BAM_CSOFT_CLIP || bam_cigar_op(cigars[i - 1]) == BAM_CHARD_CLIP)) continue;
            if (block_end == end_pos) {
                if (bam_cigar_op(cigars[n_cigar - 1]) == BAM_CREF_SKIP) break;
                int j;
                for (j = n_cigar - 2; bam_cigar_op(cigars[j]) != BAM_CREF_SKIP; --j);
                if ((bam_cigar_op(cigars[j + 1]) == BAM_CSOFT_CLIP || bam_cigar_op(cigars[j + 1]) == BAM_CHARD_CLIP)) break;
            }
        }
        /* Here an alignment block is extracted */
        if (block_end <= tr->start) continue;
        if (exon_index == -1) {
            for (exon_index = exon->idx; exon_index < exon_count && exons[exon_index]->end <= block_start; exon_index++);
            if (exon_index == exon_count) {pass = 0; break;}
        }
        exon = exons[exon_index++]; /* note exon index is plus by one here */
        if ((exon->start < block_start && block_start != pos) ||
            (exon->end > block_end && block_end != end_pos) ||
            (exon->start > block_start && exon_index != 1) ||
            (exon->end < block_end && exon_index != exon_count)){
            pass = 0;
            break;
        }
    } while (i < n_cigar && exon_index < exon_count);
    if (exon_index == -1) pass = 0; /* TODO: This is an irregular case */
    return pass;
}


int transmap_bed(bam1_t *b, bam1_t *b1, bed_t *bed, uint32_t options, uint8_t **buffer, size_t *buffer_size){
    hts_pos_t pos = b->core.pos;
    hts_pos_t end_pos = bam_endpos(b);
    uint32_t *new_cigar;
    uint32_t new_n_cigar;
    uint32_t md_clip[4] = {0, 0, 0, 0};
    if (b->core.tid != bed->tid || end_pos <= bed->start || pos >= bed->end) return TRANSMAP_UNMAPPED_NO_OVERLAP;
    if (!(options & OPTION_ALLOW_PARTIAL) && (pos < bed->start || end_pos > bed->end)) return TRANSMAP_UNMAPPED_PARTIAL;
    if (!bam_copy1(b1, b)) return -1;
    if (pos < bed->start || end_pos > bed->end || ((options & OPTION_IRREGULAR) && !(options & OPTION_NO_POLISH))){
        new_cigar = (uint32_t *) need_buffer((b->core.n_cigar << 2u) + (2u << 2u), buffer, buffer_size);
        trim_cigar(bed->start, bed->end, &pos, &end_pos, bam_get_cigar(b), b->core.n_cigar, new_cigar, &new_n_cigar, md_clip, options);
        b1->core.pos = pos;
        if (new_n_cigar == 0) return TRANSMAP_UNMAPPED_NO_OVERLAP;
        if (bam_set_cigar(b1, new_cigar, new_n_cigar) < 0) return -1;
        if (options & OPTION_FIX_MD) if (fix_MD(b1, buffer, buffer_size, md_clip, 0, options & OPTION_FIX_NM) < 0) return -1;
    }

    b1->core.tid = bed->new_tid;
    b1->core.pos = b1->core.pos - bed->start;
    if (bed->strand == '-') {
        b1->core.pos = bed->end - bed->start - bam_endpos(b1);
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
    }
    return TRANSMAP_MAPPED;
}

int transmap_gtf(bam1_t *b, bam1_t *b1, exon_t *exon, uint32_t options, uint8_t **buffer, size_t *buffer_size){
    transcript_t *tr = exon->tr;
    hts_pos_t pos = b->core.pos;
    hts_pos_t end_pos = bam_endpos(b);
    uint32_t new_n_cigar;
    int need_stitch_md;
    uint32_t md_clip[4] = {0, 0, 0, 0};
    uint32_t *new_cigar;
    if (b->core.tid != tr->tid || end_pos <= tr->start || pos >= tr->end) return TRANSMAP_UNMAPPED_NO_OVERLAP;
    if (!(options & OPTION_ALLOW_PARTIAL) && (pos < tr->start || end_pos > tr->end)) return TRANSMAP_UNMAPPED_PARTIAL;
    if (!check_exon_compatible(pos, end_pos, bam_get_cigar(b), b->core.n_cigar, exon)) return TRANSMAP_EXON_IMCOMPATIBLE;
    if (!bam_copy1(b1, b)) return -1;
    if (pos < tr->start || end_pos > tr->end || ((options & OPTION_IRREGULAR) && !(options & OPTION_NO_POLISH))){
        new_cigar = (uint32_t *) need_buffer((b->core.n_cigar << 2u) + (2u << 2u), buffer, buffer_size);
        trim_cigar(exon->start, exon->end, &pos, &end_pos, bam_get_cigar(b), b->core.n_cigar, new_cigar, &new_n_cigar, md_clip, options);
        if (new_n_cigar == 0) return TRANSMAP_UNMAPPED_NO_OVERLAP;
    } else {
        new_cigar =  bam_get_cigar(b1);
        new_n_cigar = b1->core.n_cigar;
    }
    stitch_cigar(&pos, new_cigar, new_n_cigar, &new_n_cigar, &need_stitch_md);
    b1->core.pos = pos;
    if (bam_set_cigar(b1, new_cigar, new_n_cigar) < 0) return -1;
    if (options & OPTION_FIX_MD) if (fix_MD(b1, buffer, buffer_size, md_clip, need_stitch_md, options & OPTION_FIX_NM) < 0) return -1;
    int i = exon->idx;
    exon_t **exons = tr->exons->data;
    while (exons[i]->end <= b1->core.pos) i++;
    b1->core.pos = b1->core.pos + exons[i]->tstart - exons[i]->start;
    b1->core.tid = tr->new_tid;
    if (tr->strand == '-') {
        b1->core.pos = tr->len - bam_endpos(b1);
        b1->core.flag^=16u;
        bam_rev_cigar(b1);
        bam_rev_seq(b1);
        bam_rev_qual(b1);
        uint8_t *md;
        if ((options & OPTION_FIX_MD) && (md = bam_aux_get(b1, "MD")) != NULL) {
            size_t md_len = strlen((char *)++md);
            if (need_buffer(md_len + 1, buffer, buffer_size) == NULL) return -1;
            bam_rev_aux_md(md, *buffer, md_len);
        }
    }
    return TRANSMAP_MAPPED;
}