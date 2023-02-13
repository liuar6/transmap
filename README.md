Overview
====
**transmap** takes alignments to the reference genome as input and produces alignments to genomic regions (BED mode) or transcripts (GTF mode). transmap has two major mode: BED mode and GTF mode. In BED mode, the input alignments are transformed just like they are aligned to the  genomic regions defined by the BED file. In GTF mode, the input alignments are transformed as they are aligned to  the single- or multi-exon transcripts defined by the GTF file. The main difference between BED mode and GTF mode is that the BED mode can produce intron-containing alignments while for GTF mode, the introns of the input alignments should be compatible with the exon structure of the transcripts and will be stitched in the transformed alignments. Thus GTF mode will not produce alignments with introns.

Installation
====
```sh
mkdir build && cd build 
cmake .. 
make
mv transmap ~/bin/
```
Usage
====
 | option&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | description
 ---- | ---- 
 -i / --fi | Input sam/bam file sorted (or grouped) by query name. If the input contains paired-end alignments, "HI" tag must be present to decide the paired records. Currently, transmap does not check the sanity of input sam/bam file since this requires caching the query name which could use a lot of memory. Later version might force name-sorted sam/bam file and perform sanity check.
-o / --fo | Output sam/bam file. The suffix ".sam" or ".bam" indicates the format.
-b / --bed | BED file that provides the regions on which the alignments to be generated. The name field of BED record will become the reference name of the output and should be unique.
-g / --gtf | GTF file that provides the exons of transcripts on which the alignments to be generated.
--gtf-feature | GTF feature used to define the member exons of transcripts. Default: exon. For each transcript, the member exons should be present in the same chromosome and strand and their coordinates should not be overlaped.
--gtf-attribute | GTF attribute used as the reference name of the output. Default: transcript_id. It should be noted that the transcript_name is not always unique and must be avoided.
--partial | Also process the alignment records with ranges exceed the target boundaries and these records will be trimmed from the two sides until fully contained by the targets. In GTF mode, this option allows the alignment records exceed the transcript boundaries. The records that exceed the exon ends and overlap with the introns will still be excluded no matter whether --partial is set.
--no-trim | Do not trim the marginal D or N cigars. By default, if the trimmed alignments are begined or ended with D or N cigars, these parts will be further trimmed. The additional trimming step is nessesary since tools that take sam/bam as input might not accept alignment records that begined or ended with D or N. This default behavior can be disabled with --no-trim.
--both mate |  This option indicates that for paired-end alignments, both mate should be successfully mapped to the new reference for the alignments to be reported.
--fix-NH | Fix the NH and HI tag. NH indicates number of reported alignments that contain the query in the current record and NI indicates the index of the current record of all reported alignments. Remapping  could make these information invalid, set this option to rebuild a valid NH and HI.
--fix-MD | Fix the MD tag. When an alignment record is trimmed or the target is in reverse strand, the orginal MD could become invalid, set this option to rebuild a valid MD. Fixing requires the alignments contains an original MD tag.
--fix-NM | Fix the NM tag. When an alignment record is trimmed, the original NM could become invalid, set this option to recalculate a valid NM. Fixing requires an original MD (original NM is not nessesary).

Author
====
Anrui Liu  (liuar6@gmail.com)

