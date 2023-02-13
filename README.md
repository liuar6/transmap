transmap
====
 | option&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | description
 ---- | ---- 
 -i / --fi | Input sam/bam file sorted (or grouped) by query name. If the input contains paired-end alignments, "HI" tag must be present to decide the paired records. Currently, transmap does not check the sanity of input sam/bam file since this requires caching the query name which could use a alot of memory. Later version might force name-sorted sam/bam file and perform sanity check.
-o / --fo | Output sam/bam file. The suffix ".sam" or ".bam" indicates the format.
-b / --bed | BED file defining the new targets. The name field of BED record should be unique and will become the reference name of the output.
-g / --gtf | GTF file defining the (probably) multi-exon new targets.
--gtf-feature | GTF feature used to define the member exons of transcripts. Default: exon. For each transcript, the member exons should be present in the same chromosome and strand and their coordinates should not be overlaped.
--gtf-attribute | GTF attribute used as the reference name of the output. Default: transcript_id. It should be noted that the transcript_name is not always unique and must be avoided.
--partial | Also process the records with ranges cross the boundaries of the target regions, the parts of records that excess the range of target regions will be trimmed. By default (--partial is not specified), these records will be skipped for processing.
--no-trim | Do not trim the marginal D or N cigars. By default, when --partial is specified, the alignment records will be trimmed from the two sides until fully contained by the regions or transcripts, and if the resulted alignments are begined or ended with D or N cigars, these parts will also be trimmed. This additional trimming step is nessesary since tools that take sam/bam as input might not accept alignment records that begined or ended with D or N cigars. 
--both mate |  This option indicates that for paired-end alignments, both mate should be successfully mapped to the new reference for the alignments to be reported.
--fix-NH | Fix the NH and HI tag. NH indicates number of reported alignments that contain the query in the current record and NI indicates the index of the current record of all reported alignments. Remapping  could make these information invalid, set this option to rebuild the valid NH and HI.
--fix-MD | Fix the MD tag. When an alignment record is trimmed or the target is in reverse strand, the orginal MD could become invalid, set this option to rebuild the valid MD. Fixing requires the alignments contains an original MD tag.
--fix-NM | Fix the NM tag. When an alignment record is trimmed, the original NM could become invalid, set this option to recalculate the valid NM. Fixing requires an original MD (original NM is not nessesary).
