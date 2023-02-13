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
--partial | This option informs transmap also process the bam records not fully (or partially) contained by the target regions. The alignments that spans 
--no-trim | Do not trim the marginal I or N for bam records when --partial is specified.
--both mate |  This option indicate that for paired-end alignments, both mate should be successfully mapped to the new targets for them to be report.
--fix-NH | Fix the NH and HI tag. NH indicates number of reported alignments that contain the query in the current record and NI indicates the index of the current record of all reported alignments. Remapping alignments generally will make these information invalid, set this option to rebuild the valid NH and HI.
--fix-MD | Fix the MD tag. When --partial is specified the alignments might be trimmed and make original MD invalid, set this option to rebuild the valid MD. This require the alignment contain an original MD tag.
--fix-NM | Fix the NM tag. --partial might also dirupt the validity of NM tag, set this option to rebuild the valid NM. This require the alignment contain an original MD tag (original NM is not nessesary).
