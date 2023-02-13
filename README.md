transmap
====
 | option&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | description
 ---- | ---- 
 -i / --fi | input sam/bam file sorted (or grouped) by read name. If the input bam contained alignments of paired-end reads, "HI" tag must present for transmap to identified paired alignments.
-o / --fo | output sam/bam file. The suffix ".sam" or ".bam" indicates the output format.
-b / --bed | BED file providing the coordinates of target regions. The name field should be unique will become the reference name of the output.
-g / --gtf | GTF file providing the coordinates of target transcripts.
--gtf-feature | GTF feature used to define the sub-regions of transcripts. default:exon. The coordinates of each exon should present in the same chromosome and strand and must not overlap.
--gtf-attribute | GTF attribute used as the reference name of the output. default: transcript_id. It should be noted that the transcript_name is not always unique and must be avoided.
--partial | This option informs transmap also process the bam records not fully (or partially) contained by the target regions. The alignments that spans 
--no-trim | Do not trim the marginal I or N for bam records when --partial is specified.
--both mate |  This option indicate that for paired-end alignments, both mate should be successfully mapped to the new targets for them to be report.
--fix-NH | Fix the NH and HI tag. NH indicates number of reported alignments that contain the query in the current record and NI indicates the index of the current record of all reported alignments. Remapping alignments generally will make these information invalid, set this option to rebuild the valid NH and HI.
--fix-MD | Fix the MD tag. When --partial is specified the alignments might be trimmed and make original MD invalid, set this option to rebuild the valid MD. This require the alignment contain an original MD tag.
--fix-NM | Fix the NM tag. --partial might also dirupt the validity of NM tag, set this option to rebuild the valid NM. This require the alignment contain an original MD tag (original NM is not nessesary).
