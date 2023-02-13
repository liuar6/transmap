transmap
====
 | option | description
 ---- | ---- 
 -i / --fi | input sam/bam file sorted (or grouped) by read name. If the input bam contained alignments of paired-end reads, "HI" tag must present for transmap to identified paired alignments.
-o / --fo | output sam/bam file. The suffix ".sam" or ".bam" indicates the output format.
-b / --bed | BED file providing the coordinates of target regions. The name field should be unique will become the reference name of the output.
-g / --gtf | GTF file providing the coordinates of target transcripts.
--gtf-feature | GTF feature used to define the sub-regions of transcripts. default:exon. The coordinates of each exon should present in the same chromosome and strand and must not overlap.
--gtf-attribute | GTF attribute used as the reference name of the output. default: transcript_id. It should be noted that the transcript_name is not always unique and must be avoided.
--partial | This option informs transmap also process the bam records not fully (or partially) contained by the target regions. The alignments that spans 
--no-trim | Do not trim the marginal I or N for bam records when --partial is specified.

-i/--fi                       input bam file sorted (or grouped) by read name.<br>
-o/--fo                    output bam file.<br>
-b/--bed                bed file providing the coordinates of target regions.<br>
-g/--gtf                  gtf file providing the coordinates of target transcripts.<br>
--gtf-feature        gtf feature used to define regions of transcripts. default: exon.<br>
--gtf-attribute     gtf attribute used as the name of transcripts. default: transcript_id.<br>
--partial                 also process the bam records not fully (partially) contained by the target regions.<br>
--no-trim               do not trim the marginal I or N for bam records when --partial is specified.<br>
--both-mate         require both mate of paired-end alignments to be mapped for reporting.<br>
--fix-NH                  fix the NH and HI tag.<br>
--fix-MD                  fix the MD tag if exists.<br>
--fix-NM                  fix the NM tag when --fix-MD is specified.<br>
