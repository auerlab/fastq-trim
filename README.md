# Fastq-trim
Hopefully lightening fast read trimmer

This is a back-burner project to explore the possibility of developing a
faster short read trimmer.

There is no sense of urgency since there are multiple good quality trimmers
available that are fast enough for most purposes.

However, basic testing revealed that an xz-compressed FASTQ file with 32
million reads can be scanned for adapters in 88 seconds on a ThinkPad with a
spinning disk and a 2.6 GHz i5.  This suggests that it might be possible to
perform basic trimming as is done prior to an RNA-Seq analysis in much less
time than current tools allow.

How much faster remains to be seen as we explore the detailed requirements
(in our spare time).
