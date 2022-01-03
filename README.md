# Fastq-trim
A hopefully lightening fast read trimmer.

This is a back-burner project to explore the possibility of developing a
significantly faster short read trimmer.
The ultimate goal is a trimmer that's an order of magnitude faster and good
enough that a more sophisticated trimmer won't make a meaningful difference
to the downstream analysis.

There is no sense of urgency since there are multiple highly-evolved trimmers
available that are fast enough for most purposes.

However, basic testing revealed that adapters can be removed from an
xz-compressed FASTQ file with ~32 million reads in a few minutes on a
ThinkPad with a spinning disk and a 2.6 GHz i5.
This suggests that it might be possible to perform basic trimming as is done
prior to an RNA-Seq analysis in much less time than current tools allow.

How much faster remains to be seen as we explore the detailed requirements
(in our spare time).
