Curated ribosomal sequences compiled primarily from the SILVA rRNA database by Brian Bushnell. Brian Bushnell's biostars post ( https://www.biostars.org/p/159959/ ) describes how to separate ribosomal sequences from non-ribosomal sequences using bbduk.sh:

"BBDuk works well for this purpose, if you have a large set of curated ribosomal sequences (such as Silva).  You can run it on the reads or the assemblies.

"bbduk.sh in=data.fa outm=ribo.fa outu=nonribo.fa k=31 ref=silva.fasta

"I put a link to a set of ribosomal kmers on Google drive:

"https://drive.google.com/file/d/0B3llHR93L14wS2NqRXpXakhFaEk/view?usp=sharing

"I made it mainly from Silva.  It's small (9MB) and you can use it with BBDuk like this:

"bbduk.sh in=data.fa outm=ribo.fa outu=nonribo.fa k=31 ref=ribokmers.fa.gz

"It has roughly 99.94% sensitivity against the full Silva database.

"The process was a little involved.  I started with the Silva ribosomal database, and followed this procedure:

"1) Deduplicated the sequences with dedupe.sh, since there are lots of redundant copies.

"2) Ran them through kcompress.sh, which produces a fasta file containing all the kmers of interest.  I made several different versions; one containing all 31-mers, one containing only 31-mers that occurred at least 2 times, one for 3 times, etc. up to 50 times (there is a flag for kcompress which specifies this).  That allows variable sensitivity/specificity; the 1-copy version is the biggest and most sensitive, while the 50-copy version is tiny and tends to only contain the most highly conserved ribosomal kmers (or the ones for the organisms that are most popular to sequence so they are in the database a lot).

"3) Then I generated synthetic data from Silva and tested the sensitivity of the different versions.  The point of this was to salvage kmers that were important but missing.  So, for example, I ran BBDuk with the 10-copy kmer set and kept the reads that did not match, and added in some of the most common kmers from those nonmatching reads; this keeps the file size small but increases sensitivity.  This step is not really necessary, though - you can just do 1 and 2 in general, but this is what I did.

"The file on my google drive is, I think, the version in which I kept only kmers present at least 3 times in the deduplicated Sliva database."
