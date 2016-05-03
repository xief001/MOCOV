# MOCOV
###Basic usage:

    blasr queries.fasta queries.fasta -bestn 20 -maxScore 2000 -m 4 -out mapped.m4

    filterm4.py mapped.m4 > mapped.m4.filt

    m4topre.py mapped.m4.filt mapped.m4.filt queries.fasta 24 > mapped.pre

    pbdagcon -c 1 -a mapped.pre > consensus.fasta

    MOCOV consensus.fasta [-options|--options] > corrected.fasta
###Options Description (default value)

-k/--kmer(10)

The length of k-mer. 10-mer is default.

-i/--index(3)

The num of bases you use to index k-mers.

It should be in range [1, 4].

-v/validvalue  n (5)

Every k-mer that have a value higher than n in the hash table should be defined as a valid k-mer.

And one valid k-mer is used to locate the positions of insertion/deletion error.
