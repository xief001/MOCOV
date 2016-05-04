# MOCOV
###Basic usage:

    MOCOV queries.fasta [-options|--options] > corrected.fasta
###Options Description (default value)

-k/--kmer(10)

The length of k-mer. 10-mer is default.

-i/--index(3)

The num of bases you use to index k-mers.

It should be in range [1, 4].

-v/--validvalue  n (5)

Every k-mer that have a value higher than n in the hash table should be defined as a valid k-mer.

And one valid k-mer is used to locate the positions of insertion/deletion error.
