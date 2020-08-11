## Biostrings
library(Biostrings)

# basics of DNAStrings
dna <- DNAString("TCGAGCAAT")    # define a DNAString
dna
length(dna)    # number of bases in a DNAString
DNAString("JQX")    # error - invalid bases
DNAString("NNNACGCGC-TTA-CGGGCTANN")    # valid sequence with unknowns and gaps
dna[4:6]    # extract a substring
as.character(dna)    # convert DNAString to character

# basics of DNAStringSets
set1 <- DNAStringSet(c("TCA", "AAATCG", "ACGTGCCTA", "CGCGCA", "GTT", "TCA"))    # define a DNAStringSet
set1
set1[2:3]    # extract subset of sequences
set1[[4]]    # extract one sequence as a single DNAString
length(set1)    # number of DNAstrings in set
width(set1)    # size of each DNAString
duplicated(set1)    # detect which sequences are duplicated
unique(set1)    # keep only unique sequences
sort(set1)

dna_seq <- DNAString("ATCGCGCGCGGCTCTTTTAAAAAAACGCTACTACCATGTGTGTCTATC")

# analyze DNAStrings
letterFrequency(dna_seq, "A")    # count A in sequence
letterFrequency(dna_seq, "GC")    # count G or C in sequence
dinucleotideFrequency(dna_seq)    # frequencies of all dinucleotides
trinucleotideFrequency(dna_seq)    # frequencies of all trinucleotides

# convert DNAStrings
reverseComplement(dna_seq)    # find reverse complement
translate(dna_seq)    # amino acid translation

# count and match on individual Biostrings
dna_seq <- DNAString("ATCGCGCGCGGCTCTTTTAAAAAAACGCTACTACCATGTGTGTCTATC")
dna_seq
countPattern("CG", dna_seq)    # pattern "CG" occurs 5 times
matchPattern("CG", dna_seq)    # locations of pattern "CG"
start(matchPattern("CG", dna_seq))    # start locations of the pattern
matchPattern("CTCTTTTAAAAAAACGCTACTACCATGTGT", dna_seq)    # match patterns of any length

# check for pattern and its reverse complement
countPattern("TAG", dna_seq)
countPattern(reverseComplement(DNAString("TAG")), dna_seq)

# count and match on sets of Biostrings
set2 <- DNAStringSet(c("AACCGGTTTCGA", "CATGCTGCTACA", "CGATCGCGCCGG", "TACAACCGTACA"))
set2
vcountPattern("CG", set2)    # CG counts for entire DNAStringSet
vmatchPattern("CG", set2)
vmatchPattern("CG", set2)[[1]]    # access matches for the first element of the DNAStringSet

## Exercises
eco <- DNAString("GGTTTCACCGCCGGTAATGAAAAAGGCGAACTGGTGGTGCTTGGACGCAACGGTTCCGACTACTCTGCTGCGGTGCTGGCTGCCTGTTTACGCGCCGATTGTTGCGAGATTTGGACGGACGTTGACGGGGTCTATACCTGCGACCCGCGTCAGGTGCCCGATGCGAGGTTGTTGAAGTCGA")
eco

# number of bases
length(eco)

# start codons
countPattern("ATG", eco)
matchPattern("ATG", eco)

# substrings
s <- start(matchPattern("ATG", eco))[1]
subeco <- eco[s:length(eco)]

# AAstrings
subecoAA <- translate(subeco)
length(subecoAA)

# stop codons
matchPattern("*", subecoAA)
stop <- start(matchPattern("*", subecoAA))

peptide <- subecoAA[1:(stop-1)]
as.character(peptide)

# AA frequencies and charge (at pH 7)
pos <- c("HKR")
neg <- c("DE")
letterFrequency(peptide, pos)
letterFrequency(peptide, neg)
netCharge <- letterFrequency(peptide, pos)-letterFrequency(peptide, neg)
netCharge
