#
# 'make depend' uses makedepend to automatically generate dependencies 
#               (dependencies are added to end of Makefile)
# 'make'        build executable file 'mycc'
# 'make clean'  removes all .o and executable files
#

PROGS = xtools kShuffle seqMatch PKA geneFeature geneSet geneScore TopN intersectTab example pValueLookUp alignSW mirscan4 split_file_for_cross_validation cross_validation pairwise_similarity uniqLines sequence_feature selMultiLines MatchScan polyAcluster flipStrand iFeat

# define the C++ compiler to use
CC = g++

# define any compile-time flags, use -O2
CFLAGS = -std=c++11

# define any directories containing header files other than /usr/include
#
INCLUDES = -I../include


.PHONY: depend clean

all: $(PROGS) install

xtools: utility.o xtools.o 
	$(CC) $(CFLAGS) $^ -o $@
	
cross_validation: cross_validation.o text.o utility.o
	$(CC) $(CFLAGS) $^ -o $@
	
split_file_for_cross_validation: split_file_for_cross_validation.o text.o utility.o
	$(CC) $(CFLAGS) $^ -o $@
	
mirscan4: mirscan4.o structure.o utility.o sequence.o stat.o ushuffle.o text.o container.o positional_kmer.o 
	$(CC) $(CFLAGS) $^ -o $@

kShuffle: utility.o sequence.o kShuffle.o ushuffle.o stat.o container.o positional_kmer.o
	$(CC) $(CFLAGS) $^ -o $@

seqMatch: utility.o sequence.o seqMatch.o stat.o ushuffle.o container.o positional_kmer.o
	$(CC) $(CFLAGS) $^ -o $@

alignSW: alignSW.o
	$(CC) $(CFLAGS) $^ -o $@

geneFeature: gene.o utility.o stat.o geneFeature.o
	$(CC) $(CFLAGS) $^ -o $@

geneSet: gene.o utility.o stat.o geneSet.o
	$(CC) $(CFLAGS) $^ -o $@
	
geneScore: gene.o utility.o stat.o geneScore.o
	$(CC) $(CFLAGS) $^ -o $@

PKA: sequence.o markov.o container.o utility.o text.o stat.o ushuffle.o PKA.o positional_kmer.o
	$(CC) $(CFLAGS) $^ -o $@

TopN: TopN.o text.o utility.o
	$(CC) $(CFLAGS) $^ -o $@
	
selMultiLines: selMultiLines.o text.o utility.o
	$(CC) $(CFLAGS) $^ -o $@
	
	
uniqLines: uniqLines.o text.o utility.o
	$(CC) $(CFLAGS) $^ -o $@
	
intersectTab: intersectTab.o text.o utility.o
	$(CC) $(CFLAGS) $^ -o $@

example: utility.o  example.o sequence.o stat.o ushuffle.o container.o structure.o positional_kmer.o
	$(CC) $(CFLAGS) $^ -o $@

MatchScan: utility.o  sequence.o stat.o ushuffle.o positional_kmer.o container.o text.o MatchScan.o
	$(CC) $(CFLAGS) $^ -o $@
	
polyAcluster: utility.o  container.o text.o polyAcluster.o interval.o stat.o
	$(CC) $(CFLAGS) $^ -o $@

flipStrand: utility.o interval.o flipStrand.o text.o
	$(CC) $(CFLAGS) $^ -o $@

iFeat: utility.o interval.o iFeat.o text.o
	$(CC) $(CFLAGS) $^ -o $@

pValueLookUp: utility.o pValueLookUp.o
	$(CC) $(CFLAGS) $^ -o $@

pairwise_similarity: pairwise_similarity.o sequence.o utility.o stat.o ushuffle.o container.o positional_kmer.o
	$(CC) $(CFLAGS) $^ -o $@

sequence_feature: sequence_feature.o utility.o sequence.o ushuffle.o  positional_kmer.o stat.o container.o
	$(CC) $(CFLAGS) $^ -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

# ushuffle.c
.c.o:
	gcc -O3 $(INCLUDES) -c $< -o $@	

install:
	mv $(PROGS) ../bin	

clean:
	$(RM) *.o *~ ../bin/*

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
