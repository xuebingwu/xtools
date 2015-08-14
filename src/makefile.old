#
# 'make depend' uses makedepend to automatically generate dependencies 
#               (dependencies are added to end of Makefile)
# 'make'        build executable file 'mycc'
# 'make clean'  removes all .o and executable files
#

PROGS = xtools kShuffle seqMatch PKA geneFeature geneSet geneScore TopN intersectTab example pValueLookUp alignSW mirscan4 split_file_for_cross_validation cross_validation pairwise_similarity uniqLines

# define the C++ compiler to use
CC = g++

# define any compile-time flags
CFLAGS = -std=c++11

# define any directories containing header files other than /usr/include
#
INCLUDES = -I/lab/bartel1_ata/wuxbl/lib/boost_1_57_0/  -I../include

# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
# LFLAGS = -L/home/newhall/lib  -L../lib

# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname 
#   option, something like (this will link in libmylib.so and libm.so:
# LIBS = -lmylib -lm


# define the C object files 
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
#OBJS = $(SRCS:.cpp=.o)

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: depend clean

all: $(PROGS) install

xtools: utility.o xtools.o 
	$(CC) $(CFLAGS) $^ -o $@
	
cross_validation: cross_validation.o text.o utility.o
	$(CC) $(CFLAGS) $^ -o $@
	
split_file_for_cross_validation: split_file_for_cross_validation.o text.o utility.o
	$(CC) $(CFLAGS) $^ -o $@
	
mirscan4: mirscan4.o structure.o utility.o sequence.o stat.o ushuffle.o text.o
	$(CC) $(CFLAGS) $^ -o $@

kShuffle: utility.o sequence.o kShuffle.o ushuffle.o stat.o
	$(CC) $(CFLAGS) $^ -o $@

seqMatch: utility.o sequence.o seqMatch.o stat.o ushuffle.o
	$(CC) $(CFLAGS) $^ -o $@

alignSW: alignSW.o
	$(CC) $(CFLAGS) $^ -o $@

geneFeature: gene.o utility.o stat.o geneFeature.o
	$(CC) $(CFLAGS) $^ -o $@

geneSet: gene.o utility.o stat.o geneSet.o
	$(CC) $(CFLAGS) $^ -o $@
	
geneScore: gene.o utility.o stat.o geneScore.o
	$(CC) $(CFLAGS) $^ -o $@

PKA: sequence.o markov.o container.o utility.o text.o stat.o ushuffle.o PKA.o
	$(CC) $(CFLAGS) $^ -o $@

TopN: TopN.o text.o utility.o
	$(CC) $(CFLAGS) $^ -o $@
	
uniqLines: uniqLines.o text.o utility.o
	$(CC) $(CFLAGS) $^ -o $@
	
intersectTab: intersectTab.o text.o utility.o
	$(CC) $(CFLAGS) $^ -o $@

example: utility.o  example.o sequence.o stat.o ushuffle.o
	$(CC) $(CFLAGS) $^ -o $@

pValueLookUp: utility.o pValueLookUp.o
	$(CC) $(CFLAGS) $^ -o $@

pairwise_similarity: pairwise_similarity.o sequence.o utility.o stat.o ushuffle.o
	$(CC) $(CFLAGS) $^ -o $@

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

# ushuffle.c
.c.o:
	gcc -O3 $(INCLUDES) -c $< -o $@	

install:
	mv $(PROGS) ../bin	

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
