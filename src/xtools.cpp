#include "utility.h"
#include <iostream>
using namespace std;

void help()
{
    string str =
	"\n"
    "xtools: Xuebing's tools for computational biology\n"
	"	- Xuebing Wu (wuxb07@gmail.com)\n"
    "\n"
    "Usage: xtools <program> [options]\n"
	"       - check usage of each program by omitting options, such as:\n"
	"         xtools geneFeature\n"
    "\n"
    "Programs:\n"
    "\n"
	"Gene analysis ==========\n" // 
    "    geneFeature   discover if a gene list has unique features such as short ORFs\n"
    "    geneSet       find overlaps between a gene list and annotated gene sets\n"
	"    geneScore     correlate your gene scores to gene features and / or gene sets\n"
	"\n"
	"Sequence analysis ======\n"
	"    kShuffle      Shuffle sequences preserving k-nucleotide (kmer) frequency\n"
	"    seqMatch      Find sequence match allowing certain number of mismatches\n"
	"    alignSW       Smith-waterman local sequence alignment\n"
	"    PKA           Positional kmer analysis from aligned sequences\n"
	"    PKA2          PKA with ranked or weighted sequences\n"
	"\n"
	"Interval analysis ======\n"
	"    makeBigWig    Make bigWig file from bed or bam\n"
	"    compareBigWig Compare two bigWig files to find enriched regions\n" // TODO: do binomial test or hypergenometric test using c++, also read bigwig file using c++
	"    sampleBam     Randomly sample a fraction of reads from a Bam file\n"
	"    pe2seBam      Convert paired-end Bam file to single-end\n"
	"\n"	
	"Text ===================\n"
	"    intersectTab  Intersection/subtraction/merge between two tabular files based on a shared key column\n"
	"    TopN          Keep the first N lines for each unique value in a given column\n"
    "\n";

    cerr << str;

    exit(0);
}

int main(int argc, char* argv[]) {

	if (argc < 2) help();
	
	string valid_cmds = "geneFeature, geneSet, geneScore, PKA, PKA2, intersectTab, TopN, kShuffle, seqMatch, smithWaterman, makeBigWig, compareBigWig, sampleBam, pe2seBam, ";
	
	string cmd(argv[1]);
	
	if(valid_cmds.find(cmd+",") == std::string::npos) 
	{
		cerr << "invalid command: " << argv[1] << endl;
		help();
	}
		
	for(int i=2; i < argc; i++) cmd += " "+string(argv[i]);
	
	system_run(cmd);

    return 0;
}
