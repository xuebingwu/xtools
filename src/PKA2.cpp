#include "sequence.h"
#include "markov.h"
#include "container.h"
#include "utility.h"
#include "text.h"

void print_help()
{
	string txt = "\n"
	"PKA2: Positional Kmer Analysis with ranked or weighted sequences (version 0.1)\n"  		
	"    -by Xuebing Wu (wuxb07@gmail.com)\n"  
	"    Bartel lab, Whitehead Institute\n\n"  
	
	"    Identify statistically enriched / depleted short sequences of\n"
	"	 length k (kmer) at every position in a set of aligned sequences.\n"
	"    Sequences are either weighted or ranked. Two-sample t test is\n"
	"    used to on weighted sequences, whereas Wilcoxon rank-sum test\n"
	"    (i.e. Mann-Whitney U test ) is used for ranked sequences.\n\n"  
	
	"Usage: PKA2 input [options]\n"  
	"  example usage:\n\n"  
	"    PKA2 input.fa (sequences are ranked in fasta)\n"   
	"    PKA2 input.txt -seq 1 (ranked sequences in colum 1 of the tabular input file)\n"  
	"    PKA2 input.txt -seq 1 -weight 2 (column 2 indicates sequence weight)\n\n"  
	
	"Options\n\n"  
	
	"  Input =================\n"  
	"    -alphabet <ACGT>     alphabet for generating kmers, default=ACGT, case insensitive\n"  
	"    -seq <a>             sequences are in column a. Starts at 1\n"
	"    -weight <b>          weights are in column b. Starts at 1\n"
	"    -skip <c>            skip the first c lines\n"  		 
	"  Kmer counting =========\n"  
	"    -upto <K>            consider all kmers of length 1,2,...,K. default=4 \n"  
	"    -k <K>               use fixed kmer length (K)\n"  
	"    -shift <0>           max shift (to right) allowed for kmer positions, default=0\n"  
	"    -degenerate          allow IUPAC degenerate nucleoties in kmer. Only work for DNA/RNA\n" 
	"  Statistics & output ===\n"  
	"    -o <output-prefix>   prefix for all output files, default=PKA_output\n"  
	"    -pCutoff <p>         raw p-value cut-off, default=0.01\n"  
	"    -pCutoff_B <p>       Bonferoni corrected p-value cut-off, default=0.05\n"  
	"    -startPos <n>        set position n (0,1,2,3,..) to be the start position (0) in the output\n\n";  

	cout << txt ;
	
	exit(0);
}

 
int main(int argc, char* argv[]) {

	/*******	   part 1: default parameters	   */

	string output = "PKA2"; // output prefix
	string alphabet = "ACGT";   // default DNA
	int k=0;	// k
	int max_k = 0;  // upto max_k
	int shift=0; // no shift, i.e. exact position
	bool degenerate = false;	// no degenerate bases allowed by default
	double pCutoff=0.05;	// binomial p cutoff
	double pCutoff_B=0.05;  // Bonferoni p cutoff
	int topN=1; // output the best kmer (smallest p, when p=0, use larger zscore)
	int startPos=0; // coordinates
	string left="0";
	string right="0";
	
	bool ranked = true;
	bool weighted = false;
	int skip = 0;
	int cSeq = -1; // if not specified, assume input is ranked fasta
	int cWeight = -1;
	

	string inputfile,save_to_file,str;

	/*******		part 2: get commandline arguments   */

	if (argc < 2) print_help(); // if only type PKA, print help and exit

	inputfile = argv[1]; // the first argument is input sequence file

	if (inputfile == "-h" || inputfile == "--help" ) print_help(); // PKA -h or PKA --help: print help and exit

	// other arguments
	for (int i = 2; i < argc; i++) { 
		if (i != argc) { 
			str=argv[i];
			if (str == "-o") {   // output prefix
				output = argv[i + 1];
				i=i+1;
			} else if (str == "-k") {   // fixed k
				k = atoi(argv[i + 1]);
				i=i+1;
			} else if (str == "-upto") {	// max k
				max_k = atoi(argv[i + 1]);
				i=i+1;
			} else if (str == "-shift") {
				shift = atoi(argv[i + 1]);
				i=i+1;
			} else if (str == "-alphabet") {
				alphabet = argv[i + 1];
				i=i+1;
			} else if (str == "-degenerate") {
				degenerate = true;
			} else if (str == "-pCutoff") {
				pCutoff = atof(argv[i + 1]);
				i=i+1;
			} else if (str == "-pCutoff_B") {
				pCutoff_B = atof(argv[i + 1]);
				i=i+1;
			} else if (str == "-topN") {
				topN = atoi(argv[i + 1]);
				i=i+1;
			} else if (str == "-startPos") {
				startPos = atoi(argv[i + 1]);
				i=i+1;
			} else if (str == "-highlight") {
				string s(argv[i+1]);
				vector<string> ss = string_split(s,",");
				left=ss[0];
				right=ss[1];
				i=i+1;
			} else if (str == "-save") {
				save_to_file = argv[i + 1];
				i=i+1;
			} else if (str == "-weight") {
				weighted=true;
				ranked=false; // overwrite 
				cWeight = atoi(argv[i + 1]) - 1;
				i=i+1;
			} else if (str == "-seq") {
				cSeq = atoi(argv[i + 1]) - 1;
				i=i+1;
			} else if (str == "-skip") {
				skip = atoi(argv[i + 1]);
				i=i+1;
			} else {
				message("Unknown options: "+str);
				print_help();
			}
		}
	}
	

	// determine the length of k-mer
	int min_k = 1;
	if (k == 0 && max_k == 0) // neither is specified, do upto 4-mers
	{
		max_k = 4;
	}
	else if ( k > 0 && max_k == 0) // 
	{  
		min_k = k;
		max_k = k; // start from k to k
	} // else do 1..max_k


	// print out parameters used
	message("Summary: ");
	message("   Input       :   " + inputfile);
	message("    weighted   :	" + to_string(weighted));
	message("    ranked     :	" + to_string(ranked));
	message("   Output      :   " + output +".*");  
	message("   alphabet    :   " + alphabet);
	message("   kmer from   :   " + to_string(min_k) );
	message("   kmer upto   :   " + to_string(max_k) );
	message("   shift       :   " + to_string( shift) );
	message("   degenerate  :   " + to_string(degenerate) );
	message("   pCutoff     :   " + to_string(pCutoff ));
	message("   pCutoff_B   :   " + to_string(pCutoff_B ));
	message("   top N       :   " + to_string(topN ));
	message("   start at    :   " + to_string(startPos));

	/***********	part 3: process input				 */

	vector<string> seqs;
	vector<double> weights;
	
	// load input data	
	if(weighted) load_weighted_sequences_to_vectors(inputfile,seqs,weights,skip,cSeq,cWeight);
	else seqs = load_ranked_sequences_to_vectors(inputfile,skip,cSeq);
	
	// total number of sequences
	int nSeq = seqs.size();
	// sequence lenght. assume all sequence have the same length
	int lSeq = seqs[0].size();
	
	// make sure there are sequences in the input file
	if (nSeq == 0)
	{
		message("No sequence present in input file: " + inputfile);
		exit(1);
	}
	// show the number of sequences loaded
	message(to_string(seqs.size()) + " sequences of "+to_string(lSeq)+" bases loaded from " + inputfile);
	
	// replace U with T
	if (alphabet == "ACGT" || alphabet == "ACGU")
	{
		message("Replacing U with T...");
		alphabet = "ACGT";
		for(int i=0;i<nSeq;i++)	
		{
			replace(seqs[i].begin(),seqs[i].end(),'U','T');
		}
	}

	/********   part 4: generate kmers	  */

	// total number of tests to be performed, ~ n_kmer * seq_len
	// to be used in multiple testing correction
	int nTest = 0; 
	
	// generate all exact kmers	
	vector<string> kmers = generate_kmers(min_k, alphabet);
	if (degenerate == false) nTest = kmers.size() * (lSeq - min_k + 1);
	for (k = min_k+1; k <= max_k; k++)
	{
		vector<string> tmp = generate_kmers(k, alphabet);	
		kmers += tmp;
		if (degenerate == false) nTest += tmp.size() * (lSeq - k + 1);
	}
	message(to_string(kmers.size()) +  " exact kmers to be tested" );
		
	// generate degenerate kmers, which will also include exact kmers
	if (degenerate)
	{
		kmers = degenerate_kmer(min_k);
		nTest = kmers.size() * (lSeq - min_k + 1);
		for (k = min_k+1; k <= max_k; k++)
		{
			vector<string> tmp = degenerate_kmer(k);	
			kmers += tmp;
			nTest += tmp.size() * (lSeq - k + 1);
		}		
		message(to_string(kmers.size()) + " k-mers allowing degenerate bases");
	}
	
	message(to_string(nTest) + " tests (kmer x position) will be performed");

	/******	 part 5: kmer counting and statistics		*/

	// output file
	string outfile = output+".pass.p.cutoff.txt";

	// number of significant kmers detected
	array<int,2> nSig;
	
	if (weighted) nSig = find_significant_kmer_from_weighted_sequences(seqs, weights, kmers, outfile, nTest, pCutoff, pCutoff_B, shift, startPos);
	else nSig = find_significant_kmer_from_ranked_sequences(seqs, kmers,outfile, nTest, pCutoff, pCutoff_B, shift,startPos);
	
	message( to_string (nSig[0]) +  " positional kmers with p <"+to_string(pCutoff));
	message( to_string (nSig[1]) +  " positional kmers with corrected p <"+to_string(pCutoff_B)); 
	
	// sort output
	system_run("sort -k5,5gr -k3,3n "+outfile+" > "+output+".tmp");
	system_run("mv "+output+".tmp "+outfile);
	
	// plot single nucleotide profile
	plot_PKA2_nucleotide_output( outfile,  output+".nucleotide.profile.pdf",  lSeq);
	
	// build model
	vector<positional_kmer> ranked_kmers = build_model_from_PKA2_output(outfile, -log10(pCutoff), -log10(pCutoff_B));
	message(to_string(ranked_kmers.size())+" significant positional kmers are used to build the model");
	
	// save the model to file
	save_model_to_file(ranked_kmers, output+".model.txt");	
	
	message("scoring sequences using the model...");
	//load_weighted_sequences_to_vectors(inputfile,seqs,weights,skip,cSeq,cWeight);
	for (int i=0;i<seqs.size();i++)
	{
		double score = score_sequence_using_PKA_model(ranked_kmers, seqs[i]);
		cout << seqs[i] << "\t" << weights[i] << "\t" << score << endl;		
	}
	
	// writing feature matrix
	vector<string> sig_kmers;
	vector<int> sig_positions;
	read_significant_positional_kmer_from_PKA2_output( outfile, sig_kmers, sig_positions);
	significant_feature_matrix_PKA2(seqs, weights, sig_kmers, sig_positions, output+".sig.feature.txt", shift);
		
	message("Done!");
	


	return 0;
} 

