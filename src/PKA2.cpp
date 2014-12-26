#include "sequence.h"
#include "markov.h"
#include "container.h"
#include "utility.h"
#include "text.h"

#include <boost/algorithm/string.hpp>


void print_help()
{
	string txt = "\n"
	"PKA2: Positional Kmer Analysis with ranked or weighted sequences (version 0.1)\n"  		
	"    -by Xuebing Wu (wuxb07@gmail.com)\n"  
	"    Bartel lab, Whitehead Institute\n\n"  
	
	"    Identify statistically enriched / depleted short sequences of \n"
	"    (or upto) length k (kmer) at every position in a set of aligned \n"
	"    sequences each asscociated with a weight. P-values are calculated\n"
	"    by either t-test or Wilcoxon rank-sum test (i.e. Mann-Whitney U test).\n\n"  
	
	"Usage: PKA2 input [options]\n"  
	"  example usage:\n\n"  
	"    PKA2 input.txt (sequence in column 1 and weight in column 2, t-test)\n"   
	"    PKA2 input.txt -rank_sum (rank sum test, input needs to be sorted)\n"  
	"    PKA2 input.txt -seq 2 -weight 3\n\n"  
	
	"Options\n\n"  
	
	"  Input =======================\n"  
	"    -alphabet <ACGT>     alphabet for generating kmers, default=ACGT, case insensitive\n"
	"                         note: 'dna' is equivalent to 'ACGT', \n"
	"                               'protein' is equivalent to 'ACDEFGHIJKLMNOPQRSTUVWY'\n"
	"    -seq <a>             sequences are in column a. Starts at 1\n"
	"    -weight <b>          weights are in column b. Starts at 1\n"
	"    -skip <c>            skip the first c lines, such as headers or annotations\n"  		 
	"  Kmer counting ===============\n"  
	"    -upto <K>            consider all kmers of length 1,2,...,K. default=4 \n"  
	"    -k <K>               use fixed kmer length (K)\n"  
	"    -shift_min <0>       min shift (to right) allowed for kmer positions, default=0\n" 
	"    -shift_max <1>       max shift (to right) allowed for kmer positions, default=1\n"  	 
	"    -degenerate <ACGTRYMKWSBDHVN> alphabet to use for degenerate kmers. By default \n"
	"                         all IUPAC letters (ACGTRYMKWSBDHVN) are used. One can use \n"
	"                         ACGTN to search gapped-kmers. Only work for DNA/RNA\n" 
	"    -pair                also test all possible pairs of positional monomers\n"  	 
	"  Statistics & output =========\n"  
	"    -o <output-prefix>   prefix for all output files, default=PKA_output\n"  
	"    -rank_sum            use Wilcoxon rank-sum / Mann-Whitney U test instead of t test\n"
	"                         WARNING: the input needs to be sorted \n"
	"    -pCutoff <p>         raw p-value cut-off, default=0.05\n"  
	"    -pCutoff_B <p>       Bonferoni corrected p-value cut-off, default=0.05\n"  
	"    -startPos <n>        set position n (0,1,2,3,..) to be the start position (0) in the output\n"
	"  Training & Prediction =======\n"  
	"    -predict <prefix>    use significant kmers from a previous run (-o prefix) to score input sequences\n"
	"\n"; 

	cout << txt ;
}

 
int main(int argc, char* argv[]) {

	/*******	   part 1: default parameters	   */

	string output = "PKA2"; // output prefix
	string alphabet = "ACGT";   // default DNA
	int k=0;	// k
	int max_k = 0;  // upto max_k
	int shift_min=0; // no shift, i.e. exact position
	int shift_max=1; // shift 1 bases
		
	bool degenerate = false;	// no degenerate bases allowed by default
	string degenerate_alphabet = "ACGTRYMKWSBDHVN"; // if degenerate allowed, default to use all IUPAC bases
	
	bool pair = false; // test all pairs of monomers
	
	double pCutoff=0.05;	// binomial p cutoff
	double pCutoff_B=0.05;  // Bonferoni p cutoff
	int topN=1; // output the best kmer (smallest p, when p=0, use larger zscore)
	int startPos=0; // coordinates
	string left="0";
	string right="0";
	
	bool rank_sum = false; // use rank sum test rather than t test, requie sorted input
	
	// default input format: no header, column 1 sequence and column 2 weight
	int skip = 0;
	int cSeq = 0;  // first column, 0-based
	int cWeight = 1; // second column, 0-based
	
	string prefix;
		
	string inputfile,save_to_file,str;

	/*******		part 2: get commandline arguments   */

	if (argc < 2) 
	{
		print_help(); // if only type PKA, print help and exit
		exit(0);
	}

	inputfile = argv[1]; // the first argument is input sequence file

	if (inputfile == "-h" || inputfile == "--help" ) 
	{
		print_help(); // PKA -h or PKA --help: print help and exit
		exit(0);
	}
	else if (inputfile[0] == '-')
	{
		print_help();
		message("ERROR: the first argument needs to be input file name! ");
		exit(1);
	}

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
			} else if (str == "-shift_min") {
				shift_min = atoi(argv[i + 1]);
				i=i+1;
			} else if (str == "-shift_max") {
				shift_max = atoi(argv[i + 1]);
				i=i+1;
			} else if (str == "-alphabet") {
				alphabet = argv[i + 1];
				i=i+1;
			} else if (str == "-degenerate") {
                degenerate = true;
				degenerate_alphabet = argv[i+1];
				i=i+1;
			} else if (str == "-pair") {
                pair = true;
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
			} else if (str == "-predict") {
				prefix = argv[i + 1];
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
				cWeight = atoi(argv[i + 1]) - 1;
				i=i+1;
			} else if (str == "-seq") {
				cSeq = atoi(argv[i + 1]) - 1;
				i=i+1;
			} else if (str == "-rank_sum") {
				rank_sum = true;
			} else if (str == "-skip") {
				skip = atoi(argv[i + 1]);
				i=i+1;
			} else {
				print_help();
				message("**** Unknown options: "+str);
				exit(1);
			}
		}
	}
	
	if(boost::algorithm::to_lower_copy(alphabet) == "dna") alphabet="ACGT";
	else if (boost::algorithm::to_lower_copy(alphabet) == "protein") alphabet = "ACDEFGHIJKLMNOPQRSTUVWY";


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
	message("   input       :   " + inputfile);
	message("   rank_sum    :	" + to_string(rank_sum));
	message("   output      :   " + output +".*");  
	message("   alphabet    :   " + alphabet);
	message("   kmer from   :   " + to_string(min_k) );
	message("   kmer upto   :   " + to_string(max_k) );
	message("   shift_min   :   " + to_string(shift_min) );
	message("   shift_max   :   " + to_string(shift_max) );
	if (degenerate)
    message("   degenerate  :   " + degenerate_alphabet );
    message("   pair        :   " + pair );
	message("   pCutoff     :   " + to_string(pCutoff ));
	message("   pCutoff_B   :   " + to_string(pCutoff_B ));
	message("   top N       :   " + to_string(topN ));
	message("   start at    :   " + to_string(startPos));

	/***********	part 3: process input				 */

	vector<string> seqs;
	vector<double> weights;
	
	// load input data	
	if(rank_sum) seqs = load_ranked_sequences_to_vectors(inputfile,skip,cSeq);
	else load_weighted_sequences_to_vectors(inputfile,seqs,weights,skip,cSeq,cWeight);
	
	// total number of sequences
	int nSeq = seqs.size();
	
	// make sure there are sequences in the input file
	if (nSeq == 0)
	{
		message("No sequence present in input file: " + inputfile);
		exit(1);
	}
	
	// sequence lenght. assume all sequence have the same length
	int lSeq = seqs[0].size();
	
	// show the number of sequences loaded
	message("- " + to_string(seqs.size()) + " sequences of "+to_string(lSeq)+" bases loaded from " + inputfile);
	
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
	
	if(prefix.size() > 0) // prediction mode
	{
		message("=== prediction mode ===");
		/**/
		message("building model from file: "+prefix+".pass.p.cutoff.txt");		
		vector<positional_kmer> ranked_kmers = build_model_from_PKA2_output(prefix+".pass.p.cutoff.txt", -log10(pCutoff), -log10(pCutoff_B));
		message("- " + to_string(ranked_kmers.size())+" positional kmers included in the model");
	
		// save the model to file
		// save_model_to_file(ranked_kmers, output+".model.txt");	
	
		/**/
		message("scoring input sequences using the model...");
		string scoreFile = prefix+".score";
	    ofstream out1(scoreFile.c_str());
		//load_weighted_sequences_to_vectors(inputfile,seqs,weights,skip,cSeq,cWeight);
		for (int i=0;i<seqs.size();i++)
		{
			double score = score_sequence_using_PKA_model(ranked_kmers, seqs[i]);
			out1 << seqs[i] << "\t" << weights[i] << "\t" << score << endl;		
		}
		out1.close();
		message("- done");
		
		/**/
		/**/
		/**/
		// writing feature matrix
		//vector<string> sig_kmers;
		//vector<int> sig_positions;
		//read_significant_positional_kmer_from_PKA2_output( outfile, sig_kmers, sig_positions);
		//significant_feature_matrix_PKA2(seqs, weights, sig_kmers, sig_positions, output+".sig.feature.txt", shift);
    
		//significant_feature_matrix_PKA2(seqs, weights, ranked_kmers, output+".sig.feature.txt");
	
	
		/**/
		
		if (pair == false) return 0;
		
		message("loading paired monomers from file: "+prefix+".significant.pair");
		vector<paired_kmer> paired_kmer_model = build_paired_kmer_model(prefix+".significant.pair");
		message("- " + to_string(paired_kmer_model.size())+" pairs included in the model");

		message("scoring input sequences using the model...");
		string pairedScoreFile = prefix+".pair.score";
	    ofstream out(pairedScoreFile.c_str());
		for (int i=0;i<seqs.size();i++)
		{
			double score = score_sequence_using_paired_kmer_model(paired_kmer_model, seqs[i]);
			out << seqs[i] << "\t" << weights[i] << "\t" << score << endl;		
		}	
		out.close();
		message("- done");
		
		
		return 0;
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
		kmers = degenerate_kmer(min_k,degenerate_alphabet);
		nTest = kmers.size() * (lSeq - min_k + 1);
		for (k = min_k+1; k <= max_k; k++)
		{
			vector<string> tmp = degenerate_kmer(k,degenerate_alphabet);	
			kmers += tmp;
			nTest += tmp.size() * (lSeq - k + 1);
		}		
		message(to_string(kmers.size()) + " k-mers allowing degenerate bases");
	}
	
	// multiple by shift
	nTest = nTest * (shift_max - shift_min + 1);
	message(to_string(nTest) + " tests (kmer x position x shifts) will be performed");

	/******	 part 5: kmer counting and statistics		*/

	// output file
	string outfile = output+".pass.p.cutoff.txt";

	// number of significant kmers detected
	array<int,2> nSig;
	
	if (rank_sum) nSig = find_significant_kmer_from_ranked_sequences(seqs, kmers,outfile, nTest, pCutoff, pCutoff_B, shift_max,startPos);
	else nSig = find_significant_kmer_from_weighted_sequences(seqs, weights, kmers, outfile, nTest, pCutoff, pCutoff_B, shift_min, shift_max, startPos);
	
	message( to_string (nSig[0]) +  " positional kmers with p <"+to_string(pCutoff));
	message( to_string (nSig[1]) +  " positional kmers with corrected p <"+to_string(pCutoff_B)); 
	
	// sort output
	system_run("sort -k5,5gr -k3,3n "+outfile+" > "+output+".tmp");
	system_run("mv "+output+".tmp "+outfile);
	
	message("plotting single nucleotide profile...");
	plot_PKA2_nucleotide_output( outfile,  output+".nucleotide.profile.pdf",  lSeq);
	message("- done");
	
	if (pair == false) return 0;
	
	// pairs
	int seq1_len = 1;
	int seq2_len = 1;
	int dist_min = 1;
	int dist_max = lSeq - seq1_len - seq2_len;
	vector<paired_kmer> paired_kmers = generate_paired_kmers ( alphabet, seq1_len, seq2_len, dist_max,dist_min, shift_max,shift_min);
	message(to_string(paired_kmers.size())+" paired_kmers to test in total");
	// total number of tests
	int nTest1 = 0;
	for(int i=dist_min;i<= dist_max; i++)
	{
		nTest1 += pow(alphabet.size(),seq1_len) * pow(alphabet.size(),seq2_len) * (lSeq - i - seq1_len - seq2_len +1);
	}
	message(to_string(nTest1)+" tests (kmer x position) to perform");
	array<int,2> nSig2 = find_significant_pairs_from_weighted_sequences(seqs,weights, paired_kmers, output+".significant.pair", 100, pCutoff,  pCutoff_B,startPos);
	
	message(to_string(nSig2[0])+" significant paired_kmers identified");
	

	message("Done!");
	


	return 0;
} 

