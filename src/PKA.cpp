#include "sequence.h"
#include "markov.h"
#include "container.h"
#include "utility.h"
#include "text.h"

#include <boost/algorithm/string.hpp>

void print_help()
{
    string txt = "\n"
    "PKA: Positional Kmer Analysis (version 0.3)\n"
    "\n"
    "   -by Xuebing Wu (wuxbl@wi.mit.edu), Bartel lab, Whitehead Institute\n"
    "\n"
    "   Identify statistically enriched/depleted short sequences of length k (kmer)\n"
	"   at every position in a set of aligned sequences, weighted or unweighted. \n"
	"   Degenerate nucleotides and small shift in positions can be allowed.\n"
    "\n"
    "Usage: PKA input.fa [options]\n"
    "\n"
    "       PKA input-fixed-length.fa \n" 
    "       PKA input-fixed-length.fa -bgfile background.fa\n"
    "       PKA input-variable-length.fa -first 60 \n"
    "\n"
    "Options\n"
    "\n"
	"Type of analysis\n"
	"   (default)            activated by default. Identify significant kmers using Binomial test\n"
	"                        - input can be fasta/raw/tabular format\n"
	"   -ranked              Wilcoxon rank-sum test (i.e. Mann-Whitney U test) on ranked sequences \n"
	"                        - input can be fasta/raw/tabular format, but needs to be sorted\n"
	"   -weighted            Two-sample stutent's t test on weighted sequences\n"
	"                        - input can only be tabular format, does not need to be sorted\n"
	"   -predict prefix      use significant kmers from a previous run (-o prefix) to score input sequences\n"
    "Input\n"
    "   -alphabet STR        alphabet for generating kmers, default=ACGT, case insensitive\n"
    "                        note: 'dna' is equivalent to 'ACGT', 'U' will be converted to 'T'\n"
	"                              'protein' is equivalent to 'ACDEFGHIJKLMNOPQRSTUVWY'\n"
	"   -seq INT             for tabular input: sequences are in column INT. Default 1\n"
	"   -weight INT          for tabular input: weights are in column INT. Default 2\n"
	"   -skip INT            for tabular input: skip the first INT lines, such as headers. Defaut 0\n"  
    "   -first INT           only take the first INT bases of each input sequence\n"
    "   -last INT            only take the last INT bases of each input sequence\n"
    "Kmer counting\n"
    "   -k INT               use fixed kmer length INT\n"
	"   -max_k INT           consider all kmers of length 1,2,...,INT. default=4 \n"
    "   -shift INT           max shift (to right) allowed for kmer positions\n"
	"   -max_shift INT       consider shift from 0 to INT, default=0, i.e. no shift\n"	
	"   -degenerate STR      alphabet to use for degenerate kmers. Subset of all possible IUPAC DNA \n"
	"                        letters (ACGTRYMKWSBDHVN, equivalent to 'all'). One can use \n"
	"                        ACGTN to search gapped-kmers. Only work for DNA/RNA sequences\n"
	"   -gapped              allowing gapped kmer, equivalent to '-degenerate ACGTN' \n"
	"   -pair                also test all possible pairs of positional monomers\n"  	 	
    "Statistics & output\n"
    "   -o STR               prefix for all output files, default=PKA\n"
	"   -minCount NUM        minimum number of sequences to have this kmer to include in output\n"
	"                        if smaller than 1, treat as fraction of sequences (default=5)\n"
    "   -p FLOAT             p-value cut-off, default=0.05, adjusted by Bonferroni (default) or FDR\n"
    "   -FDR                 adjust p value by FDR method instead of Bonferroni. Note FDR is slower\n"
    "   -startPos INT        re-number position INT (1,2,3,..) as 1. The position before it will be -1\n"
    "   -pseudo FLOAT        pseudocount added to background counts. default=1e-9. Ignored by -markov\n"
    "Background model for unweighted & unranked sequences (ignore if using -ranked or -weighted)\n"
	"   (default)            compare to the same kmer at other positions \n"
	"   -bgfile FILE         background sequence file\n"	
    "   -markov INT          N-th order markov model trained from input or background (with -bgfile)\n"
    "                        N=0,1,or 2. Default N=1: first order captures upto di-nucleotide bias\n"
    "   -shuffle N,M         shuffle input N times, preserving M-nucleotide frequency\n"
    "   -save FILE           save shuffled sequences or the learned markov model to a file\n"
    "   -no_bg_trim          no background sequence trimming (-first/last). valid with -markov and -bgfile\n"
    "\n\n";
	cout << txt;
}

 
int main(int argc, char* argv[]) {

    ///////////////////////////////////////////////////////////////
    //       part 1: parameters and default settings             //
	///////////////////////////////////////////////////////////////

	// type of analysis
	string analysis = "default"; // or ranked or weighted
	
	// prediction mode: -predict prefix
	string prefix;
	
	// default output prefix
    string output = "PKA"; 
	
    string alphabet = "ACGT";   // default DNA
	
	// kmer length, shift, degenerate bases
    int k=0;    // k
	int min_k = 1;
    int max_k = 0;  // upto max_k
	
    int shift=0;
	int min_shift = 0;
	int max_shift = 0;
	
    bool degenerate = false;    // no degenerate bases allowed by default
	
	string degenerate_alphabet = "ACGTRYMKWSBDHVN"; 
	// if degenerate allowed, default to use all IUPAC DNA bases. to look at gappmer only, use ACGTN
	
	/*
	string protein = "ACDEFGHIJKLMNOPQRSTUVWY";
	string dna = "ACGT";
	string DNA_gap = "ACGTN";
	string DNA_all = "ACGTRYMKWSBDHVN";
	*/
	
	// tabular input format
	int cSeq=0; // -seq, sequence in this column, first column is 0
	int cWeight = -1; // -weight, -1 means no weight
	int skip=0; // skip this number of lines at the beginning
    int first=-1; //    no 3' trim
    int last=-1; //     no 5' trim
    bool no_bg_trim = false;    // only valid when -b and -markov used
	
	// statistics
    double pseudo = 1e-9;   // pseudocounts, only used when -b or -shift used
    double pCutoff=0.05;    // binomial p cutoff
	bool Bonferroni = true;
	
	// output
    int startPos= 1; // coordinates
		
	// background 
	bool local = true;
    int shuffle_N = 0;  // shuffle N times of each input sequence
    int preserve = 2;   // preserving dinucleotide
    int markov_order = -1;  // order of markov model, -1 means not used
	
	double minCount = 5.0; // minimum number of sequences to have this kmer to be reported in output
	
	bool pair = false; // test all pairs of monomers
	
	bool build_model = false;
	
    string seqfile1,seqfile2,save_to_file,str;

	///////////////////////////////////////////////////////////////
    /*******        part 2: get commandline arguments   */
	///////////////////////////////////////////////////////////////

	// total number of arguments
    if (argc < 2) 
	{
		print_help(); // if only type PKA, print help and exit
		exit(1);
	}

	// the first argument needs to be input sequence file
    seqfile1 = argv[1]; 
    if (seqfile1 == "-h" || seqfile1 == "--help" ) 
	{
		print_help(); // PKA -h or PKA --help: print help and exit
		exit(0);
	}
	else if (seqfile1[0] == '-')
	{
		print_help();
		message("ERROR: the first argument needs to be input file name! ");
		exit(1);
	}

    // other arguments
    for (int i = 2; i < argc; i++) { 
        if (i != argc) { 
            str=argv[i];
            if (str == "-bgfile") {  // background
                seqfile2 = argv[i + 1];
				local = false;
                i=i+1;
            } else if (str == "-ranked") {
                analysis = "ranked";
				cWeight = -1;
            } else if (str == "-weighted") {
                analysis = "weighted";
				if(cWeight < 0) cWeight = 1;
			} else if (str == "-predict") {
				prefix = argv[i + 1];
				i=i+1;
            } else if (str == "-o") {   // output prefix
                output = argv[i + 1];
                i=i+1;
            } else if (str == "-k") {   // fixed k
                k = atoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-max_k") {    // max k
                max_k = atoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-shift") {
                shift = atoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-max_shift") {
                max_shift = atoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-alphabet") {
                alphabet = argv[i + 1];
				if(alphabet == "ACGU") alphabet="ACGT";
                i=i+1;
            } else if (str == "-degenerate") {
                degenerate = true;
				degenerate_alphabet = argv[i+1];
				if (degenerate_alphabet == "all") degenerate_alphabet = "ACGTRYMKWSBDHVN";
				i=i+1;
            } else if (str == "-gapped") {
                degenerate = true;
				degenerate_alphabet = "ACGTN";
			} else if (str == "-weight") {
				cWeight = atoi(argv[i + 1]) - 1;
				i=i+1;
			} else if (str == "-seq") {
				cSeq = atoi(argv[i + 1]) - 1;
				i=i+1;
			} else if (str == "-skip") {
				skip = atoi(argv[i + 1]);
				i=i+1;
            } else if (str == "-first") {
                first = atoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-last") {
                last = atoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-no_bg_trim") {
                no_bg_trim = true;
            } else if (str == "-pseudo") {
                pseudo = atof(argv[i + 1]);
                i=i+1;
            } else if (str == "-p") {
                pCutoff = atof(argv[i + 1]);
                i=i+1;
            } else if (str == "-FDR") {
                Bonferroni = false;
            } else if (str == "-startPos") {
                startPos = atoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-markov") {
                markov_order = atoi(argv[i + 1]);
				local = false;
                i=i+1;
            } else if (str == "-shuffle") {
                string s(argv[i+1]);
                vector<string> ss = string_split(s,",");
                shuffle_N = atoi(ss[0].c_str());
                preserve = atoi(ss[1].c_str());
				local = false;
                i=i+1;
            } else if (str == "-minCount") {
                minCount = atof(argv[i + 1]);
                i=i+1;
            } else if (str == "-save") {
                save_to_file = argv[i + 1];
                i=i+1;
            } else {
                message("ERROR: Unknown options: "+str);
                print_help();
				exit(1);
            }
        }
    }
	
	
	if(boost::algorithm::to_lower_copy(alphabet) == "dna") 
		alphabet="ACGT";
	else if (boost::algorithm::to_lower_copy(alphabet) == "protein") 
		alphabet = "ACDEFGHIJKLMNOPQRSTUVWY";

    // determine the length of k-mer
    if (k == 0 && max_k == 0) // neither is specified, do upto 4-mers
    {
        max_k = 4;
    }
    else if ( k > 0 && max_k == 0) // 
    {  
        min_k = k;
        max_k = k; // start from k to k
    } // else do 1..max_k
	
	// determine shift
	if(max_shift == 0 && shift > 0) // -shift but not -max_shift
	{
		min_shift = shift;
		max_shift = shift;
	}
	
	// both shift and degenerate for unranked unweighted
	if(analysis == "default" && degenerate == true && max_shift > 0 && local == false) 
	{
		message("Warnining: for unweighted sequences, only under -local mode one can use both shift and degenerate bases. Switch to -local. Ignore -markov / -bgfile / -shuffle");
		local = true;
	}

    // determine background model
    if(analysis == "default" && seqfile2.size()==0 && local == false){ // no file specified using -b
        if (markov_order > -1) // markov model
        {
            shuffle_N = 0; // ignore -shuffle
            if (markov_order>2)
            {
                message("ERROR: markov_order can only be 0, 1, or 2! ");
                exit(0);
            }
        } else if (shuffle_N == 0) // no background specified, use -markov 1
        {
            markov_order = 1;
        }
    }


    // print out parameters used
        message("Summary: ");
        message("   Input       :   " + seqfile1);
    if(first>0)
	{
    	message("                   take " + to_string(first) + " bases from 5' end");
    } else if (last >0)
    {
    	message("                   take " + to_string( last) + " bases from 3' end");
    }
    if(local)
	{
		message("   Background  :   local");
	}
	else if(seqfile2.size()>0)
	{
    	message("   Background  :   " + seqfile2);
    	if(markov_order > -1) message( "                   " + \
			to_string( markov_order) +  " order markov model");
    } else 
	{
    	if(markov_order > -1) message("   Markov model:   order " + to_string( markov_order));
        else
	    {
    	message( "   Background  :   shuffle input");
    	message( "                   preserving " + to_string(preserve) + "-nuceotide frequency");
        }
    }
        message("   Output      :   " +output +".*");  
        message("   alphabet    :   " + alphabet);
        message("   min_kmer    :   " + to_string(min_k) );
        message("   max_kmer    :   " + to_string(max_k) );
        message("   min_shift   :   " + to_string(min_shift) );
        message("   max_shift   :   " + to_string(max_shift) );
		if (degenerate)
        message("   degenerate  :   " + degenerate_alphabet );
        message("   p           :   " + to_string(pCutoff ));
        message("   start at    :   " + to_string(startPos));

///////////////////////////////////////////////////////////////
    /***********    part 3: process input                 */
///////////////////////////////////////////////////////////////
		
	vector<string> names, seqs1,seqs2;
	vector<double> weights;
	if (is_fasta(seqfile1)) 
	{
		if(analysis == "weighted")
		{
			message("ERROR: You used -weighted option which requires input to be tabular instead of fasta");
			exit(1);
		}
		message("Input file is FASTA format");
		ReadFastaToVectors(seqfile1, names, seqs1);
	}
	else 
	{
		load_sequences_from_tabular(seqfile1,seqs1,weights,skip,cSeq,cWeight);
	}

    if (seqs1.size() == 0)
    {
        message("No sequence present in file: " + seqfile1);
        exit(1);
    }
    // show the number of sequences loaded
    message(to_string(seqs1.size()) + " sequences loaded from " + seqfile1);

    // trim if -first or -last specified, note that too short sequences will be discarded
    if (first >  0) seqs1 = first_n_bases(seqs1,first);
    else if (last > 0) seqs1 = last_n_bases(seqs1,last);
	if (first >0 || last > 0) message(to_string(seqs1.size()) + 
		" sequences remain after trimming back to " + 
			to_string( seqs1[0].size() ) + " bases ");

    //debug
    //WriteFasta(vector2map(seqs1),"trimmed.fa");

    // replace U with T
    if (alphabet == "ACGT")
    {
        for (int i=0;i<seqs1.size();i++)
        {
            replace(seqs1[i].begin(),seqs1[i].end(),'U','T');
        }
    }

    markov_model markov;        // background markov model
    map<string,double> kmer_probs;  // kmer probability from markov model

	if(analysis == "default" && local == false)
	{
	    // if background not available, use markov model or shuffle the foreground
	    if (seqfile2.size() == 0) // no background file specified using -b
	    {
	        if (markov_order < 0) // shuffle
	        {
	            // shuffle sequence preserving m-let frequency, m specified by -preserve
	            seqs2 = shuffle_seqs_preserving_k_let(seqs1,shuffle_N,preserve);
			    // replace U with T
			    if (alphabet == "ACGT")
			    {
			        for (int i=0;i<seqs2.size();i++)
			        {
			            replace(seqs2[i].begin(),seqs2[i].end(),'U','T');
			        }
			    }

	            message(to_string(seqs2.size()) + " shuffled sequences generated, " + to_string( shuffle_N ) + " from each input sequence");

	            if (save_to_file.size() > 0) // save shuffled sequences
	            {        
	                WriteFasta(vector2map(seqs2),save_to_file);
	                message("Shuffled sequences saved to: "+save_to_file);
	            }
	        } else // generate markov model from input
	        {
	            markov = markov_model(markov_order,alphabet,seqs1);
	        }
	    }
	    // else load background sequence
	    else 
	    {
			if(is_fasta(seqfile2)) ReadFastaToVectors(seqfile2, names, seqs2);
			else load_sequences_from_tabular(seqfile2,seqs2,weights,skip,cSeq,cWeight);
	
	        if (seqs2.size() == 0)
	        {
	            message("No sequence present in file: " + seqfile2);
	            exit(1);
	        }
	        // show the number of sequences loaded
	        message(to_string(seqs2.size()) + " sequences loaded from " + seqfile2);
        
	        if (no_bg_trim == false || markov_order < 0)
	        {
	            // trim if -first or -last specified, note that too short sequences will be discarded
	            if (first >  0) seqs2 = first_n_bases(seqs2,first);
	            else if (last > 0) seqs2 = last_n_bases(seqs2,last);
	            if (first>0 || last>0) message( to_string(seqs2.size()) \
					+ " background sequences trimmed back to " \
						+ to_string(seqs2[0].size()) + " bases ");
	        }
    
	        if (markov_order > -1)
	        {
	            markov = markov_model(markov_order,alphabet,seqs2);
	        }
	    }

	    if(markov_order > -1 && save_to_file.size()>0) 
	    {   
	        markov.print(save_to_file);
	        message("Markov model saved to: "+save_to_file);
	    }
	}

    // total number of sequences loaded or after trimming
    // note that too short sequences after trimming will be discarded
    int nSeq1 = seqs1.size();
    int nSeq2;
    int seq_len1 = seqs1[0].size();
    int seq_len2;

    // make sure all sequences have the same size, only necessary when neither -first nor -last is used
    if (first < 0 && last < 0)
    { 
        // sequences in foreground
		vector<int> removed = filter_sequences_by_size(seqs1);
		nSeq1 = seqs1.size();
		if (removed.size() > 0) 
		{
			message(to_string(nSeq1)+" sequences left after filtering by size");
			if(analysis == "weighted")
			{
				for(int i=0;i<removed.size();i++)
				{
					weights.erase(weights.begin()+removed[i]);
				}
				// save the data
				ofstream tmpout("tmp.txt");
				for(int i=0;i<seqs1.size();i++)
					tmpout << seqs1[i] << "\t" << weights[i] << endl;
				tmpout.close();
			}
		}
        
        if (seqfile2.size()>0 && markov_order <0) // only check other sequences in background if not genrated by shuffling
        {
            filter_sequences_by_size(seqs2,seqs1[0].size());
        }  
    }


///////////////////////////////////////////////////////////////
//         part 4:  prediction mode
///////////////////////////////////////////////////////////////
	
	if(prefix.size() > 0) 
	{
		message("=== prediction mode ===");
		/**/
		message("building model from file: "+prefix+".pass.p.cutoff.txt");		
		vector<positional_kmer> ranked_kmers = build_model_from_PKA_output(prefix+".pass.p.cutoff.txt",startPos);
		message("- " + to_string(ranked_kmers.size())+" positional kmers included in the model");
		save_model_to_file(ranked_kmers, "model.txt");
		
	
		// save the model to file
		// save_model_to_file(ranked_kmers, output+".model.txt");	
	
		/**/
		message("scoring input sequences using the model...");
		string scoreFile = prefix+".score";
	    ofstream out1(scoreFile.c_str());
		//load_weighted_sequences_to_vectors(inputfile,seqs,weights,skip,cSeq,cWeight);
		for (int i=0;i<seqs1.size();i++)
		{
			//cout << i << seqs1[i] << endl;
			double score = score_sequence_using_PKA_model(ranked_kmers, seqs1[i]);
			out1 << seqs1[i] << "\t" << weights[i] << "\t" << score << endl;		
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
		for (int i=0;i<seqs1.size();i++)
		{
			double score = score_sequence_using_paired_kmer_model(paired_kmer_model, seqs1[i]);
			out << seqs1[i] << "\t" << weights[i] << "\t" << score << endl;		
		}	
		out.close();
		message("- done");
		
		
		return 0;
	}
	


/*
// debug
implant_motif(seqs1, 10, "G", 0.1);
implant_motif(seqs1, 20, "G", 0.5);
implant_motif(seqs1, 30, "CNC", 0.1);
implant_motif(seqs1, 40, "CNC", 0.5);
implant_motif(seqs1, 50, "GHG", 0.1);
WriteFasta(seqs1,"implanted.fa");
*/

    //debug
    //PrintMap(seqs1);
	

///////////////////////////////////////////////////////////////
    /********   part 5: generate kmers      */
///////////////////////////////////////////////////////////////
	
	// total number of tests to be performed, ~ n_kmer * seq_len * shift
	// to be used in multiple testing correction
	int nTest = 0; 
	
	// generate all exact kmers	
	vector<string> kmers = generate_kmers(min_k, alphabet);
	if (degenerate == false) nTest = kmers.size() * (seq_len1 - min_k + 1);
	for (k = min_k+1; k <= max_k; k++)
	{
		vector<string> tmp = generate_kmers(k, alphabet);	
		kmers += tmp;
		if (degenerate == false) nTest += tmp.size() * (seq_len1 - k + 1);
	}
	message(to_string(kmers.size()) +  " exact kmers to be tested" );
		
	vector<string> dkmers;
	// generate degenerate kmers, which will also include exact kmers
	if (degenerate)
	{
		dkmers = degenerate_kmer(min_k,degenerate_alphabet);
		nTest = dkmers.size() * (seq_len1 - min_k + 1);
		for (k = min_k+1; k <= max_k; k++)
		{
			vector<string> tmp = degenerate_kmer(k,degenerate_alphabet);	
			dkmers += tmp;
			nTest += tmp.size() * (seq_len1 - k + 1);
		}		
		message(to_string(dkmers.size()) + " k-mers allowing degenerate bases");
	}
	
	// multiple by shift
	nTest = nTest * (max_shift - min_shift + 1);
	message(to_string(nTest) + " tests (kmer x position x shifts) will be performed");


///////////////////////////////////////////////////////////////
//			part 6: search for significant kmers	
///////////////////////////////////////////////////////////////
	
	if(minCount < 1) minCount = seqs1.size() * minCount;
	
    // number of significant kmers
    int nSig = 0; // based on uncorrected and corrected p-value

	// final output file
	string out = output+".pass.p.cutoff.txt";
	
    // tmp output file for significant kmers, unsorted
    string outtmp = out+".tmp";

    // tmp output file for significnt kmers' frequency, for B significant kmers
    string output_freq = output+".Bonferroni.significant.frequency.txt";


	if(analysis != "default")
	{
		// output file
		
		if (analysis == "ranked") 
		{
			if(degenerate) nSig = find_significant_kmer_from_ranked_sequences(
				seqs1, dkmers,outtmp, nTest, pCutoff, Bonferroni, min_shift, max_shift,startPos,minCount);
			else nSig = find_significant_kmer_from_ranked_sequences(
				seqs1, kmers,outtmp, nTest, pCutoff, Bonferroni, min_shift, max_shift,startPos,minCount);
		}
		else 
		{
			if(degenerate) nSig = find_significant_kmer_from_weighted_sequences(
				seqs1, weights, dkmers, outtmp, nTest, pCutoff, Bonferroni, min_shift, max_shift, startPos,minCount);
			else nSig = find_significant_kmer_from_weighted_sequences(
				seqs1, weights, kmers, outtmp, nTest, pCutoff, Bonferroni, min_shift, max_shift, startPos,minCount);
		}
	
		message( to_string (nSig) +  " positional kmers with p <"+to_string(pCutoff));
	
		if (pair)
		{
			// pairs
			int seq1_len = 1;
			int seq2_len = 1;
			int dist_min = 1;
			int dist_max = seq_len1 - seq1_len - seq2_len;
			vector<paired_kmer> paired_kmers = generate_paired_kmers ( 
				alphabet, seq1_len, seq2_len, dist_max,dist_min, max_shift,min_shift);
			message(to_string(paired_kmers.size())+" paired_kmers to test in total");
			// total number of tests
			int nTest1 = 0;
			for(int i=dist_min;i<= dist_max; i++)
			{
				nTest1 += pow(alphabet.size(),seq1_len) * pow(alphabet.size(),seq2_len)*(seq_len1 - \
					i - seq1_len - seq2_len +1);
			}
			message(to_string(nTest1)+" tests (kmer x position) to perform");
			int nSig2 = find_significant_pairs_from_weighted_sequences(
				seqs1,weights, paired_kmers, output+".significant.pair", 100, pCutoff,Bonferroni, startPos,minCount);
	
			message(to_string(nSig2)+" significant paired_kmers identified");
		}
	}
	else
	{
	    // prepare output files
	    ofstream fout;
	    fout.open(out.c_str());
	    // data format
	    fout.precision(3);
	    // header
	    string header = "kmer\tposition\tshift\tz_score\tp\tBonferroni\tfrac1\tfrac2\tratio\tlocal_r\tFDR";
	    fout << header << endl;
	    fout.close();
    
		if(local)
		{
			if(degenerate) nSig = find_significant_degenerate_shift_kmer_from_one_set_unweighted_sequences(
				seqs1, dkmers,outtmp, nTest, pCutoff, Bonferroni, min_shift, max_shift,startPos,minCount);
			else nSig = find_significant_degenerate_shift_kmer_from_one_set_unweighted_sequences(
				seqs1, kmers,outtmp, nTest, pCutoff,Bonferroni, min_shift, max_shift,startPos,minCount);
		}
	    else if (markov_order < 0){
	        // find significant kmers by comparing two set of sequences
	        nSig = find_significant_kmer_from_two_seq_sets(
				seqs1,seqs2,kmers,dkmers,min_shift,max_shift,degenerate,
			pCutoff,Bonferroni,pseudo,startPos,nTest,outtmp,output_freq);
	    }else{
	        // find significant kmers using markov model as background
			kmer_probs = markov.probs(kmers);
	        nSig = find_significant_kmer_from_one_seq_set(
				seqs1,kmer_probs,kmers,dkmers,min_shift,max_shift,degenerate,
			pCutoff,Bonferroni, startPos,nTest,outtmp,output_freq); 
	    }

	    message( to_string (nSig) +  " significant positional kmers identified in total"); 

	}
	
	if (nSig == 0) return 0;
	
///////////////////////////////////////////////////////////////
//			part 7: post-processing	
///////////////////////////////////////////////////////////////
	
	if(Bonferroni == false)
	{
		////   calculating FDR
    	message("adjusting p-values using FDR...");
    	string script = 
        	"x=read.table('"+ outtmp +"',header=F)   \n"
			"adjusted = p.adjust(10^(-x[,5]),method='fdr',n="+ to_string(nTest) +")\n"
			"sub = adjusted < "+to_string(pCutoff)+"\n"
        	"x = cbind(x[sub,],-log10(adjusted[sub])) \n"
        	"write.table(x,file='"+outtmp+"',sep='\\t',col.names=F,row.names=F,quote=FALSE) ";
    	R_run(script);
		int n = count_lines(outtmp);
		message(to_string(n)+" significant kmers with FDR < "+to_string(pCutoff));
	}
	
	// sort output by p.value
	system_run("sort -k5,5gr "+outtmp+" > "+out);

	// most significant at each position
    system_run(" cat "+out +" | sort -k2,2g -k5,5gr  > "+outtmp);
    remove_duplicates(outtmp,output+".most.significant.each.position.txt",2,1,"");	
	
	/*
    system_run(" cat "+out +" | awk '$4>0' | sort -k2,2g -k5,5gr  > "+outtmp);
    remove_duplicates(outtmp,outtmp+".most.enriched",2,1,"");	
    system_run(" cat "+out +" | awk '$4<0' | sort -k2,2g -k5,5gr  > "+outtmp);
    remove_duplicates(outtmp,outtmp+".most.depleted",2,1,"");	
	system_run(" cat "+outtmp+".most* > "+output+".most.significant.each.position.txt");
	*/
	
    //remove intermediate files
    system_run("rm "+outtmp+"*");
	
    message("plotting the most significant kmer at each position...");

	// column to plot
	int col = 6; // Bonferroni
	if (Bonferroni == false)
	{
	    if (analysis == "ranked") col = 7;
		else if (analysis == "weighted") col = 13;
		else col = 11;
	}
	
    string plotfilename = output+".most.significant.each.position.pdf";

    plot_most_significant_kmers(output+".most.significant.each.position.txt", output+".most.significant.each.position.pdf", seq_len1, col,startPos);
	
	// if monomer is included in the analysis
	if(min_k < 2) 
	{
		message("plotting single nucleotide profile...");
		plot_nucleotide_profile( out,  output+".nucleotide.profile.pdf",  seq_len1, col, startPos);
	}
	
    message("Done!");

    return 0;
} 

