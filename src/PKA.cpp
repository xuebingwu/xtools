#include "sequence.h"
#include "markov.h"
#include "container.h"
#include "utility.h"
#include "text.h"

#include <boost/algorithm/string.hpp>

/*
to compile:

gcc -O3 -c ushuffle.c
g++ -std=c++11 -I /lab/bartel1_ata/wuxbl/lib/boost_1_57_0/ seqlib.h stat.h positional_kmer_analysis.cpp -o PKA ushuffle.o -lboost_regex-mt
    
*/

void print_help()
{
    string txt = "\n"
    "PKA: Positional Kmer Analysis (version 0.3)\n"
    "\n"
    "   -by Xuebing Wu (wuxbl@wi.mit.edu)\n"
    "       Bartel lab, Whitehead Institute\n"
    "\n"
    "   Identify statistically enriched short sequences of length k (kmer)\n"
	"   at every position in a set of aligned sequences, compared to \n"
	"   background sequences, shuffled sequences, or a background markov model. \n"
	"   Degenerate nucleotides and shift in positions can be allowed\n"
    "\n"
    "Usage: PKA input.fa [options]\n"
    "\n"
    "   example usage:\n"
    "\n"
    "       PKA input-fixed-length.fa \n" 
    "       PKA input-fixed-length.fa -b background.fa\n"
    "       PKA input-variable-length.fa -first 60 \n"
    "\n"
    "Options\n"
    "\n"
    "Input/foreground\n"
    "   -alphabet <ACGT>     alphabet for generating kmers, default=ACGT, case insensitive\n"
    "                        note: 'dna' is equivalent to 'ACGT', \n"
	"                              'protein' is equivalent to 'ACDEFGHIJKLMNOPQRSTUVWY'\n"
    "   -first <n>           take the first n bases of each input sequence\n"
    "   -last <n>            take the last n bases of each input sequence\n"
    "Background\n"
    "   -markov <N>          N-th order markov model trained from input or background (with -b)\n"
    "                        N=0,1,or 2. Default N=1: first order captures upto di-nucleotide bias\n"
    "   -b <background.fa>   background sequences in a FASTA file\n"
    "   -shuffle <N,M>       shuffle input N times, preserving M-nucleotide frequency\n"
    "   -save <filename>     save shuffled sequences or the learned markov model to a file\n"
    "   -no_bg_trim          no background sequence trimming. valid with -markov and -b\n"
    "Kmer counting\n"
    "   -upto <K>            consider all kmers of length 1,2,...,K. default=4 \n"
    "   -k <K>               use fixed kmer length (K)\n"
    "   -shift <0>           max shift (to right) allowed for kmer positions, default=0\n"
	"   -degenerate <ACGTRYMKWSBDHVN> Alphabet to use for degenerate kmers. By default \n"
	"                        all IUPAC letters (ACGTRYMKWSBDHVN) are used. One can use \n"
	"                        ACGTN to search gapped-kmers. Only work for DNA/RNA\n"
    "Statistics & output\n"
    "   -o <output-prefix>   prefix for all output files, default=PKA_output\n"
    "   -pCutoff <p>         raw p-value cut-off, default=0.01\n"
    "   -pCutoff_B <p>       Bonferoni corrected p-value cut-off, default=0.05\n"
	"	-zCutoff <z>		 z-score cutoff, default=2.0\n"
    "   -topN <N>            output top N enriched kmers at each position, default=1\n"
    "   -startPos <n>        set position n (0,1,2,3,..) to be the start position (0) in the output\n"
    "   -highlight <l,r>     highlight the region [l,r] in the plot (after adjusting -startPos) \n"
    "   -pseudo <f>          pseudocount added to background counts. default=1e-9. Ignored by -markov\n"
//    "   -plot <n>            column to plot as y-axis for top enriched kmers (default=5, see output)\n"
    "\n"
    "Outout\n"
    "\n"
    "   The following files have the same format (*: output prefix specified using -o; N: top N kmers at each position specified using -topN) \n"
    "   1. *.binomial.significant.txt: all kmer-posiiton combinations that pass raw p-value cutoff \n"
    "   2. *.binomial.significant.txt.non-overlapping: non-overlapping kmers\n"
    "   3. *.top.1.per.position.txt: the top N enriched kmers at each position \n"
    "\n"
    "   Columns are:\n"
    "   1. kmer: kmer sequence, if degenerate will be IUPAC format\n"
    "   2. regex: kmer in regex format, i.e. GHG will be G[ACT]G\n"
    "   3. position: kmer position, 1-based\n"
    "   4. frac1: fraction of foreground sequence with the kmer at this position\n"
    "   5. frac2: fraction of background sequence with the kmer at this position\n"
    "   6. ratio: frequency ratio, foreground/background\n"
    "   7. local_r: ratio of counts at this position over the mean counts at other positions, foreground only\n"
    "   8. z_score: z-score from nomal approximation of binomial distribution\n"
    "   9. p: p-value, binomial test\n"
    "   10. Bonferoni: p-value after Bonferoni correction (conservative)\n"
    "   11. FDR: q-value (FDR-adjusted p-value, less cconservative)\n"
    "\n"
    "   *.Bonferoni.significant.frequency.txt\n"
    "       kmer:position followed by frequency at each position\n"
    "       only include kmer pass Bonferoni corrected p-value cutoff\n"
    "\n"
    "   *.meme: MEME format motif file for each Bonferoni significant kmers\n"
    "   *.top.1.per.position.pdf: plot position as x and -log10(p) as y\n"
    "   *.top.1.per.position.logo.pdf: plot sequence logo of kmer extend 3 bases each side\n"
    "   *.Bonferoni.significant.frequency.pdf: plot position as x and frequency as y\n"
    "\n"
    "How does PKA work\n"
    "\n"
    "A binomial test is used to calculate p-value for each kmer at each position: \n"
    "   1. the number of trials: total number of input sequences \n"
    "   2. the number of success: number of sequences with this kmer at this position\n"
    "   3. probability of success: background probability of this kmer at this position\n"
    "The background probability can be estimated in several ways\n"
    "   1. (default) a first-order markov model capturing di-nucleotide bias in input. This is the default and no option needed. \n"
    "   2. markov model learned from background sequences: -b filename -markov 1 (or 0 or 2). \n"
    "   3. frequency in background sequences: -b filename only \n"
    "   4. frequency in shuffled input sequences preserving certain order of nucleotide bias: -shuffle N,M\n"
    "\n";
	cout << txt;
}

 
int main(int argc, char* argv[]) {

    /*******       part 1: default parameters       */

    string output = "PKA_output"; // output prefix
    string alphabet = "ACGT";   // default DNA
    int k=0;    // k
    int max_k = 0;  // upto max_k
    int shift=0; // no shift, i.e. exact position
    bool degenerate = false;    // no degenerate bases allowed by default
	string degenerate_alphabet = "ACGTRYMKWSBDHVN"; // if degenerate allowed, default to use all IUPAC bases. to look at gappmer only, use ACGTN
	/*
	string protein = "ACDEFGHIJKLMNOPQRSTUVWY";
	string dna = "ACGT";
	string DNA_gap = "ACGTN";
	string DNA_all = "ACGTRYMKWSBDHVN";
	*/
    int first=-1; //    no 3' trim
    int last=-1; //     no 5' trim
    bool no_bg_trim = false;    // only valid when -b and -markov used
    double pseudo = 1e-9;   // pseudocounts, only used when -b or -shift used
    double pCutoff=0.01;    // binomial p cutoff
    double pCutoff_B=0.05;  // Bonferoni p cutoff
    int topN=1; // output the best kmer (smallest p, when p=0, use larger zscore)
    int startPos=0; // coordinates
    int shuffle_N = 0;  // shuffle N times of each input sequence
    int preserve = 2;   // preserving dinucleotide
    int markov_order = -1;  // order of markov model
    string left="0";
    string right="0";
	
	bool build_model = false;

    string seqfile1,seqfile2,save_to_file,str;

    /*******        part 2: get commandline arguments   */

    if (argc < 2) 
	{
		print_help(); // if only type PKA, print help and exit
		exit(1);
	}

    seqfile1 = argv[1]; // the first argument is input sequence file

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
            if (str == "-b") {  // background
                seqfile2 = argv[i + 1];
                i=i+1;
            } else if (str == "-o") {   // output prefix
                output = argv[i + 1];
                i=i+1;
            } else if (str == "-k") {   // fixed k
                k = atoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-upto") {    // max k
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
				degenerate_alphabet = argv[i+1];
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
            } else if (str == "-markov") {
                markov_order = atoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-shuffle") {
                string s(argv[i+1]);
                vector<string> ss = string_split(s,",");
                shuffle_N = atoi(ss[0].c_str());
                preserve = atoi(ss[1].c_str());
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
            } else if (str == "-build_model") {
                build_model = true;
            } else {
                message("ERROR: Unknown options: "+str);
                print_help();
				exit(1);
            }
        }
    }
	
	

    // cant shift degenerate motifs
    if (degenerate && shift>0) 
    {
        message("WARNING: currently no support for shifted degenerate kmers! i.e. when -degenerate is set, shift must be 0");
        exit(0);
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

    // determine background model
    if(seqfile2.size()==0){ // no file specified using -b
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
    if(seqfile2.size()>0)
	{
    	message("   Background  :   " + seqfile2);
    	if(markov_order > -1) message( "                   " + to_string( markov_order) +  " order markov model");
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
        message("   kmer from   :   " + to_string(min_k) );
        message("   kmer upto   :   " + to_string(max_k) );
        message("   shift       :   " + to_string( shift) );
		if (degenerate)
        message("   degenerate  :   " + degenerate_alphabet );
        message("   pCutoff     :   " + to_string(pCutoff ));
        message("   pCutoff_B   :   " + to_string(pCutoff_B ));
        message("   top N       :   " + to_string(topN ));
        message("   start at    :   " + to_string(startPos));

    /***********    part 3: process input                 */

    // load foreground sequence 
    map<string,string> seqs1 = ReadFasta(seqfile1);

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
			to_string( seqs1.begin()->second.size() ) + " bases ");

    //debug
    //WriteFasta(seqs1,"trimmed.fa");

    // replace U with T
    if (alphabet == "ACGT" || alphabet == "ACGU")
    {
        alphabet = "ACGT";
        for (map<string,string>::iterator it=seqs1.begin(); it!=seqs1.end(); ++it)
        {
            replace(it->second.begin(),it->second.end(),'U','T');
        }
    }

    
    map<string,string> seqs2; // background sequences
    markov_model markov;        // background markov model
    map<string,double> kmer_probs;  // kmer probability from markov model

    // if background not available, use markov model or shuffle the foreground
    if (seqfile2.size() == 0) // no background file specified using -b
    {
        if (markov_order < 0) // shuffle
        {
            // shuffle sequence preserving m-let frequency, m specified by -preserve
            seqs2 = shuffle_seqs_preserving_k_let(seqs1,shuffle_N,preserve);

            message(to_string(seqs2.size()) + " shuffled sequences generated, " + to_string( shuffle_N ) + " from each input sequence");

            if (save_to_file.size() > 0) // save shuffled sequences
            {        
                WriteFasta(seqs2,save_to_file);
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
        seqs2 = ReadFasta(seqfile2);

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
            if (first>0 || last>0) message( to_string(seqs2.size()) + " background sequences trimmed back to " + to_string(seqs2.begin()->second.size()) + " bases ");
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

    // total number of sequences loaded or after trimming
    // note that too short sequences after trimming will be discarded
    int nSeq1 = seqs1.size();
    int nSeq2;
    int seq_len1 = seqs1.begin()->second.size();
    int seq_len2;

    // make sure all sequences have the same size, only necessary when neither -first nor -last is used
    if (first < 0 && last < 0)
    { 
        // sequences in foreground
        for (map<string,string>::iterator it=seqs1.begin(); it!=seqs1.end(); ++it) 
        {
            if(it->second.size() != seq_len1)
            {
                message(seqfile1 + ": The length of the following sequence differs from that of the first:" + it->first + ":" + it->second);
                exit(1);
            }
        }
        if (seqfile2.size()>0 && markov_order <0) // only check other sequences in background if not genrated by shuffling
        {
            seq_len2 = seqs2.begin()->second.size();
             // the first sequence
            if (seq_len1 != seq_len2)
            {
                message("sequences from the two input files should be of the same length!!");
                exit(1);
            }   
            for (map<string,string>::iterator it=seqs2.begin(); it!=seqs2.end(); ++it) 
            {
                if(it->second.size() != seq_len2)
                {  
                    message(seqfile2 + ": The length of the following sequence differs from that of the first:" + it->first + ":" + it->second);
                    exit(1);
                }
            }
        }  
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

    /********   part 4: generate kmers      */

	int nTest = 0; // total number of tests to be performed, ~ n_kmer * seq_len
	
    // generate degenerate kmers, which will also include exact kmers
    vector<string> dkmers;
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
    	message(to_string(dkmers.size()) + " k-mers in total");
	}
	
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
    
    message(to_string(nTest) + " tests (kmer x position) will be performed");

    /******     calculate kmer probability from markov model */
    kmer_probs = markov.probs(kmers);

	
	
    /******     part 5: kmer counting and statistics        */
    
    // prepare output files
    string out = output+".binomial.significant.txt";
    ofstream fout;
    fout.open(out.c_str());
    // data format
    fout.precision(3);
    // header
    string header = "kmer\tregex\tposition\tfrac1\tfrac2\tratio\tlocal_r\tz_score\tp\tBonferoni\tFDR";
    fout << header << endl;
    fout.close();
 
    // tmp output file for significant kmers, unsorted
    string outtmp = out+".tmp";

    // tmp output file for significnt kmers' frequency, for B significant kmers
    string output_freq = output+".Bonferoni.significant.frequency.txt";

    // number of significant kmers
    array<int,2> nSig= {0,0}; // based on uncorrected and corrected p-value

    if (markov_order < 0){
        // find significant kmers by comparing two set of sequences
        nSig = find_significant_kmer_from_two_seq_sets(seqs1,seqs2,kmers,dkmers,shift,degenerate,pCutoff,pCutoff_B,pseudo,startPos,nTest,outtmp,output_freq);
    }else{
        // find significant kmers using markov model as background
        nSig = find_significant_kmer_from_one_seq_set(seqs1,kmer_probs,kmers,dkmers,shift,degenerate,pCutoff,pCutoff_B,startPos,nTest,outtmp,output_freq); 
    }

    message( to_string (nSig[0]) +  " significant positional kmers identified in total"); 
    message( to_string(nSig[1]) + " remain after Bonferoni multiple testing correction" );

    /*****      part 6: sorting     */

    string cmd;

    /***    calculating q-value, i.e. FDR correction     */
    if (nSig[0]>0)
    {
        message("adjusting p-values using FDR...");

        string script = 
            "x=read.table('"+ outtmp +"',header=F)   \n"
            "x = cbind(x,p.adjust(x[,9],method='fdr',n="+ to_string(nTest) +")) \n"
            "write.table(x,file='"+outtmp+"',sep='\\t',col.names=F,row.names=F,quote=FALSE) ";
    
        R_run(script);

    }

    if (nSig[0]>0)
    {			

        message("sorting significant kmers by p-value then by z-scores...");
        // -g: allow expoential such as 1.3e-4
		system_run("cp "+out+" "+out+".enriched.and.depleted");
        system_run("sort -k9,9g -k8,8gr "+outtmp+" >> "+out+".enriched.and.depleted");

        // from now on only enriched motifs
		system_run("awk '$8>0' " + out+".enriched.and.depleted > "+out);

        // non-overlapping motifs
        // note this is not perfect, two motif could be overlapping yet distinct
        // ideally should look at raw sequence and see if highly correlated
        message("removing overlapping motifs..." );
        int total_non_overlapping = non_overlapping_sig_motifs(out,out+".non-overlapping");
        message(to_string ( total_non_overlapping) + " non-overlapping motifs found" );

        message("writing the most enriched kmer at each position...");
        // sort by position, then by p-value, then by zscore
        system_run(" cat "+outtmp+" | sort -g -k3,3g -k9,9g -k8,8gr  > "+outtmp+".tmp");
        // copy header
        system_run("head -n 1 "+out+" > "+output+".top."+to_string(topN)+".per.position.txt");
        // only keep topN
        remove_duplicates(outtmp+".tmp",outtmp+".tmp2",3,topN,"");
        // save 
        system_run("cat "+outtmp+".tmp2 >> "+output+".top."+to_string(topN)+".per.position.txt");

        //remove intermediate files
        system_run("rm "+outtmp+"*");
    }
    else
    {
        // if no significant kmers found, remove intermediate files
        system_run("rm "+output+"*");
    }

    if(nSig[0]>0)
    {
        /***    part 7: plot using R    */
        message("plotting the most enriched kmer at each position...");

        string plotfilename = output+".top."+to_string(topN)+".per.position.pdf";

        string script =
        "x=read.table('"+ output + ".top."+to_string(topN)+".per.position.txt',header=T)\n"
        "show=9 # which column to plot. 9: raw p-value, 10: corrected p-value, 11: FDR, 8: z-score\n"
        "# for p=0, fit using z-score \n"
        "x[x[,show] == 0,show] = 1e-100   \n"
        "x[,show] = -log10(x[,show])  \n"
        "x[x[,show] == 100,show] = 0.5*x[x[,show] == 100,8]^2+0.008002*x[x[,show] == 100,8]+4.492  \n"
        "pdf('"+plotfilename+"',width=10,height=6)\n"
        "miny = min(x[,show])\n"
        "maxy = max(x[,show])\n"
        "plot(x[,3],x[,show],ylim=c(miny,1.4*maxy),type='h',cex=2,lwd=2,bty='n',col='gray',main='the most enriched kmer at each position',xlab='position',ylab='-log10(p)') \n"
        "y=x[x[,11]<"+to_string(pCutoff_B)+",]\n" // FDR
        "lines(y[,3],y[,show],type='h',lwd=2,col='blue')\n"
        "y=x[x[,10]<"+to_string(pCutoff_B)+",]\n" // Bonferoni
        "lines(y[,3],y[,show],type='h',lwd=2,col='red')\n"
        "legend('topleft',lty=1,lwd=2,col=c('gray','blue','red'),legend=c('Binomial p<0.05','FDR-corrected p<0.05','Bonferoni-corrected p<0.05'),bty='n',horiz=T)\n"
        "text(x[,3],x[,show],labels=paste(x[,3],x[,1],sep=', '),pos=3,srt=90,adj=c(0,NA),offset=2,cex=0.6) \n"
        "if("+left+"<"+right+"){ rect("+left+",-10,"+right+",miny-0.01,col='green',border='green')} \n"
        "dev.off()  \n";

        R_run(script);

    }

    // make logo for top N motif passing corrected p-value cutoff
    if(nSig[1]>0 && alphabet == "ACGT")
    {
        message("creating motif logo for the top "+to_string(topN)+" kmers for each position...");
        create_logo_for_topN_sig_kmer_per_position(seqs1,output+".top."+to_string(topN)+".per.position.txt",3, pCutoff_B,startPos,output);

        // merge individual motif pdf file to a single pdf file, remove intermediate files
        system_run("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile="+output+".top."+to_string(topN)+".per.position.logo.pdf "+output+"*meme.pdf");
        system_run("rm *meme.pdf;rm *.eps");

        message("plot the frequency profile of kmers passing corrected p-value cutoff");
        plot_frequency_for_significant_kmer(output_freq,output+".Bonferoni.significant.frequency.pdf");

    }

	/*
	vector<string> sig_kmers;
	vector<int> sig_positions;
	read_significant_positional_kmer_from_file(out, sig_kmers, sig_positions);
	significant_feature_matrix(seqs1, sig_kmers, sig_positions, output + ".sig.feature.txt", "pos",false, shift);
	if (seqs2.size()>0)
		significant_feature_matrix(seqs2, sig_kmers, sig_positions, output + ".sig.feature.txt", "neg",true, shift);
*/
	//save_feature_matrix(seqs1,kmers, output + ".feature.txt", "1", false, shift);
	//if (seqs2.size()>0)
	//	save_feature_matrix(seqs2,kmers, output + ".feature.txt", "0", true, shift);
		
	// build model
	if(build_model == false) return 0;
	
	vector<positional_kmer> ranked_kmers = build_model_from_PKA_output(out);
	message(to_string(ranked_kmers.size())+" significant positional kmers are used to build the model");
	
	// save the model to file
	save_model_to_file(ranked_kmers, output+".model.txt");	
	/*
	message("scoring sequences using the model...");
	//load_weighted_sequences_to_vectors(inputfile,seqs,weights,skip,cSeq,cWeight);
	for (int i=0;i<seqs.size();i++)
	{
		double score = score_sequence_using_PKA_model(ranked_kmers, seqs[i]);
		cout << seqs[i] << "\t" << weights[i] << "\t" << score << endl;		
	}
	*/
    message("Done!");

    return 0;
} 

