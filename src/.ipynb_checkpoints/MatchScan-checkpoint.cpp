#include "utility.h"
#include "iostream"
#include "sequence.h"
#include "text.h"
#include "container.h"
#include "stat.h"

#include <filesystem>

// for vector sum
#include <boost/range/numeric.hpp>

// for getting filename from a path
//#include <boost/filesystem.hpp>

using namespace std;

void help()
{
    string str =
	"\n"
    "MatchScan\n"
	"    A program to find significant matches between a query sequence and a set of target sequences,\n"
	"    where the frequency of the matched sequences in targets correlates with target scores\n"
    "    If no query given, will test all k-mers\n"
	"\n"
	"    - Xuebing Wu (wuxb07@gmail.com)\n"
    "\n"
    "Usage\n"
	"    MatchScan -q query.fasta -t target.fasta -s score.tab -o output [ options ]\n"
    "\n"
    "Options:\n"
    "\n"
    "    -i  file   Input file (sequence in column 1 and score in column 2)\n"	
 	"    -o  file   prefix for output files (tabular data and PDF figure). Default:output\n"
    "    -q  file   fasta file containing either a long sequence or multiple kmer sequences. \n "
    "               If not given, will test all possible kmers with length(s) specified by -k \n"
    "    -t  file   fasta file containing all possible target sequences\n"
    "    -s  file   tab-delimited score file (target name in column 1, score in column 2)\n"
	"               target name should match those in target fasta file, can be a subset\n"
	"    -h  file   query sequences to highlight in green\n"
	"    -k  num(s) length of subsequence to scan at each step (i.e. k-mer)\n"
	"               e.g. to search for hexamer: -k 6; for k from 6 to 10 (default): -k 6,10\n"		
	"    -m  num    max number of target sequences to use (top and bottom)\n"		
	"               default 100000, i.e. top 50000 and bottom 50000 based on scores\n"
    "    -r         match by reverse complement. Default is identical match\n"
	"    -c  seqs   cdf plots for a list of k-mers, separated by comma, such as CAAACC,TTTACG\n"
	"               when specified, will ignore all options except -q,-t,-s, and -o\n"
    "    -d         score sequence using median difference instead of -log10(p)\n"
    "\n"
	"Example\n"
	"    MatchScan -q TERC.fa -t TERC-ChIRP.fa -s TERC-ChIRP.score.txt -o TERC\n"
	"    MatchScan -q TERC.fa -t TERC-ChIRP.fa -s TERC-ChIRP.score.txt -o TERC -k 10 \n"
	"    MatchScan -q TERC.fa -t TERC-ChIRP.fa -s TERC-ChIRP.score.txt -o TERC -h terc_probe.fa.rc \n"
    "    MatchScan -i HEK1.txt -o tmp -q ~/scripts/C++/test/MatchScan/human-mirna-seed.fa -r       \n";
    "Scoring sequences using motifs from a previous run (fixed k):\n"
    "    MatchScan -t <target_fasta> -o <finished_run>\n";
    cerr << str;

    exit(0);
}


// 
// more miR_Family_Info.txt  | grep hsa | cut -f -2 | sort | uniq | awk '{print ">"$1"\n"$2}' > human-mirna-seed.fa


int main(int argc, char* argv[]) {

    // default
    string queryFile="";
    string targetFile="";
    string scoreFile="";
	string inputFile="";
	string outputFile="output";
    string highlightFile="";
	string backgroundFile="";
	
	string kmers = "";
	
	int max_target_num = 100000;

	int k_min,k_max;
	string kmer_len="6,10";
	
    bool revcomp = false;

    bool score_with_median_dif = false;
        
	if (argc < 3) help();
	
    // parse arguments
    string str;
    for (int i = 1; i < argc; i++) { 
        if (i != argc) { 
            str=argv[i];
            if (str == "-q") {  
                queryFile = argv[i + 1];
                i=i+1;
            } else if (str == "-s") { 
                scoreFile = argv[i + 1];
                i=i+1;
            } else if (str == "-r") { 
                revcomp = true;
            } else if (str == "-d") { 
                score_with_median_dif = true;
            } else if (str == "-t") { 
                targetFile = argv[i + 1];
                i=i+1;
            } else if (str == "-i") { 
                inputFile = argv[i + 1];
                i=i+1;
            } else if (str == "-h") { 
                highlightFile = argv[i + 1];
                i=i+1;
            } else if (str == "-b") { 
                backgroundFile = argv[i + 1];
                i=i+1;
            } else if (str == "-o") { 
                outputFile = argv[i + 1];
                i=i+1;
            } else if (str == "-k") { 
                kmer_len = argv[i + 1];
                i=i+1;
            } else if (str == "-c") { 
                kmers = argv[i + 1];
                i=i+1;
            } else if (str == "-m") { 
                max_target_num = atoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-h" || str == "--help") { 
                help();
            } else {
                message("unknown option: "+str);   
                help();        
            }
        }
    }   
	
	vector<string> targetSeqs;
	vector<double> targetScores;
	vector<string> targetNames;

    // if only target sequence is given, score the sequences
    if (queryFile == "" && scoreFile == "" && inputFile == "")
    {
        // score sequences
        message("Loading kmer scores");

        // score column: 5=median difference; 6=-log10(p)
        int c = 6;
        if (score_with_median_dif) c = 5;

        map<string,double> kmer_scores = load_weight(outputFile, 2, c,true);

        message("loading target sequences...");
		message("	sequence file: "+ targetFile);
		map<string,string> targetSeqs = ReadFasta(targetFile);
        message("	"+to_string(targetSeqs.size())+" target sequences loaded");
        
        message("Scoring target sequences");

        string targetFileName = targetFile.substr(targetFile.find_last_of("/\\") + 1);
        string score_profile_file = outputFile+"-"+targetFileName+".profile.txt";
        ofstream fout(score_profile_file.c_str()); 

        int n = 0;
        for (auto i = targetSeqs.begin(); i != targetSeqs.end(); i++)
            {
                cout <<"\r["<<current_time()<<"]	" << ++n << " sequences processed             "  ;
                vector<double> scores = score_a_sequence(i->second,kmer_scores);
                //double sum = boost::accumulate(scores, 0);
                fout << i->first << "\t" << i->second << "\t" ;
                for(int j=0;j<scores.size();j++)
                    {
                        fout << "\t" << scores[j];
                    }
                fout << endl;
            }
        fout.close();
        
        cout << endl;
    
    	message("Done");
        return 0;
    }
	
	string queryName, querySeq;
    map<string, string> seq2name; // map kmer sequence to name
    
	// determine k-mer length (k)
	vector<string> kmer_lens = string_split(kmer_len,",");
	k_min = stoi(kmer_lens[0]);
	if(kmer_lens.size()>1) k_max = stoi(kmer_lens[1]);
	else k_max = k_min;
	
    // generate or load all kmers
    // if an input file is specified by -q
    //    if only a singe sequence is in the file, generate all kmers from the file
    //    else: assume each sequence is a short kmer motif
    // else
    //    generate all possible kmers given k
	vector<string> allkmers;
	vector<int> kmer_pos;
	if(queryFile.size()>0)
	{
		// load query sequence(s)
	
		message("loading query sequence(s)...");
		message("	query file: "+ queryFile);
        
        map<string,string> querySeqs = ReadFasta(queryFile);
        
        message("	" + to_string(querySeqs.size()) + " query sequences loaded" );
        
        // map from seq to name    

        for (map<string, string>::iterator i = querySeqs.begin(); i != querySeqs.end(); ++i)
        {
            queryName = i->first;
            querySeq = i->second;
            querySeq = to_upper(querySeq);
            replace(querySeq.begin(),querySeq.end(),'U','T'); // rna to dna
            if (revcomp) querySeq = reverseComplement(querySeq);
            seq2name[querySeq] = queryName + "|" + querySeq;
        }
        
        if (querySeqs.size() == 1) // a single query sequence
        {
            ifstream fin(queryFile.c_str());
            ReadOneSeqFromFasta(fin,queryName,querySeq);
            fin.close();
            querySeq = to_upper(querySeq);

            if (querySeq.size() == 0) 
                {
                    message("query file empty!");
                    exit(1);
                }
            message("	query name: "+ queryName);
            message("	query length: "+ to_string(querySeq.size()));

            // generate kmers
            for(int k=k_max;k >= k_min;k--)
            {
                for(int i = 0; i <= querySeq.size()-k;i++)
                {
                    if (revcomp) allkmers.push_back(reverseComplement(querySeq.substr(i,k))); 
                    else allkmers.push_back(querySeq.substr(i,k));
                    kmer_pos.push_back(i);
                    seq2name[allkmers[i]] = allkmers[i]; // default: name=seq
                }
            }
        } else { // loaded mulitple sequences. assum each one is a query kmer
            int i = 0;
            map<string, string>::iterator it;
            for (it = seq2name.begin(); it != seq2name.end(); ++it)
            {
                allkmers.push_back(it->first);
                kmer_pos.push_back(i);
                i = i + 1;
            }
        }
	} else // enumerate all possible kmers
	{
		message("No query sequence specified. Use all possible kmers");
		allkmers = generate_kmers(k_min, "ACGT");
		for(int k=k_min+1;k <= k_max;k++)
		{
			allkmers += generate_kmers(k, "ACGT");	
		}
		for(int i=0;i<allkmers.size();i++) 
        {
            kmer_pos.push_back(i);
            seq2name[allkmers[i]] = allkmers[i]; // default: name=seq
        }
	}
	
	if(inputFile.size() == 0) // 
	{
		message("find targets with both sequences and scores, and sort by scores...");
		srand(time(NULL));
	    string tmp = outputFile+random_string(20); // a random string for temp files
		// get all IDs in target fasta files
		system_run("cat "+targetFile + " | grep '>' | sed 's/>//g' > target."+tmp);
		// select lines in score file whose id is present in target fasta file, sort based on score
		intersectTab(scoreFile, "target."+tmp, "selected."+tmp);
		system_run("cat selected."+tmp+" | sort -grk 2 > sorted."+tmp);
		// if too many sequences, only take top and bottom
		int total_targets = count_lines("sorted."+tmp);
		message("    total target sequences sorted: "+to_string(total_targets));
	
		if(total_targets > max_target_num && kmers.size() == 0)
		{	
			int to_be_kept = int(max_target_num/2);
			message("Only keep the top "+to_string(to_be_kept)+" and the bottom "+to_string(to_be_kept)+" sequences");
			system_run("head -n "+to_string(to_be_kept)+" sorted."+tmp+" > "+scoreFile+".processed");
			system_run("tail -n "+to_string(to_be_kept)+" sorted."+tmp+" >> "+scoreFile+".processed");
		} else system_run("mv sorted."+tmp +" "+scoreFile+".processed");
	
		system_run("rm *"+tmp);

		message("loading target scores...");
		message("	score file: "+ scoreFile+".processed");
	
		int total = load_scores(scoreFile+".processed", targetNames, targetScores, 0,1);
        message("	"+to_string(total)+" scores loaded");
        
		message("loading target sequences...");
		message("	sequence file: "+ targetFile);
		map<string,string> allseqs = ReadFasta(targetFile);
        message("	"+to_string(allseqs.size())+" target sequences loaded");
		
		// match target sequences with their scores, ignore sequences without scores
        int not_found = 0;
		for(int i =0;i<targetNames.size();i++)
		{
            if (allseqs.find(targetNames[i]) == allseqs.end() ) // not found
            {
                not_found ++;
            } else {
    			targetSeqs.push_back(to_upper(allseqs[targetNames[i]]));
            }
		}
        message("	"+to_string(targetSeqs.size())+" target sequences with scores");
        if (not_found > 0)
            message("	"+to_string(not_found)+" target sequences have no scores");
		//debug
		//cout << targetSeqs[0] << endl;
		//exit(1);
	}
	else
	{
		message("loading combined sequence and score file...");
		int total = load_scores(inputFile, targetSeqs, targetScores, 0,1);
	}
	

	
	// load highlight sequences if provided
	vector<string> hltNames, hltSeqs;
	if( highlightFile.size() > 0 )
	{
		message("loading highlight sequences...");
		message("	highlight file: "+ highlightFile);
	
		ReadFastaToVectors( highlightFile, hltNames, hltSeqs );
		to_upper(hltSeqs);
		message("	number of highlight sequences: "+ to_string(hltSeqs.size()));
	}

	//background if provided
	vector<string> bgNames, bgSeqs;
	vector<double> bgScores (bgSeqs.size(),0);
	if( backgroundFile.size() > 0 )
	{
		message("loading background sequences...");
		message("	background file: "+ backgroundFile);
	
		ReadFastaToVectors( backgroundFile, bgNames, bgSeqs );
		to_upper(bgSeqs);
		message("	number of background sequences: "+ to_string(bgSeqs.size()));
	
		message("merge with background sequences...");
		targetSeqs.insert(targetSeqs.end(), bgSeqs.begin(), bgSeqs.end());
		targetScores.insert(targetScores.end(), bgScores.begin(), bgScores.end());
	}
	
	if(kmers.size()>0)
	{
		vector<string> kmers_all = string_split(kmers,",");
		for(int i=0;i<kmers_all.size();i++)
		{
			message("generating CDF plot for kmer "+kmers_all[i]);
			array<int,3> totals = kmer_cdf2(kmers_all[i], targetSeqs, targetScores, outputFile+"-"+kmers_all[i]);
            message("    "+to_string(totals[0])+" sequences with no match");
            message("    "+to_string(totals[1])+" sequences with a single match");
            message("    "+to_string(totals[2])+" sequences with multiple matchs");
		}
		return 0;
	}
	
	message("ranking sequences by scores...");
	
	vector<double> score_ranks = get_ranks(targetScores,true);
	vector<double> score_ranks_copy = score_ranks;
	
	/*
	message("generating 10 sets of shuffled scores...");
	
	vector< vector<double> > score_ranks_shuffled;
	for (int i=0;i<10;i++) 
	{
		random_shuffle(score_ranks_copy.begin(),score_ranks_copy.end());
		score_ranks_shuffled.push_back(score_ranks_copy);
	}
	*/
	
	/* 
	// normalize scores
	double max_score = max(targetScores);
	double min_score = min(targetScores);
	cout << max_score <<","<<min_score << endl;
	double range = max_score - min_score;
	for(int i=0;i<targetScores.size();i++) 
	{
		cout << targetScores[i] << "\t";
		targetScores[i] = (targetScores[i] - min_score) / range;
		cout << targetScores[i] << endl;
	}
	
    // nucleotide to position
    map<char,int> letter2pos;
    letter2pos['A'] = 0;
    letter2pos['C'] = 1;
    letter2pos['G'] = 2;
    letter2pos['T'] = 3;
	boost::numeric::ublas::matrix<double> pwm = initialize_pwm_from_one_seq(querySeq);
	print_matrix(pwm);
	pwm = update_pwm_from_seqs(targetSeqs, targetScores, pwm,  letter2pos);
	print_matrix(pwm);
	*/
		
	message("scanning query sequences and calcualting p values...");
	message("	"+to_string(allkmers.size())+" kmer to be processed...");
	
	//map<string,string> calculated_kmers; // kmer to cors

	ofstream fout(outputFile.c_str());
	
	string header = "position\tkmer_seq\tn_total_seq\tn_with_match\tmedian_diff\tlog10_p_value";
	if (hltSeqs.size()>0) header = header + "\thighlight";

	fout << header << endl;
	
	for(int i = 0; i < allkmers.size();i++)
	{
	    cout <<"\r["<<current_time()<<"] 	" << i+1 << " kmers processed"  ;
        
		string kmer_cors = to_string(kmer_pos[i])+"\t"+seq2name[allkmers[i]];
		
		//array<double,4> utest = kmer_rank_test(allkmers[i], targetSeqs, score_ranks);
        //kmer_cors = kmer_cors+"\t"+to_string(int(utest[0]))+"\t"+to_string(int(utest[1]))+"\t"+to_string(utest[2])+"\t"+to_string(utest[3]) + "\t";

        array<double,5> utest = kmer_rank_test_with_median(allkmers[i], targetSeqs, score_ranks,targetScores);
        kmer_cors = kmer_cors+"\t"+to_string(int(utest[0]))+"\t"+to_string(int(utest[1]))+"\t"+to_string(utest[4])+"\t"+to_string(utest[3]) + "\t";

		
		/*
		if(calculated_kmers.find(kmer) == calculated_kmers.end()) //  a new kmer
		{
			vector<int> counts = motif_counts_in_seqs(kmer, targetSeqs);
			vector<double> count_ranks = get_ranks(counts,true);
			kmer_cors = kmer_cors + to_string(cor(count_ranks,score_ranks));
			for (int j=0;j<10;j++) kmer_cors = kmer_cors + "\t" + to_string(cor(count_ranks,score_ranks_shuffled[j]));
			*/
			// high lights
			if (hltSeqs.size()>0)
			{
				string highlight = "none";
				for(int j=0;j<hltSeqs.size();j++)
				{
					if (hltSeqs[j].find(allkmers[i]) != std::string::npos) 
					{
						highlight = hltNames[j];
						break;
					}
				}
				kmer_cors = kmer_cors + "\t" + highlight ;
			}
			/*
			
			calculated_kmers[kmer] = kmer_cors;
		} else { // not a new kmer
			kmer_cors = calculated_kmers[kmer];
		}
		*/
		fout << kmer_cors << endl;
	}	
	fout.close();
	
	string script = 
    "pdf('"+outputFile+".pdf',width=10,height=5) \n"
	"par(cex=1.5,mar=c(5, 5, 1, 1))\n"
    "x  = read.table('"+outputFile+"', header=T) \n"
		"col=rep('blue',nrow(x))\n"
		"if(ncol(x)>6){\n"
		"col[x[,7] != 'none'] = 'green'\n"
			"}\n"
	"x = x[order(x[,6],decreasing=T),]\n"
	"plot(x[,1],x[,6],type='h',col=col, bty='n',xlab='5 prime end of match in query (nucleotides)', ylab='-log10(p)') \n"	
	"# top 3\n"
	"text(x[1,1],x[1,6],labels=x[1,2],pos=4,cex=0.5,offset=0.1)\n"
	"text(x[2,1],x[2,6],labels=x[2,2],pos=4,cex=0.5,offset=0.1)\n"
	"text(x[3,1],x[3,6],labels=x[3,2],pos=4,cex=0.5,offset=0.1)\n"
	"n = nrow(x)\n"
	"text(x[n,1],x[n,6],labels=x[n,2],pos=4,cex=0.5,offset=0.1)\n"
	"text(x[n-1,1],x[n-1,6],labels=x[n-1,2],pos=4,cex=0.5,offset=0.1)\n"
	"text(x[n-2,1],x[n-2,6],labels=x[n-2,2],pos=4,cex=0.5,offset=0.1)\n"	
	"#avg = apply(x[,8:17],1,mean,rm.nan=T) \n"
	"#err = apply(x[,8:17],1,sd) \n"
	"#plot(x[,7],type='h',col=col, bty='n',xlab='5 prime end of match in query (nucleotides)', ylab='spearman correlation') \n"
	"#arrows(1:length(avg),avg-err,1:length(avg),avg+err, code=3, angle=90, length=0,col=rgb(1,0,0,0.5),lwd=0.1) \n"
	"#for(i in 8:17){\n"
	"# lines(x[,i],type='h',col=rgb(1,0,0,0.3))\n"
	"#}\n"
    "dev.off() \n"
	"write.table(x[order(x[,6],decreasing=T),],file='"+outputFile+"',sep='\t',append=F,quote=F,row.names=F)\n"
	"\n";

	R_run(script);		
	
    cout << endl;
    
	message("Done");
		
    return 0;
}
