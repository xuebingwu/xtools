#include "utility.h"
#include "iostream"
#include "sequence.h"
#include "text.h"
#include "container.h"
#include "stat.h"



using namespace std;

void help()
{
    string str =
	"\n"
    "MatchScan\n"
	"	a program to find significant matches between a query sequence and a set of target sequences,\n"
	"   where the frequency of the matched sequences in targets correlates with target scores\n"
	"\n"
	"	- Xuebing Wu (wuxb07@gmail.com)\n"
    "\n"
    "Usage: MatchScan -q query.fasta -t target.fasta -s score.tab -o output [ options ]\n"
    "\n"
    "Options:\n"
    "\n"
    "   -q  file   fasta file containing a single query sequence\n"
    "   -t  file   fasta file containing all possible target sequences\n"
    "   -s  file   tab-delimited score file (target name in column 1, score in column 2)\n"
	"              target name should match those in target fasta file, can be a subset\n"
	"   -o  file   prefix for output files (tabular data and PDF figure)\n"
	"   -h  file   query sequences to highlight in green\n"
	"   -k  num(s) length of subsequence to scan at each step (i.e. k-mer)\n"
	"              e.g. to search for hexamer: -k 6; for k from 6 to 10 (default): -k 6,10\n"		
	"   -m  num    max number of target sequences to use (top and bottom)\n"		
	"              default 10000, i.e. top 5000 and bottom 5000 based on scores\n"		
		
    "\n"
	"Example:\n"
	" MatchScan -q TERC.fa -t TERC-ChIRP.fa -s TERC-ChIRP.score.txt -o TERC\n"
	" MatchScan -q TERC.fa -t TERC-ChIRP.fa -s TERC-ChIRP.score.txt -o TERC -k 10 \n"
	" MatchScan -q TERC.fa -t TERC-ChIRP.fa -s TERC-ChIRP.score.txt -o TERC -h terc_probe.fa.rc \n";

    cerr << str;

    exit(0);
}




int main(int argc, char* argv[]) {

    // default
    string queryFile="query";
    string targetFile="target";
    string scoreFile="score";
	string outputFile="output";
    string highlightFile="";
	string backgroundFile="";
	
	int max_target_num = 10000;

	int k_min,k_max;
	string kmer_len="6,10";
	
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
            } else if (str == "-t") { 
                targetFile = argv[i + 1];
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
	
	// determine k
	vector<string> kmer_lens = string_split(kmer_len,",");
	k_min = stoi(kmer_lens[0]);
	if(kmer_lens.size()>1) k_max = stoi(kmer_lens[1]);
	else k_max = k_min;
	
	message("kmer length from "+to_string(k_min)+" to "+to_string(k_max));
	
	
	// load query sequence, one seq only
	
	message("loading query sequence...");
	message("	query file: "+ queryFile);
	
	ifstream fin(queryFile.c_str());
	string queryName, querySeq;
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
	
	
	// intersecting score file and target file, to find targets with both sequence and score
	srand(time(NULL));
    string tmp = outputFile+random_string(20); // a random string for temp files
	// get all IDs in target fasta files
	system_run("cat "+targetFile + " | grep '>' | sed 's/>//g' > target."+tmp);
	// select lines in score file whose id is present in target fasta file, sort based on score
	intersectTab(scoreFile, "target."+tmp, "selected."+tmp);
	system_run("cat selected."+tmp+" | sort -nrk 2 > sorted."+tmp);
	// if too many sequences, only take top and bottom
	if(count_lines("sorted."+tmp) > max_target_num)
	{
		system_run("head -n "+to_string(int(max_target_num/2))+" sorted."+tmp+" > "+scoreFile+".processed");
		system_run("tail -n "+to_string(int(max_target_num/2))+" sorted."+tmp+" >> "+scoreFile+".processed");
	} else system_run("mv sorted."+tmp +" "+scoreFile+".processed");
	
	system_run("rm *"+tmp);
	
	// load score
	message("loading scores...");
	message("	score file: "+ scoreFile+".processed");
	
	vector<string> targetNames;
	vector<double> targetScores;
	int total = load_scores(scoreFile+".processed", targetNames, targetScores, 0,1);
	message("	total targets loaded: "+to_string(total));


	message("loading target sequences...");
	message("	target file: "+ targetFile);
	
	map<string,string> allseqs = ReadFasta(targetFile);
		
	// sort target sequences by their score, ignore sequences without scores
	vector<string> targetSeqs;
	for(int i =0;i<targetNames.size();i++)
	{
		targetSeqs.push_back(to_upper(allseqs[targetNames[i]]));
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
	
	message("scanning query sequences and calcualting correlation scores for each kmer...");
	
	//map<string,string> calculated_kmers; // kmer to cors
	
	for(int k=k_max;k >= k_min;k--)
	{
		string outputfilename = outputFile+"."+to_string(k);
		ofstream fout(outputfilename.c_str());
		
		message("    processing k="+to_string(k));
		for(int i = 0; i <= querySeq.size()-k;i++)
			{
				string kmer = querySeq.substr(i,k);
				string kmer_cors = to_string(i)+"\t"+kmer;
				
				// rank sum test
				array<double,4> utest = kmer_rank_test(kmer, targetSeqs, score_ranks);
				kmer_cors = kmer_cors+"\t"+to_string(int(utest[0]))+"\t"+to_string(int(utest[1]))+"\t"+to_string(utest[2])+"\t"+to_string(utest[3]) + "\t";
				
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
							if (hltSeqs[j].find(kmer) != std::string::npos) 
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
	    "pdf('"+outputfilename+".pdf',width=10,height=5) \n"
		"par(cex=1.5,mar=c(5, 5, 1, 1))\n"
	    "x  = read.table('"+outputfilename+"', header=F) \n"
			"col=rep('blue',nrow(x))\n"
			"if(ncol(x)>6){\n"
			"col[x[,7] != 'none'] = 'green'\n"
				"}\n"
		"plot(x[,6],type='h',col=col, bty='n',xlab='5 prime end of match in query (nucleotides)', ylab='-log10(p)') \n"	
		"#avg = apply(x[,8:17],1,mean,rm.nan=T) \n"
		"#err = apply(x[,8:17],1,sd) \n"
		"#plot(x[,7],type='h',col=col, bty='n',xlab='5 prime end of match in query (nucleotides)', ylab='spearman correlation') \n"
		"#arrows(1:length(avg),avg-err,1:length(avg),avg+err, code=3, angle=90, length=0,col=rgb(1,0,0,0.5),lwd=0.1) \n"
		"#for(i in 8:17){\n"
		"# lines(x[,i],type='h',col=rgb(1,0,0,0.3))\n"
		"#}\n"
	    "dev.off() \n";

		R_run(script);
		
		system_run("sort -nrk 6 "+outputfilename + " > "+tmp);
		system_run("mv "+tmp+" "+outputfilename);
		
	}
	

	
	message("Done");
		
    return 0;
}
