#include "utility.h"
#include "iostream"
#include "sequence.h"

using namespace std;

void help()
{
    string str =
	"\n"
    "proteinAlignment: score protein alignment\n"
    "\n"
    "Usage: proteinAlignment -q query.fa -t target.fa -o output.txt -s scoring_matrix -t threshold\n"
    "\n"
    "Options:\n"
    "\n"
    "   -q  <fasta>    query file, fasta format \n"
    "   -t  <fasta>    target file, fasta format \n"
	"   -m  <matrix>   scoring matrix file, alphabet in line 1\n"
	"   -o  <txt>      output file name\n"
	"   -s  <integer>  minimum alignment score for output \n"
    "\n";

    cerr << str;
	
	exit(1);
}

int main(int argc, char* argv[]) {

    // default
    string queryFile;
	string targetFile;
	string scoreMatFile;
    string outputFile;
    int t = -10000;

	if (argc < 2) help();
	
    // parse arguments
    string str;
    for (int i = 1; i < argc; i++) { 
        if (i != argc) { 
            str=argv[i];
            if (str == "-q") {  
                queryFile = argv[i + 1];
                i=i+1;
            } else if (str == "-t") { 
                targetFile = argv[i + 1];
                i=i+1;
            } else if (str == "-m") { 
                scoreMatFile = argv[i + 1];
                i=i+1;
            } else if (str == "-o") { 
                outputFile = argv[i + 1];
                i=i+1;
            } else if (str == "-s") { 
                t = atoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-h" || str == "--help") { 
                help();
            } else {
                message("unknown option: "+str);   
                help();        
            }
        }
    }   
	
	if(queryFile.size()==0 || targetFile.size() == 0 || scoreMatFile.size() == 0 || outputFile.size() == 0 || t == -10000)
	{
		help();		
		exit(1);
	}
	
	// load scoring matrix
	map<char,map<char,int>> M = LoadScoreMat(scoreMatFile);

	// load query sequences
	map<string,string> querySeqs = ReadFasta(queryFile);
	
	ofstream fout(outputFile.c_str());
	
	// for each target sequences
    ifstream fin(targetFile.c_str());
  
    string name,seq;
	int score;
	int k=0;
    while(fin.good())
    {
		cout << k++ << endl;
      ReadOneSeqFromFasta(fin,name,seq);
      for (map<string,string>::iterator it=querySeqs.begin(); it!=querySeqs.end(); ++it)
      {
		  score = seqAlignmentScore(M, it->second, seq);
		  if (score >= t) fout << it->first << "\t" << name << "\t" << score << endl;
      }
    }
    fin.close();
	fout.close();
    return 0;
}
