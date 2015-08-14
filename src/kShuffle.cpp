#include "utility.h"
#include "iostream"
#include "sequence.h"

using namespace std;

void help()
{
    string str =
	"\n"
    "kShuffle: shuffle sequences preserving k-mer frequency\n"
	"	- adapted from http://digital.cs.usu.edu/~mjiang/ushuffle/\n"
    "\n"
    "Usage: kShuffle -i input.fa -o output.fa\n"
    "\n"
    "Options:\n"
    "\n"
    "   -i  <fasta>    input file, fasta format \n"
    "   -o  <fasta>    output file, fasta format \n"
	"   -k  <integer>  preserve k-mer frequency, default=2, \n"
	"                  i.e. dinucleotide frequency is preserved \n"
    "   -n  <integer>  number of shuffled copies to generate from each \n"
	"                  input sequence, default=1. Output sequences will \n"
	"                  have names like >name-1, >name-2, ..., etc \n"
    "\n";

    cerr << str;
}

int main(int argc, char* argv[]) {

    // default
    string input;
    string output;
    int n=1;
	int k=2;

	if (argc < 2) help();
	
    // parse arguments
    string str;
    for (int i = 1; i < argc; i++) { 
        if (i != argc) { 
            str=argv[i];
            if (str == "-i") {  
                input = argv[i + 1];
                i=i+1;
            } else if (str == "-o") { 
                output = argv[i + 1];
                i=i+1;
            } else if (str == "-n") { 
                n = atoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-k") { 
                k = atoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-h" || str == "--help") { 
                help();
            } else {
                message("unknown option: "+str);   
                help();        
            }
        }
    }   

	if(output.size()==0 || input.size() == 0)
	{
		help();		
		cerr << "** ERROR: Please specify both input and output files !" << endl;
		exit(1);
	}

	/**/
	map<string,string> seqs1 = ReadFasta(input); 
    ofstream fout;
    fout.open(output.c_str());

	srand(time(NULL));

    for (map<string,string>::iterator it=seqs1.begin(); it!=seqs1.end(); ++it)
    {
		for(int i=0;i<n;i++)
		{		
			string str = shuffle_seq_preserving_k_let(it->second,k);
			fout << ">" << it->first << "-" << i+1 << endl << str << endl;
		}
    }
	fout.close();
    return 0;
}
