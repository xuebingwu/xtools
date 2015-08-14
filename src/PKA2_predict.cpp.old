#include "utility.h"
#include "iostream"
#include "sequence.h"

using namespace std;

void help()
{
    string str =
	"\n"
    "PKA2_predict: score sequences using a PKA2 model\n"
	"	- Xuebing Wu (wuxb07@gmail.com)\n"
    "\n"
    "Usage: PKA2_predict -i input -o output -m model\n"
    "\n"
    "Options:\n"
    "\n"
    "   -i  input   input sequence file, fasta format or tabular\n"
    "   -o  output  output file, tabular format\n"
    "   -m  model   model file, PKA2 output\n"
	"   -t  col     input is tabular and sequence in column col\n"
    "\n";

    cerr << str;

    exit(0);
}

int main(int argc, char* argv[]) {

    // default
    string input;
    string output;
    string model;
	int col = -1;

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
            } else if (str == "-m") { 
                model = argv[i + 1];
                i=i+1;
            } else if (str == "-t") { 
                col = stoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-h" || str == "--help") { 
                help();
            } else {
                message("unknown option: "+str);   
                help();        
            }
        }
    }   
	
	//cout << model << endl;
	
	vector<positional_kmer> ranked_kmers = load_model_from_file(model);
	message(to_string(ranked_kmers.size()) + " ranked positional kmers loaded from the model file : " + model);	
	
	if(col<0) score_fasta_using_PKA_model(input, output, ranked_kmers);
	else score_tabular_using_PKA_model(input, col,output, ranked_kmers);
	
    return 0;
}
