#include "utility.h"
#include "iostream"

#include "sequence.h"

using namespace std;

void help()
{
    string str =
	"\n"
    "Example: an example program\n"
	"	- Xuebing Wu (wuxb07@gmail.com)\n"
    "\n"
    "Usage: example -i input -o output -n number\n"
    "\n"
    "Options:\n"
    "\n"
    "   -i  input   input file\n"
    "   -o  output  output file\n"
    "   -n  number  number\n"
    "\n";

    cerr << str;

    exit(0);
}



int main(int argc, char* argv[]) {

    // default
    string input="input";
    string output="output";
    int n=0;

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
            } else if (str == "-h" || str == "--help") { 
                help();
            } else {
                message("unknown option: "+str);   
                help();        
            }
        }
    }   

	vector<string> seqs = load_ranked_sequences_to_vectors(input,0);
	pairwise_sequence_similarity_matrix( seqs,  output);
   
		

    return 0;
}
