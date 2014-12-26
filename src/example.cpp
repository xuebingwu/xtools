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
            } else if (str == "-n") { 
                n = atoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-h" || str == "--help") { 
                help();
            } else {
                message("unknown option: "+str);   
                help();        
            }
        }
    }   

	 tab_seq_to_letter_matrix( input,  output,  1,2, 1,  0);
   
		

    return 0;
}
