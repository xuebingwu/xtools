#include "utility.h"
#include "iostream"
#include "interval.h"

using namespace std;

void help()
{
    string str =
	"\n"
    "flipStrand: flip interval strand\n"
	"	- Xuebing Wu (wuxb07@gmail.com)\n"
    "\n"
    "Usage: flipStrand -i input -o output\n"
    "\n"
    "Options:\n"
    "\n"
    "   -i  input   input BED file\n"
    "   -o  output  output BED file\n"
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
	
    flipStrand(input,output);
	
    return 0;
}
