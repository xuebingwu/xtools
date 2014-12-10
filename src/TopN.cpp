#include "text.h"
#include "utility.h"
#include "iostream"

using namespace std;

void help()
{
    string str = "\n"
    "TopN: keep the first N lines for each unique value in a given column \n"
	"	- Xuebing Wu (wuxb07@gmail.com)\n"
    "\n"
    "Usage: TopN -i input -o output -c column -n N -s sortOpts \n"
    "       example: TopN -i input -o output -c 3 -n 3 -s \"-k3,3 -k4,4n\" \n"
    "\n"
    "Options:\n"
    "\n"
    "   -i  input   input file\n"
    "   -o  output  output file\n"
    "   -c  column  column number, 1-based\n"
    "   -n  number  number of lines to keep\n"
    "   -s  options sort input with given options\n"
    "\n";

    cerr << str;

    exit(0);
}

int main(int argc, char* argv[]) {

	if (argc < 2) help();

    // default
    string input="input";
    string output="output";
    int col=0;
    int topn=1;
    string sortOpts="";

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
            } else if (str == "-c") { 
                col = atoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-n") {    
                topn = atoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-s") {
                sortOpts = argv[i + 1];
                i=i+1;
            } else {
                message("unknown option: "+str);   
                help();        
            }
        }
    }   

    remove_duplicates(input,output,col,topn,sortOpts);

    return 0;
}
