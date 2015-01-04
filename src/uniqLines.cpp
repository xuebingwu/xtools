#include "text.h"
#include "utility.h"
#include "iostream"

using namespace std;

void help()
{
    string str = "\n"
    "uniqLines: keep lines with unique value in a given column \n"
	"	- Xuebing Wu (wuxb07@gmail.com)\n"
    "\n"
    "Usage: uniqLines -i input -o output -c column\n"
    "\n"
    "Options:\n"
    "\n"
    "   -i  input   input file\n"
    "   -o  output  output file\n"
    "   -c  column  column number, 1-based\n"
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
            } else {
                message("unknown option: "+str);   
                help();        
            }
        }
    }   

	int n = find_unique_lines( input, output, col);
	
	message(to_string(n)+" unique lines saved to "+output);

    return 0;
}
