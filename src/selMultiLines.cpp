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
    "   -i  file        input file\n"
	"   -c  column      IDs in this column, 1-based\n"
	"   -id_file file   file containing IDs to be selected\n"
    "   -o  file        output file\n"
	"   -n  num         each record contains n lines (e.g. 2 for fasta, 4 for fastq)\n"
	"   -id_prefix ch   prefix for ids in input file but not in id_file\n"
    "\n";

    cerr << str;

    exit(0);
}

int main(int argc, char* argv[]) {

	if (argc < 2) help();

    // default
    string input="input";
    string output="output";
	string idfile = "";
	string id_prefix = "";
    int col=1;
	int nline = 1;

    // parse arguments
    string str;
    for (int i = 1; i < argc; i++) { 
        if (i != argc) { 
            str=argv[i];
            if (str == "-i") {  
                input = argv[i + 1];
                i=i+1;
            } else if (str == "-id_file") { 
                idfile = argv[i + 1];
                i=i+1;
            } else if (str == "-o") { 
                output = argv[i + 1];
                i=i+1;
            } else if (str == "-c") { 
                col = atoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-n") { 
                nline = atoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-id_prefix") { 
                id_prefix = argv[i + 1];
                i=i+1;
            } else {
                message("unknown option: "+str);   
                help();        
            }
        }
    }   

	int selected = select_multi_lines(input, idfile,  output, nline, col, id_prefix);
	
	message(to_string(selected)+" matched IDs saved to "+output);

    return 0;
}
