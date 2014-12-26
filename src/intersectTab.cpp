#include "text.h"
#include "utility.h"
#include "iostream"
#include <string>

using namespace std;

void help()
{
    string str =
	"\n"
    "intersectTab: intersection/subtraction/merge between two tabular files based on a shared key column\n"
	"	- Xuebing Wu (wuxb07@gmail.com)\n"
    "\n"
    "Usage: intersectTab -f1 file1 -f2 file2 [options]\n"
    "\n"
    "Options:\n"
    "\n"
    "   -f1 <file>     input file 1 \n"
    "   -f2 <file>     input file 2 \n"
    "   -c1 <integer>  key column in file1, default=1 (first column)\n"
    "   -c2 <integer>  key column in file2, default=1 (first column)\n"
	"   -intersect     (default) output file1 lines with key also present in file2 \n"
	"   -subtract      output file1 lines whose key is not present in file2\n"
	"   -merge         output file1 and file2 side by side for lines with shared key\n"
	"   -fill <string> (-merge only) for keys in file 1 but not file 2, \n"
    "                  fill file 2 parts with a string in each column. \n"
	"                  To skip those lines, use -fill none \n"	
	"   -header        if both files contain a header line \n"
    "\n";

    cerr << str;

    exit(0);
}

int main(int argc, char* argv[]) {

    // default
    string f1;
    string f2;
	string outputfile;
	string job="intersect";
	unsigned c1=0;
	unsigned c2=0;
	bool header=false;
	string fill="none";

	if (argc < 5) help();
	
    // parse arguments
    string str;
    for (int i = 1; i < argc; i++) { 
        if (i != argc) { 
            str=argv[i];
            if (str == "-f1") {  
                f1 = argv[i + 1];
                i=i+1;
            } else if (str == "-f2") { 
                f2 = argv[i + 1];
                i=i+1;
            } else if (str == "-o") { 
                outputfile = argv[i + 1];
                i=i+1;
            } else if (str == "-c1") { 
                c1 = atoi(argv[i + 1]) - 1;
                i=i+1;
            } else if (str == "-c2") { 
                c2 = atoi(argv[i + 1]) - 1;
                i=i+1;
            } else if (str == "-header") { 
                header = true;
            } else if (str == "-subtract") { 
                job = "subtract";
            } else if (str == "-merge") { 
                job = "merge";
            } else if (str == "-fill") { 
                fill = argv[i + 1];
                i=i+1;
            } else if (str == "-h" || str == "--help") { 
                help();
            } else {
                message("unknown option: "+str);   
                help();        
            }
        }
    }   

	if(job == "merge") mergeTab(f1,f2,outputfile, c1,c2, header,fill);   
	else if (job == "subtract") intersectTab(f1,f2, outputfile, c1,c2, true);
	else intersectTab(f1,f2, outputfile, c1,c2, false);

    return 0;
}
