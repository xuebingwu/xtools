#include "utility.h"
#include "iostream"

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
	
	//plot_hairpin_from_file(input);
		
	
	
	/*
	
	AGAGAGAAGAUAU    ---    -GUUGC    A   CG     C   A      GUAU 
	             UGAG   GCCU      CACA ACC  UAGAU CGA CUUGUG    U
	             ||||   ||||      |||| |||  ||||| ||| ||||||    A
	             ACUC   CGGA      GUGU UGG  AUCUA GUU GAACAC    C
	--GGGGUCCAGGC    UAA    UUGUCU    A   AU     U   C      GCCU 
					 */
	//mirna_loop_definition(id,seq,fold,35);
	
	//cout << endl;
	
	//mirna_basal_stem_definition(id,seq,fold,35);
	//remove_short_stem_from_file(input, input+".no.short.stem", 20, 25,10);
	//mirna_feature_from_file(input+".no.short.stem",output,35);
	
    return 0;
}
