#include "utility.h"
#include <iostream>
using namespace std;

void help()
{
    string str =
	"\n"
    "xtools: Xuebing's tools for computational biology\n"
	"	- Xuebing Wu (wuxb07@gmail.com)\n"
    "\n"
    "Usage: xtools <command> [options]\n"
    "\n"
    "Commands:\n"
    "\n"
	"gene analysis ==========\n"
    "   geneFeature   Examine whether a list of genes have unique feaures\n"
    "   geneSet       Examine whether a list of genes have over-represented annotations\n"
    "   \n"
	"sequence analysis ======\n"
	"   PKA           Positional kmer enrichment analysis from aligned sequences\n"
	"   PKA2          PKA with ranked or weighted sequences\n"
    "\n";

    cerr << str;

    exit(0);
}

int main(int argc, char* argv[]) {

	if (argc < 2) help();
	
	string valid_cmds = "geneFeature, geneSet, PKA, PKA2, ";
	
	string cmd(argv[1]);
	
	if(valid_cmds.find(cmd+",") == std::string::npos) 
	{
		cerr << "invalid command: " << argv[1] << endl;
		help();
	}
		
	for(int i=2; i < argc; i++) cmd += " "+string(argv[i]);
	
	system_run(cmd);

    return 0;
}
