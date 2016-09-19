#include "utility.h"
#include "iostream"
#include "interval.h"

using namespace std;

void help()
{
    string str =
	"\n"
    "iFeat: assign intervals to gene features\n"
	"	- Xuebing Wu (wuxb07@gmail.com)\n"
    "\n"
    "Usage: iFeat -i input [-o output -n UTR3,ext_UTR3,CDS,intron,UTR5,uaRNA -f f1,f2,f3,f4,f5,f6 -d . ] \n"
    "\n"
    "Options:\n"
    "\n"
    "   -i  input   (required) input BED file\n"
    "   -o  output  output file name (default: output)\n"
	"   -n  list    comma-separated feature names, orderred\n"
	"   -f  list    comma-separated feature file names, same order as names\n"
	"   -d  dir     directory name for feature files\n"
	"\n"
	"Notes:\n"
	"  default feature names: UTR3,ext_UTR3,CDS,intron,UTR5,uaRNA\n"
	"  default feature directory: /lab/bartel1_ata/wuxbl/solexa_bartel/annotations/hg19 \n"
	"  default feature files: hg19.refGene.UTR3.bed,hg19.refGene.dn5kb.bed.clean,hg19.refGene.CDS.bed,hg19.refGene.intron.bed,hg19.refGene.UTR5.bed,hg19.refGene.ua5kb.bed\n"
    "\n";

    cerr << str;

    exit(0);
}


int main(int argc, char* argv[]) {

    // default
    string input;
	string output="output";
    string names="UTR3,ext_UTR3,CDS,intron,UTR5,uaRNA";
	string dir="/lab/bartel1_ata/wuxbl/solexa_bartel/annotations/hg19";
	string files="hg19.refGene.UTR3.bed,hg19.refGene.dn5kb.bed.clean,hg19.refGene.CDS.bed,hg19.refGene.intron.bed,hg19.refGene.UTR5.bed,hg19.refGene.ua5kb.bed";
	
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
                names = argv[i + 1];
                i=i+1;
            } else if (str == "-d") { 
                dir = argv[i + 1];
                i=i+1;
            } else if (str == "-f") { 
                files = argv[i + 1];
                i=i+1;
            } else if (str == "-h" || str == "--help") { 
                help();
            } else {
                message("unknown option: "+str);   
                help();        
            }
        }
    }   
	
	if(input.size()==0){
		message("Please at least specify an input file!");
		exit(1);
	}

	vector<string> all_names  = string_split(names,",");
	vector<string> all_files = string_split(files,",");
	for (int i = 0; i<all_files.size();i++){
		all_files[i] = dir + '/' + all_files[i];
	}
	interval_annotation( input,output, all_files, all_names);
	
    return 0;
}
