#include <fstream>
#include <iostream>

#include "utility.h"
#include "stat.h"
#include "gene.h"

using namespace std;

void help()
{
    string str =
	"\n"
    "geneScore: correlate your gene scores to gene features and / or gene sets\n"
	"	- Xuebing Wu (wuxb07@gmail.com)\n"
    "\n"
    "Usage: geneScore -i input -o output\n"
    "\n"
    "Options:\n"
    "\n"
    "   -i  input       input file\n"
    "   -o  output      output file\n"
	"   -d  database    gene feature (*.gft) or gene set (otherwise)\n"	
	"   -p  pcutoff     p-value cutoff, default 0.05\n"
	"   -cGene integer  which column is gene name, default 1\n"
	"   -cScore integer whcih column is gene score, default 2\n"
	"   -skip integer   skip the first N lines (i.e. header), default=0\n"
	"   -h/--help       print help message\n"
	"\n"
	"Notes: check geneFeature or geneSet for database file format\n"
    "\n";

    cerr << str;

    exit(0);
}

int main(int argc, char* argv[]) {
	
	string genefile;
	string databasefile;
	string outputfile;
	string backgroundfile;
	double p_cutoff = 2.0; // default is to output all gene features
	int skip=0; // no header
	int cGene=0; // column 1 is gene
	int cScore=1; // column 2 is score
	
	if (argc < 2) help();
	
    // parse arguments
    string str;
    for (int i = 1; i < argc; i++) { 
        if (i != argc) { 
            str=argv[i];
            if (str == "-i") {  
                genefile = argv[i + 1];
                i=i+1;
            } else if (str == "-o") { 
                outputfile = argv[i + 1];
                i=i+1;
            } else if (str == "-b") { 
                backgroundfile = argv[i + 1];
                i=i+1;
            } else if (str == "-d") { 
                databasefile = argv[i + 1];
                i=i+1;
            } else if (str == "-p") { 
                p_cutoff = atof(argv[i + 1]);
                i=i+1;
            } else if (str == "-skip") { 
                skip = atoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-cGene") { 
                cGene = atoi(argv[i + 1]) - 1;
                i=i+1;
            } else if (str == "-cScore") { 
                cScore = atoi(argv[i + 1]) - 1;
                i=i+1;
            } else if (str == "-h" || str == "--help") { 
                help();
            } else {
                message("unknown option: "+str);   
                help();        
            }
        }
    }   
	

	
	string suffix = databasefile.substr(databasefile.size()-4,4);

	if (suffix == ".gft") // gene feature
	{
		map<string,double> scores;
		load_weighted_gene_list(genefile, scores, true, skip, cGene, cScore);		
		gene_feature_correlation(scores, databasefile, outputfile);		
	}
	else // gene set
	{
		// load weighted genes
		vector<string> genes;
		vector<double> scores;
		load_weighted_gene_list(genefile, genes, scores, true, skip, cGene, cScore);

		message(to_string(scores.size()) + " genes loaded from file: "+ genefile);
		message("Testing gene sets in file: "+ databasefile);
		
		int nSig = find_sig_gene_sets_weighted(genes, scores, databasefile, outputfile, p_cutoff);
	
		message(to_string(nSig) + " significant gene features identified");
	
		// sort output file by p-value and add header
		string header = "gene.set\tdescription\tp.value\tt.statistics\tcount.f\tmean.f\tstdev.f\tcount.b\tmean.b\tstdev.b";
		sort_file_and_add_header(outputfile,header,"-k3,3g");
	}
	return 0;
}
