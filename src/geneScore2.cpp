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
    "geneScore2: find gene sets where two scores strongly correlates\n"
	"	- Xuebing Wu (wuxb07@gmail.com)\n"
    "\n"
    "Usage: geneScore2 -i input -o output\n"
    "\n"
    "Options:\n"
    "\n"
    "   -i  input       3-column tabular file:gene,score1,score2\n"
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
	int cScore1=1; // column 2 is score1
	int cScore2=2; // column 3 is score2
	
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
            } else if (str == "-d") { 
                databasefile = argv[i + 1];
                i=i+1;
            } else if (str == "-cGene") { 
                cGene = atoi(argv[i + 1]) - 1;
                i=i+1;
            } else if (str == "-cScore1") { 
                cScore1 = atoi(argv[i + 1]) - 1;
                i=i+1;
            } else if (str == "-cScore2") { 
                cScore2 = atoi(argv[i + 1]) - 1;
                i=i+1;
            } else if (str == "-p") { 
                p_cutoff = atof(argv[i + 1]);
                i=i+1;
            } else if (str == "-skip") { 
                skip = atoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-h" || str == "--help") { 
                help();
            } else {
                message("unknown option: "+str);   
                help();        
            }
        }
    }   
	
	// load gene scores
	map<string,double> scores1;
	map<string,double> scores2;
	load_weighted_gene_list(genefile, scores1, true, skip, cGene, cScore1);	
	load_weighted_gene_list(genefile, scores2, true, skip, cGene, cScore2);	
	
	message(to_string(scores1.size()) + " genes loaded from file: "+ genefile);
	message("Testing gene sets in file: "+ databasefile);
	
	string suffix = databasefile.substr(databasefile.size()-4,4);

	if (suffix == ".gft") // gene feature
	{
		// find a cut-off to see if top or bottom set of genes give better correlation
		
		//map<string,double> scores;
		//load_weighted_gene_list(genefile, scores, true, skip, cGene, cScore);		
		//gene_feature_correlation(scores, databasefile, outputfile);		
	}
	else // gene set
	{
		int nSig = find_sig_gene_sets_correlate_two_scores(scores1, scores2, databasefile, outputfile, p_cutoff);

		//int nSig = find_sig_gene_sets_weighted(genes, scores, databasefile, outputfile, p_cutoff);
	
		//message(to_string(nSig) + " significant gene features identified");
	
		// sort output file by p-value and add header
		//string header = "gene.set\tdescription\tp.value\tt.statistics\tcount.f\tmean.f\tstdev.f\tcount.b\tmean.b\tstdev.b";
		//sort_file_and_add_header(outputfile,header,"-k3,3g");
	}
	return 0;
}
