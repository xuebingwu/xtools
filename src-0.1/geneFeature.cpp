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
    "geneFeature: A program to discover what's unique about your list of genes\n"
	"	- Xuebing Wu (wuxb07@gmail.com)\n"
    "\n"
    "Usage: gene_feature -i input -o output\n"
    "\n"
    "Options:\n"
    "\n"
    "	-i  input       input file\n"
	"	-b  background  optional background gene list\n"
    "	-o  output      output file\n"
	"	-d  database    gene feature database file (*.gft)\n"	
	"	-p  pcutoff     p-value cutoff, default 0.05\n"
	"	-h/--help       print help message\n"
	"\n"
	"Notes: gene feature database file format\n"
	"	- tab-delimited\n"
	"	- row 1 is header\n"
	"	- column 1 is gene name\n"
	"	- column 2 and after: each is a numeric value list describing feature of each gene in the genome, such as transcript length, ORF length, tissue expression, etc.\n"
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
            } else if (str == "-h" || str == "--help") { 
                help();
            } else {
                message("unknown option: "+str);   
                help();        
            }
        }
    }   
	
	// load foreground genes
	set<string> genes = load_gene_list(genefile);
	message(to_string(genes.size())+" genes loaded from file: "+genefile);
	if (genes.size()<3) 
		{
			message("A minimum of 3 genes is required!");
			exit(1);
		}		
		
	// load background if specified
	set<string> backgroundgenes;
	if(backgroundfile.size()>0)
	{
		backgroundgenes = load_gene_list(backgroundfile);
		message(to_string(backgroundgenes.size())+\
			" background genes loaded from file: "+backgroundfile);
			
		// check overlap
		vector<string> tmp = set_overlap(genes,backgroundgenes);
		
		// remove background genes also in foreground
		// will add back in gene set counting
		if(tmp.size()>0) 
		{
			// remove overlapping genes from background
			backgroundgenes = set_subtract(backgroundgenes,genes);
			message(to_string(tmp.size())+\
					" genes overlap between foreground and background: "+to_string(tmp));
			message(to_string(backgroundgenes.size())+\
				" background genes remain after removing genes also present in foreground");		
		}
		
		if (backgroundgenes.size()<3) 
			{
				message("A minimum of 3 genes is required! Exit");
				exit(1);
			}				
	}	
			
	// load gene names from database, need to make sure each line has unique gene name
	vector<string> allgenes = all_genes_in_feature_table(databasefile);
	
	// remove foreground genes not present in database
	vector<bool> is_foreground_gene,is_background_gene;
			
	vector<string> removed_genes = determine_gene_position_in_feature_table(genes, \
		 allgenes, is_foreground_gene);
			
	if (removed_genes.size()>0)
		message(to_string(removed_genes.size())+\
			" genes not present in gene feature database: "+to_string(removed_genes));
			
	// remove background genes not present in database
	
	if(backgroundfile.size()>0)
	{
		removed_genes = determine_gene_position_in_feature_table (backgroundgenes, \
			 allgenes, is_background_gene);		
		if (removed_genes.size()>0)
			message(to_string(removed_genes.size()) \
				+" background genes not present in gene feature database: "\
					+to_string(removed_genes));
	}
	
	int nSig = find_sig_gene_features(is_foreground_gene, is_background_gene, \
		databasefile, outputfile, p_cutoff);
	message(to_string(nSig) + " significant gene features identified");
	
	// sort output file by p-value and add header
	string header = "feature\tp.value\tt.statistics\tcount.f\tmean.f\tstdev.f\tcount.b\tmean.b\tstdev.b";
	sort_file_and_add_header(outputfile,header,"-k2,2g");
	
	return 0;
}
