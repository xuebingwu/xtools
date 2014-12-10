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
    "geneSet: A program to discover whether some predefined gene sets are over-represented in your list of genes\n"
	"	- Xuebing Wu (wuxb07@gmail.com)\n"
    "\n"
    "Usage: geneSet -i input -o output\n"
    "\n"
    "Options:\n"
    "\n"
    "	-i  input       input file\n"
	"	-b  background  optional background gene list\n"
    "	-o  output      output file\n"
	"	-d  database    gene set annotation file (such as *.gmt)\n"	
	"	-p  pcutoff     p-value cutoff, default 0.05\n"
	"	-h/--help       print help message\n"
    "\n"
	"Notes: database file format\n"
	"	GMT: Gene Matrix Transposed file format (*.gmt)\n"
	"	http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29"
	"	- each row is a gene set\n"
	"	- columns are tab-delimited\n"
	"	- column 1: gene set name, unique\n"
	"	- column 2: a brief descripiton of the gene set\n"
	"	- column 3 and after: each column is a gene name\n";

    cerr << str;

    exit(0);
}

int main(int argc, char* argv[]) {
	
	string genefile;
	string databasefile;
	string outputfile;
	string backgroundfile;
	double p_cutoff = 0.05;
	
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
	
	// load gene names from database
	set<string> allgenes = all_genes_in_gene_sets(databasefile);
	message(to_string(allgenes.size())+" genes in gene set database: "+databasefile);
	
	// load foreground genes
	set<string> genes = load_gene_list(genefile);
	message(to_string(genes.size())+" genes loaded from file: "+genefile);
	if (genes.size()<3) 
		{
			message("A minimum of 3 genes is required!");
			exit(1);
		}		

	vector<string> removed_genes = filter_genes_not_in_gene_sets(genes, allgenes);
	if (removed_genes.size()>0)
	{
		message(to_string(removed_genes.size())+\
		" genes not present in gene set database: "+to_string(removed_genes));
		message(to_string(genes.size())+" genes left for enrichment analysis");
	}
	
	// load background if specified
	set<string> backgroundgenes;
	if(backgroundfile.size()>0)
	{
		backgroundgenes = load_gene_list(backgroundfile);
		message(to_string(backgroundgenes.size())+\
			" background genes loaded from file: "+backgroundfile);
		
		vector<string> removed_genes = filter_genes_not_in_gene_sets(backgroundgenes, allgenes);
		if (removed_genes.size()>0)
		{
			message(to_string(removed_genes.size())+\
			" background genes not present in gene set database: "+to_string(removed_genes));
		}
		// check overlap
		vector<string> tmp = set_overlap(genes,backgroundgenes);
		
		// remove background genes also in foreground
		// will add back in gene set counting
		if(tmp.size() < genes.size()) // not all foreground genes in background
		{
			message("Warning: not all foreground genes are in backgorund. \
				 Will merge foreground and background genes to make the new background");
			// remove overlapping genes from background
			message(to_string(genes.size() - tmp.size())+\
					" foreground genes added to background: "+to_string(tmp));	
		}
		backgroundgenes = set_subtract(backgroundgenes,genes);
		message(to_string(backgroundgenes.size()+genes.size())+\
			" background genes in total");	
		if (backgroundgenes.size() < genes.size()) 
			{
				message("The background gene list is too small. Should contain at least twice the number of foreground genes. Will ignore user specified background and use the entire genome with annotations!");
				backgroundfile = "";
			}				
	}	
	
	int nSig;
	if(backgroundfile.size() == 0) 
		nSig = find_sig_gene_sets(genes, databasefile, outputfile, allgenes.size(), p_cutoff);
	else
		nSig = find_sig_gene_sets_with_background(genes, backgroundgenes, databasefile, outputfile, p_cutoff);
	
	message(to_string(nSig) + " significant gene features identified");
	
	// sort output file by p-value and add header
	string header = "gene.set\tdescription\tp.value\tfold.enrichment\tobserved\texpected\ttotal.positive\ttotal.genes\toverlap.genes";
	sort_file_and_add_header(outputfile,header,"-k3,3g");
	
	return 0;
}
