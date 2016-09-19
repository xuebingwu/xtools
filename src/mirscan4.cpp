#include "structure.h"
#include "utility.h"
#include "text.h"
#include <iostream>

void help()
{
    string str =
	"\n"
    "mirscan4: de novo pri-microRNA prediction\n"
	"	- Xuebing Wu (wuxb07@gmail.com)\n"
    "\n"
    "Usage: mirscan4 -i input -o output \n"
    "\n"
    "Options:\n"
    "\n"
    "   -i  <input>              input file (*.fa, *.fasta, *.lfold, *.fold, *.hairpin.fold)\n"
    "   -o  <output>             output file prefix\n"
	"   -min_len <integer>       minimum hairpin length (all nucleotides), default=50\n"
	"   -max_loop <integer>      maximum loop length, default=100\n"
	"   -min_stem <integer>      minimum stem length, default=20\n"
	"   -max_stem <integer>      maximum stem length, default=80\n"
	"   -short_stem <integer>    maximum stem length to unpair, default=20\n"
	"   -max_bulge <integer>     maximum bulge size, default=10\n"	
	"   -mfe <double>            remove structures with minimum free energy higher than this, default=-30\n"			
	"   -profile_len <integer>   profile length, default=50\n"	
	"   -extend  <integer>       number of flanking nucleotides to consider, default=30\n"	
	"   -rnafold_opts \"options\" RNAfold options, quoted, default=\"--noPS --noClosingGU\"\n"	
/*	"training mode =============\n"
	"   -bowtieIndex <ebwt>     bowtie index for mapping and identifying flanking sequences, default hg19\n"
	"   -genomeFasta <fasta>    whole genome fasta file, default hg19\n"
	"prediction mode ===========\n"
	"   -predict  <file>        model file prefix. If not specified, run training mode\n"
*/	"\n"
	"Workflow:\n"
	"\n"
	"     sequence (*.fa, *.fasta)   -OR-  locally stable fold (*.lfold by RNALfold)    	\n"
	"                    \\                          /										\n"
	"        RNAfold -->  \\                        / <-- filter and reformat  				\n"
	"                      v                      v                          				\n"	
	"                     filtered structure (*.fold)			    						\n"
	"                               |														\n"
	"                               | force to fold as a single hairpin (by RNAduplex)		\n"
	"                               v  													\n"
	"                      single hairpin fold (*.hairpin.fold)                            	\n"
	"                               | 														\n"
	"                               | remove short stems									\n"	
	"                               v                                                       \n"
	"                   1. feature matrix for machine learning    							\n"
	"                   2. position on the genome (*.bed)	        						\n"
	"                   3. visualize structure              		    					\n"		
	"\n"
    "\n";

    cerr << str;

    exit(0);
}
/* pipeline

# for positives, take 30nt flanking cut site and fold
cd /lab/bartel1_ata/wuxbl/projects/wenwen/top56hexamer/4col_species
python convert_4col_to_fasta.py hsa_4col 30 > hsa.ext.30.fa

# fold 
RNAfold --noPS   < hsa.ext.30.fa > hsa.ext.30.fold



*/

int main(int argc, char* argv[]) {

    // default
    string input;
    string output;
	string rnafold_opts = "--noClosingGU";
	
/*	string bowtie_index;
	string bowtie_index_default="/nfs/genomes/human_gp_feb_09_no_random/bowtie/hg19";
	string genome_fasta;
	string genome_fasta_default="/nfs/genomes/human_gp_feb_09_no_random/fasta_whole_genome/hg19.fa";
	*/
	string model;
	
	int ext=30;
	int min_hairpin_length = 50;
	int max_loop_length = 100;
	int max_stem_length = 80;
	
	double mfe = -30; // minimum MFE
	
	// stem of this length if flanked by unpaired region of this length or longer, will become unpaired
	int max_short_stem_length = 20; 
	int min_pairs_left=20; // after removing short stems, at least this number of paired bases
	int max_bulge = 10; // largest bulge allowed
	int profile_length = 50;
	
	// rnalfold:
	int max_overlap_allowed = 20;
		
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
           // } else if (str == "-predict") { 
           //     model = argv[i + 1];
           //     i=i+1;
            } else if (str == "-extend") { 
                ext = stoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-min_len") { 
                min_hairpin_length = stoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-max_loop") { 
                max_loop_length = stoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-max_bulge") { 
                max_bulge = stoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-min_stem") { 
                min_pairs_left = stoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-max_stem") { 
                max_stem_length = stoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-short_stem") { 
                max_short_stem_length = stoi(argv[i + 1]);
                i=i+1;
             } else if (str == "-profile_len") { 
                profile_length = stoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-mfe") { 
                mfe = stof(argv[i + 1]);
                i=i+1;
            } else if (str == "-rnafold_opts") { 
                rnafold_opts = argv[i + 1];
                i=i+1;
 /*           } else if (str == "-bowtieIndex") { 
                bowtie_index = argv[i + 1];
                i=i+1;
            } else if (str == "-genomeFasta") { 
                genome_fasta = argv[i + 1];
                i=i+1;
  */          } else if (str == "-h" || str == "--help") { 
                help();
            } else {
                message("unknown option: "+str);   
                help();        
            }
        }
    }   
	
	bool noClosingGU = false;
	if(rnafold_opts.find("noClosingGU") != std::string::npos) noClosingGU = true;
		
	/*
	determine input type so that the program can start from intermediate outputs of a previous run
	1. *.fa/*.fasta: raw sequences
	2. *.lfold: output from RNALfold
	3. *.fold: convert *.lfold format to *.fold format
	4. *.hairpin.fold: remove branches and re-fold as a hairpin using RNAduplex
	*/
		
	//if (input.substr(input.size()-13,13) != ".hairpin.fold")
	if(ends_with(input,".hairpin.fold")){
		system_run("ln -s "+input+" "+output+".hairpin.fold");	
	} else {		
		if (ends_with(input,".fold"))
		{
			system_run("ln -s "+input+" "+output+".fold");
		} else if (ends_with(input,".lfold")) 
		{
			message("filter RNALfold output...");//except the free energy goes to sequence id
			RNALfold_to_RNAfold(input,output+".fold",min_hairpin_length,min_pairs_left,ext,mfe);
			// refold after extension
			//message("fold again including flanking sequences...");
			//system_run("RNAfold --noPS " + rnafold_opts + " < "+output+".lfold.filtered > "+output+".fold");
		} else if (ends_with(input,".fa") || ends_with(input,".fasta"))
		{
			// starting from fasta
			message("fold sequences...");
			system_run("RNAfold --noPS " + rnafold_opts + " < " + input + " > "+output+".fold");
		} 
	
		message("remove branches and re-fold hairpin using RNAduplex");
		hairpin_RNAduplex(output+".fold", output, rnafold_opts);
		RNAduplex_to_RNAfold(output+".duplex", output+".hairpin.fold");
	} 
		
	//message("generating plots...");
    //system_run("perl /lab/bartel1_ata/wuxbl/scripts/plothairpin/RNA-HairpinFigure-master/plotHairpin.pl "+output+".hairpin.fold > "+output+".hairpin.fold.plot");
	
	
	// remove short stems
	message("removing short stems...");
	remove_short_stem_from_file(output+".hairpin.fold", output+".hairpin.fold.no.short.stem", max_short_stem_length, min_pairs_left,max_bulge,noClosingGU);
	
	// get feature, output to screen
	message("generating features...");
	mirna_feature_from_file(output+".hairpin.fold.no.short.stem",output+".hairpin.basal.feat",output+".hairpin.loop.feat",35);
	
	message("generating bed file for filtered hairpins...");
	system_run("cat "+output+".hairpin.fold.no.short.stem | grep -v '((' | sed 's/U/T/g' > "+output+".hairpin.fold.no.short.stem.fa");
	system_run("/lab/bartel1_ata/wuxbl/scripts/sequence/fasta2bed -i "+output+".hairpin.fold.no.short.stem.fa -g /nfs/genomes/human_gp_feb_09_no_random/bowtie/hg19 -o "+output+".hairpin.fold.no.short.stem.bed");
	
	// plot
    system_run("perl /lab/bartel1_ata/wuxbl/scripts/plothairpin/RNA-HairpinFigure-master/plotHairpin.pl "+output+".hairpin.fold.no.short.stem  > "+output+".hairpin.fold.no.short.stem.plot");
	
	
	// make test data
	
	

    return 0;
}
















