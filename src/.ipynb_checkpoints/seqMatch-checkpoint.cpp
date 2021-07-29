#include "sequence.h"
#include "utility.h"

int main(int argc, char* argv[]) { 

  if (argc < 7)
  {
    cout << "Usage:\n";
    cout << "  seq_match -motif <motif.fa> -seq <seq.fa> -out <output_file> -mismatch <n> [-rc,-bed]" << endl;
    cout << "==========\n";
    cout << "Options:" <<endl;
    cout << "==========\n";
    cout << "  -motif    <motif.fa>	required, motif file, fasta format, could include multiple sequences\n";
    cout << "  -seq      <seq.fa>	required, sequence file, fasta format, could include multiple sequences\n";
	cout << "  -out      <output>   output file name. Default: seqmatch_output.txt\n";
    cout << "  -mismatch <n>	    max number of mismatches allowed. default 0\n";
    cout << "  -rc	                also search on reverse complement strand\n";
    cout << "  -bed	                to output BED format. Only if sequence have id like: >chr1:34197116-34197217(+), i.e. output of bedtools getfasta\n";
    cout << "==========\n";
    cout << "Output format\n";
    cout << "==========\n";
    cout << " tabular:\n";
    cout << "==========\n";
    cout << "  chr1:34197116-34197217(+)   21      5ss,T5A:CAGGAAAGT-rc    ACTTTCCTG\n";
    cout << "  column 1: sequence id\n";
    cout << "  column 2: start of the match in sequence\n";
    cout << "  column 3: motif variant id, format: motif_id,mismatch:motif_sequence[-reverse_complement]\n";
    cout << "  column 4: actual match in target sequence\n";
    cout << "==========\n";
    cout << " BED format:\n";
    cout << "==========\n";
    cout << "  column 1-3: genomic coordinates of the match, calculated using sequence ID: chr1:34197116-34197217(+) and match start in the sequence\n";
    cout << "  column 4: number of mismatches\n";
    cout << "  column 6: strand of the match sequence\n";
    cout << "  column 7: strandness of the match relative to the target sequence\n";
    return 1;
  }
  
        string motiffile, seqfile, str;
        string outfile = "seqmatch_output.txt";
        int nmismatch=0;
        bool rc = false;
        bool bed = false;
 
        for (int i = 1; i < argc; i++) { 
            if (i != argc) { 
                str=argv[i];
                if (str == "-motif") {
                   motiffile  = argv[i + 1]; 
                   i=i+1;
                } else if (str == "-seq") {
                    seqfile = argv[i + 1];
                    i=i+1;
                } else if (str == "-out") {
                    outfile = argv[i + 1];
                    i=i+1;
                } else if (str == "-mismatch") {
                    nmismatch = atoi(argv[i + 1]);
                    i=i+1;
                } else if (str == "-rc") {
                    rc = true;
                } else if (str == "-bed") {
                    bed = true;
                } else {
                    cout << "Not enough or invalid arguments, please try again.\n";
                    exit(0);
                }
            }
        }
 /*
  string motiffile = argv[1]; // motif file
  string seqfile = argv[2]; // sequence file
  string outfile = argv[3]; // output file
  int nmismatch = atoi(argv[4]);// # of mismatches allowed
  string _rc = argv[5]; // rc or norc: to search reverse complement strand or not
  string _bed = argv[6]; // to make bed output or not. only valid if sequence have id like: hg18_chr12_52642148_52644699_+
  bool rc = false;
  if (_rc == "rc") {rc=true;}
  bool bed = false;
  if (_bed == "bed") {bed=true;}
  */
    // alphabet
  string ACGT_str = "ACGT";
  set<char> alphabet;
  for (int i=0;i<4;i++)
  {
    alphabet.insert(ACGT_str[i]);
  }

  array<int,2> res = match(motiffile, seqfile, outfile, nmismatch, rc, alphabet);
  
  if ( bed && res[0]>0)
  {
    cout <<"["<<current_time()<<"] creating BED format output..."<<endl;
    string bedfile = outfile;
    bedfile += ".bed";
    tab2bed(outfile,bedfile);
    //if (tab2bed(outfile,bedfile))
    //{
    // rename(bedfile.c_str(),outfile.c_str());
    //}
  }
  cout <<"["<<current_time()<<"] Done"<<endl;
  return 0;
}

