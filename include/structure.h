#ifndef __STRUCTURE_H__
#define __STRUCTURE_H__

#include <fstream>

#include "positional_kmer.h"

using namespace std;

int plot_hairpin(string id, string seq, string fold);
int plot_hairpin_from_file(string infile);

int RNAduplex_to_RNAfold(string infile, string outfile);

int hairpin_RNAduplex(string infile, string outfile, string options="");

int trim_hairpin(string seq, string& fold,bool noClosingGU=true);

vector<int> find_loops(string structure);

int hairpin_stem_length(string fold);

int mirna_basal_stem_definition(string id, string seq, string fold, int stem_len, ofstream& out);

int mirna_loop_definition(string id,string seq, string fold, int stem_len, ofstream& out);

int mirna_feature_from_file(string infile, string outfile_basal, string outfile_loop, int stem_len);

int loop_adjacent_stem_length(string structure);

string reverse_structure(string x);

int open_large_bulge(string &structure, int max_bulge=10);

int open_large_bulge_left_side(string &structure, int max_bulge=10);

int remove_shortest_branch(string &structure);
int remove_all_branches(string &structure);

void RNALfold_to_RNAfold(string infile, string outfile, int min_hairpin_length=50,int min_pairs_left=25, int ext=20, int mfe=30);
void RNALfold_to_RNAfold_filter_overlap(string infile, string outfile, int min_hairpin_length=50, int ext=20,int max_overlap_allowed=20);

void convert_RNALfold_to_RNAfold_output(string filename);

void hairpin_scoring(string infile, string outfile, vector<positional_kmer> model_str,vector<positional_kmer> model_top,vector<positional_kmer> model_bot,vector<positional_kmer> model_loop, int model_length);

int force_unpair(string &structure, int start, int end);
int remove_short_stem(string &structure, int max_paired_length=4);
bool remove_short_stem_new(string &structure);

int remove_short_stem_from_file(string inputfile, string outputfile, int max_paired_length=4,int min_pairs_left=25,int max_bulge=10, bool noClosingGU=true);

int keep_longest_stem(string struc);

int count_pairs(string str);

int count_unpaired(string struc);

#endif
