#ifndef __SEQUENCE_H__
#define __SEQUENCE_H__

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <iomanip>
#include <array>        // std::array

//#include <boost/regex.hpp>

// for position matrix
#include <boost/numeric/ublas/matrix.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

extern "C"{
#include "ushuffle.h"
}

using namespace std;

// build a model that can be used to score sequences
// skip some lines
// specify column number (0-based) of kmer, start, end
class positional_kmer
{
public:
	string seq; // sequence, IUPAC code
	int pos;
	int size;
	double weight;
	int group; // if part of another, will have the same group number
	positional_kmer();
	positional_kmer(const positional_kmer &a);
	positional_kmer(string seq, int pos, int size, double weight, int group);
	bool equals(positional_kmer a);  
	bool is_part_of(positional_kmer a);
	const positional_kmer &operator=(const positional_kmer &a);
	string as_string(string del="_");
};


vector<positional_kmer> build_model_from_PKA2_output(string filename, double pCutoff, double pCutoff_B);
vector<positional_kmer> build_model_from_PKA_output(string filename);
void save_model_to_file(vector<positional_kmer> ranked_kmers, string filename);

vector<positional_kmer> load_model_from_file(string filename);

double score_sequence_using_PKA_model(vector<positional_kmer> ranked_kmers, string seq);



void score_fasta_using_PKA_model(string seqfile, string outputfile, vector<positional_kmer> ranked_kmers);
void score_tabular_using_PKA_model(string tabfile, int col, string outputfile, vector<positional_kmer> ranked_kmers);


class paired_kmer
{
public:
	string seq1;
	string seq2;
	int dist; // distance between seq1 start and seq2 start
	int pos; // position of seq1
	int shift;
	double weight;
	paired_kmer();
	paired_kmer(const paired_kmer &a);
	paired_kmer(string seq1, string seq2, int dist=0, int shift=0, int pos=0, double weight=0);
	const paired_kmer &operator=(const paired_kmer &a);
	string as_string(string del="_", bool add_shift=true,bool add_pos=false, bool add_weight=false);
	bool equals(paired_kmer a, bool identical_pos=true, bool identical_shift=true);
	int len(); // total length = seq1 + seq2 + gap
	int gap();
};

vector<paired_kmer> build_paired_kmer_model(string filename);

double score_sequence_using_paired_kmer_model(vector<paired_kmer> model,  string seq);

vector<paired_kmer> generate_paired_kmers (
	string alphabet,
int seq1_len,
int seq2_len,
int max_dist,
int min_dist=1,
int max_shift=0,
int min_shift=0	
);

array<int,2> find_significant_pairs_from_weighted_sequences(
	vector<string> seqs,
	vector<double> weights, 
	vector<paired_kmer> paired_kmers, 
	string outfile, 
	int nTest, 
	double pCutoff=0.05, 
	double pCutoff_B=0.05,
	int startPos=0) ;

	
void plot_PKA2_nucleotide_output(string infile, string outfile, int lSeq);




void load_weighted_sequences_to_vectors(string filename, vector<string> &seqs, vector<double> &weights,int skip=0, int cSeq=1,int cWeight=2);

vector<string> load_ranked_sequences_to_vectors(string filename, int skip=0, int cSeq=1);


array<int,2> find_significant_kmer_from_ranked_sequences(vector<string> seqs, vector<string> kmers, string outfile, int nTest, double pCutoff=0.05, double pCutoff_B=0.05, int shift=0, int startPos=0 );

array<int,2> find_significant_kmer_from_weighted_sequences(vector<string> seqs,vector<double> weights, vector<string> kmers, string outfile, int nTest, double pCutoff=0.05, double pCutoff_B=0.05,int shift_min=0, int shift_max=2, int startPos=0);

void save_feature_matrix(map<string,string> seqs, vector<string> kmers, string outfile, string label, bool append=false, int shift=0);

void read_significant_positional_kmer_from_file(string inputfile, vector<string> &kmers, vector<int> &positions);

map<string, string> seq_vector2map(vector<string> seqs);

void read_significant_positional_kmer_from_PKA2_output(string inputfile, vector<string> &kmers, vector<int> &positions);

void significant_feature_matrix(map<string,string> seqs, vector<string> kmers, vector<int> positions, string outfile, string label, bool append=false, int shift=0);

void significant_feature_matrix_PKA2(vector<string> seqs, vector<double> weights, vector<positional_kmer> ranked_kmers, string outfile);

// convert fasta file to a letter matrix
void fasta_to_letter_matrix(string input, string output);

// convert fasta file to a letter matrix, no header
void tab_seq_to_letter_matrix(string input, string output, int k_min, int k_max, int col, int skip);

// PKA: remove overlapping motifs
// input: all significant motifs
int non_overlapping_sig_motifs(string inputfile, string outputfile);

// count substring in a map of sequences, i.e. from fasta
// used in generating markov model
int countSubstringInSeqs(map<string,string>seqs, string sub);

// of not of identical length, return -1
int sequence_similarity(string a, string b);

// calculate and write pairwise similarity matrix of
void pairwise_sequence_similarity_matrix(vector<string> seqs, string filename);

	
set<int> findall(string seq, string motif);

// dict for IUPAC degenerate nucleotides
map<char,string> define_IUPAC();


// generate all possible combinations of letters in alphabet of fixed length k, i.e. kmer
vector<string> generate_kmers(int k, string alphabet);

// degenerate kmers based on IUPAC, DNA only
// remove those with terminal N, which is not a kmer, but (k-i)mer
vector<string> degenerate_kmer(int k, string alphabet="ACGTRYMKWSBDHVN");


// convert one degenerate kmer to regular expression, e.g. ARG -> A[AG]G
string degenerate_kmer_to_regex(string kmer,map<char,string> iupac);

// expand a degenerate kmer to all possible exact kmers
vector<string> expand_degenerate_kmer(string seq, map<char,string> iupac);

// implant a motif to a set of sequences
// motif can contain degenerate nucleotides
void implant_motif(map<string,string> &seqs, int position, string motif, double fraction);


// write pwm in meme format
void write_pwm_in_meme_format(boost::numeric::ublas::matrix<double> pwm, string motifname, string filename);


// build position weight matrix for a positional kmer, can include flanking sequence
boost::numeric::ublas::matrix<double> create_position_weight_matrix_from_kmer(map<string,string> seqs, string kmer, int position, map<char,string> iupac, int d, int startPos);

//PKA : create logo for a single kmer
void create_logo_for_kmer(map<string,string> seqs, string kmer, int position, map<char,string> iupac, int d, int startPos,string output);

// create logos for top 1 kmer passing pCutoff_B
void create_logo_for_topN_sig_kmer_per_position(map<string,string> seqs, string filename,int d, double pCutoff_B,int startPos,string output);


//vector<int> count_one_kmer_in_all_seqs_regex(map<string,string> seqs, boost::regex pattern, int k, int seq_len_minus_k_plus_1, int k_plus_shift);



// doesn't work with shift+degenerate
// for each kmer, count its freq at each position(start), allowing shifts (to the right)
map<string,vector<int> > count_all_kmer_in_seqs(vector<string> kmers, map<string,string> seqs, int shift);  

//PKA
void print_kmer_positional_profile(map<string,vector<int> > data);


map<string,vector<int> > degenerate_kmer_counts(vector<string> dkmers,map<string,vector<int> > data, map<char,string> define_iupac);

array<int,2> find_significant_kmer_from_one_seq_set(map<string,string>seqs1, map<string,double> probs_kmer,vector<string>kmers, vector<string> dkmers, int shift,bool degenerate,double pCutoff, double pCutoff_B, int startPos,int nTest, string outfile, string output_count_file);

// two file comparison, not allow shift and degenerate at the same time
array<int,2> find_significant_kmer_from_two_seq_sets(map<string,string>seqs1, map<string,string>seqs2, vector<string>kmers, vector<string> dkmers, int shift,bool degenerate,double pCutoff,double pCutoff_B, double pseudo,int startPos,int nTest, string outfile,string output_count_file);


void plot_frequency_for_significant_kmer(string inputfile, string outputfile);
   

char complement(char ch);

string reverseComplement(string seq);

void ReadOneSeqFromFasta(ifstream& infile, string& name, string& seq);

map<string,string> ReadFasta(string filename);

void WriteFasta(map<string,string> seqs, string filename);


//uShuffle
//http://digital.cs.usu.edu/~mjiang/ushuffle/
// shuffle sequence preserving k-let
string shuffle_seq_preserving_k_let(string str,int k);

//
map<string,string> shuffle_seqs_preserving_k_let(map<string,string> seqs, int N, int k);

map<string,string> first_n_bases(map<string,string> seqs,int n);

map<string,string> last_n_bases(map<string,string> seqs,int n);


void mismatches(map<string,string>& mutant,map<string,int>& dist, string motif, int n, set<char> alphabet);

map<string,string> ExpandMotifs(map<string,string>& motifs, int nmismatch, bool rc, set<char> alphabet) ;


array<int,2> match(string motiffile, string seqfile, string outfile, int nmismatch, bool rc, set<char> alphabet);



int tab2bed_galaxy(string infile, string outfile);


int tab2bed_bedtools(string infile, string outfile);

#endif
