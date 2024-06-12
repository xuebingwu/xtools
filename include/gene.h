#ifndef __GENE_H__
#define __GENE_H__

#include <vector>
#include <set>
#include <string>
#include <map>

using namespace std;


// load gene names 
set<string> load_gene_list(string filename, bool toupper=true, int skip=0, int col=0);

int find_sig_gene_sets_correlate_two_scores(map<string,double> scores1, map<string,double> scores2, string geneSetFile, string outputfile, double r_cutoff);

void load_weighted_gene_list(string filename, vector<string> &genes, vector<double> &scores, bool toupper=true, int skip=0, int col1=0, int col2=1);

void load_weighted_gene_list(string filename, map<string, double> &scores, bool toupper=true, int skip=0, int col1=0, int col2=1);

void gene_feature_correlation(map<string,double> scores, string geneFeatureFile, string outputfile);

// return the list of all gene names in the database
set<string> all_genes_in_gene_sets(string filename);

// return the list of all gene names in feature table
// assume all gene ids are unique
vector<string> all_genes_in_feature_table(string filename);

// remove gene ids not recognized in gene sets
vector<string> filter_genes_not_in_gene_sets(set<string> &genes, set<string> allgenes);

// find position/row number and also remove genes not found
vector<string> determine_gene_position_in_feature_table(set<string> genes, vector<string> allgenes, vector<bool> &is_foreground_gene);

int find_sig_gene_features(vector<bool> is_foreground_gene, vector<bool> is_background_gene, string geneFeatureFile, string outputfile, double p_cutoff=0.05);

// find significnat gene sets
// need total number of genes as input
// at least two genes to do the test
int find_sig_gene_sets(set<string> genes, string geneSetFile, string outputfile, unsigned N, double p_cutoff=0.05);
int find_sig_gene_sets_with_background(set<string> genes, set<string> backgroundgenes, string geneSetFile, string outputfile, double p_cutoff=0.05);

int find_sig_gene_sets_weighted(vector<string> genes, vector<double> scores, string geneSetFile, string outputfile, double p_cutoff=0.05);


#endif
