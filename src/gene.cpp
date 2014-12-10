#include "utility.h"
#include "stat.h"
#include "gene.h"
#include "container.h"

#include <fstream>
#include <vector>
#include <set>
#include <array>
#include <map>
#include <string>
#include <iostream>
#include <algorithm>
#include <ctime>



using namespace std;

// load gene names 
set<string> load_gene_list(string filename, bool toupper/*=true*/, int skip/*=0*/, int col/*=0*/)
{
	set<string> genes;
	
	ifstream fin;
	fin.open(filename.c_str());

	string line;
	vector<string> flds;
	
	int i=0;
	while(fin && i<skip) getline(fin,line);
	
	while(fin)
	{
		getline(fin,line);
		if(line.length() == 0)
			continue;
		flds = string_split(line,"\t");
		if (toupper) genes.insert(to_upper(flds[col]));
		else genes.insert(flds[col]);
	}
	fin.close();
	
	return genes;	
}

// load gene names and their scores
void load_weighted_gene_list(string filename, vector<string> &genes, vector<double> &scores, bool toupper/*=true*/, int skip/*=0*/, int col1/*=0*/, int col2/*=1*/)
{	
	ifstream fin;
	fin.open(filename.c_str());

	string line;
	vector<string> flds;
	
	int i=0;
	while(fin && i++ < skip) getline(fin,line);
	
	while(fin)
	{
		getline(fin,line);
		if(line.length() == 0)
			continue;
		flds = string_split(line,"\t");
		//cout << flds[col1] << ","<< to_upper(flds[col1]) << endl;
		if (toupper) genes.push_back(to_upper(flds[col1]));
		else genes.push_back(flds[col1]);
		scores.push_back(stof(flds[col2]));
	}
	fin.close();
}

// load gene names and their scores, to a map
// for geneScore with features
void load_weighted_gene_list(string filename, map<string, double> &scores, bool toupper/*=true*/, int skip/*=0*/, int col1/*=0*/, int col2/*=1*/)
{	
	ifstream fin;
	fin.open(filename.c_str());

	string line;
	vector<string> flds;
	
	int i=0;
	while(fin && i++ < skip) getline(fin,line);
	
	while(fin)
	{
		getline(fin,line);
		if(line.length() == 0)
			continue;
		flds = string_split(line,"\t");
		if (toupper) scores[to_upper(flds[col1])] = stof(flds[col2]);
		else scores[flds[col1]] = stof(flds[col2]);
	}
	fin.close();
}

// return the list of all gene names in the database
set<string> all_genes_in_gene_sets(string filename)
{
	set<string> allgenes;
	
	ifstream fin;
	fin.open(filename.c_str());
	
	string line;
	vector<string> flds;

	while(fin)
	{
		getline(fin,line);
		if(line.length() == 0)
			continue;
		flds = string_split(line,"\t");	
		for(unsigned i=2;i<flds.size();i++) allgenes.insert(flds[i]);
	}
	fin.close();	
	return allgenes;
}

// return the list of all gene names in feature table
// assume all gene ids are unique
vector<string> all_genes_in_feature_table(string filename)
{
	vector<string> allgenes;
	
	ifstream fin;
	fin.open(filename.c_str());
	
	string line;
	vector<string> flds;
	
	// skip header
	getline(fin,line);
		
	while(fin)
	{
		getline(fin,line);
		if(line.length() == 0)
			continue;
		flds = string_split(line,"\t");	
		allgenes.push_back(flds[0]);
	}
	
	fin.close();
		
	return allgenes;
}

// remove gene ids not recognized in gene sets
vector<string> filter_genes_not_in_gene_sets(set<string> &genes, set<string> allgenes)
{
	vector<string> genes_removed;
	
	// need reverse iterator otherwise will miss genes next to gene to be deleted
	set<string>::reverse_iterator it; 
	for (it=genes.rbegin();it!=genes.rend();it++)
	{
		if (allgenes.find(*it) == allgenes.end())
		{
			genes_removed.push_back(*it);
			genes.erase(*it);
		}
	}
	return genes_removed;
}

// find position/row number and also remove genes not found
vector<string> determine_gene_position_in_feature_table(set<string> genes, vector<string> allgenes, vector<bool> &has_gene)
{
	vector<string> genes_removed;
	
	// initialize has_gene = false at every position
	has_gene.clear();
	for(unsigned i=0;i<allgenes.size();i++) has_gene.push_back(false);
	// need reverse iterator otherwise will miss genes next to gene to be deleted
	set<string>::reverse_iterator it;
	for (it=genes.rbegin();it!=genes.rend();it++)
	{
		int found = -1;
		for(unsigned i=0;i<allgenes.size();i++)
		{
			//cout << *it << "\t" << allgenes[i] << endl;
			if (*it == allgenes[i]) found = i;
		}
		if (found == -1) // not found
		{
			genes_removed.push_back(*it);
			//cout <<"not found " << *it << endl;
			//genes.erase(*it);
		} else {
			has_gene[found] = true;
			//cout << "found "<< *it << endl;
		}
	}
	return genes_removed;
}

int find_sig_gene_features(vector<bool> is_foreground_gene, vector<bool> is_background_gene, string geneFeatureFile, string outputfile, double p_cutoff)
{	
	ifstream fin;
	fin.open(geneFeatureFile.c_str());
	
	ofstream fout;
	fout.open(outputfile.c_str());

	string line;
	vector<string> flds;
	
	int nSig = 0;
	
	// header
	getline(fin,line);
	vector<string> header = string_split(line,"\t");
		
	for (unsigned i=1;i< header.size();i++)// for each feature
	{
		fin.clear(); // clear bad state after eof
		fin.seekg( 0 ); // back to line 1
		getline(fin,line);// skip header
		
		int line_number = 0;
	
		vector<double> foreground;
		vector<double> background;
		
		// processing feature i
		while(fin)
		{
			getline(fin,line);
			if(line.length() == 0)
				continue;
			flds = string_split(line,"\t");
			if(is_foreground_gene[line_number]) foreground.push_back(stof(flds[i]));
			else if (is_background_gene.size() == 0) // no background genes provided, then all other genes are background
				{
					background.push_back(stof(flds[i]));
				}
			else if (is_background_gene[line_number]) // if background genes provided, check if in background
				background.push_back(stof(flds[i]));
			
			line_number ++;
		}
		//cout << foreground.size() << endl;
		//cout << background.size() << endl;
		
		if (foreground.size()<3) 
			{
				message("A minimum of 3 genes is required ("+to_string(foreground.size())+" given)");
				exit(1);
			}		
		if (background.size()<3) 
			{
				message("A minimum of 3 background genes is required ("+to_string(background.size())+" given)");
				exit(1);
			}		
		// test
		
		array<double,6> ttest = t_test(foreground, background);
		if(ttest[1] > p_cutoff) continue;
		fout << header[i] << "\t" << ttest[1] <<  "\t" << ttest[0] << "\t" << foreground.size() << "\t" << ttest[2] << "\t" << ttest[3] << "\t" << background.size() << "\t" << ttest[4] << "\t"  << ttest[5] << endl;
		nSig ++;
	}
	
	fin.close();
	fout.close();
	
	return nSig;
}

// given a list of genes with weight
// compare to each feature in the database
// shuffle input gene list 100 times to compute p-value
void gene_feature_correlation(map<string,double> scores, string geneFeatureFile, string outputfile)
{
	// determine a common set of genes 
	ifstream fin;
	fin.open(geneFeatureFile.c_str());
	
	ofstream fout;
	fout.open(outputfile.c_str());

	string line;
	vector<string> flds;
		
	// header
	getline(fin,line);
	vector<string> header = string_split(line,"\t");
			
    vector< vector<double> > all_scores; // all scores for shared genes
	for (int i=0; i< header.size(); i++) 
	{
		vector<double> tmp;
		all_scores.push_back(tmp);
	}
	
	map<string,double>::iterator it;
		
	while(fin)	
	{
		getline(fin,line);// skip header
		if(line.length()==0) continue;
		flds = string_split(line,"\t");
		it = scores.find(flds[0]);
		if (it != scores.end())// shared genes
		{
			all_scores[0].push_back(it->second);
			for (int i=1; i< header.size(); i++)
			{
				all_scores[i].push_back(stof(flds[i]));
			}
		}		
	}
	message(to_string(all_scores[0].size())+" genes in common between query genes and database");
	
	vector<double> random_score = all_scores[0];
	for (int i=1;i<header.size();i++)
		{
			double r = cor(all_scores[0],all_scores[i]);
			message("correlation with "+header[i]+": r="+to_string(r));
			/*
			// shuffle, not very useful
			int n = 0;
			for (int j=0;j<100;j++)
			{
				srand (unsigned ( std::time(0) ) * j);
				random_shuffle(random_score.begin(),random_score.end());
				cout << cor(random_score,all_scores[i]) << endl;
				if( r < cor(random_score,all_scores[i])) n ++;
			}
			message("permutaiton test p value = "+to_string(double(n)/100));
			*/
			
		}	
}


// find significnat gene sets
// need total number of genes as input
// at least two genes to do the test
int find_sig_gene_sets(set<string> genes, string geneSetFile, string outputfile, unsigned N, double p_cutoff)
{
	unsigned n = genes.size(); // total samples
	
	ifstream fin;
	fin.open(geneSetFile.c_str());
	
	ofstream fout;
	fout.open(outputfile.c_str());

	string line;
	vector<string> flds;
	
	int nSig = 0;
	
	while(fin)
	{
		getline(fin,line);
		if(line.length() == 0)
			continue;
		// replace space with underscore
		replace(line.begin(),line.end(),' ','_');
		replace(line.begin(),line.end(),',',';');
		flds = string_split(line,"\t");
		//message("processing gene set "+to_string(i++)+": "+flds[0],"\r");
		//cout << flds[0] << "\t" << flds[1] << endl;	
		set<string> geneSet(flds.begin()+3,flds.end());
		//cout << geneSet.size() <<" genes in gene set" << endl;
		vector<string> overlap = set_overlap(genes,geneSet);
		//cout << overlap.size() << " genes overlap" << endl;
		// statistics
		unsigned k = overlap.size(); // observd
		if (k < 2) continue;
		unsigned r = geneSet.size(); // total positive
		double p = hypergeometric_test( k,  r,  n,  N);
		if (p> p_cutoff) continue;
		nSig ++;
		double expected = double(n) * r / N;
		double fold_enrichment = k / expected;
		// geneset_name, geneset_annotation, p_value, enrichment, observed, expected, total_pos, total, all_overlap_genes
		fout << flds[0] << "\t" << flds[1] << "\t" << p << "\t" << fold_enrichment << "\t" << k << "\t" << expected << "\t" << r << "\t" << N  << "\t" << to_string(overlap,",") << endl;
		geneSet.clear();
		overlap.clear();
	}
	fin.close();
	fout.close();
	
	return nSig;
}

// treat each gene set as foreground do two sample t test
int find_sig_gene_sets_weighted(vector<string> genes, vector<double> scores, string geneSetFile, string outputfile, double p_cutoff)
{	
	ifstream fin;
	fin.open(geneSetFile.c_str());
	
	ofstream fout;
	fout.open(outputfile.c_str());

	string line;
	vector<string> flds;
	
	int nSig = 0;
	
	int i = 1;
	while(fin)
	{
		getline(fin,line);
		if(line.length() == 0)
			continue;
		// replace space with underscore
		replace(line.begin(),line.end(),' ','_');
		replace(line.begin(),line.end(),',',';');
		flds = string_split(line,"\t");
		//message("processing gene set "+to_string(i++)+": "+flds[0]);
		//cout << flds[0] << "\t" << flds[1] << endl;	
		set<string> geneSet(flds.begin()+3,flds.end());
		//cout << geneSet.size() <<" genes in gene set" << endl;
		
		vector<double> foreground;
		vector<double> background;
		
		set<string>::iterator it;
		for (int i =0; i< genes.size();i++)
		{
			it = geneSet.find(genes[i]);
			if(it == geneSet.end()) background.push_back(scores[i]);
			//if(find(geneSet.begin(),geneSet.end(),genes[i]) == geneSet.end()) background.push_back(scores[i]);
			else foreground.push_back(scores[i]);
		}
		if (foreground.size()<3 || background.size() < 3) continue;	
		
		// test
		array<double,6> ttest = t_test(foreground, background);
		if(ttest[1] > p_cutoff) continue;
		fout << flds[0] << "\t" << flds[1] << "\t" << ttest[1] <<  "\t" << ttest[0] << "\t" << foreground.size() << "\t" << ttest[2] << "\t" << ttest[3] << "\t" << background.size() << "\t" << ttest[4] << "\t"  << ttest[5] << endl;
		nSig ++;
		//geneSet.clear();
	}
	fin.close();
	fout.close();
	
	return nSig;
}


int find_sig_gene_sets_with_background(set<string> genes, set<string> backgroundgenes, string geneSetFile, string outputfile, double p_cutoff)
{
	unsigned n = genes.size(); // total samples
	
	unsigned N = backgroundgenes.size(); // total background
	
	ifstream fin;
	fin.open(geneSetFile.c_str());
	
	ofstream fout;
	fout.open(outputfile.c_str());

	string line;
	vector<string> flds;
	
	int nSig = 0;
	
	while(fin)
	{
		getline(fin,line);
		if(line.length() == 0)
			continue;
		// replace space with underscore
		replace(line.begin(),line.end(),' ','_');
		replace(line.begin(),line.end(),',',';');
		flds = string_split(line,"\t");
		//message("processing gene set "+to_string(i++)+": "+flds[0],"\r");
		//cout << flds[0] << "\t" << flds[1] << endl;	
		set<string> geneSet(flds.begin()+3,flds.end());
		//cout << geneSet.size() <<" genes in gene set" << endl;
		vector<string> overlap1 = set_overlap(genes,geneSet);
		vector<string> overlap2 = set_overlap(backgroundgenes,geneSet);
		
		//cout << overlap.size() << " genes overlap" << endl;
		// statistics
		unsigned k = overlap1.size(); // observd
		if (k < 2) continue;
		unsigned r = overlap2.size(); // total positive
		double p = hypergeometric_test( k,  r,  n,  N);
		if (p> p_cutoff) continue;
		nSig ++;
		double expected = double(n) * r / N;
		double fold_enrichment = k / expected;
		// geneset_name, geneset_annotation, p_value, enrichment, observed, expected, total_pos, total, all_overlap_genes
		fout << flds[0] << "\t" << flds[1] << "\t" << p << "\t" << fold_enrichment << "\t" << k << "\t" << expected << "\t" << r << "\t" << N  << "\t" << to_string(overlap1,",") << endl;
		geneSet.clear();
		overlap1.clear();
		overlap2.clear();
	}
	fin.close();
	fout.close();
	
	return nSig;
}
