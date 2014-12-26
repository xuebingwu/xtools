#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <array>        // std::array

//#include <boost/regex.hpp>

// for position matrix
#include <boost/numeric/ublas/matrix.hpp>

	
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>       /* log10 */
#include <algorithm>
#include "stat.h"
#include "utility.h"
#include "sequence.h"
#include "container.h"
#include "markov.h"
	
extern "C"{
#include "ushuffle.h"
}

using namespace std;

// number of identical bases
int sequence_similarity(string a, string b)
{
	int s = 0;
	if(a.size() != b.size()) return -1;
	for(int i=0;i<a.size();i++)
	{
		if(a[i] == b[i]) s++;
	}
	return s;
}

// calculate and write pairwise similarity matrix of
void pairwise_sequence_similarity_matrix(vector<string> seqs, string filename)
{
	ofstream out(filename.c_str());
	map<string,int> data;
	for(int i =0;i<seqs.size();i++)
	{
		for(int j=0;j<seqs.size();j++)
		{
			if (j<i) out << data[to_string(j)+","+to_string(i)] << "\t";
			else if (i == j) out << seqs[i].size() << "\t";
			else
			{
				data[to_string(i)+","+to_string(j)] = sequence_similarity(seqs[i],seqs[j]);
				out << data[to_string(i)+","+to_string(j)] << "\t";
			}
		}
		out << endl;
	}
	out.close();
}

// load weighted sequence file into two vectors: seqs and weights
// the first three columns should be: id, seq, weight (tab-delimited)
void load_weighted_sequences_to_vectors(string filename, vector<string> &seqs, vector<double> &weights,int skip/*=0*/, int cSeq/*=1*/, int cWeight/*=2*/) {
	// minimum number of columns in input
	int nCol = max(cSeq,cWeight)+1;
	
	ifstream fin;
	fin.open(filename.c_str());

    string line;
    vector<string> flds;

	// skip header lines
	int i=0;
	while(i<skip) getline(fin,line);
	
    while(fin)
    {
        getline(fin,line);
        if (line.length() == 0)
            continue;

        flds = string_split(line,"\t");
		if (flds.size() < nCol) message("too few columns in line: "+line);
		else
		{
			seqs.push_back(flds[cSeq]);
			weights.push_back(stof(flds[cWeight]));
		}
	}
	fin.close();
}

// load ranked sequences to vectors
// default: no header, first column is id, second column is sequence
// r: number of header lines to skip.
// c: column for sequence, 0-based, i.e. first column is 0
vector<string> load_ranked_sequences_to_vectors(string filename, int skip/*=0*/, int cSeq/*=1*/){
	vector<string> seqs;
	
	ifstream fin;
	fin.open(filename.c_str());

    string line;
    vector<string> flds;
	
	if (cSeq<0) // input is fasta
	{
	    while(fin)
	    {
	        getline(fin,line);
	        if (line.length() == 0)
	            continue;
			if(line[0] != '>') seqs.push_back(line);
		}
		fin.close();
		return seqs;	
	}
	
	// skip header lines
	int i=0;
	while(i<skip) getline(fin,line);

	// read sequences into a vector
    while(fin)
    {
        getline(fin,line);
        if (line.length() == 0)
            continue;

        flds = string_split(line,"\t");
		if (flds.size() < cSeq+1) message("skip lines with not enough columns: "+line);
		else seqs.push_back(flds[cSeq]);
	}
	fin.close();
	
	return seqs;
}


// find significant positional kmer from ranked sequences
// use non-parametric U test
// input: a vector of ranked sequence with the same length, only ACGT 
array<int,2> find_significant_kmer_from_ranked_sequences(
	vector<string> seqs, 
	vector<string> kmers, 
	string outfile, 
	int nTest, 
	double pCutoff/*=0.05*/, 
	double pCutoff_B/*=0.05*/, 
	int shift/*=0*/, 
	int startPos/*=0*/ )
{
    int nSeq = seqs.size();		// total number of sequences
	int lSeq = seqs[0].size(); // length of the first sequence, assume all have the same length
	
	//message(to_string(nSeq)+" sequences");
	//message(to_string(lSeq)+" length");
	
	// the number of significant positional kmers found
	array<int,2> nSig = {0,0};
		
	// a map defining IUPAC degenerate nucletodes
	map<char,string> define_iupac = define_IUPAC();

	// will output found significant positional kmers
	ofstream outstream;
	outstream.open(outfile.c_str());

	// start of kmer counting and test
	for( int i=0; i<kmers.size();i++) // for each kmer
	{
		// expand a degenerate kmer to all possible element/exact kmers
		//cout << kmers[i] << endl;
        vector<string> exp_kmers = expand_degenerate_kmer(kmers[i],define_iupac);
		int k = kmers[i].size();
		for( int pos=0; pos < lSeq-k+1; pos ++) // at each position
		{
			// ranks of sequences with this kmer at this position
			vector<int> ranks;
			// for each sequence, 
			for( int j=0;j<nSeq;j++) 
			{
				//find if any of the expanded kmer is present at position pos
				for( int n=0;n<exp_kmers.size();n++)
				{
					/* speed not affected by shift */
					size_t found = seqs[j].substr(pos,k+shift).find(exp_kmers[n]); // if found any kmer allowing shift 
					if (found!=std::string::npos)
					{
						// add this sequence's rank to sample 1
						ranks.push_back(j);
						// stop looking for the next expanded kmer in the same sequence, continue to the next sequence
						break;
					} /**/
					// another implementation, slower, linear with shift
					/*
					bool found = false;
					for(int s=0;s<=shift;s++)
					{
						if(seqs[j].substr(pos+s,k) == exp_kmers[n]) 
						{
							ranks.push_back(j);
							//cout << exp_kmers[n] <<","<< seqs[j].substr(pos+s,k) << endl;
							found = true;
							break;
						}
					}
					if (found) break;
					*/
				}
			}
			// now ranks includes ranks of all sequences containing kmer i at position pos
			// if less than 3 sequence contain this kmer, or less than 3 sequence don't have this kmer,
			// just go on to the next position
			if (ranks.size() < 3 || ranks.size() > nSeq-3 ) 
			{
				continue;
			}
			array<double,2> utest = Mann_Whitney_U_test(ranks, nSeq);
			if (utest[1] < pCutoff)
			{
				nSig[0] ++;
				double pB = min(1.0,utest[1]*nTest);
				if (pB < pCutoff_B) nSig[1]++;
	            outstream << kmers[i] << "\t" << pos-startPos << "\t" << k+shift << "\t" << utest[0] << "\t" << -log10(utest[1]) << "\t" << -log10(pB) << endl;
			}
		}
	}

    return nSig;
} // end of function

// find significant positional kmer from weighted sequences
// use t-test
// input: a vector of weighted sequence with the same length, only ACGT 
// 
array<int,2> find_significant_kmer_from_weighted_sequences(
	vector<string> seqs,
	vector<double> weights, 
	vector<string> kmers, 
	string outfile, 
	int nTest, 
	double pCutoff/*=0.05*/, 
	double pCutoff_B/*=0.05*/,
	int shift_min/*=0*/, 
	int shift_max/*=2*/,
	int startPos/*=0*/) 
{
    int nSeq = seqs.size();		// total number of sequences
	int lSeq = seqs[0].size(); // length of the first sequence, assume all have the same length
	
	//message(to_string(nSeq)+" sequences");
	//message(to_string(lSeq)+" length");
	
	// the number of significant positional kmers found
	array<int,2> nSig = {0,0};
		
	// a map defining IUPAC degenerate nucletodes
	map<char,string> define_iupac = define_IUPAC();

	// will output found significant positional kmers
	ofstream outstream;
	outstream.open(outfile.c_str());

	// start of kmer counting and test
	for( int i=0; i<kmers.size();i++) // for each kmer
	{
		// expand a degenerate kmer to all possible element/exact kmers		
        vector<string> exp_kmers = expand_degenerate_kmer(kmers[i],define_iupac);
		int k = kmers[i].size();
		for (int shift = shift_min; shift <= shift_max; shift ++)
		{
			for( int pos=0; pos < lSeq-k+1; pos ++) // at each position
			{
				//debug cout << kmers[i] << "@" << pos << endl;
			
				// whether each sequence is positive
				vector<bool> positive;
					// for each sequence, 
				for( int j=0;j<nSeq;j++) 
				{
					//debug cout << "seq " << j << endl;
				
					// initialize
					bool present = false;
					//find if any of the expanded kmer is present at position pos
					for( int n=0;n<exp_kmers.size();n++)
					{
						//debug cout << "exact kmer: " << exp_kmers[n] << endl;
						/* speed not affected by shift */
						size_t found = seqs[j].substr(pos,k+shift).find(exp_kmers[n]); // if found any kmer allowing shift 
						if (found!=std::string::npos)
						{
							// add this sequence's rank to sample 1
							//weights1.push_back(weights[j]);
							// stop looking for the next expanded kmer in the same sequence, continue to the next sequence
							present = true; 
							break;
						} /**/
					}
				
					//debug cout << "present = " << present << endl;
				
					positive.push_back(present);
				}
				vector<double> weights1,weights2;
				for( int j=0;j<nSeq;j++)
				{
					if(positive[j]) weights1.push_back(weights[j]);
					else weights2.push_back(weights[j]);
				}
				// now weights includes weights of all sequences containing kmer i at position pos
				// if less than 3 sequences with/without this kmer, just go on to the next position
				// in such cases the kmer is unlikely to be significant
				// also the t.test will not work well. it requires at least 2 data points in each sample
				if (weights1.size() < 3 || weights2.size() < 3 ) 
				{
					continue;
				}
				array<double,6> ttest = t_test(weights1,weights2,false);
				if (ttest[1] < pCutoff)
				{
					nSig[0]++;
					// Bonferoni correction
					double pB = min(1.0,ttest[1]*nTest);
					if(pB < pCutoff_B) nSig[1]++;
		            outstream << kmers[i] << "\t" << pos-startPos << "\t" << shift+k << "\t" << ttest[0] << "\t" << -log10(ttest[1]) << "\t" << -log10(pB) << "\t" << weights1.size() << "\t" << ttest[2] << "\t" << ttest[3] << "\t" << weights2.size() << "\t" << ttest[4] << "\t"  << ttest[5] << endl;
				}
			}
		}
	}

    return nSig;
} // end of function

// generate paired kmer,
vector<paired_kmer> generate_paired_kmers (
	string alphabet,
int seq1_len,
int seq2_len,
int max_dist,
int min_dist,
int max_shift,
int min_shift
){
	vector<paired_kmer> paired_kmers;
	vector<string> kmers1 = generate_kmers(seq1_len,alphabet);	
	vector<string> kmers2 = generate_kmers(seq2_len,alphabet);	
	for (int i=0;i<kmers1.size();i++)
	{
		for (int j=0;j<kmers2.size();j++)
		{
			for (int d = min_dist; d <= max_dist; d++)
			{
				for (int s= min_shift; s <= max_shift; s++)
				{
					paired_kmer a(kmers1[i],kmers2[j],d, s, 0 , 0);
					paired_kmers.push_back(a);
				}
			}
		}
	}
	return paired_kmers;
}
	

// pairwise
array<int,2> find_significant_pairs_from_weighted_sequences(
	vector<string> seqs,
	vector<double> weights, 
	vector<paired_kmer> paired_kmers, 
	string outfile, 
	int nTest, 
	double pCutoff/*=0.05*/, 
	double pCutoff_B/*=0.05*/,
	int startPos/*=0*/) 
{
    int nSeq = seqs.size();		// total number of sequences
	int lSeq = seqs[0].size(); // length of the first sequence, assume all have the same length
	
	//message(to_string(nSeq)+" sequences");
	//message(to_string(lSeq)+" length");
	
	// the number of significant positional kmers found
	array<int,2> nSig = {0,0};
		
	// will output found significant positional kmers
	ofstream outstream;
	outstream.open(outfile.c_str());

	// start of kmer counting and test
	for( int i=0; i<paired_kmers.size();i++) // for each kmer
	{
		int k = paired_kmers[i].len();
		for( int pos=0; pos < lSeq-k+1; pos ++) // at each position
		{
			vector<bool> positive;
			// whether each sequence is positive
			// for each sequence, 
			for( int j=0;j<nSeq;j++) 
			{	
				bool present = false;	
				for (int offset = 0; offset <= paired_kmers[i].shift; offset ++)
				{
					if (pos + k + offset > lSeq - 1) break;
					if (seqs[j].substr(pos+offset,paired_kmers[i].seq1.size()) == paired_kmers[i].seq1 && seqs[j].substr(pos+offset+paired_kmers[i].dist,paired_kmers[i].seq2.size()) == paired_kmers[i].seq2) 
					{
						present = true;
						break;
					}
				}
				positive.push_back(present);
			}
			vector<double> weights1,weights2;
			for( int j=0;j<nSeq;j++)
			{
				if(positive[j]) weights1.push_back(weights[j]);
				else weights2.push_back(weights[j]);
			}
			// now weights includes weights of all sequences containing kmer i at position pos
			// if less than 3 sequences with/without this kmer, just go on to the next position
			// in such cases the kmer is unlikely to be significant
			// also the t.test will not work well. it requires at least 2 data points in each sample
			if (weights1.size() < 3 || weights2.size() < 3 ) 
			{
				continue;
			}
			array<double,6> ttest = t_test(weights1,weights2,false);
			if (ttest[1] < pCutoff)
			{
				nSig[0]++;
				// Bonferoni correction
				double pB = min(1.0,ttest[1]*nTest);
				if(pB < pCutoff_B) nSig[1]++;
	            outstream << paired_kmers[i].as_string("_") << "\t" << pos-startPos << "\t" << paired_kmers[i].shift+k << "\t" << ttest[0] << "\t" << -log10(ttest[1]) << "\t" << -log10(pB) << "\t" << weights1.size() << "\t" << ttest[2] << "\t" << ttest[3] << "\t" << weights2.size() << "\t" << ttest[4] << "\t"  << ttest[5] << endl;
			}
		}
	}
    return nSig;
} // end of function



// nucleotide plot from PKA2 weighted output
void plot_PKA2_nucleotide_output(string infile, string outfile, int lSeq){
	
	string script = "# R script \n"
	"lSeq="+to_string(lSeq)+" \n"
	"data = numeric(4*lSeq)\n"
	"dict = new.env()\n"
	"dict[[\"A\"]] = 1\n"
	"dict[[\"C\"]] = 2\n"
	"dict[[\"G\"]] = 3\n"
	"dict[[\"T\"]] = 4\n"

	"x=read.table('"+infile+"',header=F)\n"
	"x = x[x[,3]==1,]\n"
	"pv = x[,5] * ((x[,4]>0)*2-1)\n"

	"pos=numeric(nrow(x))\n"
	"for (i in 1:nrow(x)){\n"
	"	pos[i] = x[i,2]*4+dict[[as.character(x[i,1])]]\n"
	"}\n"
	"data[pos] = pv\n"

	"maxy = max(abs(data))\n"
	"color=c('red','green','blue','yellow') \n"
	"label=c('A','C','G','T')\n"
	"pdf('"+outfile+"',width=10,height=5)\n"
	"bp = barplot(data,col=rep(color,4),ylab='disfavored <=== log10(p) ===> favored',ylim=c(-maxy,maxy)*1.3,border=NA)\n"
		
	"text(bp,data,labels=rep(label,lSeq),pos= (data>0)*2+1,offset=0.2,col='gray') \n"

	"for(i in 1:lSeq){\n"
	"	abline(v=bp[i*4+1]/2+bp[i*4]/2,col='gray')\n"
	"	text(bp[i*4-2]/2+bp[i*4-1]/2,-maxy*1.2,labels=i)\n"
	"}\n"

	"legend('topright',legend=label,col=color,pch=15,bty='n')\n"

	"dev.off() \n";
	
	R_run(script);
}

positional_kmer::positional_kmer(string seq1, int pos1, int size1, double weight1, int group1){
	seq = seq1;
	pos = pos1;
	size = size1;
	weight = weight1;
	group = group1;
}

positional_kmer::positional_kmer() 
{
  // allocate variables
	seq = "";
	pos = -1;
	size = -1;
	weight = -1.0;
	group = -1;
}

positional_kmer::positional_kmer(const positional_kmer &a) 
{
  // allocate variables
  positional_kmer();
  // copy values
  operator = (a);
}

bool positional_kmer::equals(positional_kmer a)
{
	return seq == a.seq && pos == a.pos && size == a.size;
}

bool positional_kmer::is_part_of(positional_kmer a)
{
	// only if a is part of this at the same position
	if (pos >= a.pos && ( (pos + seq.size()) <= (a.pos + a.seq.size()) ) && seq.size() != a.seq.size())
		if (a.seq.substr(pos - a.pos, seq.size()) == seq)
			return true;
	return false;
}

const positional_kmer &positional_kmer::operator = (const positional_kmer &a)
{
	seq = a.seq;
	pos = a.pos;
	size = a.size;
	weight = a.weight;
	group = a.group;
	return *this;
}

string positional_kmer::as_string(string del/*="_"*/)
{
	return seq+del+to_string(pos)+del+to_string(size)+del+to_string(weight)+del+to_string(group);
}


paired_kmer::paired_kmer(string seqa, string seqb, int distance, int offset, int position, double score)
{
	seq1 = seqa;
	seq2 = seqb;
	dist = distance;
	pos = position;
	shift=offset;
	weight=score;
}

paired_kmer::paired_kmer() 
{
  // allocate variables
	seq1 = "A";
	seq2 = "A";
	dist = 1;
	shift=0;
	pos=0;
	weight=0.0;
}

paired_kmer::paired_kmer(const paired_kmer &a) 
{
  // allocate variables
  paired_kmer();
  // copy values
  operator = (a);
}

// everything equal except weight
bool paired_kmer::equals(paired_kmer a, bool identical_pos, bool identical_shift)
{
	if (identical_pos && pos != a.pos) return false;
	if (identical_shift && shift != a.shift) return false;
	return seq1 == a.seq1 && dist == a.dist && seq2 == a.seq2;
}


const paired_kmer &paired_kmer::operator = (const paired_kmer &a)
{
	seq1 = a.seq1;
	seq2 = a.seq2;
	dist = a.dist;
	pos = a.pos;
	shift = a.shift;
	weight = a. weight;
	return *this;
}

string paired_kmer::as_string(string del/*="_"*/, bool add_shift/*=false*/,bool add_pos/*=false*/, bool add_weight/*=false*/)
{
	string res = seq1+del+to_string(dist)+del+seq2;
	if (add_shift) res += del + to_string(shift);
	if (add_pos) res += del + to_string(pos);
	if (add_weight) res += del + to_string(weight);
	return res;
}

int paired_kmer::len()
{
	return dist+int(seq2.size());
}

int paired_kmer::gap()
{
	return dist - int(seq1.size());
}

// note pCutoff and pCutoff_B here are -log10 transformed
// note input should be sorted by weight then by size
vector<positional_kmer> build_model_from_PKA2_output(string filename, double pCutoff, double pCutoff_B)
{
	int skip = 0;
	int cKmer = 0;
	int cStart = 1;
	int cSize = 2;
	int cStat = 3;
	int cWeight = 4; // p, log10 
	int cpB = 5; // corrected p, -log10
		
	vector<string> kmers;
	vector<int> starts, sizes;
	vector<double> weights;
	string kmer;
		
	map<char,string> iupac = define_IUPAC();
	
	set <string> ranked_kmer_ids;
	vector<positional_kmer> ranked_kmers;
	
	ifstream fin;
	fin.open(filename.c_str());
	
	string line;
	vector<string> flds;
	// skip header lines
	int i=0;
	while(fin && i++ < skip) getline(fin,line);

	int group = 0;
    while(fin)
    {
        getline(fin,line);
        if (line.length() == 0)
            continue;

        flds = string_split(line,"\t");
		
		double weight = stof(flds[cWeight]);
		double logpb = stof(flds[cpB]);
		
		if (weight < pCutoff || logpb < pCutoff_B) continue;
		
		if(stof(flds[cStat])<0) weight = -weight;
		
		vector<string> tmp = expand_degenerate_kmer(flds[cKmer],iupac);
		for(int i=0;i<tmp.size();i++)
		{
			string id = tmp[i]+","+flds[cStart]+","+flds[cSize];
			//debug 
			//cout << id << endl;
			// if new kmer instead of the same kmer with lower weight (can happen when degenerate base is used)
			if(ranked_kmer_ids.find(id) == ranked_kmer_ids.end()) 
			{
				ranked_kmer_ids.insert(id);
				positional_kmer x(tmp[i],stoi(flds[cStart]),stoi(flds[cSize]),weight,0);
				/*
				// tried to remove longer but weaker motif, or ignore shorter when stronger present, worse in training
				// try test with crossvalidation
				bool found = false;
				for(int j=0;j<ranked_kmers.size();j++)
				{
					if (ranked_kmers[j].is_part_of(x)) 
					{
						found = true;
						// debug cout << "ignore " << x.as_string() << " given " << ranked_kmers[j].as_string() << endl;
						break;
					}
				}
				if (found) continue;
				// if shorter but weaker, assign the same group number
				found = false;
				for(int j=0;j<ranked_kmers.size();j++)
				{					
					if (x.is_part_of(ranked_kmers[j])) 
					{
						found = true;
						x.group = ranked_kmers[j].group;
						// debug cout << "add " << x.as_string() << " to group " << to_string(x.group) << endl;
						break;
					}
				}
				if (!found) 
				{
					group++;
					x.group = group;
				}
				*/
				ranked_kmers.push_back(x);
			}
		}
	}
	return ranked_kmers;
}

// no degenerate nucleotides
vector<paired_kmer> build_paired_kmer_model(string filename)
{
	int skip = 0; // no header
	int cKmer = 0; // first column is kmer
	int cStart = 1; // second column is start 
	int cStat = 3; // statistics, tells you sign of signficance
	int cWeight = 4; // -log10(p)

	vector<paired_kmer> paired_kmers;
	
	ifstream fin;
	fin.open(filename.c_str());
	
	string line;
	vector<string> flds;
	
	// skip header lines
	int i=0;
	while(fin && i++ < skip) getline(fin,line);

    while(fin)
    {
        getline(fin,line);
        if (line.length() == 0)
            continue;

        flds = string_split(line,"\t");
		
		string kmer = flds[cKmer];
		int pos = stoi(flds[cStart]);
		double score = stof(flds[cWeight]);
		if(stof(flds[cStat]) < 0 ) score = -score;
		
		vector<string> tmp = string_split(kmer,"_");
		string seq1 = tmp[0];
		string seq2 = tmp[2];
		int dist = stoi(tmp[1]);
		int shift = stoi(tmp[3]);
		
		paired_kmer x(seq1, seq2, dist, shift, pos, score);
		paired_kmers.push_back(x);
	}
	fin.close();
	return paired_kmers;
}

// no degenerate nucleotides
vector<positional_kmer> build_model_from_PKA_output(string filename)
{
	int skip = 1;
	int cKmer = 0;
	int cStart = 2;
	int cWeight = 7; // z-score

	vector<positional_kmer> ranked_kmers;
	
	ifstream fin;
	fin.open(filename.c_str());
	
	string line;
	vector<string> flds;
	// skip header lines
	int i=0;
	while(fin && i++ < skip) getline(fin,line);

    while(fin)
    {
        getline(fin,line);
        if (line.length() == 0)
            continue;

        flds = string_split(line,"\t");
							
		positional_kmer x(flds[cKmer],stoi(flds[cStart]),flds[cKmer].size(),stof(flds[cWeight]),0);
		ranked_kmers.push_back(x);
	}
	fin.close();
	return ranked_kmers;
}


void save_model_to_file(vector<positional_kmer> ranked_kmers, string filename)
{
	ofstream out;
	out.open(filename.c_str());
	
	for( int i=0;i<ranked_kmers.size();i++)
		out << ranked_kmers[i].as_string("\t") << endl;
	
	out.close();
}

vector<positional_kmer> load_model_from_file(string filename)
{
	vector<positional_kmer> ranked_kmers;
	
	ifstream fin;
	fin.open(filename.c_str());
	
	string line;
	vector<string> flds;
	
    while(fin)
    {
        getline(fin,line);
        if (line.length() == 0)
            continue;

        flds = string_split(line,"\t");
		positional_kmer x(flds[0],stoi(flds[1]),stoi(flds[2]),stof(flds[3]),stoi(flds[4]));
		ranked_kmers.push_back(x);
	}
	return ranked_kmers;
}


// ignore group information
double score_sequence_using_PKA_model(vector<positional_kmer> ranked_kmers,  string seq)
{
	double score = 0;
	for( int i=0;i<ranked_kmers.size();i++)
	{
		//cout << ranked_kmers[i].as_string() ;
		//cout << "\t" << seq.substr(ranked_kmers[i].pos,ranked_kmers[i].size) ;
		size_t found = seq.substr(ranked_kmers[i].pos,ranked_kmers[i].size).find(ranked_kmers[i].seq);
		if (found!=std::string::npos) 
		{
			score += ranked_kmers[i].weight;
			//cout << "\t" << score;
		}
		//cout << endl;
	}
	return score;
}


double score_sequence_using_paired_kmer_model(vector<paired_kmer> model,  string seq)
{
	double score = 0;
	int lSeq = seq.size();
	for( int i=0;i<model.size();i++)
	{
		bool present = false;
		int l = model[i].len();	
		
		for (int offset = 0; offset <= model[i].shift; offset ++)
		{
			if (model[i].pos + l + offset > lSeq - 1) break;
			if (seq.substr(model[i].pos+offset,model[i].seq1.size()) == model[i].seq1 && seq.substr(model[i].pos+offset+model[i].dist,model[i].seq2.size()) == model[i].seq2) 
			{
				present = true;
				break;
			}
		}
		if (present) 
		{
			score += model[i].weight;
			//cout << "\t" << score;
		}
		//cout << endl;
	}
	return score;
}

// use group information, worse
double score_sequence_using_PKA_model_use_group(vector<positional_kmer> ranked_kmers,  string seq)
{
	double score = 0;
	set<int> scored;
	for( int i=0;i<ranked_kmers.size();i++)
	{
		//cout << ranked_kmers[i].as_string() ;
		//cout << "\t" << seq.substr(ranked_kmers[i].pos,ranked_kmers[i].size) ;
		size_t found = seq.substr(ranked_kmers[i].pos,ranked_kmers[i].size).find(ranked_kmers[i].seq);
		if (found!=std::string::npos) 
		{
			if(scored.find(ranked_kmers[i].group) == scored.end())  // if the group was not used before
			{
				score += ranked_kmers[i].weight;
				scored.insert(ranked_kmers[i].group);
				//cout << "\t" << score;
			}
		}
		//cout << endl;
	}
	return score;
}

// for both PKA and PKA2
void score_fasta_using_PKA_model(string seqfile, string outputfile, vector<positional_kmer> ranked_kmers)
{
    ifstream fin(seqfile.c_str());
	ofstream fout(outputfile.c_str());
  
    string name,seq;
    while(fin.good())
    {
		ReadOneSeqFromFasta(fin,name,seq);
     	double score = score_sequence_using_PKA_model(ranked_kmers, seq);
  		fout << name << "\t" << seq << "\t" << score << endl;		
    }
    fin.close();
	fout.close();
}

// sequence in column col, 1 based
void score_tabular_using_PKA_model(string tabfile, int col, string outputfile, vector<positional_kmer> ranked_kmers)
{
	col = col - 1;
	
    ifstream fin(tabfile.c_str());
	ofstream fout(outputfile.c_str());
  
    string line;
	vector<string> flds;
	
    while(fin.good())
    {
		getline(fin,line);
		if(line.length() == 0) continue;
		flds = string_split(line,"\t");
     	double score = score_sequence_using_PKA_model(ranked_kmers, flds[col]);
  		fout << line << "\t" << score << endl;		
    }
    fin.close();
	fout.close();
}


// count and write feature matrix, format
// seqid, label (i.e. 1/0)
void save_feature_matrix(map<string,string> seqs, vector<string> kmers, string outfile, string label, bool append/*=false*/, int shift/*=0*/)
{
    int nSeq = seqs.size();		// total number of sequences
	int lSeq = seqs.begin()->second.size(); // length of the first sequence, assume all have the same length
	
	//message(to_string(nSeq)+" sequences");
	//message(to_string(lSeq)+" length");
	
	// a map defining IUPAC degenerate nucletodes
	map<char,string> define_iupac = define_IUPAC();

	ofstream outstream;
	if (append) 
	{
		outstream.open(outfile.c_str(),ios::app);
	} 
	else
	{		
		outstream.open(outfile.c_str());
		// output header
		outstream << "SeqID\t" << label;
		for( int i =0; i< kmers.size();i++)
		{
			for(int j=0;j< lSeq-kmers[i].size()+1; j++)
			{
				outstream << "\t" << kmers[i] << "_" << j << "_" << shift ;
			}
		}
		outstream << endl;
	}


	// start of kmer counting and test
	for(map<string,string>::iterator it=seqs.begin();it!=seqs.end();it++) // each sequence is a line
	{ 
		outstream << it->first << "\t" << label; 
		for( int i=0; i<kmers.size();i++) // for each kmer
		{
			// expand a degenerate kmer to all possible element/exact kmers		
	        vector<string> exp_kmers = expand_degenerate_kmer(kmers[i],define_iupac);
			int k = kmers[i].size();
			for( int pos=0; pos < lSeq-k+1; pos ++) // at each position
			{
				//debug cout << kmers[i] << "@" << pos << endl;
				// initialize
				int present = 0;
				//find if any of the expanded kmer is present at position pos
				for( int n=0;n<exp_kmers.size();n++)
				{
					//debug cout << "exact kmer: " << exp_kmers[n] << endl;
					/* speed not affected by shift */
					size_t found = it->second.substr(pos,k+shift).find(exp_kmers[n]); // if found any kmer allowing shift 
					if (found!=std::string::npos)
					{
						// add this sequence's rank to sample 1
						//weights1.push_back(weights[j]);
						// stop looking for the next expanded kmer in the same sequence, continue to the next sequence
						present = 1; 
						break;
					} /**/
				}
				outstream << "\t" << present;
			}
		}
		outstream << endl;
	}
} // end of function

// read PKA output into two vectors
void read_significant_positional_kmer_from_file(string inputfile, vector<string> &kmers, vector<int> &positions)
{
	ifstream fin;
	fin.open(inputfile.c_str());
	
	string line;
	vector<string> flds;
	// skip header lines
	getline(fin,line);

    while(fin)
    {
        getline(fin,line);
        if (line.length() == 0)
            continue;

        flds = string_split(line,"\t");
		
		kmers.push_back(flds[0]);
		positions.push_back(stoi(flds[2]));
	}
	fin.close();
}

// read PKA2 output into two vectors
void read_significant_positional_kmer_from_PKA2_output(string inputfile, vector<string> &kmers, vector<int> &positions)
{
	ifstream fin;
	fin.open(inputfile.c_str());
	
	string line;
	vector<string> flds;
	// skip header lines
	getline(fin,line);

    while(fin)
    {
        getline(fin,line);
        if (line.length() == 0)
            continue;

        flds = string_split(line,"\t");
		
		kmers.push_back(flds[0]);
		positions.push_back(stoi(flds[1]));
	}
	fin.close();
}

// vector<string> to map<string,string>
map<string, string> seq_vector2map(vector<string> seqs)
{
	map<string,string> res;
	for(int i=0;i<seqs.size();i++)
		res[to_string(i)] = seqs[i];
	return res; 
}

void significant_feature_matrix(map<string,string> seqs, vector<string> kmers, vector<int> positions, string outfile, string label, bool append/*=false*/, int shift/*=0*/)
{
    int nSeq = seqs.size();		// total number of sequences
	int lSeq = seqs.begin()->second.size(); // length of the first sequence, assume all have the same length
	
	//message(to_string(nSeq)+" sequences");
	//message(to_string(lSeq)+" length");
	
	// a map defining IUPAC degenerate nucletodes
	map<char,string> define_iupac = define_IUPAC();

	ofstream outstream;
	if (append) 
	{
		outstream.open(outfile.c_str(),ios::app);
	} 
	else
	{		
		outstream.open(outfile.c_str());
		// output header
		outstream << "SeqID\tLabel";
		for( int i =0; i< kmers.size();i++)
		{
			outstream << "\t" << kmers[i] << "_" << positions[i] << "_" << shift ;
		}
		outstream << endl;
	}


	// start of kmer counting and test
	for(map<string,string>::iterator it=seqs.begin();it!=seqs.end();it++) // each sequence is a line
	{ 
		outstream << it->first << "\t" << label; 
		for( int i=0; i<kmers.size();i++) // for each kmer
		{
			// expand a degenerate kmer to all possible element/exact kmers		
	        vector<string> exp_kmers = expand_degenerate_kmer(kmers[i],define_iupac);
			int k = kmers[i].size();
			
			int present = 0;
			//find if any of the expanded kmer is present at position pos
			for( int n=0;n<exp_kmers.size();n++)
			{
				//debug cout << "exact kmer: " << exp_kmers[n] << endl;
				/* speed not affected by shift */
				size_t found = it->second.substr(positions[i],k+shift).find(exp_kmers[n]); // if found any kmer allowing shift 
				if (found!=std::string::npos)
				{
					// add this sequence's rank to sample 1
					//weights1.push_back(weights[j]);
					// stop looking for the next expanded kmer in the same sequence, continue to the next sequence
					present = 1; 
					break;
				} /**/
			}
			outstream << "\t" << present;
		}
		outstream << endl;
	}
	outstream.close();
} // end of function


// should work for ranked kmer with degenerate nucleotides
void significant_feature_matrix_PKA2(vector<string> seqs, vector<double> weights, vector<positional_kmer> ranked_kmers, string outfile)
{
    int nSeq = seqs.size();		// total number of sequences
	int lSeq = seqs[0].size(); // length of the first sequence, assume all have the same length
	
	//message(to_string(nSeq)+" sequences");
	//message(to_string(lSeq)+" length");
	
	// a map defining IUPAC degenerate nucletodes
	map<char,string> define_iupac = define_IUPAC();

	ofstream outstream;
	outstream.open(outfile.c_str());
	// output header
	outstream << "SeqID\tScore";
	for( int i =0; i< ranked_kmers.size();i++)
	{
		outstream << "\t" << ranked_kmers[i].as_string() ;
	}
	outstream << endl;


	// start of kmer counting and test
	for(int j=0;j<seqs.size();j++) // each sequence is a line
	{ 
		outstream << "Seq-"<< j << "\t" << weights[j]; 
		for( int i=0; i<ranked_kmers.size();i++) // for each kmer
		{
			// expand a degenerate kmer to all possible element/exact kmers		
	        vector<string> exp_kmers = expand_degenerate_kmer(ranked_kmers[i].seq,define_iupac);
			int k = ranked_kmers[i].seq.size();
			
			int present = 0;
			//find if any of the expanded kmer is present at position pos
			for( int n=0;n<exp_kmers.size();n++)
			{
				//debug cout << "exact kmer: " << exp_kmers[n] << endl;
				/* speed not affected by shift */
				size_t found = seqs[j].substr(ranked_kmers[i].pos,ranked_kmers[i].size).find(exp_kmers[n]); // if found any kmer allowing shift 
				if (found!=std::string::npos)
				{
					// add this sequence's rank to sample 1
					//weights1.push_back(weights[j]);
					// stop looking for the next expanded kmer in the same sequence, continue to the next sequence
					present = 1; 
					break;
				} /**/
			}
			outstream << "\t" << present;
		}
		outstream << endl;
	}
	outstream.close();
} // end of function

// PKA: remove overlapping motifs
// input: all significant motifs
int non_overlapping_sig_motifs(string inputfile, string outputfile)
{
    // load file, store pos, kmer, z-score, line number
    // from top
    vector<int> position;
    vector<int> length;

    vector<bool> removed; // indicator whether removed or not

    ifstream fin;
    fin.open(inputfile.c_str());

    string line;
    vector<string> flds;

    // skip header;
    getline(fin,line);

    while(fin)
    {  
        getline(fin,line);
        if (line.length() == 0)
            continue;

        flds = string_split(line,"\t");
        length.push_back(flds[0].size());
        position.push_back(stoi(flds[2]));
        removed.push_back(false);
    }
    fin.close();
    
    //
    int current_best = 0; // position of currently best kmer
    int next_best = 1;
    while(next_best>0) // if next best can be found
    {
        // keep the current best motif
        removed[current_best] = false;
        next_best = -1;
        // remove from remaining set thsoe overlapping with the top motif
        int start = position[current_best];
        int end = start + length[current_best] - 1;
        for(int i = current_best+1;i<removed.size();i++)
        {
            if (removed[i] == false)
            {
                int start2 = position[i];
                int end2 = start2+ length[i] - 1;
                if (start <= start2 && end >= start2 || start <= end2 && end >= end2) // overlap
                {
                    removed[i] = true;        
                }
                else if (next_best < 0) next_best = i; 
            }
        }
        current_best = next_best;
    }
    
    // output selected lines
    ofstream fout;
    fout.open(outputfile.c_str());
    fin.open(inputfile.c_str());
    getline(fin,line);
    int n = -1;
    int total_selected = 0;
    while(fin)
    {
        getline(fin,line);
        if(line.length() == 0) continue;
        n++;
        if (removed[n] == false) 
        {
            fout << line << endl;
            total_selected++;
        }
    }
    fin.close();
    fout.close();
    return total_selected;
}


// count substring in a map of sequences, i.e. from fasta
// used in generating markov model
int countSubstringInSeqs(map<string,string>seqs, string sub)
{
    int count = 0;
    for(map<string,string>::iterator it=seqs.begin();it!=seqs.end();it++)
    {
        count += findall(it->second,sub).size();
    } 
    return count;
}



set<int> findall(string seq, string motif)
{
    // find all match positions of motif in seq
    // allowing overlap 
    set<int> allpos; // to store all found positions
    size_t pos = 0; //start at pos
    pos = seq.find(motif); // the first match
    while(pos != string::npos)
    {
        allpos.insert(pos);
        pos = seq.find(motif,pos+1); // start search for the next match from pos + 1
    }
    return allpos;
}


// dict for IUPAC degenerate nucleotides
map<char,string> define_IUPAC()
{
    map<char,string> m;

    m['A']="A";
    m['C']="C";
    m['G']="G";
    m['T']="T";
    m['U']="T";
    m['R']="AG";
    m['Y']="CT";
    m['M']="AC";
    m['K']="GT";
    m['W']="AT";
    m['S']="CG";
    m['B']="CGT";
    m['D']="AGT";
    m['H']="ACT";
    m['V']="ACG";
    m['N']="ACGT";

    return m;
}

/*
// for output only
map<char,string> interpret_IUPAC()
{
    map<char,string> m;

    m['A']="A";
    m['C']="C";
    m['G']="G";
    m['T']="T";
    m['U']="T";
    m['R']="AG";
    m['Y']="CT";
    m['M']="AC";
    m['K']="GT";
    m['W']="AT";
    m['S']="CG";
    m['B']="a";
    m['D']="c";
    m['H']="g";
    m['V']="t";
    m['N']="N";

    return m;
}
*/

// generate all possible combinations of letters in alphabet of fixed length k, i.e. kmer
vector<string> generate_kmers(int k, string alphabet)
{
    vector<string> kmers;
    
    // first position, each letter in alphabet
    for(int i = 0; i< alphabet.size(); i++)
    {
        string s(1,alphabet[i]);
        kmers.push_back(s);
    }

    // for the rest k-1 positions, fill all possible letters 
    for(int i = 0; i< k-1;i++)
    {
        int n = kmers.size();
        // for each current kmer, append all possible letters
        for(int j=0; j< n; j++)
        {
            string tmp = kmers[j];
            kmers[j] = kmers[j]+alphabet[0];
            for( int m =0;m< alphabet.size()-1;m++)
            {
                string s(1,alphabet[m+1]);
                kmers.push_back(tmp+s);
            }
        }
    }
    return kmers;
}

// degenerate kmers based on IUPAC, DNA only
// remove those with terminal N, which is not a kmer, but (k-i)mer
vector<string> degenerate_kmer(int k, string alphabet/*ACGTRYMKWSBDHVN*/)
{
    vector<string> tmp = generate_kmers(k,alphabet);
    // remove those with terminal N
    vector<string> dkmers;
    for( int i = 0; i < tmp.size(); i++)
    {  
        if (tmp[i][0] != 'N' && tmp[i][k-1] != 'N') dkmers.push_back(tmp[i]);
    }
    return dkmers;
}

// convert one degenerate kmer to regular expression, e.g. ARG -> A[AG]G
string degenerate_kmer_to_regex(string kmer,map<char,string> iupac)
{
    string str;
    for( int i  = 0; i< kmer.size();i++)
    {
        if (iupac[kmer[i]].size() > 1) str = str + "["+iupac[kmer[i]]+"]";
        else str = str + iupac[kmer[i]];
    }    
    return str;
}

// expand a degenerate kmer to all possible exact kmers
vector<string> expand_degenerate_kmer(string seq, map<char,string> iupac)
{
    vector<string> res (1,"");
    vector<string> tmp;

    for( int i=0;i<seq.size();i++)
    {  
        string exp = iupac[seq[i]];
        int L = res.size();
        tmp.clear();
        for( int j=0;j<L;j++)
        {  
            for( int k=0;k<exp.size();k++)
            {  
                tmp.push_back(res[j]+exp[k]);
            }
        }
        res = tmp;
    }
    return res;
}

// implant a motif to a set of sequences
// motif can contain degenerate nucleotides
void implant_motif(map<string,string> &seqs, int position, string motif, double fraction)
{
    int k = motif.size();
    int nSeq = seqs.size() * fraction; // number of sequences to be implanted
    // expand degenerate nucleotides
    vector<string> kmers = expand_degenerate_kmer(motif, define_IUPAC());

    // implant into the first nSeq
    for(map<string,string>::iterator it=seqs.begin();it!=seqs.end();it++)
    {
        it->second.replace(position,k,kmers[rand()%kmers.size()]);
        nSeq--;
        if (nSeq==0) return ;     
    }
}


// write pwm in meme format
void write_pwm_in_meme_format(boost::numeric::ublas::matrix<double> pwm, string motifname, string filename)
{
    ofstream out;
    out.open(filename.c_str());
    out << "MEME version 4.4\n\n";
    out << "ALPHABET= ACGT\n\n";
    out << "strands: + -\n\n";
    out << "Background letter frequencies\n";
    out << "A 0.25 C 0.25 G 0.25 T 0.25" << "\n\n";
    out << "MOTIF "<<motifname<<"\n\n";
    out << "letter-probability matrix: alength= 4 w= "<<pwm.size2()<<" nsites= 1000 E= 0"<<endl;
    for (int i = 0; i < pwm.size2(); ++ i)
    {
        for (int j = 0; j < pwm.size1(); ++ j)
        {
            out << pwm(j,i) << "\t";
        }
        out << endl;
    }
    out << endl;
    out.close();
}


// build position weight matrix for a positional kmer, can include flanking sequence
boost::numeric::ublas::matrix<double> create_position_weight_matrix_from_kmer(map<string,string> seqs, string kmer, int position, map<char,string> iupac, int d, int startPos)
{
    //position: 0-based
    // need to back-calculate position = position + startPos
    position = position + startPos; 
    
    //cout << kmer << "\t" << position;     //debug

    // length of sequences
    int seqLen = seqs.begin()->second.size();
    int k = kmer.size();

    // begin and end of region in the sequence
    int start = max(0,position-d);
    int end = min(seqLen-1,position + k + d - 1);

    // initialize an empty matrix
    int L = end - start + 1;
    boost::numeric::ublas::matrix<double> pwm (4,L);
    for (int i = 0; i < pwm.size1(); ++ i)
        for (int j = 0; j < pwm.size2(); ++ j)
            pwm(i,j) = 0.0;
  
    //cout << "matrix initiaze ok" << endl;

    // nucleotide to position
    map<char,int> letter2pos;
    letter2pos['A'] = 0;
    letter2pos['C'] = 1;
    letter2pos['G'] = 2;
    letter2pos['T'] = 3;
 
    // expand kmer
    vector<string> dkmersexp = expand_degenerate_kmer(kmer,iupac); 

    // find sequence with motif match and update matrix       
    int total_positive_seq = 0; // total number of sequences contain the motif
    for(map<string,string>::iterator it=seqs.begin();it!=seqs.end();it++)
    {
        for( int i =0;i<dkmersexp.size();i++)
        {
            if (it->second.substr(position,k) == dkmersexp[i])
            {
                // update matrix
                for( int j=start;j<=end;j++)
                {
                    pwm(letter2pos[it->second[j]],j-start) += 1;
                }                
                total_positive_seq ++;
                break;
            }
        }
    }

    //cout << "matrix done" << endl;


    // divide
    pwm /= total_positive_seq;

    //debug
    //cout << "after division" << endl;
    //cout << kmer << "," << total_positive_seq << endl;
    //print_matrix(pwm);
    //    cout << kmer << "\t" << total_positive_seq << endl;

    return pwm;
}

//PKA : create logo for a single kmer
void create_logo_for_kmer(map<string,string> seqs, string kmer, int position, map<char,string> iupac, int d, int startPos,string output)
{
    boost::numeric::ublas::matrix<double> pwm;
    pwm = create_position_weight_matrix_from_kmer(seqs, kmer, position,iupac,d,startPos);
    
    string meme_filename = output+"_"+kmer+"_"+to_string(position)+".meme"; // file name use 1-based position
    write_pwm_in_meme_format(pwm,kmer,meme_filename);

    string cmd = "ceqlogo -i "+meme_filename+" -o "+meme_filename+".eps -t "+meme_filename+" -x '' -b 2 -c 1 -w 10 -h 5";  
    system(cmd.c_str());
    cmd = "ps2pdf -dEPSCrop "+meme_filename+".eps "+meme_filename+".pdf";
    system(cmd.c_str());
}

// create logos for top 1 kmer passing pCutoff_B
void create_logo_for_topN_sig_kmer_per_position(map<string,string> seqs, string filename,int d, double pCutoff_B,int startPos,string output)
{
    ifstream fin;
    fin.open(filename.c_str());
    string line;
    vector<string> flds;

    map<char,string> iupac = define_IUPAC();

    // skip header;
    getline(fin,line);

    while(fin)
    {  
        getline(fin,line);
        if (line.length() == 0)
            continue;
        flds = string_split(line,"\t");
        //cout << flds[9] <<","<< pCutoff_B << endl;
        if(stod(flds[9]) < pCutoff_B)
        {
            create_logo_for_kmer(seqs,flds[0],stoi(flds[2]),iupac,d,startPos,output); // change back to 0-based coordinates
        }
    }    
    fin.close();
}

/*
vector<int> count_one_kmer_in_all_seqs_regex(map<string,string> seqs, boost::regex pattern, int k, int seq_len_minus_k_plus_1, int k_plus_shift)
{
    vector<int> counts;
    for(int j=0;j<seq_len_minus_k_plus_1;j++) counts.push_back(0);
    for(map<string,string>::iterator it=seqs.begin();it!=seqs.end();it++)
    {   
        string seq = (*it).second+"@@@@@@@@@@@@@@";
         // for each position + shift window
        for( int m=0;m<seq_len_minus_k_plus_1;m++)
        {   
            if (boost::regex_search(seq.substr(m,k_plus_shift),pattern)) counts[m]++;
        }
    }
    return counts;
}
*/

/*
// find significant degenerate kmers with shift
array<int,2> find_significant_degenerate_kmer_with_shift(vector<string> kmers, map<string,string> seqs1, map<string,string> seqs2, int shift, double pCutoff, double pCutoff_B,float pseudo, int startPos, int nTest, string output)
{
    map<char,string> define_iupac = define_IUPAC();
    //map<char,string> interpret_iupac = interpret_IUPAC();

    int seq_len = seqs1.begin()->second.size();

    int nSeq1 = seqs1.size();
    int nSeq2 = seqs2.size();

    // significant kmers
    array<int,2> nSig  = {0,0};

    ofstream fout;
    fout.open(output.c_str());
    fout.precision(3);

    for( int i=0;i<kmers.size();i++)
    { 
        // each kmer could be of different size
        int k = kmers[i].size();
        // convert each kmer to pattern
        string regex_format = degenerate_kmer_to_regex(kmers[i],define_iupac);
        boost::regex pattern(regex_format); 
        vector<int> counts1 = count_one_kmer_in_all_seqs_regex(seqs1,pattern, k, seq_len - k + 1, k + shift);
        vector<int> counts2 = count_one_kmer_in_all_seqs_regex(seqs2,pattern, k, seq_len - k + 1, k + shift);

        //
        int total_counts1 = sum(counts1);

        for( int m=0; m < seq_len - k + 1;m++)
        {   
            //if (counts1[m] == 0) continue;
            double f2 = float(counts2[m]+pseudo)/nSeq2;
            double p = binom_test(nSeq1,counts1[m],f2);
            if (p < pCutoff)
            {
                nSig[0]++;
                double f1 = float(counts1[m])/nSeq1;
                double expected = nSeq1 * f2;
                double z = (counts1[m] - expected) / sqrt(expected*(1-f2));
                double local_r = double(counts1[m]) / (total_counts1 - counts1[m]) * (counts1.size()-1);
                fout  << kmers[i] << "\t" << regex_format << "\t"  << m-startPos <<"\t"<< p << "\t"<< min(1.0,p*nTest)<< "\t"<< z << "\t" << f1/f2 << "\t" << f1 << "\t" << f2 << "\t" << local_r << endl;
                if (p*nTest<pCutoff) nSig[1]++;
            }
        }
    }

    return nSig;
}
*/

// doesn't work with shift+degenerate
// for each kmer, count its freq at each position(start), allowing shifts (to the right)
map<string,vector<int> > count_all_kmer_in_seqs(vector<string> kmers, map<string,string> seqs, int shift)
{
    // kmers: vector of all kmers
    // 
    map<string, vector<int> > data;

    int i,j,k;
	int m;
    string name,seq;
    map<string,string>::iterator it;
    vector<int> counts;
    set<int> found;

    int seq_len = seqs.begin()->second.size();

    for (i=0;i<kmers.size();i++)
    {  
        // kmers could be of different size
        k = kmers[i].size();
        // initialize counts
        for(j=0;j<seq_len-k+1;j++) counts.push_back(0);
        // loop through each sequencs, count this kmer
        for(it=seqs.begin();it!=seqs.end();it++)
        {   
            name = (*it).first;
            seq = (*it).second;

            // find all matches
            found = findall(seq,kmers[i]);
            // add shifts
            vector<int> to_be_inserted;
            for (set<int>::iterator it=found.begin(); it!=found.end(); ++it)
            {   
                for (m=1;m<shift+1;m++)
                {   
                    if( *it > m) to_be_inserted.push_back(*it-m);
                }
            }
            for (m=0;m<to_be_inserted.size();m++)  found.insert(to_be_inserted[m]);
            for (set<int>::iterator it=found.begin(); it!=found.end(); ++it) counts[*it]++;
        }
        data[kmers[i]] = counts;
        counts.clear();
    }
    return data;
}   

//PKA
void print_kmer_positional_profile(map<string,vector<int> > data)
{   
    for (map<string,vector<int> >::iterator it=data.begin();it!=data.end();it++)
    {
        cout << (*it).first;
        for( int i=0;i<(*it).second.size();i++) cout << ',' << (*it).second[i];
        cout << endl;
    }
}


map<string,vector<int> > degenerate_kmer_counts(vector<string> dkmers,map<string,vector<int> > data, map<char,string> define_iupac)
{
    map<string,vector<int> > new_data;
    vector<int> tmp;
   // for each degenerate kmer, combine element kmer counts
    for( int i=0;i<dkmers.size();i++)
    {  
        vector<string> dkmersexp = expand_degenerate_kmer(dkmers[i],define_iupac);
        new_data[dkmers[i]] = data[dkmersexp[0]];
//        cout << i << "\t" << dkmers[i] << "\t" << dkmersexp.size() << endl;
        for( int j=1;j<dkmersexp.size();j++)
        {
  //          cout << i << "," << j << endl;
            tmp = sum(new_data[dkmers[i]],data[dkmersexp[j]]);
            new_data[dkmers[i]] = tmp;
    //        cout << "done" << endl;
        }
    }
    return new_data;
}

array<int,2> find_significant_kmer_from_one_seq_set(
	map<string,string>seqs1, 
	map<string,double> probs_kmer,vector<string>kmers, 
	vector<string> dkmers, 
	int shift,
	bool degenerate,
	double pCutoff, 
	double pCutoff_B, 
	int startPos,
	int nTest, 
	string outfile, 
	string output_count_file)
{
    int nSeq1 = seqs1.size();
    array<int,2> nSig = {0,0};

    // kmer counts in foreground
    // note that each count vector could be of different length due to kmer length difference
    message("counting exact kmers in foreground sequences...");
    map<string, vector<int> > data1 = count_all_kmer_in_seqs(kmers, seqs1, shift);

    ofstream outstream;
    outstream.open(outfile.c_str());

    ofstream outcounts;
    outcounts.open(output_count_file.c_str());

    // compute p-value for exact kmers
    message("computing p-values for exact kmers...");
    for (map<string,vector<int> >::iterator it=data1.begin(); it!=data1.end(); ++it) 
    {
        // counts in foreground
        vector<int> counts1 = it->second;

        // total counts 
        int total_counts1 = sum(counts1);

        double f2 = min(1.0,probs_kmer[it->first]*(1+shift));


        // for each position
        for( int m=0;m<counts1.size();m++)
        { 
            if (counts1[m] == 0 && f2 ==0) continue; 
            // compare background p estimate from shuffling and markov model
            double p = binom_test(nSeq1,counts1[m],f2);
            if (p < pCutoff)
            {
                nSig[0]++;
                double f1 = float(counts1[m])/nSeq1;
                double expected = nSeq1 * f2;
                double z = (counts1[m] - expected) / sqrt(expected*(1-f2));
                // local enrichment
                double local_r = double(counts1[m]) / (total_counts1 - counts1[m]) * (counts1.size()-1);
                double corrected_p = min(1.0,p*nTest);
                outstream << it->first << "\t" << it->first << "\t"  << m-startPos << "\t" << f1 << "\t" << f2 << "\t" << f1/f2 << "\t"  << local_r << "\t" << z << "\t" << p << "\t"<< corrected_p << endl;
                if (corrected_p < pCutoff_B)
                {
                    nSig[1]++;
                    outcounts << it->first << ":" << m-startPos;
                    for( int x=0;x<counts1.size();x++) outcounts << "\t" << double(counts1[x])/nSeq1;
                    outcounts << endl;
                }
            }
        }
    }

    // compute p-value for non-exact kmers, only work without shift
    if (degenerate)
    {   

        message( to_string( nSig[0])+ " significant exact kmers identified");
        message( to_string( nSig[1] ) + " remain after Bonferoni multiple testing correction");
        message( "computing p-values for degenerate kmers...");

        // map from degeneate bases to all allowed bases, e.g. R => AG
        map<char,string> define_iupac = define_IUPAC();
        // for output purpurse only, e.g. will use g to replace H = [ACT] = non-G
        // map<char,string> interpret_iupac = interpret_IUPAC();

        // for each degenerate kmer, combine element kmer counts
        for( int i=0;i<dkmers.size();i++)
        {
            // if already in exact kmers, i.e. contain no degenerate bases, skip
            if ( find(kmers.begin(), kmers.end(), dkmers[i])!=kmers.end() ) continue;

            // expand a degenerate kmer to all possible element/exact kmers
            vector<string> dkmersexp = expand_degenerate_kmer(dkmers[i],define_iupac);

            // for output only, also generate regex version of the degenerate kmer
            string dkmersexp_readable = degenerate_kmer_to_regex(dkmers[i],define_iupac);

            // initialize the count at each position with the first element kmer
            vector<int> counts1  = data1[dkmersexp[0]];
            // sum over the rest element kmers
            for( int j=1;j<dkmersexp.size();j++) counts1 = sum(counts1,data1[dkmersexp[j]]);

            // calculate probs from markov model: sum of element prob
            double f2 = probs_kmer[dkmersexp[0]];
            for( int j=1;j<dkmersexp.size();j++) f2 += probs_kmer[dkmersexp[j]]; 
            //f2 = f2 * (1+shift); // can't have shift on degenerate kmers

            // for output only
            string regex_format = degenerate_kmer_to_regex(dkmers[i],define_iupac);

            // total counts
            int total_counts1 = sum(counts1);

            // calculate  p-value 
            for( int m=0;m<counts1.size();m++)
            { 
                if (counts1[m] == 0 && f2==0) continue;
                double p = binom_test(nSeq1,counts1[m],f2);
                if (p < pCutoff)
                {
                    nSig[0]++;
                    double f1 = float(counts1[m])/nSeq1;
                    double expected = nSeq1 * f2;
                    double z = (counts1[m] - expected) / sqrt(expected*(1-f2)); 
                    double local_r = double(counts1[m]) / (total_counts1 - counts1[m]) * (counts1.size()-1);
                    double corrected_p = min(1.0,p*nTest);
                    outstream  << dkmers[i] << "\t" << regex_format << "\t"  << m-startPos << "\t" << f1 << "\t" << f2 << "\t" << f1/f2 << "\t"  << local_r << "\t" << z << "\t" << p << "\t"<< corrected_p << endl;

                    if (corrected_p < pCutoff_B) 
                    { 
                        nSig[1]++;
                        outcounts << dkmers[i] << ":" << m-startPos;
                        for( int x=0;x<counts1.size();x++) outcounts << "\t" << double(counts1[x])/nSeq1;
                        outcounts << endl;
                    }

                }
            }
        }
    }

    outstream.close();

    outcounts.close();

    return nSig;
} // end of function


// two file comparison, not allow shift and degenerate at the same time
array<int,2> find_significant_kmer_from_two_seq_sets(map<string,string>seqs1, map<string,string>seqs2, vector<string>kmers, vector<string> dkmers, int shift,bool degenerate,double pCutoff,double pCutoff_B, double pseudo,int startPos,int nTest, string outfile,string output_count_file)
{
    int nSeq1 = seqs1.size();
    int nSeq2 = seqs2.size();
    array<int,2> nSig = {0,0};

    // kmer counts in foreground
    // note that each count vector could be of different length due to kmer length difference
    message("counting exact kmers in foreground sequences...");
    map<string, vector<int> > data1 = count_all_kmer_in_seqs(kmers, seqs1, shift);

    // kmer counts in foreground
    message("counting exact kmers in background sequences...");
    map<string, vector<int> > data2 = count_all_kmer_in_seqs(kmers, seqs2, shift);

    ofstream outstream;
    outstream.open(outfile.c_str());

    ofstream outcounts;
    outcounts.open(output_count_file.c_str());

    // compute p-value for exact kmers
    message("computing p-values for exact kmers...");
    for (map<string,vector<int> >::iterator it=data1.begin(); it!=data1.end(); ++it) 
    {
        // counts in foreground
        vector<int> counts1 = it->second;
        // counts in background
        vector<int> counts2 = data2[it->first];

        // total counts 
        int total_counts1 = sum(counts1);

        // for each position
        for( int m=0;m<counts1.size();m++)
        { 
            if (counts1[m] == 0 && counts2[m] == 0) continue; 
            double f2 = float(counts2[m]+pseudo)/nSeq2;
            // compare background p estimate from shuffling and markov model
            double p = binom_test(nSeq1,counts1[m],f2);
            if (p < pCutoff)
            {
                nSig[0]++;
                double f1 = float(counts1[m])/nSeq1;
                double expected = nSeq1 * f2;
                double z = (counts1[m] - expected) / sqrt(expected*(1-f2));
                // local enrichment
                double local_r = double(counts1[m]) / (total_counts1 - counts1[m]) * (counts1.size()-1);
                double corrected_p = min(1.0,p*nTest);
                outstream << it->first << "\t" << it->first << "\t"  << m-startPos << "\t" << f1 << "\t" << f2 << "\t" << f1/f2 << "\t"  << local_r << "\t" << z << "\t" << p << "\t"<< corrected_p << endl;
                if (corrected_p < pCutoff_B)
                {
                    nSig[1]++;
                    outcounts << it->first << ":" << m-startPos;
                    for( int x=0;x<counts1.size();x++) outcounts << "\t" << double(counts1[x])/nSeq1;
                    outcounts << endl;
                }

            }
        }
    }

    // compute p-value for non-exact kmers, only work without shift
    if (degenerate)
    {   

        message(to_string(nSig[0]) + " significant exact kmers identified");
        message(to_string( nSig[1] ) +" remain after Bonferoni multiple testing correction");
        message("computing p-values for degenerate kmers...");

        // map from degeneate bases to all allowed bases, e.g. R => AG
        map<char,string> define_iupac = define_IUPAC();
        // for output purpurse only, e.g. will use g to replace H = [ACT] = non-G
        // map<char,string> interpret_iupac = interpret_IUPAC();

        // for each degenerate kmer, combine element kmer counts
        for( int i=0;i<dkmers.size();i++)
        {
            // if already in exact kmers, i.e. contain no degenerate bases, skip
            if ( find(kmers.begin(), kmers.end(), dkmers[i])!=kmers.end() ) continue;

            // expand a degenerate kmer to all possible element/exact kmers
            vector<string> dkmersexp = expand_degenerate_kmer(dkmers[i],define_iupac);

            // for output only, also generate regex version of the degenerate kmer
            string dkmersexp_readable = degenerate_kmer_to_regex(dkmers[i],define_iupac);

            // initialize the count at each position with the first element kmer
            vector<int> counts1  = data1[dkmersexp[0]];
            vector<int> counts2  = data2[dkmersexp[0]];

            // sum over the rest element kmers
            for( int j=1;j<dkmersexp.size();j++)
            {
                counts1 = sum(counts1,data1[dkmersexp[j]]);
                counts2 = sum(counts2,data2[dkmersexp[j]]);
            }

            // for output only
            string regex_format = degenerate_kmer_to_regex(dkmers[i],define_iupac);

            // total counts
            int total_counts1 = sum(counts1);

            // calculate  p-value 
            for( int m=0;m<counts1.size();m++)
            { 
                if (counts1[m] == 0 && counts2[m]==0) continue;
                double f2 = float(counts2[m]+pseudo)/nSeq2;
                double p = binom_test(nSeq1,counts1[m],f2);
                if (p < pCutoff)
                {
                    nSig[0]++;
                    double f1 = float(counts1[m])/nSeq1;
                    double expected = nSeq1 * f2;
                    double z = (counts1[m] - expected) / sqrt(expected*(1-f2)); 
                    double local_r = double(counts1[m]) / (total_counts1 - counts1[m]) * (counts1.size()-1);
                    double corrected_p = min(1.0,p*nTest); 
                    outstream  << dkmers[i] << "\t" << regex_format << "\t"  << m-startPos << "\t" << f1 << "\t" << f2 << "\t" << f1/f2 << "\t"  << local_r << "\t" << z << "\t" << p << "\t"<< corrected_p << endl;

                    if (corrected_p < pCutoff_B)
                    {
                        nSig[1]++;
                        outcounts << dkmers[i] << ":" << m-startPos;
                        for( int x=0;x<counts1.size();x++) outcounts << "\t" << double(counts1[x])/nSeq1;
                        outcounts << endl;
                    }

                }
            }
        }
    }

    outstream.close();
    outcounts.close();

    return nSig;
} // end of function



void plot_frequency_for_significant_kmer(string inputfile, string outputfile)
{
    string script = 
    "pdf('"+outputfile+"',width=10,height=5) \n"
    "con  <- file('"+inputfile+"', open = 'r') \n"
    "while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) { \n"
    "   myVector <- (strsplit(oneLine, '\\t')) \n"
    "   x=myVector[[1]] \n"
    "   name=x[1] \n"
    "   freq=as.numeric(x[2:length(x)])\n"
    "   plot(freq,type='h',xlab='position',ylab='frequency',main=name) \n"
    "} \n"
    "close(con) \n"
    "dev.off() \n";

    R_run(script);
    
} 
   



char complement(char ch){
        switch (ch)
        {
        case 'A':return 'T';
        case 'C':return 'G';
        case 'G':return 'C';
        case 'T':return 'A';
        default:return 'N';
        }
}

string reverseComplement(string seq){
        int L = seq.length();
        string rc (L,'0');
        for( int i=0;i<L; i++)
        {
                rc[L-i-1] = complement(seq[i]);
        }
        return rc;
}

void ReadOneSeqFromFasta(ifstream& infile, string& name, string& seq){
  // read one sequence from fasta file
  getline(infile,name); // read identifier line
  name = name.substr(1);// remove leading '>'
  seq = "";		// initialize sequence
  string str;		
  while(infile.peek() != '>' && infile.good())
  {// before next '>' and before hitting the end of the file
        getline(infile,str);
        seq.append(str);
    }
  seq = to_upper(seq);
  // trim
  seq.erase(seq.find_last_not_of(" \n\r\t")+1);
} 

//read all sequences in a fasta file
map<string,string> ReadFasta(string filename){
  ifstream fin(filename.c_str());
  
  map<string,string> seqs;
  string name,seq;
  while(fin.good())
  {
    ReadOneSeqFromFasta(fin,name,seq);
    seqs[name] = seq;
  }
  fin.close();
  return seqs;
}

void WriteFasta(map<string,string> seqs, string filename){
    ofstream fout;
    fout.open(filename.c_str());
    for (map<string,string>::iterator it=seqs.begin(); it!=seqs.end(); ++it)
    {
        fout << ">" << it->first << endl << it->second << endl;
    }
}

map<string,string> first_n_bases(map<string,string> seqs,int n){
    // take the first n bases of each sequences
    // if the sequences is shorter than n, discard it
    map<string,string> res;
    for (map<string,string>::iterator it=seqs.begin(); it!=seqs.end(); ++it)
    {   
        if(it->second.size() >= n)
        {
            res[it->first] = it->second.substr(0,n);
        }
    }
    return res;
}

map<string,string> last_n_bases(map<string,string> seqs,int n){
    // take the last n bases of each sequences
    // if the sequences is shorter than n, discard it
    map<string,string> res;
    for (map<string,string>::iterator it=seqs.begin(); it!=seqs.end(); ++it)
    {   
        if(it->second.size() >= n)
        {  
            res[it->first] = it->second.substr(it->second.size()-n,n);
        }
    }
    return res;
}

// convert fasta file to a letter matrix
void fasta_to_letter_matrix(string input, string output){
	ofstream out(output.c_str());
	map<string,string> seqs = ReadFasta(input);
	
	// sequence length, assume fixed
	int L = int(seqs.begin()->second.size());
	
	//header
	out << "SeqID";
	for(int i=0;i<L;i++) out << "\t" << "k1:p" << i+1 ;
	out << endl;
	
    for (map<string,string>::iterator it=seqs.begin(); it!=seqs.end(); ++it)
	{
	    out << it->first;
		for(int i=0;i<L;i++) out << "\t" << it->second[i] ;
		out << endl;
	}
	out.close();
}

// convert fasta file to a letter matrix, no header
void tab_seq_to_letter_matrix(string input, string output, int k_min, int k_max, int col, int skip){
	col = col - 1;
		
	ifstream in(input.c_str());
	ofstream out(output.c_str());

	string line;
	vector<string> flds;

	int i = 0;
	while(in.good() && i++ < skip) getline(in,line);
	
	while(in.good())
	{
		getline(in,line);
		if(line.length()==0) continue;
		flds = string_split(line,"\t");
		// output other columns first
		for (int i=0;i<flds.size();i++) 
		{
			if(i != col) out << flds[i] << "\t";
		}
		// split sequence as letters
		for(int k = k_min;k<=k_max;k++)
			for(int i=0;i<=flds[col].size()-k;i++) 
				out << flds[col].substr(i,k) << "\t" ;
		out << endl;
	}
	in.close();
	out.close();
}

/**************** ushuffle *************/

//uShuffle
//http://digital.cs.usu.edu/~mjiang/ushuffle/
// shuffle sequence preserving k-let
string shuffle_seq_preserving_k_let(string str,int k){
    const char * c = str.c_str();
    char *t = new char[str.size() + 1];
    shuffle(c,t, str.size(), k);
    t[str.size()] = '\0';
    string res(t);
    return res;
}

//
map<string,string> shuffle_seqs_preserving_k_let(map<string,string> seqs, int N, int k){
    // obtain a time-based seed
    srand(time(NULL));
    map<string,string> seqs2;
    for (map<string,string>::iterator it=seqs.begin(); it!=seqs.end(); ++it)
    {
        for( int i =0;i<N;i++) seqs2[it->first + "-" + to_string(i)] = shuffle_seq_preserving_k_let(it->second,k);
    }
    return seqs2;
}



/*******  seq_match   *******/


void mismatches(map<string,string>& mutant,map<string,int>& dist, string motif, int n, set<char> alphabet){
  set<char>::iterator it;
  if (mutant.count(motif) == 0)
  {
    mutant[motif] = "";
    dist[motif]=n;
  }
  if(n==0){return;}
  for( int i=0;i<motif.length();i++)
  {
      string str=motif;
      set<char> ab = alphabet;
      ab.erase(str[i]);
      for (it = ab.begin(); it!=ab.end(); it++)
      {
         str[i] = *it;
         //cout << "mutate "<<motif<<" to "<<str<<endl;
         if (mutant.count(str) >0)
         {
           if(dist[str] >= n)
           {
             //cout << mutant[str] <<endl;
             continue;
           }
         }
         
         //mutated to a new sequence
           //cout <<"new mutation"<<endl;
           mutant[str] = mutant[motif];
           mutant[str].push_back(','); 
           mutant[str].push_back(motif[i]);
           mutant[str].append(to_string(i+1));
           mutant[str].push_back(str[i]);
           dist[str]=n;
           //cout << "tag="<<mutant[str]<<" dist="<<n<<endl;
   
         if (n>1)
         {
           //cout << "subproc" <<endl;
           mismatches(mutant,dist,str,n-1,alphabet);
         }
      }

  }
}

map<string,string> ExpandMotifs(map<string,string>& motifs, int nmismatch, bool rc, set<char> alphabet) { 
  map<string,string> expandedmotifs;
  // generate mismatched motifs
  map<string,string> mutants;
  map<string,int> dist;
  map<string,string>::iterator it;   // iterator for motifs
  map<string,int>::iterator it2; // iterator for mutants
  string name,seq,tmp;
  //cout<<"input motifs"<<endl;
  //PrintMap(motifs);
  for(it=motifs.begin();it!=motifs.end();it++)
  {
    name = (*it).first;
    seq = (*it).second;
   
    mismatches(mutants,dist,seq,nmismatch,alphabet);
    //cout << mutants.size()<<" mutants identified" <<endl;
    //PrintMap(mutants);
    // add mutants to motifs
    for(it2=dist.begin();it2!=dist.end();it2++)
    {
      string tmp = name;
      tmp.append(",").append((*it2).first);
      tmp.append(mutants[(*it2).first]);
      expandedmotifs[tmp] = (*it2).first;
      //cout << name <<","<<tmp<<","<<expandedmotifs[tmp]<<endl;
    }
    // clear the mutants list
    mutants.clear();
    dist.clear();
  }
  //PrintMap(expandedmotifs);
  //cout << expandedmotifs.size() <<" expanded motifs"<<endl;
  //cout <<"add reverse complement"<<endl;
  map<string,string> expandedmotifs_rc = expandedmotifs;
  if (rc)
  {
    for(it=expandedmotifs.begin();it!=expandedmotifs.end();it++)
    {
      name = (*it).first;
      expandedmotifs_rc[name.append(",rc")] = reverseComplement((*it).second);
    }
  }   
 
  return expandedmotifs_rc;
 
}



array<int,2> match(string motiffile, string seqfile, string outfile, int nmismatch, bool rc, set<char> alphabet) {
  int nsite = 0; // total number of matches to report
  int nseq = 0;
  ifstream fmotif, fseq;
  ofstream fout;
  
  // load motifs
  map<string,string> motifs = ReadFasta(motiffile);
  cout <<"["<<current_time()<<"] "<<motifs.size()<< " motifs loaded from "<<motiffile<<endl;
  
  // expand motifs
  map<string,string> expandedmotifs = ExpandMotifs(motifs,nmismatch,rc,alphabet);
  cout <<"["<<current_time()<<"] "<<expandedmotifs.size()<< " motifs after expanding (mismatch/reverse-complement)"<<endl;
  
  //PrintMap(expandedmotifs);
  
  // searching motifs in each sequence
  fseq.open(seqfile.c_str());
  fout.open(outfile.c_str());
  
  string seqname,seq,motifname,motif;
  while(fseq.good())
  {
    // read one sequence
    ReadOneSeqFromFasta(fseq,seqname,seq);
    nseq = nseq + 1;

    cout.flush();
    // iterate over motifs
    map<string,string>::iterator it;
    for(it=expandedmotifs.begin();it!=expandedmotifs.end();it++)
    {
      motifname = (*it).first;
      motif = (*it).second;
      //cout << "searching for "<<motifname<<":"<< motif <<endl;
      set<int> found = findall(seq,motif);
      for (set<int>::iterator it=found.begin(); it!=found.end(); ++it)
      {
        fout <<seqname<<"\t"<< *it << "\t"<< motifname <<"\t"<<motif<<endl; 
      }
      nsite = nsite + found.size();
    }
    cout <<"\r["<<current_time()<<"] " << nsite << " sites found in "<< nseq << " sequences             "  ;
  }
  
  cout << endl; 
  fseq.close();
  fout.close();
  
  array<int,2> res = {{nsite, nseq}};
  
  return res;
}


int tab2bed_galaxy(string infile, string outfile){
  //hg18_chr6_122208322_122209078_+     635     5ss,A7C,G8T-rc  AAGTACCTG
  //hg18_chr6_122208322_122209078_+     553     5ss,C1G,G3A     GAAGTAAGT
  ifstream fin;
  ofstream fout;
  fin.open(infile.c_str());
  fout.open(outfile.c_str());
  string line;
  vector<string> flds;
  vector<string> pos;
  vector<string> nm;
  while(fin)
  {
    getline(fin,line);
    if (line.length() == 0)
      continue;
    flds = string_split(line,"\t");
    pos = string_split(flds[0],"_");
    if (pos.size() < 5)
    {
      cout << "\n!! incorrect sequence name format!\n make sure sequence name looks like: hg18_chr6_122208322_122209078_+\n\n";
      return 0;
    }
    if (pos.size() > 5)
    {// something like chr1_random, skip
      continue;
    }
    string chr = pos[1];
    int start = atoi(pos[2].c_str());
    int end = atoi(pos[3].c_str());
    int match_start = atoi(flds[1].c_str());
    int motifLen = flds[3].length();
    // check if match on the other strand
    string strandness = "sense";
    if (flds[2].find("rc") !=string::npos)
      strandness = "antisense";
    string strand = pos[4];
    if (strand== "+")
    {
      start = start + match_start;
      if (strandness == "antisense")
      {
        strand = "-";
      }
    }
    else//sequence on the - strand of the genome
    {
      start = end - match_start - motifLen;
      if (strandness == "antisense")
      {
        strand = "+";
      }
    }
    end = start + motifLen;
    // number of mismatches
    nm = string_split(flds[2],",");
    int score = nm.size()-2;
    if (strandness == "antisense") {score--;}

    fout << chr <<"\t"<<start<<"\t" <<end<<"\t"<<flds[3]<< "\t"<<score <<"\t"<<strand<< "\t"<<strandness<<"\t"<<flds[2]<<"\t"<<flds[0]<<"\t"<<flds[1]<<endl;

  }
  fin.close();
  fout.close();
  return 1;//return 1 if successfully created bed file
}


int tab2bed_bedtools(string infile, string outfile){
  //chrX:20597309-20645164(-)       7533    u1,CAGGTAAGT    CAGGTAAGT
  //chr17:70312707-70951085(+)      486494  u1,CAGGTAAGT,rc ACTTACCTG
  //can be without strand

  ifstream fin;
  ofstream fout;
  fin.open(infile.c_str());
  fout.open(outfile.c_str());
  string line,strand,chr;
  vector<string> flds,flds0,flds1,flds2,pos,nm;
  while(fin)
  {
    getline(fin,line);
    if (line.length() == 0)
      continue;
    flds = string_split(line,"\t");
    flds0 = string_split(flds[0],":");
    chr = flds0[0];
    strand = "+";
    // seee if there is ( 
    if (flds0[1].find("(") == string::npos)
    { 
      pos = string_split(flds0[1],"-");
    }
    else
    {
      flds1 = string_split(flds0[1],")");
      flds2 = string_split(flds1[0],"(");
      strand = flds2[1];
      pos = string_split(flds2[0],"-");
    }
    int start = atoi(pos[0].c_str());
    int end = atoi(pos[1].c_str());
    int match_start = atoi(flds[1].c_str());
    int motifLen = flds[3].length();
    // check if match on the other strand
    string strandness = "sense";
    if (flds[2].find("rc") !=string::npos)
      strandness = "antisense";
    
    if (strand== "+")
    {
      start = start + match_start;
      if (strandness == "antisense") 
      {
        strand = "-";
      }    
    }
    else//sequence on the - strand of the genome
    {
      start = end - match_start - motifLen;
      if (strandness == "antisense") 
      {
        strand = "+";
      }
    }
    end = start + motifLen;
    // number of mismatches
    nm = string_split(flds[2],",");
    int score = nm.size()-2;
    if (strandness == "antisense") {score--;} 
   
    fout << chr <<"\t"<<start<<"\t" <<end<<"\t"<<flds[3]<< "\t"<<score <<"\t"<<strand<< "\t"<<strandness<<"\t"<<flds[2]<<"\t"<<flds[0]<<"\t"<<flds[1]<<endl;
    
  }
  fin.close();
  fout.close();
  return 1;//return 1 if successfully created bed file
}


