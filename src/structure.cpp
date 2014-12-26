#include <fstream>

#include <stdio.h>
#include <stdlib.h>

#include "utility.h"
#include "structure.h"

#include <map>
#include <iostream>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
using namespace boost;


// open short stems next to loop
// short stem: more unpaired on both sides than paired, .(.  ..((..,
int remove_short_stem(string &structure, int max_paired_length/*=3*/)
{
	// check >mir302locus:1476:-44.10： (((.(((......)))....)))
	// also handle this (((...(((......)))..)))
	
	size_t found_left,found_right;
    for (int L = 1; L <= max_paired_length; L++)
	{
		string pattern_left,pattern_right;
		string allunpair=".......................";
		for (int i =0;i<L;i++) 
			{
				pattern_left += "(";
				pattern_right += ")";
			}
		for (int i=0;i<L;i++) 
			{
				pattern_left = "."+pattern_left+".";
				pattern_right = "."+pattern_right+".";
			}
			
		found_right = structure.find(pattern_right);
		
		// if found right ...))...
		if(found_right != std::string::npos)
		{
			structure.replace(found_right,pattern_right.size(),allunpair.substr(0,pattern_right.size()));
			int right = 0; // how many ) 
			int left = L;
			for(int i = 1;i<=found_right;i++)
			{
				if(structure[found_right-i] == ')') right++;
				else if(structure[found_right-i] == '(')
				{
					if (right>=0) right --;
					if(right<0)
					{
						structure.replace(found_right-i,1,".");
						//cout << structure << "\t" << right << "\t" << left << endl;
						left --;
					}
					if (left == 0) break;
				} 
			}
		}
		
		found_left = structure.find(pattern_left);
		
		if(found_left != std::string::npos)
		{
			structure.replace(found_left,pattern_left.size(),allunpair.substr(0,pattern_left.size()));
			int right = L; // how many ) 
			int left = 0;
			for(int i = found_left+pattern_left.size(); i<structure.size();i++)
			{
				if(structure[i] == '(') left++;
				else if(structure[i] == ')')
				{
					if (left>=0) left --;
					if(left<0)
					{
						structure.replace(i,1,".");
						right --;
					}
					if (right == 0) break;
				} 
			}
		}

		if (found_left == std::string::npos && found_right == std::string::npos) continue;
		// if some changes made, do it again
		remove_short_stem(structure, max_paired_length);
	}
	if (found_left == std::string::npos && found_right == std::string::npos) return 0;
}

// count number of pairs in a structure
int count_pairs(string struc)
{
	int count = 0;
	for (int i =0;i<struc.size();i++)
		if (struc[i] != '.') 
			count ++ ;
	return count;
}

// find the longest stem in a structure
// does not count unpaired positions 
// input  (((((.......((((((((......)))).))))....((((((.....((((((((((((......))))))))))))))))))...)))))
// output .......................................((((((.....((((((((((((......))))))))))))))))))...)))))
/*
	when hit ), start counting paired positions (s), stop when hit ( or the end of the struc (e)
    go backwards, unpair before max.e
*/

int keep_longest_stem(string struc)
{
	vector<int> open_start, open_end, close_start, close_end, open_length, close_length;
	bool at_open = false;
	
	 
	for(int i=0;i<struc.size();i++)
	{
		if(struc[i] == ')')
		{
			if(at_open)
			{
				at_open = false;
				close_start.push_back(i);
				close_length.push_back(1);
				open_end.push_back(i-1);
			}
			else
			{
				close_length[close_length.size()-1] ++;
			}
		}
		else if (struc[i] == '(' )
		{
			if (at_open == false) // hit another stem
				{
					at_open = true;
					if (close_start.size()>0) close_end.push_back(i-1); // the end of last stem
					open_start.push_back(i); // start of a new stem
					open_length.push_back(1);
				} 
				else // still in opn
				{
					open_length[open_length.size()-1] ++;
				}
		}
	}
	close_end.push_back(struc.size()-1);
	
	// find the longest stem
	int longest = 0;
	int longest_len = -1;
	for (int i=0;i<open_start.size();i++)
	{
		cout << "open: " << open_start[i] << "-" << open_end[i] << "," << open_length[i] << endl;
		cout << "close: " << close_start[i] << "-" << close_end[i] << "," << close_length[i] << endl;
		
		int length = open_length[i];
		if (close_length[i] < length ) length = close_length[i];
		
		if (length > longest_len)
		{
			longest_len = length;
			longest = i;
		}
	}
	cout << "longest = " << longest << ", " << longest_len << endl;
	
	/*	
	// unpair other stems
	string all_unpair = "......................................................................";
	struc.replace(longest_end,struc.size()-longest_end-1,all_unpair.substr(0,struc.size()-longest_end-1));
	int i = longest_start;
	while(i>=0)
	{
		if (struc[i] == '(' && longest_len > 0) 
		{
			longest_len --;
		}
		else struc.replace(i,1,".");
		i--;
	}

    cout << struc << endl;
	*/
	return 0;
}

// read a RNAfold file and remove short stem from each
// remove the fold if less than N pairs found
int remove_short_stem_from_file(string infile, string outfile, int max_paired_length/*=4*/, int min_paired_bases_left/*=50*/)
{
	ifstream fin;
	fin.open(infile.c_str());
	
	ofstream fout;
	fout.open(outfile.c_str());

	string line, id_and_seq;
	vector<string> flds;
	
	while(fin)
	{
		getline(fin,line);
		if(line.length()==0) continue;
		if (line[0] == '.' || line[0] == '(' ) //structure line 
		{
			remove_short_stem(line,max_paired_length);
			if (count_pairs(line) >= min_paired_bases_left) fout << id_and_seq << endl << line << endl;
			id_and_seq.clear();
		}
		else if (line[0] == '>') id_and_seq = line;
		else id_and_seq = id_and_seq + "\n" + line;
	}
}

// from RNALfold make fasta file for each hairpin, optionally extend some bases
void RNALfold_to_RNAfold(string infile, string outfile, int min_hairpin_length/*=50*/, int ext/*=20*/)
{
	ifstream fin;
	fin.open(infile.c_str());
	
	ofstream fout;
	fout.open(outfile.c_str());
	
	string line;
	vector<string> flds;
	
	string seqID;
	vector<string> structures;
	string extend_structure;
	for(int i=0;i<ext;i++) extend_structure += ".";
	
	while(fin)
	{
		getline(fin,line);
		if(line.length()==0) continue;
		if (line[0] == '>') // read id
		{
			seqID = line;
			structures.clear();
		}
		else if(line[0] == '.'  || line[0] == '(' ) structures.push_back(line); // read structure
		else // read sequence
		{
			for(int i =0;i<structures.size();i++) // for each structure
			{
				// find first space, left is structure, right is score
				size_t found = structures[i].find(" "); // ((..)). (-10.60) 2799 z= -2.699
				// cout << structures[i].substr(0,found) << "," << structures[i].substr(found) << endl;
				string struc = structures[i].substr(0,found);
				if (struc.size() < min_hairpin_length) continue; // filter by hairpin size
				string tmp = structures[i].substr(found+2); // -10.60) 2799 z= -2.699
				split( flds, tmp, is_any_of(")z"), token_compress_on );
				int pos = stoi(trim_copy(flds[1])); //
				if (pos<ext || pos + struc.size() + ext > line.size() ) continue; // skip
				string name = seqID + ":" + to_string(pos) + ":" + flds[0];
				erase_all(name," ");
				//cout << struc << endl;
				//remove_short_stem(struc, 3,3);
				//cout << struc << endl;
				struc = extend_structure + struc + extend_structure;
				fout << name << endl << line.substr(pos - ext - 1 , struc.size()) << endl << struc << endl;
			}
			getline(fin,line); // skip total score line
		}
	}
	fin.close();
	fout.close();
}

// do not use ! will remove reall miRNA since this is prior to filtering
// from RNALfold make fasta file for each hairpin, optionally extend some bases
void RNALfold_to_RNAfold_filter_overlap(string infile, string outfile, int min_hairpin_length/*=50*/, int ext/*=20*/,int max_overlap_allowed/*=20*/)
{
	ifstream fin;
	fin.open(infile.c_str());
	
	ofstream fout;
	fout.open(outfile.c_str());
	
	string line;
	vector<string> flds;
	
	string seqID;
	vector<string> structures;
	string extend_structure;
	for(int i=0;i<ext;i++) extend_structure += ".";
	
	while(fin)
	{
		getline(fin,line);
		if(line.length()==0) continue;
		if (line[0] == '>') // read id
		{
			seqID = line;
			structures.clear();
		}
		else if(line[0] == '.'  || line[0] == '(' ) structures.push_back(line); // read structure
		else // read sequence
		{
			int previous_pos=1000000000;
			for(int i =0;i<structures.size();i++) // for each structure
			{
				// find first space, left is structure, right is score
				size_t found = structures[i].find(" "); // ((..)). (-10.60) 2799 z= -2.699
				// cout << structures[i].substr(0,found) << "," << structures[i].substr(found) << endl;
				string struc = structures[i].substr(0,found);
				if (struc.size() < min_hairpin_length) continue; // filter by hairpin size
				string tmp = structures[i].substr(found+2); // -10.60) 2799 z= -2.699
				split( flds, tmp, is_any_of(")z"), token_compress_on );
				int pos = stoi(trim_copy(flds[1])); //
				if (pos<ext || pos + struc.size() + ext > line.size() ) continue; // skip if too close to boundary
				// skip if overlap with previous. note that rnalfold reports 3' structures first 
				if(pos + int(struc.size()) - previous_pos > max_overlap_allowed) continue;
				previous_pos = pos;
				string name = seqID + ":" + to_string(pos) + ":" + flds[0];
				erase_all(name," ");
				//cout << struc << endl;
				//remove_short_stem(struc, 3,3);
				//cout << struc << endl;
				struc = extend_structure + struc + extend_structure;
				fout << name << endl << line.substr(pos - ext - 1 , struc.size()) << endl << struc << endl;
			}
			getline(fin,line); // skip total score line
		}
	}
	fin.close();
	fout.close();
}

void hairpin_scoring(string infile, string outfile, vector<positional_kmer> model_str,vector<positional_kmer> model_top,vector<positional_kmer> model_bot,vector<positional_kmer> model_loop, int profile_length)
{
	ifstream fin;
	fin.open(infile.c_str());
	
	ofstream fout;
	fout.open(outfile.c_str());
	
	string line;
	vector<string> flds;
	
    while(fin)
    {
        getline(fin,line);
        if (line.length() == 0)
            continue;

        flds = string_split(line,"\t");
		string id = flds[0];
		string top = to_upper(flds[1]);
		string bot = to_upper(flds[2]);
		string str = to_upper(flds[3]);
		string loop = to_upper(flds[5]);
		
		double best_score = -1000000;
		int best_pos = -1;
		// for each position
		for(int i=0; i <= top.size() - profile_length;i++)
		{
			double score_str = score_sequence_using_PKA_model(model_str, str.substr(i,profile_length));
			double score_top = score_sequence_using_PKA_model(model_top, top.substr(i,profile_length));
			double score_bot = score_sequence_using_PKA_model(model_bot, bot.substr(i,profile_length));
			/*
			double score_total = score_str + score_top + score_bot;
			if (score_total > best_score)
			{
				best_score = score_total;
				best_pos = i;
			}
			*/
			fout << id << "\t" << i << "\t" << str.substr(i,profile_length) << "\t" << score_str << "\t" << top.substr(i,profile_length) << "\t" << score_top << "\t" << bot.substr(i,profile_length) << "\t" << score_bot << endl;
		}
		//double score_loop = score_sequence_using_PKA_model(model_loop, loop);
		//fout << id << "\t" << best_score << "\t" << best_pos << "\t" << score_loop << endl;
	}
	fin.close();
	fout.close();
}
