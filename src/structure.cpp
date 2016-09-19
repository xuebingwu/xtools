#include <fstream>

#include <stdio.h>
#include <stdlib.h>

#include "utility.h"
#include "structure.h"
#include "container.h"
#include "sequence.h"
#include "stat.h"

#include <array>
#include <map>
#include <iostream>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>

#include <algorithm>

using namespace boost;
using namespace std;





// only keep the longest stem, return number of base pairs in this stem
int remove_all_branches(string &structure)
{
	int res = 2;
	while(res > 1) res = remove_shortest_branch(structure);
    return -res;
}

// remove the shortest branch, return the number of branches before removal
// after removal, the number of branches can be less than n-1
// if only one branch found, return -1 * npairs
// if no branch found, return 0
int remove_shortest_branch(string &structure)
{
	// first remove the second one
	// ..((((...))))...(((...)))..(((((((...)))))))...
	
	int n_stem = 0;
	int start=0;
	int start_shortest=0;
	int npairs=0;
	int npairs_shortest=10000000;
	
	bool new_stem = true;
	int n_left = 0; // how many ( before hitting a new )
	for (int i=0;i<structure.size();i++)
	{
		if(structure[i] == ')')
		{
			if (new_stem)
			{
				n_stem ++;
				new_stem = false;
				npairs = 0;
				start = i;
				//cout << "new stem: " << i << endl;
			}
			npairs ++;
		}
		else if (structure[i] == '(')
		{
			if (new_stem == false) // the first time hit a ( after )
			{
				new_stem = true;
				npairs = min(npairs, n_left);
				if (npairs < npairs_shortest)
				{
					npairs_shortest = npairs;
					start_shortest = start;
					//cout << start << "," << npairs << endl;
				}
				n_left = 0;
			}
			n_left ++;
		}
	}
	// the last stem
	
	npairs = min(npairs, n_left);
	if (npairs < npairs_shortest)
	{
		npairs_shortest = npairs;
		start_shortest = start;
		//cout << start << "," << npairs << endl;
	}
	
	//cout << "shortest stem: " << start_shortest << "," << npairs_shortest << endl;
	//cout << "before " << structure << endl;

	if (n_stem == 1)
	{
		return -1 * npairs_shortest;
	}
	else if (n_stem == 0) 
	{
		return 0;
	}	

	// remove the shortest branch
	int L = npairs_shortest;
	for(int i = start_shortest;i < structure.size();i++)
	{
		if (structure[i] == ')' && L > 0)
		{
			L --;
			structure.replace(i,1,".");
		}
	}
	// cout << "after  " << structure << endl;
	
    L = npairs_shortest;
   	for(int i = start_shortest;i >= 0;i--)
   	{
   		if (structure[i] == '(' && L > 0)
   		{
   			L --;
   			structure.replace(i,1,".");
   		}
   	}
	// cout << "after  " << structure << endl;

	return n_stem;
}

//  force a region and its paired region to unpair
int force_unpair(string &structure, int start, int end)
{
	//cout << structure << endl;
	//cout << "unpair " << start << "-" << end << endl;
	
	// how many ( after subtracting )
	int direction = 0;
	for(int i=start;i<=end;i++){
		if (structure[i] == '(') direction++;
		else if (structure[i] == ')')direction--; 
	}
	
	//cout << direction << endl;
	
	string allunpair = ".......................................";
	int L = end - start + 1;
	
	structure.replace(start,L,allunpair.substr(0,L));
	
	if (direction == 0) return 0;
	else if (direction > 0) // (( removed, now search and remove ))
	{
		int newleft = 0;
		int right_to_delete = direction;
		for(int i = end+1;i<structure.size();i++)
		{
			//cout << i << ", newleft=" << newleft << ",right_to_delete=" << right_to_delete << endl;
			if (structure[i] == '(') // if found (
			{
				newleft++; // 
			} else if (structure[i] == ')'){ // if found )
				if(newleft > 0)
				{
					newleft--;
				} else {
					structure.replace(i,1,".");
					right_to_delete--;
					if(right_to_delete == 0)
					{
						break;
					}
				}
			}
		}
	} else// )) removed, now search and remove ((
	{
		int newright = 0;
		int left_to_delete = -direction;
		for(int i = start-1;i>=0;i--)
		{
			if (structure[i] == ')') // if found )
			{
				newright++; // 
			} else if (structure[i] == '('){ // if found (
				if(newright > 0)
				{
					newright--;
				} else {
					structure.replace(i,1,".");
					left_to_delete--;
					if(left_to_delete ==0)
					{
						break;
					}
				}
			}
		}
	}
	//cout << structure << endl;
	return direction;
}


// try to handle interrupted stems like ......)).)).....

bool remove_short_stem_new(string &structure)
{
	vector<int> unpaired_block_start;
	vector<int> unpaired_block_end;
	vector<int> unpaired_block_length;

	int L = structure.size();
	
	bool previous_position_paired;
	int cur_start,cur_len;
	
	if(structure[0] == '.') {
		previous_position_paired = false;
		cur_start = 0;
		cur_len = 1;
		unpaired_block_start.push_back(0);
	} else {
		previous_position_paired = true;
		cur_start = -1;
		cur_len = 0;
	}
	
	for(int i=1;i<L;i++)
	{
		if(structure[i] == '.')
		{
			if (previous_position_paired == false) // ...
			{
				cur_len ++;
			} else { // (.
				// from paired to unpaired
				unpaired_block_start.push_back(i);
				cur_len = 1;
			}
			previous_position_paired = false;
		} else {
			if (previous_position_paired == false) { // .)
				// from unpaired to paired
				unpaired_block_length.push_back(cur_len);
				unpaired_block_end.push_back(i-1);
				cur_len = 0;
			}
			previous_position_paired = true;
		}
	}
	if(structure[L-1] == '.')
	{
		unpaired_block_length.push_back(cur_len);
	}

	// find loops, and change its length to very large number
	for(int i =0;i< unpaired_block_start.size();i++)
	{
		//cout << unpaired_block_start[i] << "-" <<unpaired_block_end[i] << "," << unpaired_block_length[i]  << endl ;
		if(unpaired_block_start[i] != 0 && unpaired_block_end[i] != L-1) // if not terminal unpaired region
		{
			if( structure[unpaired_block_start[i]-1] == '(' && structure[unpaired_block_end[i]+1] == ')') // (...)
			{
				//found a loop
				unpaired_block_length[i] = 10;
			}
		}
		//cout << unpaired_block_start[i] << "-" <<unpaired_block_end[i] << "," << unpaired_block_length[i]  << endl ;
	}
	
	
	// find 
	for(int i =0;i< unpaired_block_start.size();i++)
	{
		for(int j = i+1; j< unpaired_block_start.size();j++)
		{
			// number of bases in between 
			int dist = unpaired_block_start[j] - unpaired_block_end[i] - 1;
			// calculate the effective distance: paired - unpaired
			// dist -= 2*count_unpaired(structure.substr(unpaired_block_end[i]+1,dist));
			if(unpaired_block_length[i] < dist) // next unpaired block is too far away
				break;
			else if (unpaired_block_length[j] >= dist ) {
				//cout << "force unpair " << structure.substr(unpaired_block_end[i]+1,dist) << " "<< unpaired_block_end[i]+1 << "-" << unpaired_block_start[j]-1 << endl;
				force_unpair(structure,unpaired_block_end[i]+1,unpaired_block_start[j]-1);
				//cout << structure << endl;
				remove_short_stem_new(structure);
				return true;
			}
		}
	}
	return false;
}

// open short stems next to loop
// short stem: more or the same unpaired on both sides than paired, .(.  ..((..,
// does not handle interruppted regions such as ......)).)).....
// loop is considered as very large unpaired region
int remove_short_stem(string &structure, int max_paired_length/*=3*/)
{
	// check >mir302locus:1476:-44.10： (((.(((......)))....)))
	// also handle this (((...(((......)))..)))
	


	size_t found_left,found_right;
    for (int L = 1; L <= max_paired_length; L++)
	{
		// cout << "L=" << L << endl;
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
			
		// only right
		string str_right = structure;
		replace(str_right.begin(),str_right.end(),'(','.');
	
		found_right = str_right.find(pattern_right);
		
		// if found right ...))...
		if(found_right != std::string::npos)
		{
			// cout << "found right " << pattern_right << "\t" << found_right <<  endl;
			
			found_right = found_right + L;
			
			structure.replace(found_right,L,allunpair.substr(0,L));
			// cout << structure << endl;
			
			int right = 0; // how many ) 
			int left = L;
			// cout << " start removing pattern " << endl;
			for(int i = 1;i<=found_right ;i++)
			{
				if(structure[found_right-i] == ')') right++;
				else if(structure[found_right-i] == '(')
				{
					if (right>=0) right --;
					if(right<0)
					{
						structure.replace(found_right-i,1,".");
						// cout << structure << "\t" << right << "\t" << left << endl;
						left --;
					}
					if (left == 0) break;
				} 
			}
			// cout << " after removing pattern" << endl;
			//  cout << structure << endl;
		}
		
		// only left 
		string str_left = structure;
		replace(str_left.begin(),str_left.end(),')','.');
		
		found_left = str_left.find(pattern_left);
		
		// ...(((...
		if(found_left != std::string::npos)
		{
		 	// cout << "found left " << pattern_left << "\t" << found_left <<  endl;
			
			found_left = found_left + L;
			structure.replace(found_left,L,allunpair.substr(0,L));
			int right = L; // how many ) 
			int left = 0;
			for(int i = found_left+L+1; i<structure.size();i++)
			{
				if(structure[i] == '(') left++;
				else if(structure[i] == ')')
				{
					if (left>=0) left --;
					if(left<0)
					{
						structure.replace(i,1,".");
						// cout << structure << endl;
						right --;
					}
					if (right == 0) break;
				} 
			}
			//cout << pattern_left << endl;
			// cout << structure << endl;
		}

		if (found_left == std::string::npos && found_right == std::string::npos) continue;
		// if some changes made, do it again
		remove_short_stem(structure, max_paired_length);
	}
	if (found_left == std::string::npos && found_right == std::string::npos) return 0;
}

// find loop positions in a structure
// return: loop1_start, loop1_end, loop2_start, loop2_end
// ...((..))...((..))..
vector<int> find_loops(string structure)
{
	vector<int> res;
	int start = -1;
	for(int i=0;i<structure.size();i++)
	{
		string di = structure.substr(i,2);
		if(di == "(.")
		{
			start = i+1;
		} else if (di == ".)")
		{
			if (start > 0) // first time hit .) after (.
			{
				res.push_back(start);
				res.push_back(i);
				start = -1; // reset start 
			}
		}
	}
	return res;
}

// removing branches, find the loop, use RNAcofold to fold the hairpin without consindering structure within each arm
int hairpin_RNAduplex(string infile, string outfile, string options)
{
	ifstream fin;
	fin.open(infile.c_str());
	
	ofstream fout;
	fout.open((outfile+".tmp").c_str());

	string line, id, seq, fold;
	vector<string> flds;
	
	while(fin)
	{
		getline(fin,line);
		if(line.length()==0) continue;
		if (line[0] == '.' || line[0] == '(' ) //structure line 
		{
			//remove newline
			flds = string_split(line," ");
			fold = flds[0];
			fold.erase(fold.find_last_not_of(" \n\r\t")+1);
			remove_all_branches(fold);
			
			// find the loop
			vector<int> loops = find_loops(fold);
			if(loops.size() != 2) 
			{
				if(loops.size() == 0) message("No loops found! "+id);
				else message("More than one loops found! "+id);
				continue;
			} 
	
			// determine the middle of the loop to separate 5p and 3p arms
			int loop_size = loops[1] - loops[0] + 1;
			double loop_mid_point = loops[0] + loop_size/2;	
			
			string arm5 = seq.substr(0,loop_mid_point);
			string arm3 = seq.substr(loop_mid_point,string::npos) ;
			fout << id << endl << arm5 << endl << ">seq;" << arm5 << ";" << arm3 << endl << arm3 << endl;
		}
		else if (line[0] == '>') id = line;
		else seq = line;//.erase(line.find_last_not_of(" \n\r\t")+1);
	}
	fin.close();
	fout.close();
	system_run("RNAduplex "+options+ " < " + outfile + ".tmp > " + outfile+".duplex");
	// clean up
	system_run("rm "+outfile+".tmp");
	return 1;
}

// RNAduplex output to RNAfold output

/*
>chr14[shredderedfragment]:51722:-45.00
>seq;UAAAUUCAGUUAGGUUGGCCACUAACUAUUAGGUUGGUGCAAAAGUAAUUGUGGUUUUUGCCAUUGAA;AGUAACAGUAAAAGCCAUGAUUACUUUUGCACCAACCAAAUAGUUUGUGUAAUUAGAUUGGGCAGCUU
.(((((((.(((((..(((.(((((((.((((((((((((((((((((((((((((((.&.)))))))))))))))))))))))))))))).))))))).)))))))).))))))).   4,62  :   7,63  (-62.40)
	*/
	
int RNAduplex_to_RNAfold(string infile, string outfile)
{
	ifstream fin;
	fin.open(infile.c_str());

	ofstream fout;
	fout.open(outfile.c_str());

	string line, id, seq1,seq2,seq, fold1,fold2;
	vector<string> flds;

	int n_removed = 0;

	while(fin)
	{
		getline(fin,line);
		if(line.length()==0) continue;
		if(line.substr(0,5) == ">seq;") // sequence line
		{
			flds = string_split(line,";");
			seq1 = flds[1];
			seq2 = flds[2];
		}
		else if (line[0] == '.' || line[0] == '(' ) //structure line 
		{
			// convert : to spaces
			replace(line.begin(),line.end(),':',' ');
			//replace(line.begin(),line.end(),'\t',' ');
			// combine spaces
			line = combine_spaces(line,' ');
			//remove newline
			flds = string_split(line," &,");
			//cout << to_string(flds,"_") << endl;
			
			fold1 = flds[0];
			fold2 = flds[1];
			
			int start1 = stoi(flds[2]);
			int end1 = stoi(flds[3]);
			int start2 = stoi(flds[4]);
			int end2 = stoi(flds[5]);
			
						//if loop size 0
			if (end1 == seq1.size() && start2 == 1 && fold1[fold1.size()-1] == '(' && fold2[0] == ')')
			{
				fold1.replace(fold1.size()-1,1,1,'.');
				fold2.replace(0,1,1,'.');
			}
			 
			fout << id << endl << seq1 << seq2 << endl;
			while(start1>1)
			{
				start1--;
				fold1 = "."+fold1;
			}
			while(end1 < seq1.size())
			{
				end1++;
				fold1 = fold1 + ".";
			}
			while(start2>1)
			{
				start2--;
				fold2 = "."+fold2;
			}
			while(end2 < seq2.size())
			{
				end2++;
				fold2 = fold2 + ".";
			}
			fout << fold1 << fold2 << endl;
			
		}
		else if (line[0] == '>') id = line;
		else message("ERROR: "+line);
	}

	return 1;
}

/*
basal stem definition
starting from the basal stem, go up 35 base pairs, ignore and remove bulges

>hsa-let-7g
CGCUCCGUUUCCUUUUGCCUGAUUCCAGGCUGAGGUAGUAGUUUGUACAGUUUGAGGGUCUAUGAUACCACCCGGUACAGGAGAUAACUGUACAGGCCACUGCCUUGCCAGGAACAGCGCGCCAGCUGCCAAGUGGGG
...............((((((.((((.(((..(((((((.((((((((((((.................................)))))))))))).))))))).))).))))))).))).................
--CGCUCCGUUUCCUUU   -   A    A   UG       A            UGAGGGUCUAUGAUA 
                 UGC CUG UUCC GGC  AGGUAGU GUUUGUACAGUU               C
                 ||| ||| |||| |||  ||||||| ||||||||||||               C
                 GCG GAC AAGG CCG  UCCGUCA CGGACAUGUCAA               C
GGGGUGAACCGUCGACC   C   -    A   -U       C            UAGAGGACAUGGCCC 

output columns / features:
1	stem_5p
2	stem_3p 
3	flank_5p 
4	flank_3p
5	loop_seq
6	actual_stem_length
7	loop_size
8	bulge_count
9	bulge_imbalance_count
10	bulge_imbalance_length 
11	mismatch_count 
12	mismatch_longest
13	pairing_score 
14	to_string(pairing) 
15	to_string(bulge_size)

*/


// assume there is only one loop
int mirna_basal_stem_definition(string id, string seq, string fold, int stem_len, ofstream& out)
{
	// find the loop
	vector<int> loops = find_loops(fold);
	if(loops.size() != 2) 
	{
		//if(loops.size() == 0) message("No loops found! "+id);
		//else message("More than one loops found! "+id);
		return -1;
	} 
	
	// determine the middle of the loop to separate 5p and 3p arms
	int loop_size = loops[1] - loops[0] + 1;
	double loop_mid_point = loops[0] + loop_size/2;	

	//  find the start of the stem 5p
	int left = 0;
	while(fold[left] == '.') left++;
	string flank_5p = seq.substr(0,left);
		
	// find the end of the stem 3p
	int right = seq.size()-1;
	while(fold[right] == '.') right--;
	string flank_3p = seq.substr(right+1);	
	
	// clean up structure by removing asymmetric bulge nucleotides until 35 bps from the basal stem, not counting bulges	
	int length = 0; // length from loop not counting bulges

	string stem_5p,stem_3p;

	
	vector<int> bulge_size (stem_len,0); // bulges along the stem, for each position, startign from the loop, top + bottom -	 
	int bulge_count = 0;
	int bulge_imbalance_count = 0;
	int bulge_imbalance_length = 0;
	int bulge_largest = 0;
	int mismatch_longest;
	vector<int> pairing(stem_len,0); // mismatch:0,GU:1,AU:2,GC:3
	int mismatch_count = 0;
	string a,b,ab;
	int l=0;
	while(length < stem_len)
	{
		if(fold[left] == '(' && fold[right] == ')' ) // match
		{
			a = seq.substr(left,1);
			b = seq.substr(right,1);
			ab = a+b;
			stem_5p += a;
			stem_3p += b;
			bulge_size[length] = l;
			l=0;
			if (ab == "AU" || ab == "UA") pairing[length] = 2;
			else if (ab == "CG" || ab == "GC") pairing[length] = 3;
			else pairing[length] = 1; //UG/GU
			left++;
			right--;
			length++;
		} else if (fold[left] == '.' && fold[right] == '.') // mismatch
		{
			stem_5p += seq.substr(left,1);
			stem_3p += seq.substr(right,1);
			bulge_size[length] = l;
			l=0;
			pairing[length] = 0; 
			left++;
			right--;
			length++;
			mismatch_count++;
		} else if (fold[left] == '.') // insertion on the left
		{
			// from: ((.(...)))
			// to:   (((...)))
			l = 0;
			while(fold[left] == '.')
			{
				left++;
				l++;
			}
			//bulge_size[length] = l;
			bulge_count++;
			bulge_imbalance_count++;
			bulge_imbalance_length += l;
			if (l > bulge_largest) bulge_largest = l;
		} else if(fold[right] == '.')
		{
			// from: (((...)..))
			// to:   (((...).))
			l = 0;
			while(fold[right] == '.') 
			{
				l--;
				right--;
			}
			//bulge_size[length] = l;
			bulge_count ++;
			bulge_imbalance_count--;
			bulge_imbalance_length += l;
			if (-l > bulge_largest) bulge_largest = -l;
		}
		else
		{
			//message("ERROR: left= "+fold[left],+" right= "+fold[right]);
			return -3;
		}
		if(left > loop_mid_point || right < loop_mid_point+1) 
		{
			//cout << seq[left] << "," << loop_mid_point << "," << seq[right] << endl;
			//message("ERROR not enough sequence:"+fold);
			//return -4; // not enough sequence
		}
	}
	
	// mismatch_longest
	string pairing_str = to_string(pairing,"");
	std::array<int,2> tmp = find_longest_run(pairing_str, "0");
	mismatch_longest = tmp[0];
	
	// get stem length
	int actual_stem_length = stem_len;
	int i = stem_len - 1;
	while(pairing[i] == 0)
	{
		actual_stem_length--;
		i--;
	}
	
	// get actual loop sequence
	string loop_seq = seq.substr(left,right-left+1);
	loop_size = loop_seq.size();
	
	// pairing_score
	int pairing_score = sum(pairing);
	
	int total_stem_length = hairpin_stem_length(fold);
	
	// sequence features
	/*
	1. basal_U
	2. basal_G
	3. basal_UG
	4. CNNC1-10
	5. CNNC
	7. apical_GUG
	8. apical_UGU
	9. apical_UGU_GUG
	*/
	
	// basal UG
	int basal_U = 0;
	if(flank_5p[flank_5p.size()-1] == 'U') basal_U=1;
	int basal_G = 0;
	if(stem_5p[0] == 'G') basal_G = 1;
	int basal_UG = 0;
	if (basal_U == 1 && basal_G == 1) basal_UG = 1;
	
	// CNNC
	vector<int> CNNC(11,0);
	int l_tmp = 10;
	if (flank_3p.size() < l_tmp) l_tmp = flank_3p.size();
	for(int i=0;i<l_tmp;i++)
	{
		if(flank_3p[i] == 'C' && flank_3p[i+3] =='C') 
		{
			CNNC[i] = 1;
			CNNC[10]++;
		}
	}
	
	// apical GUG
	vector<int> apical_GUG(3,0); // loop position 1,2,1 or 2
	vector<int> apical_UGU(3,0); // loop position 1,2,1 or 2
	
	if(loop_seq.size() > 2){
		if(loop_seq.substr(0,3) == "GUG") apical_GUG[0] = 1;
		else if(loop_seq.substr(0,3) == "UGU") apical_UGU[0] = 1;
	}
	if(loop_seq.size() > 3)
	{
		if (loop_seq.substr(1,3) == "GUG") apical_GUG[1] = 1;
		else if (loop_seq.substr(1,3) == "UGU") apical_UGU[1] = 1;
	}
	apical_GUG[2] = apical_GUG[0] + apical_GUG[1];
	apical_UGU[2] = apical_UGU[0] + apical_UGU[1];
	
	
	// GHG
	int GHG7 = 0;
	int G7 = 0;
	int H8 =0;
	int G9 = 0;
	int pos = 0;
	i = 0;
	while(pos < 6)
	{
		if(stem_5p[i] != '-' && stem_3p[i] != '-') pos++;
		i++;
	}
	if(stem_3p[i] == 'G') G7 = 1;
	if(stem_3p[i+1] != 'G')H8=1;
	if(stem_3p[i+2] == 'G') G9 = 1;
	GHG7 = G7+H8+G9;
	
	out << id << "\t" << total_stem_length << "\t" << actual_stem_length << "\t" << loop_size << "\t" << bulge_count << "\t" << bulge_imbalance_count  << "\t" << bulge_imbalance_length << "\t" << mismatch_count << "\t" << mismatch_longest << "\t" << pairing_score  << "\t" << basal_U << "\t" << basal_G << "\t" << basal_UG << "\t" << to_string(CNNC) << "\t" << to_string(apical_UGU) << "\t" << to_string(apical_GUG) << "\t" << G7<< "\t" << H8 << "\t" << G9 << "\t" << GHG7 << "\t" << to_string(pairing) << "\t" << to_string(bulge_size) << endl;
	
	
	/*
	// fill sequence with X
	string fillstr = "XXXXXXXXXXXXXXXX";
	
	cout << "flank 5p " << flank_5p ;	
	flank_5p = fillstr + flank_5p;
	flank_5p = flank_5p.substr(flank_5p.size()-2,2);
	cout << " => " << flank_5p << endl;
	
	cout << "flank 3p " << flank_3p ;
	flank_3p = flank_3p + fillstr;
	flank_3p = flank_3p.substr(0,12);
	cout << " => " << flank_3p << endl;
	
	cout << "loop_seq " << loop_seq << endl ;	
	string loop_seq_left = loop_seq + fillstr;	
	loop_seq_left = loop_seq_left.substr(0,4);
	string loop_seq_right = fillstr + loop_seq;
	loop_seq_right = fillstr + loop_seq;	
	loop_seq_right = loop_seq_right.substr(loop_seq_right.size()-4,4);
	cout << loop_seq_left << "," << loop_seq_right << endl;	
	
	
	
	// output
	//cout << id << "\t" <<  stem_5p << "\t" << stem_3p << "\t" << flank_5p << "\t" << flank_3p << "\t" << loop_seq << "\t" << total_stem_length << "\t" << actual_stem_length << "\t" << loop_size << "\t" << bulge_count << "\t" << bulge_imbalance_count  << "\t" << bulge_imbalance_length << "\t" << mismatch_count << "\t" << mismatch_longest << "\t" << pairing_score << "\t" << to_string(pairing,",") << "\t" << to_string(bulge_size,",") << endl ;
	
	out << id << "\t" << total_stem_length << "\t" << actual_stem_length << "\t" << loop_size << "\t" << bulge_count << "\t" << bulge_imbalance_count  << "\t" << bulge_imbalance_length << "\t" << mismatch_count << "\t" << mismatch_longest << "\t" << pairing_score << "\t"  << basal_U << "\t" << basal_G << "\t" << basal_UG << "\t" << to_string(CNNC) << "\t" << to_string(apical_UGU) << "\t" << to_string(apical_GUG) << "\t" << to_string(GHG) << "\t"<< to_string(pairing) << "\t" << to_string(bulge_size) << endl;
	*/
	return 0;
	
}


/*
loop definition
starting from the loop, go back 35 base pairs, ignore and remove bulges

>hsa-let-7g
CGCUCCGUUUCCUUUUGCCUGAUUCCAGGCUGAGGUAGUAGUUUGUACAGUUUGAGGGUCUAUGAUACCACCCGGUACAGGAGAUAACUGUACAGGCCACUGCCUUGCCAGGAACAGCGCGCCAGCUGCCAAGUGGGG
...............((((((.((((.(((..(((((((.((((((((((((.................................)))))))))))).))))))).))).))))))).))).................
--CGCUCCGUUUCCUUU   -   A    A   UG       A            UGAGGGUCUAUGAUA 
                 UGC CUG UUCC GGC  AGGUAGU GUUUGUACAGUU               C
                 ||| ||| |||| |||  ||||||| ||||||||||||               C
                 GCG GAC AAGG CCG  UCCGUCA CGGACAUGUCAA               C
GGGGUGAACCGUCGACC   C   -    A   -U       C            UAGAGGACAUGGCCC 

output columns / features:
1	stem_5p
2	stem_3p 
3	flank_5p 
4	flank_3p
5	loop_seq
6	actual_stem_length
7	loop_size
8	bulge_count
9	bulge_imbalance_count
10	bulge_imbalance_length 
11	mismatch_count 
12	mismatch_longest
13	pairing_score 
14	to_string(pairing) 
15	to_string(bulge_size)

*/

int mirna_loop_definition(string id, string seq, string fold, int stem_len, ofstream& out)
{
	// find the loop
	vector<int> loops = find_loops(fold);
	if(loops.size() != 2) 
	{
		//if(loops.size() == 0) message("No loops found! "+id);
		//else message("More than one loops found! "+id);
		return -1;
	} 
	
	// if the loop starts too early, there is not enough sequences to make a long stem
	if(loops[0] < stem_len || loops[1] + stem_len > seq.size()) 
	{
		//message("hairpin not long enough! "+ id);
		return -2;
	}
		
	int loop_size = loops[1] - loops[0] + 1;
	string loop_seq = seq.substr(loops[0],loop_size);
	
	//cout << "loop seq:" << loop_seq << ", size " << loop_size << endl;

	// clean up structure by removing asymmetric bulge nucleotides until 35 bps from the loop, not counting bulges	
	int left = loops[0]-1; // point to (
	int right = loops[1]+1; // point to )

	string stem_5p;
	string stem_3p;
	string paired;
	
	vector<int> bulge_size (stem_len,0); // bulges along the stem, for each position, startign from the loop, top + bottom -	 
	int bulge_count = 0;
	int bulge_imbalance_count = 0;
	int bulge_imbalance_length = 0;
	int bulge_largest = 0;
	int mismatch_longest;
	vector<int> pairing(stem_len,0); // mismatch:0,GU:1,AU:2,GC:3
	int mismatch_count = 0;
	string a,b,ab;
	int l = 0;
	int length = stem_len-1; // length from loop not counting bulges
	while(length >= 0)
	{
		if(fold[left] == '(' && fold[right] == ')' ) // match
		{
			a = seq.substr(left,1);
			b = seq.substr(right,1);
			ab = a+b;
			stem_5p = a + stem_5p;
			stem_3p = b + stem_3p;
			paired = "|" + paired;
			bulge_size[length] = l;
			l=0;
			if (ab == "AU" || ab == "UA") pairing[length] = 2;
			else if (ab == "CG" || ab == "GC") pairing[length] = 3;
			else pairing[length] = 1; //UG/GU
			left--;
			right++;
			length--;
		} else if (fold[left] == '.' && fold[right] == '.') // mismatch
		{
			stem_5p = seq.substr(left,1) + stem_5p;
			stem_3p = seq.substr(right,1) + stem_3p;
			paired = " " + paired;
			bulge_size[length] = l;
			l=0;
			pairing[length] = 0; 
			left--;
			right++;
			length--;
			mismatch_count++;
		} else if (fold[left] == '.') // insertion on the left
		{
			// from: ((.(...)))
			// to:   (((...)))
			l = 0;
			while(fold[left] == '.')
			{
				stem_3p = "-" + stem_3p;
				stem_5p = seq.substr(left,1) + stem_5p;
				paired = " " + paired;
				left--;
				l++;
			}
			bulge_count++;
			bulge_imbalance_count++;
			bulge_imbalance_length += l;
			if (l > bulge_largest) bulge_largest = l;
		} else if(fold[right] == '.')
		{
			// from: (((...)..))
			// to:   (((...).))
			l = 0;
			while(fold[right] == '.') 
			{
				stem_5p = "-" + stem_5p;
				stem_3p = seq.substr(right,1) + stem_3p;
				paired = " " + paired;
				l--;
				right++;
			}
			bulge_count ++;
			bulge_imbalance_count--;
			bulge_imbalance_length += l;
			if (-l > bulge_largest) bulge_largest = -l;
		}
		else
		{
			//message("ERROR: left= "+fold.substr(left,1)+" right= "+fold.substr(right,1) + ": id="+id);
			return -3;
		}
		if(left < 0 || right >= fold.size()) 
		{
			//message("not enough sequence: "+id);
			return -4; // not enough sequence
		}
	}
	

	
	// mismatch_longest
	string pairing_str = to_string(pairing,"");
	std::array<int,2> tmp = find_longest_run(pairing_str, "0");
	mismatch_longest = tmp[0];
	
	// get stem length
	int actual_stem_length = stem_len;
	int i = 0;
	while(pairing[i] == 0)
	{
		actual_stem_length--;
		i++;
	}
	
	// pairing_score
	int pairing_score = sum(pairing);
	
	string flank_5p = seq.substr(0,left+1);
	string flank_3p = seq.substr(right);
		
	int total_stem_length = hairpin_stem_length(fold);
	
	// sequence features
	// basal UG
	int basal_U = 0;
	if(flank_5p[flank_5p.size()-1] == 'U') basal_U=1;
	int basal_G = 0;
	if(stem_5p[0] == 'G') basal_G = 1;
	int basal_UG = 0;
	if (basal_U == 1 && basal_G == 1) basal_UG = 1;
	
	// CNNC
	vector<int> CNNC(11,0);
	int l_tmp = 10;
	if (flank_3p.size() < l_tmp) l_tmp = flank_3p.size();
	for(int i=0;i<l_tmp;i++)
	{
		if(flank_3p[i] == 'C' && flank_3p[i+3] =='C') 
		{
			CNNC[i] = 1;
			CNNC[10]++;
		}
	}
	
	// apical GUG
	vector<int> apical_GUG(3,0); // loop position 1,2,1 or 2
	vector<int> apical_UGU(3,0); // loop position 1,2,1 or 2
	
	if(loop_seq.size() > 2){
		if(loop_seq.substr(0,3) == "GUG") apical_GUG[0] = 1;
		else if(loop_seq.substr(0,3) == "UGU") apical_UGU[0] = 1;
	}
	if(loop_seq.size() > 3)
	{
		if (loop_seq.substr(1,3) == "GUG") apical_GUG[1] = 1;
		else if (loop_seq.substr(1,3) == "UGU") apical_UGU[1] = 1;
	}
	apical_GUG[2] = apical_GUG[0] + apical_GUG[1];
	apical_UGU[2] = apical_UGU[0] + apical_UGU[1];
	
	// GHG
	int GHG7 = 0;
	int G7 = 0;
	int H8 =0;
	int G9 = 0;
	int pos = 0;
	i = 0;
	while(pos < 6)
	{
		if(stem_5p[i] != '-' && stem_3p[i] != '-') pos++;
		i++;
	}
	if(stem_3p[i] == 'G') G7 = 1;
	if(stem_3p[i+1] != 'G')H8=1;
	if(stem_3p[i+2] == 'G') G9 = 1;
	GHG7 = G7+H8+G9;
	
	out << id << "\t" << total_stem_length << "\t" << actual_stem_length << "\t" << loop_size << "\t" << bulge_count << "\t" << bulge_imbalance_count  << "\t" << bulge_imbalance_length << "\t" << mismatch_count << "\t" << mismatch_longest << "\t" << pairing_score  << "\t" << basal_U << "\t" << basal_G << "\t" << basal_UG << "\t" << to_string(CNNC) << "\t" << to_string(apical_UGU) << "\t" << to_string(apical_GUG) << "\t" << G7<< "\t" << H8 << "\t" << G9 << "\t" << GHG7 << "\t" << to_string(pairing) << "\t" << to_string(bulge_size) << endl;

	/**
	// plot flanking sequences
	stem_5p = flank_5p + stem_5p;
	stem_3p = reverse(flank_3p) + stem_3p;
	if(flank_5p.size () > flank_3p.size())
	{
		paired = string(flank_5p.size(),' ') + paired;
		for(int i=flank_3p.size();i<flank_5p.size();i++) 
		{
			stem_3p = " "+stem_3p; 
		}
	} else {
		paired = string(flank_3p.size(), ' ') + paired;
		for(int i=flank_5p.size();i<flank_3p.size();i++) 
		{
			stem_5p = " "+stem_5p; 
		}
	}
	int loop_mid_point = loop_size/2;
	stem_5p += loop_seq.substr(0,loop_mid_point);
	stem_3p += reverse(loop_seq.substr(loop_mid_point,string::npos)); 	
	
	// indicate cleavate sites
	string cleavage5(stem_5p.size(),' ');
	string cleavage3(stem_3p.size(),' ');

	// find 5' cleavage site
	int i1,pos1;
	i = 0;
	pos = 0;
	while(pos < 30)
	{
		if(stem_5p[i] != '-' && stem_5p[i] != ' ') pos++;
		i++;
	}
	cleavage5.replace(i,1,1,'-');
	// find -14U, -13G
	i1 = i-1;
	pos1 = 1;
	while(pos1 < 13)
	{
		if(stem_5p[i1] != '-' && stem_3p[i1] != '-') pos1++;
		i1--;
	}
	while(stem_5p[i1] == '-' || stem_3p[i1] == '-') i1--;
	cleavage5.replace(i1,1,1,'G');
	while(pos1 < 14)
	{
		if(stem_5p[i1] != '-' && stem_3p[i1] != '-') pos1++;
		i1--;
	}
	cleavage5.replace(i1,1,1,'U');
	
	
	// find the end of 22 bp stem
	i++;
	pos++;
	while(pos < 51)
	{
		if(stem_5p[i] != '-' && stem_3p[i] != '-') pos++;
		cleavage5.replace(i,1,1,'-');
		i++;
	}
	while(stem_5p[i] == '-' || stem_3p[i] == '-') i++;
	cleavage5.replace(i,4,"UGUG");
	
	i = 0;
	pos = 0;
	while(pos < 30)
	{
		if(stem_3p[i] != '-' && stem_3p[i] != ' ') pos++;
		i++;
	}
	cleavage3.replace(i,1,1,'-');
	// find GHG
	i1 = i-1;
	pos1 = 1;
	while(pos1 < 3)
	{
		if(stem_5p[i1] != '-' && stem_3p[i1] != '-') pos1++;
		i1--;
	}
	while(stem_5p[i1] == '-' || stem_3p[i1] == '-') i1--;
	cleavage3.replace(i1,1,1,'G');
	while(pos1 < 4)
	{
		if(stem_5p[i1] != '-' && stem_3p[i1] != '-') pos1++;
		i1--;
	}
	while(stem_5p[i1] == '-' || stem_3p[i1] == '-') i1--;
	cleavage3.replace(i1,1,1,'H');
	while(pos1 < 5)
	{
		if(stem_5p[i1] != '-' && stem_3p[i1] != '-') pos1++;
		i1--;
	}
	while(stem_5p[i1] == '-' || stem_3p[i1] == '-') i1--;
	cleavage3.replace(i1,1,1,'G');
	while(pos1 < 17)
	{
		if(stem_5p[i1] != '-' && stem_3p[i1] != '-') pos1++;
		i1--;
	}
	cleavage3.replace(i1-4,5,"CCNCC");
	i++;
	pos++;
	while(pos < 51)
	{
		if(stem_3p[i] != '-' && stem_5p[i] != '-') pos++;
		cleavage3.replace(i,1,1,'-');
		i++;
	}
	cleavage3.replace(i,1,1,'-');

	
	/**
	cout << id << endl;//<< seq << endl << fold << endl;
	cout << cleavage5 << endl << stem_5p << endl << paired << string(loop_mid_point, ' ') << ")"  << endl << stem_3p << endl << cleavage3 << endl << endl;
	/**/

	return 0;
	
}

// plot hairpin
int plot_hairpin(string id, string seq, string fold)
{
	// find the loop
	vector<int> loops = find_loops(fold);
	if(loops.size() != 2) 
	{
		//if(loops.size() == 0) message("No loops found! "+id);
		//else message("More than one loops found! "+id);
		return -1;
	} 
	
	// determine the middle of the loop to separate 5p and 3p arms
	int loop_size = loops[1] - loops[0] + 1;
	double loop_mid_point = loops[0] + loop_size/2;	
	string loop_seq = seq.substr(loops[0],loop_size);

	//  find the start of the stem 5p
	int left = 0;
	while(fold[left] == '.') left++;
	string flank_5p = seq.substr(0,left);
		
	// find the end of the stem 3p
	int right = seq.size()-1;
	while(fold[right] == '.') right--;
	string flank_3p = seq.substr(right+1);	
	
	// clean up structure by removing asymmetric bulge nucleotides until 35 bps from the basal stem, not counting bulges	
	int length = 0; // length from loop not counting bulges

	string stem_5p,stem_3p;
	string paired;
	
	//vector<int> bulge_size (stem_len,0); // bulges along the stem, for each position, startign from the loop, top + bottom -	 
	int bulge_count = 0;
	int bulge_imbalance_count = 0;
	int bulge_imbalance_length = 0;
	int bulge_largest = 0;
	int mismatch_longest;
	//vector<int> pairing(stem_len,0); // mismatch:0,GU:1,AU:2,GC:3
	int mismatch_count = 0;
	string a,b,ab;
	int l=0;
	while(left < loops[0] && right > loops[1])
	{
		if(fold[left] == '(' && fold[right] == ')' ) // match
		{
			a = seq.substr(left,1);
			b = seq.substr(right,1);
			ab = a+b;
			stem_5p += a;
			stem_3p += b;
			
			//bulge_size[length] = l;
			l=0;
			if (ab == "AU" || ab == "UA") 
				{
					paired += "|"; 
					//pairing[length] = 2;
				}
			else if (ab == "CG" || ab == "GC")
				{
					paired += "|"; 
					//pairing[length] = 3;
				} 
			else 
			{
				//pairing[length] = 1; //UG/GU
				paired += "."; 
			}
			left++;
			right--;
			length++;
		} else if (fold[left] == '.' && fold[right] == '.') // mismatch
		{
			stem_5p += seq.substr(left,1);
			stem_3p += seq.substr(right,1);
			paired += " ";
			//bulge_size[length] = l;
			l=0;
			//pairing[length] = 0; 
			left++;
			right--;
			length++;
			mismatch_count++;
		} else if (fold[left] == '.') // insertion on the left
		{
			// from: ((.(...)))
			// to:   (((...)))
			l = 0;
			while(fold[left] == '.')
			{
				stem_3p += "-";
				stem_5p += seq.substr(left,1);
				paired += " ";
				left++;
				l++;
			}
			//bulge_size[length] = l;
			bulge_count++;
			bulge_imbalance_count++;
			bulge_imbalance_length += l;
			if (l > bulge_largest) bulge_largest = l;
		} else if(fold[right] == '.')
		{
			// from: (((...)..))
			// to:   (((...).))
			l = 0;
			while(fold[right] == '.') 
			{
				stem_5p += "-";
				stem_3p += seq.substr(right,1);
				paired += " ";
				l--;
				right--;
			}
			//bulge_size[length] = l;
			bulge_count ++;
			bulge_imbalance_count--;
			bulge_imbalance_length += l;
			if (-l > bulge_largest) bulge_largest = -l;
		}
		else
		{
			//message("ERROR: left= "+fold[left],+" right= "+fold[right]);
			return -3;
		}
		if(left > loop_mid_point || right < loop_mid_point+1) 
		{
			//cout << seq[left] << "," << loop_mid_point << "," << seq[right] << endl;
			//message("ERROR not enough sequence:"+fold);
			//return -4; // not enough sequence
		}
	}

	//cout << id << endl << seq << endl << fold << endl;
	//cout << stem_5p << endl << paired << endl << stem_3p << endl;
	
	
	/**/
	// plot flanking sequences
	stem_5p = flank_5p + stem_5p;
	stem_3p = reverse(flank_3p) + stem_3p;
	if(flank_5p.size () > flank_3p.size())
	{
		paired = string(flank_5p.size(),' ') + paired;
		for(int i=flank_3p.size();i<flank_5p.size();i++) 
		{
			stem_3p = " "+stem_3p; 
		}
	} else {
		paired = string(flank_3p.size(), ' ') + paired;
		for(int i=flank_5p.size();i<flank_3p.size();i++) 
		{
			stem_5p = " "+stem_5p; 
		}
	}

	int half_loop_sie = loop_size/2;
	stem_5p += loop_seq.substr(0,half_loop_sie);
	stem_3p += reverse(loop_seq.substr(half_loop_sie,string::npos)); 	
	
	// indicate cleavate sites
	string cleavage5(stem_5p.size(),' ');
	string cleavage3(stem_3p.size(),' ');

	// find 5' cleavage site
	int i,pos,i1,pos1,p5;
	i = 0;
	pos = 0;
	while(pos < 30)
	{
		if(stem_5p[i] != '-' && stem_5p[i] != ' ') pos++;
		i++;
	}
	while(stem_5p[i] == '-') i++;
	p5 = i;// position of 5p cleavage
	cleavage5.replace(i,1,1,'-');
	// find -14U, -13G
	i1 = i-1;
	pos1 = 1;
	while(pos1 < 13)
	{
		if(stem_5p[i1] != '-' && stem_3p[i1] != '-') pos1++;
		i1--;
	}
	while(stem_5p[i1] == '-' || stem_3p[i1] == '-') i1--;
	cleavage5.replace(i1,1,1,'G');
	while(pos1 < 14)
	{
		if(stem_5p[i1] != '-' && stem_3p[i1] != '-') pos1++;
		i1--;
	}
	while(stem_5p[i1] == '-') i1--;
	cleavage5.replace(i1,1,1,'U');
	
	
	// find the end of 22 bp stem
	i++;
	pos++;
	while(pos < 51)
	{
		if(stem_5p[i] != '-' && stem_3p[i] != '-') pos++;
		cleavage5.replace(i,1,1,'-');
		i++;
	}
	while(stem_5p[i] == '-' || stem_3p[i] == '-') 
	{
		cleavage5.replace(i,1,1,'-');
		i++;
	}
	cleavage5.replace(i,4,"UGUG");
	
	// find 3p cleavage
	i = p5-1;
	pos = 0;
	while(pos < 1) // count 2 base pairs to the left of 5p cleavage is the 3p cleavage
	{
		if(stem_5p[i] != '-' && stem_3p[i] != '-') pos++;
		i--;
	}
	while(stem_5p[i] == '-' || stem_3p[i] == '-') i--;
	cleavage3.replace(i,1,1,'-');
	// find GHG
	i1 = i-1;
	pos1 = 1;
	while(pos1 < 3)
	{
		if(stem_5p[i1] != '-' && stem_3p[i1] != '-') pos1++;
		i1--;
	}
	while(stem_5p[i1] == '-' || stem_3p[i1] == '-') i1--;
	cleavage3.replace(i1,1,1,'G');
	while(pos1 < 4)
	{
		if(stem_5p[i1] != '-' && stem_3p[i1] != '-') pos1++;
		i1--;
	}
	while(stem_5p[i1] == '-' || stem_3p[i1] == '-') i1--;
	cleavage3.replace(i1,1,1,'H');
	while(pos1 < 5)
	{
		if(stem_5p[i1] != '-' && stem_3p[i1] != '-') pos1++;
		i1--;
	}
	while(stem_5p[i1] == '-' || stem_3p[i1] == '-') i1--;
	cleavage3.replace(i1,1,1,'G');
	while(pos1 < 17)
	{
		if(stem_5p[i1] != '-' && stem_3p[i1] != '-') pos1++;
		i1--;
	}
	cleavage3.replace(i1-4,5,"CCNCC");
	i++;
	pos++;
	while(pos < 22)
	{
		if(stem_3p[i] != '-' && stem_5p[i] != '-') pos++;
		cleavage3.replace(i,1,1,'-');
		i++;
	}
	cleavage3.replace(i,1,1,'-');

	paired = paired +  string(half_loop_sie, ' ') ;

	/**/
	cout << id << "\tL=" << length << endl << endl;//<< seq << endl << fold << endl;
	cout << cleavage5 << endl << stem_5p << endl << paired << ")"  << endl << stem_3p << endl << cleavage3 << endl << endl;
	
	// repair
	for (i=0;i<paired.size();i++)
	{
		string tmp = stem_5p.substr(i,1) + stem_3p.substr(i,1);
		if(tmp == "AU" || tmp == "UA" ||tmp == "GC" || tmp == "CG" ) paired.replace(i,1,1,'|');
		else if (tmp == "GU" || tmp == "UG" )paired.replace(i,1,1,'.');
	}
	
	cout << cleavage5 << endl << stem_5p << endl << paired << ")"  << endl << stem_3p << endl << cleavage3 << endl << endl;
	
	// anchor at 5' cleavage, 25 bp left, 25 bp right, remove bulge
	pos = 0;
	i=p5;
	string stem_top,stem_bot,stem_mid;
	while(pos < 25)
	{
		if(stem_3p[i] != '-' && stem_5p[i] != '-')
		{
			stem_top += stem_5p[i];
			stem_bot += stem_3p[i];
			stem_mid += paired[i];
			pos++;	
		} 
		i++;
	}
	pos = 0;
	i=p5-1;
	while(pos < 25)
	{
		if(stem_3p[i] != '-' && stem_5p[i] != '-')
		{
			stem_top = stem_5p[i] + stem_top;
			stem_bot = stem_3p[i] + stem_bot;
			stem_mid = paired[i] + stem_mid;
			pos++;	
		} 
		i--;
	}
	
	cout << "TOP\t" << stem_top << endl << "MID\t" << stem_mid << endl << "BOT\t"<< stem_bot << endl;
	
	return 0;
	
}


int plot_hairpin_from_file(string infile)
{
	ifstream fin;
	fin.open(infile.c_str());
	

	string line, id, seq;
	vector<string> flds;
		
	while(fin)
	{
		getline(fin,line);
		if(line.length()==0) continue;
		line.erase(line.find_last_not_of(" \n\r\t")+1);
		if (line[0] == '.' || line[0] == '(' ) //structure line 
		{
			flds = string_split(line," ");
			plot_hairpin( id,  seq, flds[0]);
		}
		else if (line[0] == '>') id = line;
		else seq = line;
	}
	return 1;
}




int mirna_feature_from_file(string infile, string outfile_basal, string outfile_loop, int stem_len)
{
	ifstream fin;
	fin.open(infile.c_str());
	
	ofstream fout_basal;
	fout_basal.open(outfile_basal.c_str());
	
	ofstream fout_loop;
	fout_loop.open(outfile_loop.c_str());
	
	// header 
	string header = "ID\ttotal_stem_length\tactual_stem_length\tloop_size\tbulge_count\tbulge_imbalance_count\tbulge_imbalance_length\tmismatch_count\tmismatch_longest\tpairing_score";

	// basal UG motif
	header += "\tbasal_U\tbasal_G\tbasal_UG";
	// CNNC
	for(int i=0;i<10;i++) header += "\tCNNC_" + to_string(i);
	header += "\tCNNC_count";
	header += "\tapical_UGU_0\tapical_UGU_1\tapical_UGU\tapical_GUG_0\tapical_GUG_1\tapical_GUG" ;
	
	header += "\tG7\tH8\tG9\tGHG7" ;
	
	// pairing status at each position
	for(int i=0;i<stem_len;i++) header += "\tpairing_"+to_string(i);
	// bulge size at each position
	for(int i=0;i<stem_len;i++) header += "\tbulge_" + to_string(i);
	
	fout_basal << header << endl;
	fout_loop << header << endl;

	string line, id, seq;
	vector<string> flds;
		
	while(fin)
	{
		getline(fin,line);
		if(line.length()==0) continue;
		line.erase(line.find_last_not_of(" \n\r\t")+1);
		if (line[0] == '.' || line[0] == '(' ) //structure line 
		{
			flds = string_split(line," ");
			mirna_loop_definition(id,seq, flds[0], stem_len,fout_loop);
			mirna_basal_stem_definition(id,seq, flds[0], stem_len,fout_basal);
		}
		else if (line[0] == '>') id = line;
		else seq = line;
	}
	return 1;
}


// open terminal GU
bool open_GU(string seq, string& fold, int left, int right)
{
	if(seq[left] == 'G' && seq[right] == 'U' || seq[left] == 'U' && seq[right] == 'G' )
	{
		fold.replace(left,1,1,'.');
		fold.replace(right,1,1,'.');
		return true;
	}
	return false;
}


/*
trim hairpins from basal stem and loop

value:
-1: incorrect number of loops (not 1)
-2: incorrect structure
*/
// assume there is only one loop
int trim_hairpin(string seq, string& fold, bool noClosingGU)
{
	
	// find the loop
	vector<int> loops = find_loops(fold);
	if(loops.size() != 2) 
	{
		//if(loops.size() == 0) message("No loops found! "+id);
		//else message("More than one loops found! "+id);
		return -1;
	} 
	
	// determine the middle of the loop to separate 5p and 3p arms
	int loop_size = loops[1] - loops[0] + 1;
	double loop_mid_point = loops[0] + loop_size/2;	

	//  find the start of the stem 5p
	int left = 0;
	while(fold[left] == '.') left++;
	string flank_5p = seq.substr(0,left);
		
	// find the end of the stem 3p
	int right = seq.size()-1;
	while(fold[right] == '.') right--;
	string flank_3p = seq.substr(right+1);	

	if(noClosingGU)
	{
		if(open_GU( seq,  fold, left, right))
		{
			left++;
			right--;
		}
	}

	string stem_5p;
	string stem_3p;
	string pairing;
	string a,b,ab;
	
	int score = 0;
	int left0 = left;
	int right0 = right;
	
	while(left < loops[0]) // before hit the loop
	{
		if(fold[left] == '(' && fold[right] == ')' ) // match
		{
			a = seq.substr(left,1);
			b = seq.substr(right,1);
			ab = a+b;
			stem_5p += a;
			stem_3p += b;
			if (score < 0) score = 0; // if not started calculating score, now start
			if (ab == "AU" || ab == "UA") { pairing += "2"; score += 2;}
			else if (ab == "CG" || ab == "GC") { pairing += "3"; score += 3;}
			else {pairing += "1"; score+=1;} //UG/GU
			left++;
			right--;
		} 
		else // non-match
		{
			if (fold[left] == '.' && fold[right] == '.') // mismatch
			{
				stem_5p += seq.substr(left,1);
				stem_3p += seq.substr(right,1);
				pairing += "0"; 
				if(score >= 3) score -= 3;
				left++;
				right--;
			} else if (fold[left] == '.') // insertion on the left
			{
				// from: ((.(...)))
				// to:   (((...)))
				while(fold[left] == '.')
				{
					stem_5p += seq.substr(left,1);
					stem_3p += "-";
					pairing += " "; 
					left++;
					if(score >= 3) score -= 3;
				}
			} else if(fold[right] == '.')
			{
				// from: (((...)..))
				// to:   (((...).))
				while(fold[right] == '.') 
				{
					stem_3p += seq.substr(right,1);
					stem_5p += "-";
					pairing += " "; 
					right--;
					if(score >= 3) score -= 3;
				}
			}
			else
			{
				//message("ERROR: left= "+fold[left],+" right= "+fold[right]);
				return -2;
			}
			
			//
			if (score < 3 && score >= 0)
			{

				int len = left - left0;
				fold.replace(left0,len,len,'.'); 
				len = right0 - right;
				fold.replace(right+1,len,len,'.');
				left0 = left;
				right0 = right;
				//cout << fold << endl;
				score = -10;
			}
			//cout << left << "\t" << score << endl;
		}
	}
	
	//cout << stem_5p << endl;
	//cout << pairing << endl;
	//cout << stem_3p << endl;
	
	// starting from loop
	// first find the start of the stem
	while(fold[left0] == '.') left0++;
	while(fold[right0] == '.') right0--;
	// open terminal GU
	if(noClosingGU)
	{
		if(open_GU( seq,  fold, left0, right0))
		{
			left0++;
			right0--;
		}
	}
	int basal_start = left0;
	left = loops[0]-1; // point to (
	right = loops[1]+1; // point to )
	
	if(noClosingGU)
	{
		if(open_GU(seq,  fold, left, right))
		{
			left--;
			right++;
		}
	}
	
	
	left0 = left;
	right0 = right;
	score = 0;
	//cout << "basal start: "<<basal_start << endl;
	while(left >= basal_start)
	{
		if(fold[left] == '(' && fold[right] == ')' ) // match
		{
			a = seq.substr(left,1);
			b = seq.substr(right,1);
			ab = a+b;
			//stem_5p += a;
			//stem_3p += b;
			if (score < 0) score = 0; // if not started calculating score, now start
			if (ab == "AU" || ab == "UA") score += 2;
			else if (ab == "CG" || ab == "GC") score += 3;
			else score += 1; //UG/GU
			left--;
			right++;
		} 
		else
		{
			if (fold[left] == '.' && fold[right] == '.') // mismatch
			{
				//stem_5p += seq.substr(left,1);
				//stem_3p += seq.substr(right,1); 
				left--;
				right++;
				if(score >= 3) score -= 3;
			} else if (fold[left] == '.') // insertion on the left
			{
				// from: ((.(...)))
				// to:   (((...)))
				while(fold[left] == '.')
				{
					left--;
					if(score >= 3) score -= 3;
				}
			} else if(fold[right] == '.')
			{
				// from: (((...)..))
				// to:   (((...).))
				while(fold[right] == '.') 
				{
					if(score >= 3) score -= 3;
					right++;
				}
			}
			else
			{
				//message("ERROR: left= "+fold.substr(left,1)+" right= "+fold.substr(right,1) + ": id="+id);
				return -2;
			}
			if (score < 3 && score >= 0) // score > 0 means 
			{
				//cout << "left=" << left << endl;
				//cout << "left0=" << left0 << endl;
				//cout << "right=" << right << endl;
				//cout << "right0=" << right0 << endl;
				int len = left0 - left;
				fold.replace(left+1,len,len,'.'); 
				len = right - right0;
				fold.replace(right0,len,len,'.');
				left0 = left;
				right0 = right;
				score = -10; //
				//cout << "loop," << left << "," << score << endl;
				//cout << fold << endl;
			}
		}
	}
	
	while(fold[left0] == '.') left0--;
	while(fold[right0] == '.') right0++;
	// open terminal GU
	if(noClosingGU) open_GU( seq,  fold, left0, right0);

	return score/3;
	
}


/*
	
	
*/


/*

get stem length, including mismatches but not bulges

>hsa-let-7g
CGCUCCGUUUCCUUUUGCCUGAUUCCAGGCUGAGGUAGUAGUUUGUACAGUUUGAGGGUCUAUGAUACCACCCGGUACAGGAGAUAACUGUACAGGCCACUGCCUUGCCAGGAACAGCGCGCCAGCUGCCAAGUGGGG
...............((((((.((((.(((..(((((((.((((((((((((.................................)))))))))))).))))))).))).))))))).))).................
--CGCUCCGUUUCCUUU   -   A    A   UG       A            UGAGGGUCUAUGAUA 
                 UGC CUG UUCC GGC  AGGUAGU GUUUGUACAGUU               C
                 ||| ||| |||| |||  ||||||| ||||||||||||               C
                 GCG GAC AAGG CCG  UCCGUCA CGGACAUGUCAA               C
GGGGUGAACCGUCGACC   C   -    A   -U       C            UAGAGGACAUGGCCC 
*/
int hairpin_stem_length(string fold)
{
	// find the loop
	vector<int> loops = find_loops(fold);
	if(loops.size() != 2) 
	{
		if(loops.size() == 0) return -1; // no loop
		else return -2; // more than 1 loop
	} 
	
	// clean up structure by removing asymmetric bulge nucleotides until 35 bps from the loop, not counting bulges	
	int left = loops[0]-1; // point to (
	int right = loops[1]+1; // point to )
	int length = 0; // length from loop not counting bulges
	int last_mismatch_length = 0;
	while(1)
	{
		if(fold[left] == '(' && fold[right] == ')' ) // match
		{
			left--;
			right++;
			length++;
			last_mismatch_length = 0;
		} else if (fold[left] == '.' && fold[right] == '.') // mismatch
		{
			left--;
			right++;
			length++;
			last_mismatch_length++;
		} else if (fold[left] == '.') // insertion on the left
		{
			// ((.(...)))
			while(fold[left] == '.')
			{
				left--;
			}
			last_mismatch_length = 0;
		} else if(fold[right] == '.')
		{
			//(((...)..))
			while(fold[right] == '.') 
			{
				right++;
			}
			last_mismatch_length = 0;
		}
		else
		{
			message("ERROR:"+fold);
			return -3;
		}
		if(left < 0 || right >= fold.size()) break;
	}
	return length - last_mismatch_length;
}

string reverse_structure(string x)
{
    string y = x;
	int L = x.size();
	for(int i=0;i<L;i++)
		{
			int pos = L - i - 1;
			if(x[i] == '(') y[pos] = ')';
			else if (x[i] == ')') y[pos] = '(';
			else y[pos] = '.';
		}	
	return y;
}

int open_large_bulge(string &structure, int max_bulge)
{
	open_large_bulge_left_side(structure, max_bulge);
	string r = reverse_structure(structure);
	open_large_bulge_left_side(r, max_bulge);
	structure = reverse_structure(r);
	return count_pairs(structure);
}

// for a hairpin after removing branches and short stems
// if found internal loop/buldge longer than 9, open the basal part
// only look for 
int open_large_bulge_left_side(string &structure, int max_bulge)
{
	// ........((.((.............((.((.....))))))))....
	// ....((((((((.....)).)).............)).))........
	int start = -1;
	int npair = 1; // stem length
	int bulge_len = 0;
	// find the first (
	int k;
	for (k=0;k<structure.size();k++)
	{
		if (structure[k] == '(') break;	
	}
	start = k;
	for (int i=k+1;i<structure.size();i++)
	{
		if (structure[i] == ')') return 0;
	    if (structure[i] == '(')	
		{
			if (bulge_len > 0) // the end of a bulge
			{
				if(bulge_len <= max_bulge) 
				{
					bulge_len = 0;
					continue;
				}
				// remove stem left
				for(int j=0;j<i;j++)
				{
					if(structure[j] == '(') 
					{
						structure[j] = '.';
						//cout << structure << endl;
						npair ++;
					}
				}
				int left = 1;
				for(int j=i+1;j<structure.size();j++)
				{
					if (structure[j] == '(') left ++;
					else if (structure[j] == ')') 
					{
						if(left >= 0) left --;
						if(left <0) 
						{
							structure[j] = '.';
							//cout << structure << endl;
							npair--;
						}
						if (npair == 0) break;
					}
				}
				// start a new stem
				start = i;
				bulge_len = 0;
			}
		}
		else if (structure[i] == '.')
		{
			if (start > 0)
			{
				bulge_len ++;
			}
		}
	}
}


// count number of pairs in a structure
int count_pairs(string struc)
{
	int count = 0;
	for (int i =0;i<struc.size();i++)
		if (struc[i] == '(') 
			count ++ ;
	return count;
}

// count number of unpaired positions in a structure
int count_unpaired(string struc)
{
	int count = 0;
	for (int i =0;i<struc.size();i++)
		if (struc[i] == '.') 
			count ++ ;
	return count;
}

// NOT USED
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

// length of the stem next to loop, no mismatch allowed
int loop_adjacent_stem_length(string structure)
{
	int i = 0;
	while(structure[i] != ')') i++;
	int j = 0; // consecutive right stem length
	while(structure[i+j] == ')') j++;
	
	// left (
	int l=1;
	while(structure[i-l] != '(') l++;
	int k=0;
	while(structure[i-l-k] == '(')k++;
	
	return min(j,k);
}

// read a RNAfold file and remove short stem from each
// remove the fold if less than N pairs found
int remove_short_stem_from_file(string infile, string outfile, int max_paired_length/*=4*/, int min_pairs_left/*=25*/, int max_bulge/*=10*/,bool noClosingGU)
{
	ifstream fin;
	fin.open(infile.c_str());
	
	ofstream fout;
	fout.open(outfile.c_str());

	string line, id, seq, fold;
	vector<string> flds;
	
	int n_removed = 0;
	
	while(fin)
	{
		getline(fin,line);
		if(line.length()==0) continue;
		if (line[0] == '.' || line[0] == '(' ) //structure line 
		{
			//remove newline
			flds = string_split(line," ");
			fold = flds[0];
			fold.erase(fold.find_last_not_of(" \n\r\t")+1);
			//remove_all_branches(fold);
			//cout << seq << "|" << endl;
			//cout << fold << "|" << endl;
			//cout << id << endl;
			//cout << seq << endl;
			//cout << fold << endl;
			
			// to handle terminal short stems, add ..... to each side
		    fold = "...................." + fold + "...................."; 
			remove_short_stem_new(fold);
			fold = fold.substr(20,fold.size()-40);
			
			//cout << seq.size() << "," << fold.size() << endl;
			int score = trim_hairpin(seq,fold,noClosingGU);
			//cout << fold << "|" << endl;
			
			//open_large_bulge(line, max_bulge);
			//cout << score << endl;
			int npair = count_pairs(fold);
			//cout << npair << endl;
			if ( npair >= min_pairs_left) 
			{
				//if(loop_adjacent_stem_length(line) == 35 && npair == 35)
				//{
					fout << id << endl << seq << endl << fold << endl;
				//}
			}
			else
			{
				n_removed ++;
			}
		}
		else if (line[0] == '>') id = line;
		else seq = line.erase(line.find_last_not_of(" \n\r\t")+1);
	}
	message(" - discard "+to_string(n_removed) + " too short longest stems after filtering");
	return n_removed;
}

// from RNALfold make fasta file for each hairpin, optionally extend some bases
// can filter by MFE
void RNALfold_to_RNAfold(string infile, string outfile, int min_hairpin_length/*=50*/,int min_pairs_left/*=25*/, int ext/*=20*/, int mfe/*=-30*/)
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
	
	int n_invalid = 0;
	int n_too_short = 0;
	int n_no_ext = 0;
	int n_not_stable = 0;
	int n_total = 0;
	int n_main_stem_too_short = 0;
	int n_valid = 0;
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
			n_total += structures.size();
			for(int i =0;i<structures.size();i++) // for each structure
			{
				// find first space, left is structure, right is score
				size_t found = structures[i].find(" "); // ((..)). (-10.60) 2799 z= -2.699
				// cout << structures[i].substr(0,found) << "," << structures[i].substr(found) << endl;
				//cout << "1" << endl;
				string struc = structures[i].substr(0,found);
				if (struc.size() < min_hairpin_length) 
				{
					n_too_short ++;
					continue; // filter by hairpin size
				}
				//cout << "2" << endl;
				
				string tmp = structures[i].substr(found+2); // -10.60) 2799 z= -2.699
				split( flds, tmp, is_any_of(")z"), token_compress_on );
				
				// mfe filter
				if (stof(flds[0]) > mfe)
				{
					n_not_stable ++;
					continue;
				}
				
				int pos = stoi(trim_copy(flds[1])) - 1; // 1-based to 0-based
				if (pos<ext || pos + struc.size() + ext >= line.size() ) 
				{
					n_no_ext++;
					continue; // skip
				}
				
				// remove branches
				remove_all_branches(struc);
				// to handle terminal short stems, add ..... to each side
				struc = extend_structure + struc + extend_structure;
				remove_short_stem_new(struc);
				int npair = count_pairs(struc);
				if ( npair < min_pairs_left) 
				{
					n_main_stem_too_short++;
					continue;
				}
						
				string name = seqID + ":" + to_string(n_valid) + ":" + flds[0];
				erase_all(name," ");
				string seq = to_upper(line.substr(pos - ext , struc.size()));
				if (valid_sequence(seq,"ACGU")) 
				{
					fout << name << endl << seq << endl << struc << endl;
					n_valid ++;
				}
				else n_invalid ++;
			
			}
			getline(fin,line); // skip total score line
		}
	}
	fin.close();
	fout.close();
	message(" - loaded "+to_string(n_total) + " structures in total");
	message(" - discard " +to_string(n_too_short) + " hairpins shorter than "+to_string(min_hairpin_length));
	message(" - discard " +to_string(n_not_stable) + " hairpins with MFE > "+to_string(mfe));
	message(" - discard " +to_string(n_no_ext) + " structures can not be extended");
	message(" - discard " +to_string(n_main_stem_too_short) + " structures whose longest stem is too short");
	message(" - discard " +to_string(n_invalid) + " structures containing non ACGU letters");
	message(" - output "+to_string(n_total-n_too_short-n_not_stable-n_no_ext-n_invalid-n_main_stem_too_short) + " structures");
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
