#include "utility.h"
#include "iostream"

using namespace std;

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
		
		// if found right ...((...
		if(found_right != std::string::npos)
		{
			cout << "right\t"<< L << "\t" << structure  << endl;
			
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
						cout << "right\t"<< L << "\t" << structure << "\t" << found_right-i << "\t" << found_right << endl;
						left --;
					}
					if (left == 0) break;
				} 
			}
		}
		
		found_left = structure.find(pattern_left);
		
		if(found_left != std::string::npos)
		{
			cout << "left\t"<< L << "\t" << structure  << endl;
			 
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
						cout << "left\t"<< L << "\t" << structure << "\t" << i << "\t" << found_left << endl;
						
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


void help()
{
    string str =
	"\n"
    "Example: an example program\n"
	"	- Xuebing Wu (wuxb07@gmail.com)\n"
    "\n"
    "Usage: example -i input -o output -n number\n"
    "\n"
    "Options:\n"
    "\n"
    "   -i  input   input file\n"
    "   -o  output  output file\n"
    "   -n  number  number\n"
    "\n";

    cerr << str;

    exit(0);
}

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

int main(int argc, char* argv[]) {

	//string structure="......(((((.(((......)))....))))).....(((...)))..((((...))))..";
	//string structure = ".((((((((((.......(((.........)))....(....).((..))...)))))))))).";
	//string structure ="..(((...(((((....).))))...)))...";
	string structure ="(((((.......((((((((......)))).))))....((((((.....((((((((((((......))))))))))))))))))...)))))..((...))...";
	cout << structure << endl;
    keep_longest_stem(structure);	

    return 0;
}
