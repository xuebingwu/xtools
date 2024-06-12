#include <string>
#include <vector>

#include <iostream>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>

#include "utility.h"


using namespace std;

// load the top fraction of a list of ranked values, return the total number of data
// col: start 0
int load_sorted_table(string filename, int col, vector<double> &table, double fraction)
{
	//bool descending = true;
	
	ifstream fin;
	fin.open(filename.c_str());
	
	string line;
	vector<string> flds;
		
	// save line 1 data
	getline(fin,line);
	flds = string_split(line,"\t");
	//double first_data = stof(flds[col]);
	// save line 2 data
	getline(fin,line);
	flds = string_split(line,"\t");
	//double second_data = stof(flds[col]);
	// cout << "first two data points = " << first_data << "," << second_data << endl;
	//if( first_data < second_data ) descending = false;
	
	int total = 1;
	
	// count total number of lines
	while(fin)
	{
		getline(fin,line);
		if(line.length()==0) continue;
		total ++;
	}
	fin.close();
	
	// cout << total << " data points in total" << endl;

	// the line to stop
	int line_to_stop = total * fraction;

	// cout << "keep the top fraction " << fraction << endl;
	// cout << "keep the top lines " << line_to_stop << endl;

	int n = 0;
	fin.open(filename.c_str());
	while(fin)
	{
		getline(fin,line);
		if(line.length()==0) continue;
		flds = string_split(line,"\t");
		double v = stof(flds[col]);
		table.push_back(v);
		n++;
		if(n == line_to_stop) break;
	}
	fin.close();
	// cout << "finished" << endl;
	return total;
}

// 
double lookup_table(double x, vector<double> table, int total)
{
	int L = table.size();
	
	bool descending = table[0] > table[1];
	
	if ((descending == true && x < table[L-1]) || (descending == false && x > table[L-1])) return 1.0;

	int i=0;

	if (descending) 
	{
		while(x<=table[i]) i++;
	}
	else
	{
		while(x>=table[i]) i++;
	}
	return (i+1.0)/total;
}

void calculate_p_from_file(string infile, string outfile, int col, vector<double> table, int total, double fraction)
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
		if(line.length() == 0) continue;
		flds = string_split(line,"\t");
		double p = lookup_table(stof(flds[col]),table,total) / fraction;
		if(p<1) fout << line << "\t" << p << endl; 
	}
	fin.close();
	fout.close();
}

int main(int argc, char* argv[])
{
	string inputfile = argv[1];
	int col1 = atoi(argv[2]);
	string tablefile = argv[3];
	int col2 = atoi(argv[4]);
	double fraction = atof(argv[5]);
	string outputfile = argv[6];

	vector<double> table;
	int total = load_sorted_table(tablefile, col2, table, fraction);
	
	// cout << table.size() << endl;
	// cout << total << " values in table" << endl;
	
	calculate_p_from_file(inputfile, outputfile, col1, table, total, fraction);

	return 0;
}