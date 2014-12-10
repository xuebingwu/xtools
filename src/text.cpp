#include <fstream>

#include <stdio.h>
#include <stdlib.h>

#include "utility.h"
#include "text.h"

#include <map>
#include <iostream>

// subtract file 2 from file 1, based on shared key column
void intersectTab(string file1, string file2, string outputfile, unsigned col1/*=0*/, unsigned col2/*=0*/, bool subtract/*=false*/)
{
    // read file2 first 
    set<string> data2;
    vector<string> flds;
	string line;
	
	ifstream f2;
	f2.open(file2.c_str());
    
	while(f2)
	{
		getline(f2,line);
		if(line.length() == 0) continue;
		flds = string_split(line,"\t");
		data2.insert(flds[col2]);
	}
	f2.close();
	
	ifstream f1;
	f1.open(file1.c_str());
	
	ofstream out;
	out.open(outputfile.c_str());
    
	unsigned total = 0;
	unsigned comm = 0;
	unsigned diff = 0;
	
	while(f1)
	{
		getline(f1,line);
		if(line.length() == 0) continue;
		flds = string_split(line,"\t");
		total ++;
		//cout << flds[col1] << endl;
		if (data2.find(flds[col1]) != data2.end()) // shared
		{
		    if(subtract == false) 
			{
				out << line << endl; // output common lines
				comm ++;
			}
			
		}
		else if (subtract) // not shared and output difference
		{
			out << line << endl;
			diff ++;
		}
	}          
	
	if (subtract) cout << diff << " of " << total << " lines unique to " << file1 << endl;
	else cout << comm << " of " << total << " lines shared" << endl;
	
    f1.close();  
	out.close();         
}

void mergeTab(string file1, string file2, string outputfile, unsigned col1/*=0*/, unsigned col2/*=0*/, bool header/*=false*/, string fill/*="None"*/)
{
	/*
    // add file2 to file1 side by side, using column col1 in file1 and column col2 in file2 to match rows
    // col1 and col2 are 0-based, but the program may be 1-based
    // combine header directly if header = true
    // for lines in file1 but not file2, fill with 'fill' if fill != "none"
	*/
	
    // read file2 first 
    map<string,string> data2;
    vector<string> flds;
	string line;
	
	ifstream f2;
	f2.open(file2.c_str());

	string header2;
	if(header) getline(f2,header2);
    
	while(f2)
	{
		getline(f2,line);
		if(line.length() == 0) continue;
		flds = string_split(line,"\t");
		data2[flds[col2]] = line;
	}
	f2.close();
	
    unsigned n2 = data2.size();
    unsigned nc = flds.size();
	
	//cout << n2 << "," << nc << endl;
	
    // match and add to file1
	ifstream f1;
	f1.open(file1.c_str());
	
	ofstream out;
	out.open(outputfile.c_str());
	
	string header1;
	if(header) 
	{
		getline(f1,header1);
		header1 = header1 + "\t" +header2;
	}
    //cout << header1 << endl;
	
    unsigned n1 = 0; // lines in file 1
    unsigned n3 = 0; // common lines
    
    string fillline = fill;
	if (fill != "none") 
	{
		for (int i=1;i<nc;i++) fillline += "\t" + fill;
	}

	while(f1)
	{
		getline(f1,line);
		if(line.length() == 0) continue;
		flds = string_split(line,"\t");
		n1 ++;
		//cout << flds[col1] << endl;
		if (data2.find(flds[col1]) != data2.end())
		{
			n3 ++;
			out << line << "\t" << data2[flds[col1]] << endl;
		}
		else if (fill != "none")
		{
			out << line+ "\t" + fillline << endl;
		}
	}          
    f1.close();      
	out.close();     
}
		
void remove_duplicates(string input, string output, int col, int max, string sort_opts="")
{
    if (sort_opts.size()>0)
    {
        string cmd = "sort "+sort_opts+" "+input+" > xxx.tmp";
        system(cmd.c_str());
        input="xxx.tmp";
    }


    col = col - 1;
    ifstream fin;
    ofstream fout;
    fin.open(input.c_str());
    fout.open(output.c_str());
    string line;
    vector<string> flds;

    string curr = "";
    int k = max;

    while(fin)
    {
        getline(fin,line);
        if (line.length() == 0)
            continue;
        flds = string_split(line,"\t");
        if (curr != flds[col])
        {
            fout << line << endl;;
            k = max - 1;
            curr = flds[col];
        }
        else if (k>0)
        {
            fout << line << endl;
            k--;
        }
    }

    if(sort_opts.size()>0) system("rm xxx.tmp");

}

