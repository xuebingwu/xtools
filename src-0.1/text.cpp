#include <fstream>

#include <stdio.h>
#include <stdlib.h>

#include "utility.h"
#include "text.h"

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

