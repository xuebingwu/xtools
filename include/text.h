
#ifndef __TEXT_H__
#define __TEXT_H__

#include <string>
#include <vector>

using namespace std;

void intersectTab(string file1, string file2, string outputfile, unsigned col1=0, unsigned col2=0, bool subtract=false);

void mergeTab(string file1, string file2, string outputfile, unsigned col1=0, unsigned col2=0, bool header=false, string fill="none");

void remove_duplicates(string input, string output, int col, int max, string sort_opts);

#endif

