#ifndef __INTERVAL_H__
#define __INTERVAL_H__

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <iomanip>
#include <array>        // std::array

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>

// for position matrix
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;

int single_cleavage_sites(string inputfile, string output, int min_count, int window, double max_frac, double min_freq);
int cluster_sites(string inputfile, string outputfile,  bool start, int d);
void flipStrand(string inputfile, string outputfile);
void interval_annotation(string inputfile, string outputfile, vector<string> annotationfiles, vector<string> names);

#endif
