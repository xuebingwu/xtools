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

#include "container.h"

// for position matrix
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;


// print a matrix
void print_matrix(boost::numeric::ublas::matrix<double> m)
{
    for (unsigned i = 0; i < m.size1(); ++ i)
    {
        for (unsigned j = 0; j < m.size2(); ++ j)
        {
            cout << m(i,j) << "\t";
        }
        cout << endl;
    }
}


