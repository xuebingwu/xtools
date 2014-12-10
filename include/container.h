#ifndef __CONTAINER_H__
#define __CONTAINER_H__

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


// merge two vectors
template <typename T>
std::vector<T> operator+(const std::vector<T> &A, const std::vector<T> &B)
{
    std::vector<T> AB;
    AB.reserve( A.size() + B.size() );                // preallocate memory
    AB.insert( AB.end(), A.begin(), A.end() );        // add A;
    AB.insert( AB.end(), B.begin(), B.end() );        // add B;
    return AB;
}

template <typename T>
std::vector<T> &operator+=(std::vector<T> &A, const std::vector<T> &B)
{
    A.reserve( A.size() + B.size() );                // preallocate memory without erase original data
    A.insert( A.end(), B.begin(), B.end() );         // add B;
    return A;                                        // here A could be named AB
}

// print a matrix
void print_matrix(boost::numeric::ublas::matrix<double> m);

// sum of two vectors
template<typename T>
vector<T> sum(vector<T> a, vector<T> b)
{
    vector<T> res;
    for (T i =0;i< a.size();i++) res.push_back(a[i]+b[i]);
    return res;
}



template<class key,class val>
void PrintMap(map<key,val> m)
{// print out a map
  for (typename map<key,val>::iterator it = m.begin();it!=m.end(); it++)
  {
    cout << (*it).first << "\t=>\t" <<(*it).second<<endl;
  }
}

#endif

