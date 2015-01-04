#ifndef __STAT_H__
#define __STAT_H__

#include <vector>
#include <array>
#include <algorithm>


using namespace std;

double p_value_for_z_score(double z);

// only enrichment
double binom_test(int trials,int success,double success_fraction);

double hypergeometric_test(unsigned k, unsigned r, unsigned n, unsigned N);

array<double,2> Mann_Whitney_U_test(vector<int> ranks, int N);

array<double,2> two_samples_t_test_equal_sd(double Sm1, double Sd1, unsigned Sn1, double Sm2, double Sd2, unsigned Sn2);

array<double,2> two_samples_t_test_unequal_sd(double Sm1, double Sd1, unsigned Sn1, double Sm2, double Sd2, unsigned Sn2);

array<double,6> t_test(vector<double> x, vector<double> y, bool equal_var=false);

template<typename T>
T sum(vector<T> v)
{
	return accumulate(v.begin(), v.end(), 0.0);
}

template<typename T>
double mean(vector<T> v)
{
    return double(sum(v))/v.size();
}

template<typename T>
double var(vector<T> v)
{
	vector<double> diff(v.size());
	transform(v.begin(), v.end(), diff.begin(),
    bind2nd(minus<double>(), mean(v)));
	double sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
	double var = sq_sum / v.size();
	return var;
}

// correlation
template<typename T>
double cor(vector<T> x, vector<T> y)
{
	double mean_x = mean(x);
	double mean_y = mean(y);
	
	double cov = 0;
	for (int i=0;i<x.size();i++)
		cov += (x[i] - mean_x)*(y[i]- mean_y);
	
	return cov / x.size() / sqrt(var(x)*var(y));
}


template<typename T>
double sd(vector<T> v)
{
	return sqrt(var(v));
}

#endif
