#include "utility.h"
#include "iostream"
#include <cmath>

using namespace std;

void help()
{
    string str =
	"\n"
    "mRNA_loop: simulating mRNA closed-loop formation\n"
	"	- Xuebing Wu (wuxb07@gmail.com)\n"
    "\n"
    "Usage: mRNA_loop -L mRNA_length -TL polyA_tail_length -T 1000000\n"
    "\n"
    "Options:\n"
    "\n"
    "   -L  NUM   mRNA length (default=1000 nucleotides)\n"
    "   -TL  NUM  polyA tail length (default=100 nucleotides)\n"
	"   -T NUM    longest simulation time (1 step = 1 nucleotide, default 1000000)\n"
	"   -R NUM    integer seed for random generator\n"
    "\n";

    cerr << str;

    exit(0);
}

int loop2d(double L, double minD, double P, int max_t, int seed)
{
	// initialization
	srand (time(NULL)+seed);
	double x=0;
	double y=sqrt(2*P*L);
	int t = 0;
	double x1 = 0;
	double y1 = 0;
	double min_d2 = minD*minD;
	double L2 = L*L;
	//cout << y << endl;
	//cout << min_d2 << endl;
	double pi180 = 3.14159265/180;
	while (x*x+y*y > min_d2 && t < max_t)
	{
	    t = t + 1;
	    x1 = x;
		y1 = y;
		double r = rand()%360 * pi180;
		x1 = x + cos(r);
		y1 = y + sin(r); 
	    if (x1*x1+y1*y1 < L2) {x = x1;y = y1;}
		//cout << t <<"\t" << x << "\t" << y << endl;
	}
	return t;
}

int loop3d(double L, double minD, double P, int max_t, int seed)
{
	// initialization
	srand (time(NULL)+seed);
	double x=0;
	double y=sqrt(2*P*L);
	double z = 0;
	int t = 0;
	double x1 = 0;
	double y1 = 0;
	double z1 = 0;
	double min_d2 = minD*minD;
	double L2 = L*L;
	double pi180 = 3.14159265/180;
	while (x*x+y*y+z*z > min_d2 && t < max_t)
	{
	    t = t + 1;
	    x1 = x;
	    y1 = y;
		z1 = z;
		double r1 = rand()%360 * pi180;
		double r2 = rand()%360 * pi180;
		x1 = x + sin(r1)*cos(r2);
		y1 = y + sin(r1)*sin(r2);
		z1 = z + cos(r1);
	    if (x1*x1+y1*y1+z1*z1 < L2) {x = x1;y = y1;z=z1;}
	}
	return t;
}




int main(int argc, char* argv[]) {

    // default
	double P = 0.79; // persistence length of mRNA
	double nm_per_nt = 0.3; // nm per nucleotide
    double L=1000;
	double taillen=100;
	double radius_cap = 4.5;
	int max_t = 1000000;
	int seed = 1;

	if (argc < 2) help();
	
    // parse arguments
    string str;
    for (int i = 1; i < argc; i++) { 
        if (i != argc) { 
            str=argv[i];
            if (str == "-L") {  
                L = stof(argv[i + 1]);
                i=i+1;
            } else if (str == "-TL") { 
                taillen = stof(argv[i + 1]);
                i=i+1;
            } else if (str == "-T") { 
                max_t = stoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-R") { 
                seed = seed + stoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-h" || str == "--help") { 
                help();
            } else {
                cout << "unknown option: "+str << endl;   
                help();        
            }
        }
    }   
	
	double radius_pA = taillen / 250 * 5.9;
	double minD = (radius_pA + radius_cap) / nm_per_nt;
	
	//cout << "minD = " << minD << endl;
	
	cout << loop2d(L,minD,P,max_t,seed) << "\t" << loop3d(L,minD,P,max_t,seed) << endl;
	
	
    return 0;
}
