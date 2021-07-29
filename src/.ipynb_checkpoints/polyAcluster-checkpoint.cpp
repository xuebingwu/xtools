#include "utility.h"
#include "iostream"
#include "interval.h" 
#include "stat.h"  


using namespace std;

void help()
{
    string str =
	"\n"
    "polyAcluster: clustering polyA sites\n"
	"	- Xuebing Wu (wuxb07@gmail.com)\n"
    "\n"
    "Usage: polyAcluster -i input -o output [options]\n"
    "\n"
    "Options:\n"
    "\n"
    "   -i  input   input file, 6 column bed format \n"
    "   -o  output  prefix for output files\n"
    "   -d  number  max offset to peak site in a cluster.\n"
	"               default is to determine d automatically\n"
	"   -s          sort input\n"
	"   -5          use 5' end of each site in the input\n"
	"Parameters for determining d:\n"
	"   -r  number  to search d within (0,r). Default=50\n"
	"   -c  number  min read count for a site to be used for learning d (default=1000)\n"
	"   -f  number  max read count ratio for a site in (0,r) w.r.t the peak (default=1.0)\n"
	"   -p  number  min relative frequency to be considered as offset cleavage (default=0.01)\n"
    "\n";
	
    cerr << str;
	

    exit(0);
}



int main(int argc, char* argv[]) {

    // default parameters
    string input="input";
    string output="output";
	
	// is input sorted 
	bool sorted = true;
	
	// 
	bool start = false;
    int d= -1;
	string pos = "3";
	int min_count = 1000;
	int radius=50;
	double max_frac = 1.0;
	double min_freq = 0.01;

	if (argc < 2) help();
	
    // parse arguments
    string str;
    for (int i = 1; i < argc; i++) { 
        if (i != argc) { 
            str=argv[i];
            if (str == "-i") {  
                input = argv[i + 1];
                i=i+1;
            } else if (str == "-o") { 
                output = argv[i + 1];
                i=i+1;
            } else if (str == "-d") { 
                d = atoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-r") { 
                radius = atoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-c") { 
                min_count = atoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-f") { 
                max_frac = atof(argv[i + 1]);
                i=i+1;
            } else if (str == "-p") { 
                min_freq = atof(argv[i + 1]);
                i=i+1;
            } else if (str == "-5") {  
                start = true;
				pos = "2";
            } else if (str == "-s") { 
                sorted = false;
            } else if (str == "-h" || str == "--help") { 
                help();
            } else {
                message("unknown option: "+str);   
                help();        
            }
        }
    }   
	
	// learn d automatically
	if(d < 0){
		message("learning cluster size from abundant sites...");
		d = single_cleavage_sites(input, output+".cluster.size", min_count, radius,max_frac,min_freq); 
		message("cluster size determined to be "+to_string(d));
	}
		
	message("split input by strand");
	// split by strand 
	system_run("awk '$6==\"+\"' "+input + " > " + output+".input+");
	system_run("awk '$6==\"-\"' "+input + " > " + output+".input-");
	
	if (sorted == false)
	{
		message("sort plus strand");
		system_run(" sort -k1,1 -k2,2n "+output+".input+ > "+output+".input+.sorted");
		system_run(" mv "+output+".input+.sorted "+ output+".input+");
		message("sort minus strand"); 
		system_run(" sort -k1,1 -k2,2n "+output+".input- > "+output+".input-.sorted");
		system_run(" mv "+output+".input-.sorted "+ output+".input-");
	}
	
	message("find peaks on plus strand");
	int total_plus = cluster_sites(output+".input+", output+".cluster+", start, d);
	message("find peaks on minus strand"); 
	int total_minus = cluster_sites(output+".input-", output+".cluster-", start, d);
	system_run(" cat "+output+".cluster? > "+output+".cluster"); 	
	message(to_string(total_plus+total_minus)+" sites left");
    return 0;
}
