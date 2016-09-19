#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <algorithm>

#include "utility.h" 
#include "stat.h"

#include "text.h"
#include "time.h"
#include <map>
#include <iostream>

#include "interval.h" 
#include "text.h" 




// flip strand
void flipStrand(string inputfile, string outputfile){
	ifstream fin(inputfile.c_str());
    ofstream fout(outputfile.c_str());

    string line;
    vector<string> flds;

    while(fin)
    {
        getline(fin,line);
        if (line.length() == 0)
            continue;

        line.erase(line.find_last_not_of(" \n\r\t")+1);

        if(line[0] == '#') continue;

        flds = string_split(line);
		if (flds[5]=="+") flds[5] = "-";
		else flds[5] = "+";
		fout << to_string(flds,"\t") << endl;
	}
	fin.close();
	fout.close();
}

// find sites with >= min_count reads, accounting for >min_fraction of reads in a window
// not start or end
int single_cleavage_sites(string inputfile, string output,int min_count, int radius, double max_frac, double min_freq){
	// strand??
	// find all sites with min_count reads 
	system_run(" awk '$5 >= " + to_string(min_count) + " ' " + inputfile+ " > " + inputfile+".abundant" );
	// generate a window around those sites
	system_run(" awk '{print $1\"\t\"$2-"+to_string(radius)+"\"\t\"$3+"+to_string(radius)+"\"\t\"$4\"\t\"$5\"\t\"$6}' "+inputfile+".abundant > "+inputfile+".abundant.window");
	// find all other sites in those windows, calculate fraction 
	system_run(" bedtools intersect -a "+inputfile+".abundant.window -b " + inputfile + " -s -wo > "+inputfile+".nearby  " );
	// load pos fracs
	map<string ,vector<double> > fractions; 
	string filename=inputfile+".nearby";
	ifstream fin(filename.c_str());
    string line;
    vector<string> flds;
	vector<string> cur_pos_rel;
	vector<double> cur_frac;	
	// read first line
    getline(fin,line);
	line.erase(line.find_last_not_of(" \n\r\t")+1);
    flds = string_split(line);
	string cur_chr = flds[0];
	string cur_pos = flds[1];
	string rel_pos = to_string(stoi(flds[7])-stoi(flds[1])-radius);
	if (rel_pos != "0") {
		cur_pos_rel.push_back(rel_pos);
		cur_frac.push_back(stof(flds[10])/stof(flds[4]));
	}
	
	int total_eligible = 0;
    while(fin)
    {
        getline(fin,line);
        if (line.length() == 0)
            continue;
		
		line.erase(line.find_last_not_of(" \n\r\t")+1);
		
        flds = string_split(line);
		
		string rel_pos = to_string(stoi(flds[7])-stoi(flds[1])-radius);
		
		if (cur_chr != flds[0] || cur_pos != flds[1]) {
			// a new window 
			// first check if the current window is eligible
			if (max(cur_frac) < max_frac )
			{	
				total_eligible ++;
				for(int i=0;i<cur_frac.size();i++)
				{
					fractions[cur_pos_rel[i]].push_back(cur_frac[i]);
				}
			}
			// set as current window
			cur_chr = flds[0];
			cur_pos = flds[1];
			cur_pos_rel.clear();
			cur_frac.clear();
		} else {
			if (rel_pos != "0") {
				cur_pos_rel.push_back(rel_pos);
				cur_frac.push_back(stof(flds[10])/stof(flds[4]));
			}
		}
	}
	fin.close();
	// calcualte mean
	vector<double> avgs(radius*2+1,0);
	vector<int> n(radius*2+1,0);
	for (int i= -radius;i<= radius;i++)
	{
		n[i+radius] = fractions[to_string(i)].size();
	    avgs[i+radius] = sum(fractions[to_string(i)])/total_eligible;
	}
	n[radius] = total_eligible;
	avgs[radius] = 1.0;
	
	// save and plot
	filename = output + ".summary.txt";
	ofstream out(filename.c_str());
	out << "position" << '\t' << "sites" << "\t" << "relFreq" << endl;
	for(int i=0;i<n.size();i++)
		out << i-radius << "\t" << n[i] << "\t" << avgs[i]  << endl;
	out.close();
	
	// plot 
	string rscript = ""
	"x=read.table('"+filename+"',header=T)\n"
	"pdf('"+output+".pdf')\n"
	"par(cex=1.5)\n"
	"plot("+to_string(-radius)+":"+to_string(radius)+",x[,3],type='h',lwd=2,xlab='distance to major site',ylab='relative frequency',col='blue')\n"
	"abline(h=0.01,col='gray',lty=2)\n"
	"plot("+to_string(-radius)+":"+to_string(radius)+",x[,2],type='h',lwd=2,xlab='distance to major site',ylab='number of sites',col='blue')\n"
	"dev.off()\n";
	R_run(rscript);
	
	// determine d
	int d=0;
	while(radius > d  && avgs[radius+d] > min_freq && avgs[radius-d] > min_freq) d++;
	return d-1;
}


// cluster sites
// input file should be sorted
int cluster_sites(string inputfile, string outputfile,  bool start, int d){
	int total = 0;
	
	ifstream fin(inputfile.c_str());
	ofstream fout(outputfile.c_str());
	

    string line;
    vector<string> flds;
		
	string cur_chr=""; // currnt chromosome
	int cur_pos = -1;
	int cur_count = 0;
	string cur_line = "";
	
	int pc = 2; // default use column 3 as position, i.e. 3' end
	if(start) pc = 1;
	
    while(fin)
    {
        getline(fin,line);
        if (line.length() == 0)
            continue;
		
		line.erase(line.find_last_not_of(" \n\r\t")+1);
		
		if(line[0] == '#') continue;

        flds = string_split(line);
		int pos = stoi(flds[pc]); 
		int count = stoi(flds[4]);
		if (flds.size() < 6) message("too few columns in line: "+line);
		else
		{
			if (cur_chr == flds[0] && pos < cur_pos + d ){
				// the same chr within distance
				if (cur_count < count){
					cur_pos = pos;
					cur_count = count;
					cur_line = line;
				}
			} 
			else { // the next site is too far away 
					// output current peak
				fout << cur_line << endl;
				total ++;
				// goes to the next site
				cur_chr = flds[0];
				cur_pos = pos;
				cur_count = count;
				cur_line = line;
			}
		}
	}
	fin.close();
	fout.close();
	return total;
}






// assign intervals to annotations
// inputfile: BED file, a list of intervals
// annotationfiles: a list of annotation BED file names separated by ,
// names: a list of names for the annotation files
// min_frac: min fraction overlap
void interval_annotation(string inputfile, string outputfile, vector<string> annotationfiles, vector<string> names){
		
	// make sure the number of annotation files match the number of names
	if (annotationfiles.size() != names.size()) {
		message("the number of annotation files should be the same as the number of names");
		exit(1);
	}

	// outputfile will be rewrite
	system_run("if [ -f "+outputfile+" ]; then rm "+outputfile +"; fi ");
	
	// initialization
	system_run("cp "+inputfile+" OTHER");
	
	message("intersecting with...");
	for (int i=0;i < annotationfiles.size();i++ ){
		system_run("mv OTHER TOTAL");
		message(" - "+names[i]);
		system_run(" bedtools intersect -u -a TOTAL -b " + annotationfiles[i] + " -s | awk '{print $0\"\t"+names[i]+"\"}' >> "+outputfile);
		system_run(" bedtools intersect -v -a TOTAL -b " + annotationfiles[i] + " -s > OTHER");
	}
	system_run(" cat OTHER | awk '{print $0\"\tOther\"}' >> "+outputfile); 
	system_run(" rm OTHER TOTAL");
	
	// total number of lines
	int total = count_lines(outputfile);
	
	// calculate fraction
	system_run(" awk '{print $NF}' "+outputfile+" | uniq -c | awk '{print $2\"\t\"$1\"\t\"$1/"+to_string(total)+"}' > "+outputfile+".summary");
	
	message("Summary (also see "+outputfile+".summary)");
	system_run("cat "+outputfile+".summary");
}



