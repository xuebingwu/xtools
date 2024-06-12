#include "utility.h"
#include "iostream"
#include "text.h"
#include <boost/algorithm/string.hpp>
using namespace boost;

using namespace std;

void help()
{
    string str =
	"\n"
    "cross_validation: cross-validadtion for kpLogo2 models\n"
	"	- Xuebing Wu (wuxb07@gmail.com)\n"
    "\n"
    "Usage: cross_validation -data data -train_cmd cmd -test_cmd cmd -nfold n\n"
	"       example : cross_validation -data kpLogo2_test_Doench2014.txt  -train_cmd \"kpLogo2 data.train.@ -seq 1 -weight 2  -o output.@ -upto $k -shift_max $shift -pair 2>> log \" -test_cmd \"kpLogo2 data.test.@ -seq 1 -weight 2 -pCutoff_B 2 -predict output.@ -pair 2>> log\" -nfold 10 \n\n"
    "Options:\n"
    "\n"
    "   -data      data   data file\n"
    "   -train_cmd cmd    cmd for training, quoted, use @ for i\n"
    "   -test_cmd  cmd    cmd for testing, quoted, use @ for i\n"
    "   -nfold     n      n fold cross-validation\n"
    "\n";

    cerr << str;

    exit(0);
}

int main(int argc, char* argv[]) {

    // default
    string data, train_cmd, test_cmd;
    int nfold=5;

	if (argc < 2) help();
	
    // parse arguments
    string str;
    for (int i = 1; i < argc; i++) { 
        if (i != argc) { 
            str=argv[i];
            if (str == "-data") {  
                data = argv[i + 1];
                i=i+1;
            } else if (str == "-train_cmd") {  
                train_cmd = argv[i + 1];
                i=i+1;
            } else if (str == "-test_cmd") { 
                test_cmd = argv[i + 1];
                i=i+1;
            } else if (str == "-nfold") { 
                nfold = atoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-h" || str == "--help") { 
                help();
            } else {
                message("unknown option: "+str);   
                help();        
            }
        }
    }   

	message("spliting data to training and test sets...");
	nfold = split_file_for_cross_validation( data,  "data", nfold);

	string cmd;
	for (int i=0;i<nfold;i++)
	{
		message("cross-validation: "+to_string(i+1)+" of "+to_string(nfold));
		cmd = replace_all_copy(train_cmd, "@", to_string(i));
		message(cmd);
		system_run(cmd);
		cmd = replace_all_copy(test_cmd, "@", to_string(i));
		message(cmd);
		system_run(cmd);
	}
	
	/* 
	
	for k in 1 2 3 4
	do 
	    for shift in 0 1 2
	    do
	    echo $k, $shift
	    rm log; cross_validation -data kpLogo2_test_Doench2014.txt  -train_cmd "kpLogo data.train.@ -o output.@ -max_k $k -max_shift $shift -weighted 2>> log " -test_cmd "kpLogo data.test.@ -predict output.@ -weighted 2>> log" -nfold 10
	    rm tmp
	    cat output.*.score | grep -v inf > tmp
	    Rscript cor.r
	    done 
	done

#cor.r
x=read.table('tmp',header=F)
s1 = summary(lm(x[,2]~x[,3]))
s2 = summary(lm(x[,2]~x[,4]))
s3 = summary(lm(x[,2]~x[,3]*x[,4]))
c(s1$r.squared,s2$r.squared,s3$r.squared)
	
	
	*/
	

    return 0;
}
