#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <sys/time.h>
#include <array>

using namespace std;

int char2int(char ch){
        switch (ch) {
                case 'A' : return 0;
                case 'C' : return 1;
                case 'G' : return 2;
                case 'T' : return 3;
                default  : return 4;
                }
}

char complement(char ch){
	switch (ch)
	{
	case 'A':return 'T';
	case 'C':return 'G';
	case 'G':return 'C';
	case 'T':return 'A';
	default:return 'N';
	}
}

string reverseComplement(string seq){
	int L = seq.length();
	string rc (L,'0');
	for (int i=0;i<L; i++)
	{
		rc[L-i-1] = complement(seq[i]);
	}
	return rc;
}


int **initializeMatrix(int ROWS,int COLUMNS){
  int **dynamicArray = 0;
  //memory allocated for elements of rows.
  dynamicArray = new int *[ROWS] ;
  //memory allocated for  elements of each column.
  for( int i = 0 ; i < ROWS ; i++ ){
    dynamicArray[i] = new int[COLUMNS];
 }
  return dynamicArray;
}

void deleteMatrix(int **x, int ROWS, int COLUMNS){
    for (int i=0;i<ROWS;i++)
    {
	delete [] x[i];
    }
    delete [] x; 
}

int **readScoreMatrix(char filename[]){
	FILE *f = fopen(filename,"r");
	int ab_size = 5;
	int temp;
	int **scoremat = initializeMatrix(ab_size,ab_size);
        for (int i=0;i<ab_size;i++)
                for (int j=0;j<ab_size;j++) {
                        if (fscanf(f, "%d ", &temp) == EOF) break;
                        scoremat[i][j]=temp;
                }
        //printf("read finished\n");
        fclose(f);
        /*
        // show the score matrix
        for (int i=0;i<ab_size;i++) {
            for (int j=0;j<ab_size;j++) {
                printf("%d ",scoremat[i][j]);
            }
            printf("\n");
        }
        */
	return scoremat;
}

// remove character such as newline
int remove_char(const char *src,char *dest,char c)
{
    int removed=0;
    
    while (*src) {
        if (*src!=c) {
            *dest++=*src;
        } else {                
            ++removed;
        }
        ++src;
    }
    *dest=0;
    return removed;
}

int getSeqfromFasta(FILE *in, char descr[], char seq[])
{
        int d, s, length = 0;

        char tmp[100000];

        if(in == NULL)
                return -1;

        d = fscanf(in, "%[^\n]%*c", tmp); // read until '\n' is found

        remove_char(tmp,descr,'>');

        s = fscanf(in, "%[^>]%n%*c", tmp, &length); // read until '>' is found

        if(d == EOF || s == EOF)
                return -1;
        remove_char(tmp,seq,'\n');
        // to upper case
        string str = seq;
        for (int i=0;i<str.length();i++)
            seq[i] = toupper(seq[i]);
        return string(seq).length();
}


void findMax(int **H,int N_a, int N_b, int& i_max, int& j_max, int& H_max){
  // search H for the maximal score
  H_max = 0;
  i_max=0,j_max=0;
  for(int i=1;i<=N_a;i++){
    for(int j=1;j<=N_b;j++){
      if(H[i][j]>H_max){
        H_max = H[i][j];
        i_max = i;
        j_max = j;
      }
    }
  }
}



void printMatrix(int **H, int N_a, int N_b){
  // Print the matrix H to the console
  cout<<"**********************************************"<<endl;
  cout<<"The scoring matrix is given by  "<<endl<<endl;
  for(int i=1;i<=N_a;i++){
    for(int j=1;j<=N_b;j++){
      cout<<H[i][j]<<" ";
    }
    cout<<endl;
    }
}

void backtrack(int **I_i, int **I_j, int i_max, int j_max, string seq_a, string seq_b,int **scoremat){
     // Backtracking from H_max

  int N_a = seq_a.length();
  int N_b = seq_b.length();
  int current_i=i_max,current_j=j_max;
  int next_i=I_i[current_i][current_j];
  int next_j=I_j[current_i][current_j];
  int tick=0;
  char consensus_a[N_a+N_b+2],consensus_b[N_a+N_b+2];

  while(((current_i!=next_i) || (current_j!=next_j)) && (next_j != 0) && (next_i != 0)){

    if(next_i==current_i)  consensus_a[tick] = '-';                  // deletion in A
    else                   consensus_a[tick] = seq_a[current_i-1];   // match/mismatch in A

    if(next_j==current_j)  consensus_b[tick] = '-';                  // deletion in B
    else                   consensus_b[tick] = seq_b[current_j-1];   // match/mismatch in B

    current_i = next_i;
    current_j = next_j;
    next_i = I_i[current_i][current_j];
    next_j = I_j[current_i][current_j];
    tick++;
    }
	 // Output of the consensus motif to the console
  for(int i=tick-1;i>=0;i--) cout<<consensus_a[i];
  cout<<endl;
  for(int i=tick-1;i>=0;i--) {
      if (scoremat[char2int(consensus_a[i])][char2int(consensus_b[i])]>0) cout<<'|';
      else cout <<' ';
  }
  cout<<endl;
  for(int j=tick-1;j>=0;j--) cout<<consensus_b[j];
  cout<<endl;
}


array<int,2> find_array_max(array<int,4> a,int length){

  int max = a[0];            // start with max = first element
  int ind = 0;

  for(int i = 1; i<length; i++){
      if(a[i] > max){
	max = a[i];
	ind = i; 
      }
  }
  array<int,2> res = {max,ind};
  return res;                    // return highest value in array
}


void computeHmatrix(int **H,int **I_i, int **I_j, string seq_a, string seq_b, int **scoremat,int gap_penalty) {
  int N_a = seq_a.length();                     // get the actual lengths of the sequences
  int N_b = seq_b.length();
  array<int,4> temp;
  array <int,2> maxp;
  for(int i=1;i<=N_a;i++){
    for(int j=1;j<=N_b;j++){
      temp[0] = H[i-1][j-1]+scoremat[char2int(seq_a[i-1])][char2int(seq_b[j-1])];
      temp[1] = H[i-1][j]-gap_penalty;
      temp[2] = H[i][j-1]-gap_penalty;
      temp[3] = 0;
	  maxp = find_array_max(temp,4);
	  H[i][j] = maxp[0];
      switch(maxp[1]){
      case 0:                                  // score in (i,j) stems from a match/mismatch
        I_i[i][j] = i-1;
        I_j[i][j] = j-1;
        break;
      case 1:                                  // score in (i,j) stems from a deletion in sequence A
        I_i[i][j] = i-1;
        I_j[i][j] = j;
        break;
      case 2:                                  // score in (i,j) stems from a deletion in sequence B
        I_i[i][j] = i;
        I_j[i][j] = j-1;
        break;
      case 3:                                  // (i,j) is the beginning of a subsequence
        I_i[i][j] = i;
        I_j[i][j] = j;
        break;
      }
    }
  }
}

void localalignment(string seq_a,string seq_b,int **scoremat,int gap_penalty){
	int N_a = seq_a.length();                     // get the actual lengths of the sequences
  	int N_b = seq_b.length();
  	int **H = initializeMatrix(N_a+1,N_b+1);

  for(int i=0;i<=N_a;i++){
    for(int j=0;j<=N_b;j++){
      H[i][j]=0;
    }
  } 

  int **I_i = initializeMatrix(N_a+1,N_b+1);
  int **I_j = initializeMatrix(N_a+1,N_b+1);     // Index matrices to remember the 'path' for backtracking
 

  computeHmatrix(H,I_i,I_j,seq_a,seq_b,scoremat,gap_penalty); 

   int i_max, j_max, H_max = 0;
   findMax(H,N_a,N_b,i_max,j_max,H_max);
   cout<<"score\t"<<H_max<<"\ti\t"<<i_max<<"\tj\t"<<j_max<<endl;
   backtrack(I_i,I_j,i_max,j_max,seq_a,seq_b,scoremat);
   
   // free memory
   deleteMatrix(H,N_a+1,N_b+1);
   deleteMatrix(I_i,N_a+1,N_b+1);
   deleteMatrix(I_j,N_a+1,N_b+1);
}


void show_help(){
    cout << endl << "alignSW: local sequence alignment using the smith-waterman algorithm" << endl;
    cout << "	- Xuebing Wu (wuxb07@gmail.com)" << endl;
    cout<< endl << "Usage:"<<endl;
    cout<<"     alignSW -query query.fa -target target.fa -scoring_matrix identity.mat -gap_penalty 6" << endl;
    cout<< endl << "Input/Options:" <<endl;
    cout<<" -q | -query   <query.fa>        query sequence, fasta format" << endl;
    cout<<" -t | -tartet  <target.fa>       target sequence, fasta format" << endl;
    cout<<" -s | -scoring_matrix  <identity.mat>    scoring matrix, see notes" << endl;
    cout<<" -g | -gap_penalty    <integer>    gap penalty, positive integer, default 6" << endl;
    cout<<" -r | -reverse_complement	reverse complement query sequence" << endl;
    cout<< endl << "Notes"<<endl;
    cout<<" scoring matrix format (only the numeric part)" <<endl;
    cout<<"          ____target____ " << endl;
    cout<<"          _A__C__G__T__N " <<endl;
    cout<<"      |A|  1 -2 -2 -2 -2 " << endl;
    cout<<"      |C| -2  1 -2 -2 -2 " << endl;
    cout<<" query|G| -2 -2  1 -2 -2 " << endl;
    cout<<"      |T| -2 -2 -2  1 -2 " << endl;
    cout<<"      |N| -2 -2 -2 -2 -2 " << endl;
}

int main(int argc, char** argv){

  
  
  // read arguments
  if(argc < 8){
    show_help();
    exit(1);
  }
 

  char *nameof_seq_a;
  char *nameof_seq_b;
  char *nameof_scoremat;
  int gap_penalty = 6;
  bool rc = false;
  string str;
  for(int i=1;i<argc;i++){
    str = argv[i];
    if (str == "-query" || str == "-q"){nameof_seq_a=argv[i+1];i=i+1;}
    else if (str == "-target" || str == "-t"){nameof_seq_b=argv[i+1];i=i+1;}
    else if (str == "-scoring_matrix" || str == "-s"){nameof_scoremat=argv[i+1];i=i+1;}
    else if (str == "-gap_penalty" || str == "-g"){gap_penalty=atof(argv[i+1]);i=i+1;}
    else if (str == "-reverse_complement" || str == "-r"){rc=true;}
    else { cout << "invalid option: " << str << endl;show_help();exit(1);}
}

//cout << nameof_seq_a << endl << nameof_seq_b << endl << nameof_scoremat <<endl<<gap_penalty << endl;

  // load score matrix
  int **scoremat = readScoreMatrix(nameof_scoremat);

  // load sequences
  string seq_a,seq_b; 

  char desc1[100000];
  char desc2[100000];
  char seq[100000];
  
  FILE *fq,*fs;
  fq = fopen(nameof_seq_a,"r");
 
  while(getSeqfromFasta(fq,desc1,seq)>0)
  {
      seq_a = seq;
      //cout <<"seq_a:" <<desc1<<",len="<<seq_a.length()<<endl;

      if (rc)
      {
	  string seq_a_rc = reverseComplement(seq_a);
	  //cout <<"reverse complement:"<<endl<<seq_a_rc<<endl;
	  seq_a = seq_a_rc;
      }

      fs = fopen(nameof_seq_b,"r");
      while(getSeqfromFasta(fs,desc2,seq) > 0){
          seq_b = seq;
          cout <<desc1<<","<<seq_a.length()<<","<<desc2<<","<<seq_b.length()<<endl;  
          //cout <<seq_b<<endl;
          localalignment(seq_a,seq_b,scoremat,gap_penalty);
      }
      fclose(fs);
  }
  fclose(fq);
} // END of main

