/** dna.cpp
 *Class that reads in files containing DNA strings, calculates statistics
 *and then write those statistics and generated DNA string to an outfile.
 *@author Mary Kate Crittenden
 *2278514
 *critt102@mail.chapman.edu
 *CPSC 350-01
 *Assignment 1
 */
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;
//The main method
int main(int argc, char** argv){

  int repeat=1;
  //while loop for allowing users to process multiple files while program is running
  while(repeat>=1){
    //prompt user for command line argument
    cout<<"File name: ";
    //take in arguement from command line
    cin>>argv[0];
    //current line being read
    string line;
    //reponse from user on whether they want to process another file
    char again;
    //total length of all lines in file, used to calculate mean
    double length_of_lines=0;
    //total number of dna strings in the file
    double num_of_lines=0;
    //total number of nucleotides in the dna strings
    double num_of_nucleos;
    //place to store the variance calculation final value
    double variance=0;
    //store value of number of A nucleotides
    int num_of_a=0;
    //store value of number of A nucleotides
    int num_of_c=0;
    //store value of number of A nucleotides
    int num_of_g=0;
    //store value of number of A nucleotides
    int num_of_t=0;
    //number of aa bigrams
    double aa=0;
    //number of ac bigrams
    double ac=0;
    //number of at bigrams
    double at=0;
    //number of ag bigrams
    double ag=0;
    //number of ca bigrams
    double ca=0;
    //number of cc bigrams
    double cc=0;
    //number of ct bigrams
    double ct=0;
    //number of cg bigrams
    double cg=0;
    //number of tt bigrams
    double tt=0;
    //number of ta bigrams
    double ta=0;
    //number of tc bigrams
    double tc=0;
    //number of tg bigrams
    double tg=0;
    //number of gg bigrams
    double gg=0;
    //number of gs bigrams
    double ga=0;
    //number of gc bigrams
    double gc=0;
    //number of gt bigrams
    double gt=0;
    //total number of bigrams
    double num_of_bigrams=0;
    //placeholder of blank string, use to concatenate strnigs
    string blank="";
    //final mean value
    double mean;
    //final standard deviation value
    double sd;
    //a variable for distribution function
    double a=0.00;
    //a variable for distribution function
    double b=0.00;
    //c variable for distribution function
    double c=0.00;
    //final integer calculated using distribution function
    int d;
    //Probability of nucleotide a
    double prob_a;
    //Probability of nucleotide c
    double prob_c;
    //Probability of nucleotide t
    double prob_t;
    //Probability of nucleotide g
    double prob_g;
    //input and output streams
    ifstream dna_file;
    ofstream out_file;
    //open input file
    dna_file.open(argv[0], ios::in);
    //for processing first file
    if(repeat==1){
      out_file.open("marykatecrittenden.out", ios::out);
    }
    //for processing 2+ files (to append to end of existing output)
    else if(repeat>1){
          out_file.open("marykatecrittenden.out", ios::out | ios::app);
    }
    //read each string of dna from files and calculate statistics
    if(dna_file.is_open()){
      cout<<"File found."<<endl;
      while(getline(dna_file, line)){
        num_of_lines++;
        length_of_lines+=line.size();
        //count the number of each nucleotide and bigram type
        for(int i=0; i<line.size();++i){
          char nucleo=tolower(line[i]);
          if(nucleo=='a')
      			num_of_a++;
      		else if(nucleo=='c')
      			num_of_c++;
      		else if(nucleo=='t')
      			num_of_t++;
      		else if(nucleo=='g')
      			num_of_g++;
          //store the adjacent nucleotide for counting bigrams
          char pair=tolower(line[i+1]);
          if(pair=='a')
        		num_of_a++;
        	else if(pair=='c')
      			num_of_c++;
      		else if(pair=='t')
      			num_of_t++;
          else if(pair=='g')
        		num_of_g++;
          //create a string to represent the current bigram pair
          string bigram=blank+nucleo+pair;
          if(bigram.compare("aa")==0)
            aa++;
          else if(bigram.compare("ac")==0)
            ac++;
          else if(bigram.compare("at")==0)
            at++;
          else if(bigram.compare("ag")==0)
            ag++;
          else if(bigram.compare("ca")==0)
            ca++;
          else if(bigram.compare("cc")==0)
            cc++;
          else if(bigram.compare("ct")==0)
            ct++;
          else if(bigram.compare("cg")==0)
            cg++;
          else if(bigram.compare("tt")==0)
            tt++;
          else if(bigram.compare("ta")==0)
            ta++;
          else if(bigram.compare("tc")==0)
            tc++;
          else if(bigram.compare("tg")==0)
            tg++;
          else if(bigram.compare("gg")==0)
            gg++;
          else if(bigram.compare("ga")==0)
            ga++;
          else if(bigram.compare("gc")==0)
            gc++;
          else if(bigram.compare("gt")==0)
            gt++;
          num_of_bigrams++;
          //additional i to prevent counting overlapping bigrams
          i++;

        }
      }
      num_of_nucleos=num_of_a+num_of_c+num_of_g+num_of_t;
      dna_file.close();
    }
    dna_file.open(argv[0], ios::in);
    //calculate variance
    if(dna_file.is_open()){
      //variable for storing part of variance calculation
      double vari2=0;
      while(getline(dna_file, line)){
        //initial part of variance calculation
        double vari=(line.size())-(length_of_lines/num_of_lines);
        vari2+=(vari*vari);
      }
      variance=vari2/num_of_lines;
      dna_file.close();
    }
    //write header to file only once
    if(out_file.is_open() && repeat==1){
      //Header
      out_file<<"Mary Kate Crittenden \n"<<"2278514 \n"<<"critt102@mail.chapman.edu \n"<<"CPSC 350-01 \n"<<"Assignment 1 \n"<<endl;
    }
    //writing statistics to file and generation 1000 dna strings
    if(out_file.is_open()){
      //initialize rand function to be different each time program is run
      srand (time(NULL));
      prob_a=(num_of_a/num_of_nucleos)*100;
      prob_c=(num_of_c/num_of_nucleos)*100;
      prob_t=(num_of_t/num_of_nucleos)*100;
      prob_g=(num_of_g/num_of_nucleos)*100;
      mean=length_of_lines/num_of_lines;
      sd=sqrt(variance);

      //summary statistics
      out_file<<"Sum: "<<length_of_lines<<endl;
      out_file<<"Mean: "<<mean<<endl;
      out_file<<"Variance: "<<variance<<endl;
      out_file<<"Standard Deviation: "<<sqrt(variance)<<endl;
      out_file<<"\n"<<"Nucleotide Probabilities: "<<endl;
      out_file<<"Probability of A: "<<prob_a<<"%"<<endl;
      out_file<<"Probability of C: "<<prob_c<<"%"<<endl;
      out_file<<"Probability of T: "<<prob_t<<"%"<<endl;
      out_file<<"Probability of G: "<<prob_g<<"%"<<endl;
      out_file<<"\n"<<"Bigram Probabilities: "<<endl;
      out_file<<"Probability of aa bigram: "<<(aa/num_of_bigrams)*100<<"%"<<endl;
      out_file<<"Probability of ac bigram: "<<(ac/num_of_bigrams)*100<<"%"<<endl;
      out_file<<"Probability of at bigram: "<<(at/num_of_bigrams)*100<<"%"<<endl;
      out_file<<"Probability of ag bigram: "<<(ag/num_of_bigrams)*100<<"%"<<endl;
      out_file<<"Probability of ca bigram: "<<(ca/num_of_bigrams)*100<<"%"<<endl;
      out_file<<"Probability of cc bigram: "<<(cc/num_of_bigrams)*100<<"%"<<endl;
      out_file<<"Probability of ct bigram: "<<(ct/num_of_bigrams)*100<<"%"<<endl;
      out_file<<"Probability of cg bigram: "<<(cg/num_of_bigrams)*100<<"%"<<endl;
      out_file<<"Probability of tt bigram: "<<(tt/num_of_bigrams)*100<<"%"<<endl;
      out_file<<"Probability of ta bigram: "<<(ta/num_of_bigrams)*100<<"%"<<endl;
      out_file<<"Probability of tc bigram: "<<(tc/num_of_bigrams)*100<<"%"<<endl;
      out_file<<"Probability of tg bigram: "<<(tg/num_of_bigrams)*100<<"%"<<endl;
      out_file<<"Probability of gg bigram: "<<(gg/num_of_bigrams)*100<<"%"<<endl;
      out_file<<"Probability of ga bigram: "<<(ga/num_of_bigrams)*100<<"%"<<endl;
      out_file<<"Probability of gc bigram: "<<(gc/num_of_bigrams)*100<<"%"<<endl;
      out_file<<"Probability of gt bigram: "<<(gt/num_of_bigrams)*100<<"%"<<endl;
      out_file<<"\n"<<"Generated DNA strings: "<<endl;
      //dna string generation
      for(int j=0; j<1000; ++j){
        //placeholder for current dna string
        string dna_string="";
        a=((double)rand())/(RAND_MAX);
        b=((double)rand())/(RAND_MAX);
        c=(sqrt((-2)*log(a)))*(cos(2*M_PI*b));
        d=(sd*c)+mean;
        //if distribution calculation rounds down to 0, make it round up to 1 nucleotide
        if(d<1){
          d=1;
        }
        /* generate nucleotides using length from distribution calculation and
        probability of each nucleotide */
        for(int k=1;k<=d;++k){
          //represents the random number used to determine the type of current nucelotide
          double this_nucleo=rand()%100;
          if(this_nucleo<prob_a){
            dna_string+="a";
          }
          else if(prob_a<this_nucleo && this_nucleo<(prob_a+prob_c)){
            dna_string+="c";
          }
          else if((prob_a+prob_c)<this_nucleo && this_nucleo<(prob_a+prob_c+prob_t)){
            dna_string+="t";
          }
          else if((prob_a+prob_c+prob_t)<this_nucleo && this_nucleo<(prob_a+prob_c+prob_t+prob_g)){
            dna_string+="g";
          }
        }
        out_file<<dna_string<<endl;
      }
      out_file<<endl;
      out_file.close();
    }
    cout<<"Results printed to output file."<<endl;
    //ask users if they want to process another file
    cout<<"Do you want to process another list? y or n?"<<endl;
    cin>>again;
    again=tolower(again);
    if(again=='y'){
      repeat+=1;
    }
    else{
      repeat=0;
    }
  }
  return 0;
}
