#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <Rcpp.h>
using namespace Rcpp;

// sourceCpp("/media/30deser/user/hadlich/Austausch/Quails/RefSeq/read_file.cpp")

// library(Rcpp)
// setwd("/disk1/R/reviseArrayAssignment/TEST_PIG/in/ensembl_Sscrofa11.1/")
// sourceCpp("/disk1/R/reviseArrayAssignment/SKRIPTE/get-fasta.cpp")


//' @export
void sample_defaults(NumericVector x = NumericVector::create(), //Size 0 vector
                     bool bias =true, //Set to true
                     std::string method = "rcpp rules!"){ //Default string
  Rcout << "x size: " << x.size() << ", ";
  Rcout << "bias value: " << bias << ", ";
  Rcout << "method value: " << method << std::endl;
}

// [[Rcpp::export]]
CharacterVector read_file_cpp(CharacterVector path) {
  std::string fname = as<std::string>(path);
  std::ifstream t(fname.c_str());
  std::stringstream buffer;
  buffer << t.rdbuf();
  return buffer.str();
}

// [[Rcpp::export]]
CharacterVector read_file_cpp2(std::string path) {
  std::ifstream in(path.c_str());
  std::string contents;
  in.seekg(0, std::ios::end);
  contents.resize(in.tellg());
  in.seekg(0, std::ios::beg);
  in.read(&contents[0], contents.size());
  in.close();
  return(contents);
}


// [[Rcpp::export]]
CharacterVector reversecomp(CharacterVector x) {
  CharacterVector y(x.length());
  for (int i=0; i<=x.length(); i++) {
    std::string str = as<std::string>(x[i]);
    int n = str.length();
    for (int j = 0; j < n / 2; j++)
      std::swap(str[j], str[n - j - 1]);
    y[i] = str;
  }
  return(y);
}


// [[Rcpp::export]]
void prepare_vcf(std::string file_in, std::string folder_out, std::string fs = "/") {
//std::string fname = as<std::string>(file_in);
  std::ifstream in(file_in.c_str());
  std::string rec;
  std::string chromosome("");
  std::ofstream outfile;
  while(getline(in, rec, '\n')){ // get text line by line
    if (rec.length()==0 || rec.at(0) == '#') continue;
    CharacterVector val(4);
    std::string item;
    std::stringstream ss(rec);
    size_t i=0;
    for (i=0; i<4; i++) {
      getline(ss, item, '\t');
      if (i==2) getline(ss, item, '\t');
      val[i] = item;
      if (i>1 && item.length()>1) break;
    }
    if (i<4) continue;

    if (chromosome.compare(val[0]) != 0) {
      if (chromosome.length() > 0) {
        outfile.close();
      }
      chromosome = val[0];
      ss.str("");
      ss << folder_out << fs << chromosome << ".vcf.cut";
      if (chromosome.length() < 10)  Rcout << "\t\tchromosome" << chromosome << std::endl;
      outfile.open(ss.str().c_str());
    }
    //Rcout << i << " " << rec << std::endl;
    outfile << val[1] << "\t" << val[2] << "\t" << val[3] << std::endl;
  }
  outfile.close();
}


// https://lemire.me/blog/2012/06/20/do-not-waste-time-with-stl-vectors/
// [[Rcpp::export]]
List get_fasta(CharacterVector path, List chrlist) {
  List records = List::create();
  std::string fname = as<std::string>(path);
  std::ifstream in(fname.c_str());
//Rcout << "l" << chrlist.length() << (chrlist["a"]) << std::endl;
  in.get(); // remove first '>'
  std::string rec;
  CharacterVector chrs = chrlist.names();
  while(getline(in,rec,'>')){ // get text until next '>'

    rec.erase(std::remove(rec.begin(),rec.end(),'\r'),rec.end()); // remove carriage return (CR) symbols [important for windows]
    int newLineLoc = rec.find('\n');
    std::string header = rec.substr(0,newLineLoc);
    int pos = rec.find(' ');
    header = header.substr(0, pos);
    if (std::find(chrs.begin(), chrs.end(), header) == chrs.end())
      continue;
    std::string sequence = rec.substr(newLineLoc+1, rec.length()-newLineLoc-2);
    sequence.erase(std::remove(sequence.begin(),sequence.end(),'\n'),sequence.end()); // remove new line symbols
//records[header] = sequence;

    List df = chrlist[header];
    NumericVector start = df["start"];
    NumericVector end = df["end"];
//    Rcout << header << ":" << start.length() << std::endl;
    CharacterVector seqs(start.length());
    for (int i=0; i<start.length(); i++) {
      //    Rcout << "method value: " << start[i] << " : " << stop[i] << " len=" << stop[i]-start[i]+1 << std::endl;
      seqs[i] = sequence.substr(start[i]-1, end[i]-start[i]+1);
    }
    records[header] = seqs;
  }
  return(records);
}

// [[Rcpp::export]]
CharacterVector read_fasta(CharacterVector path) {
  CharacterVector records;
  std::string fname = as<std::string>(path);
  std::ifstream in(fname.c_str());
  in.get(); // remove first '>'
  std::string rec;
  while(getline(in,rec,'>')){ // get text until next '>'
    int newLineLoc = rec.find('\n');
    std::string header = rec.substr(0,newLineLoc);
    std::string sequence = rec.substr(newLineLoc+1, rec.length()-newLineLoc-2);
    sequence.erase(std::remove(sequence.begin(),sequence.end(),'\n'),sequence.end()); // remove new line symbols
    records[header]=sequence;
  }
  return(records);
}
