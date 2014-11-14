#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#define DIM 2

using namespace std;

#if 0
template <class T>
ostream &operator<<(ostream &stream, vector <T> vec){
  for(int i = 0; i < vec.size(); ++i){
    stream << vec.at(i) << " ";
  }
  stream << endl;
  return stream;
}
#endif

class data{
  vector <string> head;
  vector <vector <double> > body;
public:
  data(int dim){ head.resize(dim); }
  int dim(){ return head.size(); }
  friend data read_data(istream &);
  friend ostream &operator<<(ostream &, data);
};

ostream &operator<<(ostream &stream, data dat){
  /* dump the contents of the struct data */
  for(int i = 0; i < dat.dim(); ++i){ /* header */

    stream << dat.head.at(i) << "\t";
  }
  stream << endl;
  for(int n = 0; n < dat.body.size(); ++n){ /* body */
    for(int i = 0; i < dat.dim(); ++i){
      stream << dat.body.at(n).at(i) << "\t";
    }
    stream << endl;
  }
  return stream;
}

inline data read_data(istream &ifs){
  /* read data from a stream 
   * - the first line of the data file is header line
   * - label of the variables are written in the header line.
   * - the following line is the data.
   */

  data input_data(DIM);
  for(int i = 0; i < input_data.dim(); i++){
    ifs >> input_data.head.at(i); /* read header line */
  }

  vector <double> one_item(input_data.dim());
  while(!ifs.eof()){
    /* read all data from file */
    for(int i = 0; i < input_data.dim(); i++){
      ifs >> one_item.at(i); 
    }
    input_data.body.push_back(one_item);
  }
  cout << input_data;
  return input_data;
}

data read_data(char *filename){
  /* read data from a file */
  ifstream ifs(filename);
  if ( ifs.fail() ){
    cerr << "cannot open data file" << endl; exit(1);
  }
  data input_data = read_data(ifs);
  ifs.close(); /* file close */
  return input_data;
}

int main(){
  read_data((char *)"input/sample_mix.txt");
  return 0;
}
