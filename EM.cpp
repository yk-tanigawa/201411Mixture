#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random> /* require C++11 or later */
#include <cmath>
#include "my_vector.hpp" /* vector operation */

#define DIM 2

using namespace std;

class data;
class params;
class em;

/* ここからデータの読み込み */

class data{
  vector <string> head;
  vector <vector <double> > body;
public:
  data(){}
  data(int dim){ head.resize(dim); }
  int dim(){ return (int)head.size(); }
  int size(){ return (int)body.size(); }
  vector<double> at(int i){ return body.at(i); }
  friend data read_data(istream &);
  friend ostream &operator<<(ostream &, data);
  friend vector<double> normal_distribution(data, params);
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
  /* 
   * read data from a stream 
   * - the first line of the data file is header line
   * - label of the variables are written in the header line.
   * - the following lines are the body of data.
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

/*
 *
 * ここからparameters
 *
 */

class parameter{
  /* parameters of a Mixture Gaussian model */
  vector<double> pi_;
  vector<vector <double> > mu_;
  vector<double> sigma_;
public:
  parameter(){}
  parameter(int k, int dim = 2){ /* default: dim = 2 */
    /* initialize parameters 
     *  k   = #{ components }
     *  dim = #{ dimention } */
    pi_.resize(k, 1.0 / k); /* pi_k  = 1.0 / dim */
    sigma_.resize(k, 0.3);    /* sigma^2_k = 0.3 */
    mu_.resize(k); /* assign random values to mu */
    {
      std::random_device rd;
      std::mt19937 mt(0); /* DEBUG用にseedを固定 */
      //std::mt19937 mt(rd());
      std::uniform_real_distribution<double> rand(0.0, 1.0);
      for(int i = 0; i < mu_.size(); ++i){
	mu_.at(i).resize(dim);
	for(int j = 0; j < dim; ++j){ mu_.at(i).at(j) = rand(mt); }
      }
    }
  }
  int k(){ return (int)pi_.size(); } /* returns #{ components } */
  int dim(){ return (int)mu_.at(0).size(); } /* returns dimention */
  vector<double> mu(int i){ return mu_.at(i); }
  void set_mu(int i, vector<double> mu_new){
    mu_.at(i) =  mu_new;
  }
  double sigma(int i){ return sigma_.at(i); }
  void set_sigma(int i, double sigma_new){
    sigma_.at(i) = sigma_new;
  }
  double pi(int i){ return pi_.at(i); }
  void set_pi(int i, double pi_new){
    pi_.at(i) = pi_new;
  }
  double det_sigma(){
    double determinant = 1;
    for(int i = 0; i < (int)sigma_.size(); ++i){
      determinant *= sigma_.at(i);
    }
    return determinant;
  }
  ostream &export_params_head(ostream &stream){
    for(int i = 0; i < k(); ++i){
      stream << i << "_" << "p" << "\t" /*pi*/
	     << i << "_" << "s" << "\t"; /* sigma */
      for(int d = 0; d < dim(); ++d){
	stream << i << "_" << "u" << d << "\t"; /* mu */
      }
    }
    stream << endl;
    return stream;
  }
  ostream &export_params(ostream &stream){
    /* pi sigma mu*/
    for(int i = 0; i < k(); ++i){
      stream << pi_.at(i) << "\t"
	     << sigma_.at(i) << "\t";
      for(int d = 0; d < dim(); ++d){
	stream << mu_.at(i).at(d) << "\t";
      }
    }
    stream << endl;
    return stream;
  }
  friend ostream &operator<<(ostream &, parameter);
  friend vector<double> normal_distribution(data, params);
};

ostream &operator<<(ostream &stream, parameter params){
  /* dump the contents of the struct parameter */
  stream << "number of components : " << params.k() << endl;
  stream << "i|pi\tsigma\t<mu_>" << endl;
  for(int i = 0; i < params.k(); ++i){
    stream << i << "|" << params.pi_.at(i) << "\t"
	   << params.sigma_.at(i) << "\t<";
    for(int d = 0; d < params.dim() - 1; ++d){
      stream << params.mu_.at(i).at(d) << ", ";
    }
    stream << params.mu_.at(i).at(params.dim() - 1)
	   << ">" << endl;
  }
  return stream;
}

/*
 *
 * EM algorithm
 *
 */

double NormalDistribution(vector<double> x, vector<double> mu, double sigma){
  int dim = (int)x.size() < (int)mu.size() ? (int)x.size() : (int)mu.size();
  x.resize(dim); mu.resize(dim);
  double expterm = - (1.0 / 2.0) * norm(x - mu) / (sigma * sigma);
  return (1.0 / (pow(2 * M_PI, dim / 2.0 ) * sigma)) * exp(expterm);
}

class em{
  data data;
  parameter params;
  vector< vector<double> > gamma;
  vector< vector<double> > w;
public:
  em(class data &d, class parameter &p){
    data = d; params = p;
    w.resize(p.k());
    for(int z = 0; z < params.k(); ++z){
      w.at(z).resize(data.size());
    }    
  }
  void estep(){
    for(int i = 0; i < data.size(); ++i){
      double partition_i = 0;
      for(int z = 0; z < params.k(); ++z){	
	w.at(z).at(i)
	  = NormalDistribution(data.at(i),
			       params.mu(z),
			       params.sigma(z))
	  * params.pi(z);
	partition_i += w.at(z).at(i);
      }
      for(int z = 0; z < params.k(); ++z){	
	w.at(z).at(i) /= partition_i;
      }
    }
  }
  void mstep(){
    vector<double> N(params.k(), 0);
    for(int z = 0; z < params.k(); ++z){
      for(int i = 0; i < data.size(); ++i){
	N.at(z) += w.at(z).at(i);
      }
      {
	vector<double> mu_new(data.dim(), 0);
	for(int i = 0; i < data.size(); ++i){
	  mu_new += w.at(z).at(i) * data.at(i);
	}
	mu_new *= (1.0 / N.at(z));
	params.set_mu(z, mu_new);
      }
      {
	double sigma_new = 0;
	for(int i = 0; i < data.size(); ++i){
	  sigma_new +=
	    norm(data.at(i) - params.mu(z))
	    * w.at(z).at(i);
	}
	sigma_new /= (params.dim() * N.at(z));
	params.set_sigma(z, sqrt(sigma_new));
      }
      params.set_pi(z, N.at(z) / data.size());
    }
  }
  void show_params_head(){
    params.export_params_head(cout);
  }
  void show_params(){
    params.export_params(cout);
  }
  double loglikelihood(){
    double logsum = 0;
    for(int i = 0; i < data.size(); ++i){
      double sum = 0;
      for(int z = 0; z < params.k(); ++z){
	sum += params.pi(z)
	  * NormalDistribution(data.at(i),
			       params.mu(z),
			       params.sigma(z));
      }
      logsum += log(sum);
    }
    return logsum;
  }
};


int main(){
  data d = read_data((char *)"input/sample_mix.txt");
  parameter params(3);
  em em(d, params);

#if 1
  cout << "n\t";
  //em.show_params_head();
  for(int i = 0; i < 10; ++i)
  {
    cout << i << "\t";
    cout << em.loglikelihood() << endl;
    em.estep();
    em.mstep();
    //em.show_params();
  }
  #endif
  return 0;
}
