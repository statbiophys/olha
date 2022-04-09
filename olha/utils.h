#ifndef UTILS
#define UTILS

#include <chrono>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <map>
#include <tuple>
#include <functional>
#include <array>
#include <random>
#include <stdexcept>
#include <time.h>

// maximal distance between sequences, used as a default return value
// when the distances are above the threshold.
#define MAX_DISTANCE 10000

// size of the amino alphabet
#define ALPHABET_SIZE 21


/* Nucleotides: A, C, G, T = 0, 1, 2, 3 */
/* Amino-acids: 
* , A , C , D , E , F , G , H , I , K , L , M , N , P , Q , R , S , T , V , W , Y
0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20
*/
  const std::array<std::size_t, 64> DNA_TO_AMINO = {
    9, 12, 9, 12, 17, 17, 17, 17, 15, 16, 15, 16, 8, 8, 11, 8, 14, 7, 14, 7, 13, 13, 13, 13, 15, 15, 15, 15, 10, 10, 10, 10, 4, 3, 4, 3, 1, 1, 1, 1, 6, 6, 6, 6, 18, 18, 18, 18, 0, 20, 0, 20, 16, 16, 16, 16, 0, 2, 19, 2, 10, 5, 10, 5};

  const std::map<char, int> AMINO_ACID_TO_INT{{'*',0},{'A',1},{'C',2},{'D',3},{'E',4},{'F',5},{'G',6},{'H',7},{'I',8},{'K',9},{'L',10},{'M',11},{'N',12},{'P',13},{'Q',14},{'R',15}, {'S',16},{'T',17},{'V',18},{'W',19},{'Y',20}};


  using amino_acid=short unsigned int;
  using dna=std::vector<short unsigned int>;
  using amino=std::vector<short unsigned int>;
 
  /* Return a dna sequence from std::string */
  inline dna string_to_dna(const std::string& x){
    dna ns;
    for(char l:x){
      if(l == 'A')
	ns.push_back(0);
      else if(l == 'C')
	ns.push_back(1);
      else if(l == 'G')
	ns.push_back(2);
      else if(l == 'T')
	ns.push_back(3);
    }
    return ns;
  }

    /* Return a dna sequence from std::string */
  inline dna string_to_amino(const std::string& x){
    amino as;
    for(char l:x){
      as.push_back(AMINO_ACID_TO_INT.at(l));
    }
    return as;
  }



/* Return a std::string of the dna sequence */
inline std::string dna_to_string(const dna& x){
  const std::string dna_letter = "ACGT";
  std::string st = "";
  for(auto l: x){
    st += dna_letter[l];
  }
  return st;
}

/* Return a std::string of the amino sequence */
inline std::string amino_to_string(const amino& x){
  const std::string amino_letter = "*ACDEFGHIKLMNPQRSTVWY";
  std::string st = "";
  for(auto l: x){
    st += amino_letter[l];
  }
  return st;
}


/* Transcribe a dna sequence to amino acids */
inline amino translate(const dna& x){
  amino y;
  y.reserve(x.size()/3);
  if(x.size()%3 != 0){
    std::cout << "Translation not possible, invalid length.";
    return amino();
  }
  for(uint ii=0; ii < x.size(); ii += 3){
    y.push_back(DNA_TO_AMINO[16*x[ii] + 4*x[ii+1] + x[ii+2]]);
    }
  return y;
}


/* Add a uniform error on the dna sequence. 
Each nucleotide has a probability "error_rate" of switching  */
inline void add_error(dna& x, double error_rate, std::mt19937_64& e){
  std::uniform_real_distribution<> uniform(0, 1);
  std::uniform_int_distribution<> random_nucleotide(0, 3);
  for(std::size_t ii=0; ii < x.size(); ++ii){
    if(uniform(e) < error_rate){
      x[ii] = random_nucleotide(e);
    }
  }
}


/* Create a non-uniform discrete distribution 
   Methods:
   @ DiscreteDistribution(const std::vector<double>& weights, int seed)
   Create a DiscreteDistribution, take for argument the vector of each elements
   probability (not necessarily normed), and an integer seed.
   @ generate
   Generate one element
*/
class DiscreteDistribution{
 private:
  std::discrete_distribution<> dis;
 public:
  DiscreteDistribution(){
    dis = std::discrete_distribution<>();
  }
  template< class InputIt >
    DiscreteDistribution(InputIt first, InputIt last){
    dis = std::discrete_distribution<>(first, last);
  }

  inline std::size_t generate(std::mt19937_64& gen){
    return dis(gen);
  }
};



/* Create an object that can generate chains of
   nucleotides of arbitrary sizes 
   Methods:
   @ MarkovDNA(initial, matrix)
   takes as argument an array of size 16, each block of 4 is a distribution.
   @ generate(std::size_t nn)
   Generate a nucleotide sequence of size nn
*/
class MarkovDNA{
 private:
  DiscreteDistribution dd_initial;
  std::array<DiscreteDistribution, 4> dd_matrix;
 public:
  MarkovDNA(){
  }
  MarkovDNA(const std::array<double, 4>& initial,
	   const std::array<double, 16>& matrix){
   dd_initial = DiscreteDistribution(initial.begin(), initial.end());
   for(std::size_t ii=0; ii < 4; ++ii)
     dd_matrix[ii] = DiscreteDistribution(matrix.begin() + 4*ii, matrix.begin() + 4*(ii+1));
 }

  inline dna generate(std::size_t nn, std::mt19937_64& random_generator){
    dna x;
    if(nn == 0)
      return x;
    std::size_t nuc = dd_initial.generate(random_generator);
    x.push_back(nuc);
    for(std::size_t ii=0; ii < nn-1; ++ii){
      nuc = dd_matrix[nuc].generate(random_generator);
      x.push_back(nuc);
    }
    return x;
  }
};




#endif //UTILS
