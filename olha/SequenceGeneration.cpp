#include "SequenceGeneration.h"


SequenceGenerationVDJ::SequenceGenerationVDJ(const std::string& file_gen, uint64_t seed , bool verbose_arg){
  if(seed)
    random_generator = std::mt19937_64(seed);
  else
    random_generator = std::mt19937_64(std::random_device()());

  verbose = verbose_arg;
  load_file(file_gen);
}


void SequenceGenerationVDJ::load_file(const std::string& file_gen){
  std::ifstream file(file_gen);
  std::istringstream iss;
  std::string line;
  if (file.is_open()){
      /* Load genes */
      std::getline(file, line); // # V
      std::getline(file, line); // nb of V
      iss = std::istringstream(line);
      if(not (iss >> nb_V)) throw std::runtime_error("Invalid file");
      for(std::size_t ii=0; ii < nb_V; ++ii){
	std::getline(file, line);
	seq_V_CDR3.push_back(string_to_dna(line));
      }

      std::getline(file, line); // # D
      std::getline(file, line); // nb of D
      iss = std::istringstream(line);
      if(not (iss >> nb_D)) throw std::runtime_error("Invalid file");
      for(std::size_t ii=0; ii < nb_D; ++ii){
	std::getline(file, line);
	seq_D_CDR3.push_back(string_to_dna(line));
      }

      std::getline(file, line); // # J
      std::getline(file, line); // nb of J
      iss = std::istringstream(line);
      if(not (iss >> nb_J)) throw std::runtime_error("Invalid file");
      for(std::size_t ii=0; ii < nb_J; ++ii){
	std::getline(file, line);
	seq_J_CDR3.push_back(string_to_dna(line));
      }

      /* Marginals */
      std::getline(file, line); // # Marginals
      std::getline(file, line); // # P(V)
      std::getline(file, line); // nbV
      std::getline(file, line);
      iss = std::istringstream(line);
      std::vector<double> probas_V;
      for(std::size_t ii = 0; ii < nb_V; ++ii){
	double p;
	if(not (iss >> p)) throw std::runtime_error("Invalid file");
	probas_V.push_back(p);
      }
      d_V = DiscreteDistribution(probas_V.begin(), probas_V.end());

      std::getline(file, line); // # P(D, J)
      std::getline(file, line); // nbD nbJ
      std::getline(file, line);
      iss = std::istringstream(line);
      std::vector<double> probas_DJ;
      for(std::size_t ii = 0; ii < nb_D * nb_J; ++ii){
	double p;
	if(not (iss >> p)) throw std::runtime_error("Invalid file");
	probas_DJ.push_back(p);
      }
      d_DJ = DiscreteDistribution(probas_DJ.begin(), probas_DJ.end());

      std::getline(file, line); // # P(delV|V)
      std::getline(file, line); // nb_V  nb_delV
      iss = std::istringstream(line);
      if(not (iss >> nb_V >> nb_delV)) throw std::runtime_error("Invalid file");
      for(std::size_t ii=0; ii < nb_V; ++ii){
	std::getline(file, line);
	iss = std::istringstream(line);
	std::vector<double> delV;
	delV.clear();
	for(std::size_t jj = 0; jj < nb_delV; ++jj){
	  double p;
	  if(not (iss >> p)) throw std::runtime_error("Invalid file");
	  delV.push_back(p);
	}
	d_delV.push_back(DiscreteDistribution(delV.begin(), delV.end()));
      }


      std::getline(file, line); // # P(delJ|J)
      std::getline(file, line); // nb_J  nb_delJ
      iss = std::istringstream(line);
      if(not (iss >> nb_J >> nb_delJ)) throw std::runtime_error("Invalid file");
      for(std::size_t ii=0; ii < nb_J; ++ii){
	std::getline(file, line);
	iss = std::istringstream(line);
	std::vector<double> delJ;
	delJ.clear();
	for(std::size_t jj = 0; jj < nb_delJ; ++jj){
	  double p;
	  if(not (iss >> p)) throw std::runtime_error("Invalid file");
	  delJ.push_back(p);
	}
	d_delJ.push_back(DiscreteDistribution(delJ.begin(), delJ.end()));
      }


      std::getline(file, line); // # P(delDl, delDr | D)
      std::getline(file, line); // nb_delDr nb_delDl nb_D
      iss = std::istringstream(line);
      if(not (iss >> nb_delDr >> nb_delDl >> nb_D)) throw std::runtime_error("Invalid file");
      for(std::size_t ii=0; ii < nb_D; ++ii){
	std::getline(file, line);
	iss = std::istringstream(line);
	std::vector<double> delD;
	for(std::size_t jj = 0; jj < nb_delDl * nb_delDr; ++jj){
	  double p;
	  if(not (iss >> p)) throw std::runtime_error("Invalid file");
	  delD.push_back(p);
	}
	d_delD3_delD5.push_back(DiscreteDistribution(delD.begin(), delD.end()));
      }



      std::getline(file, line); // # P(insVD)
      std::getline(file, line); // nb_insVD
      iss = std::istringstream(line);
      if(not (iss >> nb_insVD)) throw std::runtime_error("Invalid file");
      std::getline(file, line);
      iss = std::istringstream(line);
      std::vector<double> proba_insVDs;
      for(std::size_t ii=0; ii < nb_insVD; ++ii){
	double p;
	if(not (iss >> p)) throw std::runtime_error("Invalid file");
	proba_insVDs.push_back(p);
      }

      d_insVD = DiscreteDistribution(proba_insVDs.begin(),
				      proba_insVDs.end());


      std::getline(file, line); // # Markov VD
      std::getline(file, line);
      iss = std::istringstream(line);
      std::array<double, 16> markov_VD_mat;
      for(std::size_t ii=0; ii < 16; ++ii){
	double p;
	if(not (iss >> p)) throw std::runtime_error("Invalid file");
	markov_VD_mat[ii] = p;
      }

      std::getline(file, line); // # First nucleotide bias VD
      std::getline(file, line);
      iss = std::istringstream(line);
      std::array<double, 4> first_nuc_VD;
      for(std::size_t ii=0; ii < 4; ++ii){
	double p;
	if(not (iss >> p)) throw std::runtime_error("Invalid file");
	first_nuc_VD[ii] = p;
      }

      markov_VD = MarkovDNA(first_nuc_VD, markov_VD_mat);

      std::getline(file, line); // # P(insDJ)
      std::getline(file, line); // nb_insDJ
      iss = std::istringstream(line);
      if(not (iss >> nb_insDJ)) throw std::runtime_error("Invalid file");
      std::getline(file, line);
      iss = std::istringstream(line);
      std::vector<double> proba_insDJs;
      for(std::size_t ii=0; ii < nb_insDJ; ++ii){
	double p;
	if(not (iss >> p)) throw std::runtime_error("Invalid file");
	proba_insDJs.push_back(p);
      }

      d_insDJ = DiscreteDistribution(proba_insDJs.begin(),
				      proba_insDJs.end());

      std::getline(file, line); // # Markov DJ
      std::getline(file, line);
      iss = std::istringstream(line);
      std::array<double, 16> vec_markov_DJ;
      for(std::size_t ii=0; ii < 16; ++ii){
	double p;
	if(not (iss >> p)) throw std::runtime_error("Invalid file");
	vec_markov_DJ[ii] = p;
      }

      std::getline(file, line); // # First nucleotide bias DJ
      std::getline(file, line);
      iss = std::istringstream(line);
      std::array<double, 4> first_nuc_DJ;
      for(std::size_t ii=0; ii < 4; ++ii){
	double p;
	if(not (iss >> p)) throw std::runtime_error("Invalid file");
	first_nuc_DJ[ii] = p;
      }

      markov_DJ = MarkovDNA(first_nuc_DJ, vec_markov_DJ);

      std::getline(file, line); // # Error rate
      std::getline(file, line);
      iss = std::istringstream(line);
      if(not (iss >> error_rate)) throw std::runtime_error("Invalid file");

      std::getline(file, line); // # Thymus selection parameter
      std::getline(file, line);
      iss = std::istringstream(line);
      if(not (iss >> thymic_Q)) throw std::runtime_error("Invalid file");


      // std::getline(file, line); // # Conserved J residues
      // std::getline(file, line);
      // iss = std::istringstream(line);
      // std::copy(std::istream_iterator<std::string>(iss),
      // 		std::istream_iterator<std::string>(),
      // 		std::back_inserter(conserved_J_residues));

      if(verbose){
	std::cout << "# File " << file_gen << " loaded. Summary:\n";
	std::cout << "# Number of V genes: " << nb_V << "\n";
	std::cout << "# Number of D genes: " << nb_D << "\n";
	std::cout << "# Number of J genes: " << nb_J << "\n";
	std::cout << "# Number of delV: " << nb_delV << "\n";
	std::cout << "# Number of delJ: " << nb_delJ << "\n";
	std::cout << "# Number of delDl: " << nb_delDl << "\n";
	std::cout << "# Number of delDr: " << nb_delDr << "\n";
	std::cout << "# Error rate: " << error_rate << "\n";
	std::cout << "# Thymic selection parameter: " << thymic_Q << "\n";
      }
    }
  else{
    std::cout << "Invalid file\n";
  }
}

std::tuple<std::string, std::string, std::size_t, std::size_t> SequenceGenerationVDJ::generate(bool functional){
  while(true){

    std::size_t v_index = d_V.generate(random_generator);
    std::size_t dj_index = d_DJ.generate(random_generator);
    std::size_t d_index = dj_index / nb_J;
    std::size_t j_index = dj_index % nb_J;

    dna seqV = seq_V_CDR3[v_index];
    std::size_t lenV = seqV.size();
    dna seqD = seq_D_CDR3[ d_index ];
    dna seqJ = seq_J_CDR3[ j_index ];

    std::size_t lenD = seqD.size();
    std::size_t lenJ = seqJ.size();
    std::size_t delV = d_delV[v_index].generate(random_generator);

    if (lenV <= delV) // maybe remove
      continue;

    std::size_t delD = d_delD3_delD5[d_index].generate(random_generator);
    std::size_t delD5 = delD / nb_delDr;
    std::size_t delD3 = delD % nb_delDr;
    std::size_t delJ = d_delJ[j_index].generate(random_generator);

    if(lenD < (delD3 + delD5) or lenJ < delJ) //same
      continue;

    std::size_t insVD = d_insVD.generate(random_generator);
    std::size_t insDJ = d_insDJ.generate(random_generator);

    // Check if the sequence is in-frame
    if(functional and (lenV - delV + lenD - delD5 - delD3 + lenJ - delJ + insVD + insDJ)%3 != 0)
      continue;

    // Generate the insertions
    dna ins_seq_VD = markov_VD.generate(insVD, random_generator);
    dna ins_seq_DJ = markov_DJ.generate(insDJ, random_generator);

    // Generate the sequence
    dna seq;
    seq.insert(seq.end(), seqV.begin(), seqV.end() - delV);
    seq.insert(seq.end(), ins_seq_VD.begin(), ins_seq_VD.end());
    seq.insert(seq.end(), seqD.begin() + delD5, seqD.end() - delD3);
    seq.insert(seq.end(), ins_seq_DJ.rbegin(), ins_seq_DJ.rend()); // Need to reverse DJ
    seq.insert(seq.end(), seqJ.begin() + delJ, seqJ.end());

    // translate
    amino seq_aa = translate(seq);

    // Check for stop codon
    if(functional and std::find(seq_aa.begin(), seq_aa.end(), 0) != seq_aa.end())
      continue;

    // removed to match olga/sonia behavior
    // Check for conserved extremities
    // if(functional and seq_aa[0] != 2)
    //   continue;
    // amino_acid last = seq_aa.back();
    // if(functional and last != 5 and last != 18 and last != 19)
    //   continue;

    // apply error rate, on the nucleotide sequence
    add_error(seq, error_rate, random_generator);

    return std::make_tuple(dna_to_string(seq),
			   amino_to_string(translate(seq)),
			   v_index,
			   j_index);
  }
}




SequenceGenerationVJ::SequenceGenerationVJ(const std::string& file_gen, uint64_t seed, bool verbose_arg){
  if(seed)
    random_generator = std::mt19937_64(seed);
  else
    random_generator = std::mt19937_64(std::random_device()());
  verbose = false; // x_verbose;
  load_file(file_gen);
}


void SequenceGenerationVJ::load_file(const std::string& file_gen){
  std::ifstream file(file_gen);
  std::istringstream iss;
  std::string line;
  if (file.is_open()){
      /* Load genes */
      std::getline(file, line); // # V
      std::getline(file, line); // nb of V
      iss = std::istringstream(line);
      if(not (iss >> nb_V)) throw std::runtime_error("Invalid file");
      for(std::size_t ii=0; ii < nb_V; ++ii){
	std::getline(file, line);
	seq_V_CDR3.push_back(string_to_dna(line));
      }

      std::getline(file, line); // # J
      std::getline(file, line); // nb of J
      iss = std::istringstream(line);
      if(not (iss >> nb_J)) throw std::runtime_error("Invalid file");
      for(std::size_t ii=0; ii < nb_J; ++ii){
	std::getline(file, line);
	seq_J_CDR3.push_back(string_to_dna(line));
      }

      /* Marginals */
      std::getline(file, line); // # P(V, J)
      std::getline(file, line); // nbV nbJ
      std::getline(file, line);
      iss = std::istringstream(line);
      std::vector<double> probas_VJ;
      for(std::size_t ii = 0; ii < nb_V * nb_J; ++ii){
	double p;
	if(not (iss >> p)) throw std::runtime_error("Invalid file");
	probas_VJ.push_back(p);
      }
      d_VJ = DiscreteDistribution(probas_VJ.begin(), probas_VJ.end());

      std::getline(file, line); // # P(delV|V)
      std::getline(file, line); // nb_V  nb_delV
      iss = std::istringstream(line);
      if(not (iss >> nb_V >> nb_delV)) throw std::runtime_error("Invalid file");
      for(std::size_t ii=0; ii < nb_V; ++ii){
	std::getline(file, line);
	iss = std::istringstream(line);
	std::vector<double> delV;
	delV.clear();
	for(std::size_t jj = 0; jj < nb_delV; ++jj){
	  double p;
	  if(not (iss >> p)) throw std::runtime_error("Invalid file");
	  delV.push_back(p);
	}
	d_delV.push_back(DiscreteDistribution(delV.begin(), delV.end()));
      }


      std::getline(file, line); // # P(delJ|J)
      std::getline(file, line); // nb_J  nb_delJ
      iss = std::istringstream(line);
      if(not (iss >> nb_J >> nb_delJ)) throw std::runtime_error("Invalid file");
      for(std::size_t ii=0; ii < nb_J; ++ii){
	std::getline(file, line);
	iss = std::istringstream(line);
	std::vector<double> delJ;
	delJ.clear();
	for(std::size_t jj = 0; jj < nb_delJ; ++jj){
	  double p;
	  if(not (iss >> p)) throw std::runtime_error("Invalid file");
	  delJ.push_back(p);
	}
	d_delJ.push_back(DiscreteDistribution(delJ.begin(), delJ.end()));
      }

      std::getline(file, line); // # P(insVJ)
      std::getline(file, line); // nb_insVJ
      iss = std::istringstream(line);
      if(not (iss >> nb_insVJ)) throw std::runtime_error("Invalid file");
      std::getline(file, line);
      iss = std::istringstream(line);
      std::vector<double> proba_insVJs;
      for(std::size_t ii=0; ii < nb_insVJ; ++ii){
	double p;
	if(not (iss >> p)) throw std::runtime_error("Invalid file");
	proba_insVJs.push_back(p);
      }

      d_insVJ = DiscreteDistribution(proba_insVJs.begin(),
				      proba_insVJs.end());


      std::getline(file, line); // # Markov VJ
      std::getline(file, line);
      iss = std::istringstream(line);
      std::array<double, 16> markov_VJ_mat;
      for(std::size_t ii=0; ii < 16; ++ii){
	double p;
	if(not (iss >> p)) throw std::runtime_error("Invalid file");
	markov_VJ_mat[ii] = p;
      }

      std::getline(file, line); // # First nucleotide bias VJ
      std::getline(file, line);
      iss = std::istringstream(line);
      std::array<double, 4> first_nuc_VJ;
      for(std::size_t ii=0; ii < 4; ++ii){
	double p;
	if(not (iss >> p)) throw std::runtime_error("Invalid file");
	first_nuc_VJ[ii] = p;
      }

      markov_VJ = MarkovDNA(first_nuc_VJ, markov_VJ_mat);

      std::getline(file, line); // # Error rate
      std::getline(file, line);
      iss = std::istringstream(line);
      if(not (iss >> error_rate)) throw std::runtime_error("Invalid file");

      std::getline(file, line); // # Thymus selection parameter
      std::getline(file, line);
      iss = std::istringstream(line);
      if(not (iss >> thymic_Q)) throw std::runtime_error("Invalid file");


      // TODO: should reimplement that at some point
      // std::getline(file, line); // # Conserved J residues
      // std::getline(file, line);
      // iss = std::istringstream(line);
      // std::copy(std::istream_iterator<std::string>(iss),
      // 		std::istream_iterator<std::string>(),
      // 		std::back_inserter(conserved_J_residues));

      if(verbose){
	std::cout << "# File " << file_gen << " loaded. Summary:\n";
	std::cout << "# Number of V genes: " << nb_V << "\n";
	std::cout << "# Number of J genes: " << nb_J << "\n";
	std::cout << "# Number of delV: " << nb_delV << "\n";
	std::cout << "# Number of delJ: " << nb_delJ << "\n";
	std::cout << "# Error rate: " << error_rate << "\n";
	std::cout << "# Thymic selection parameter: " << thymic_Q << "\n";
      }
    }
  else{
    std::cout << "Invalid file\n";
  }
}


std::tuple<std::string, std::string, std::size_t, std::size_t> SequenceGenerationVJ::generate(bool functional){
  while(true){

    std::size_t vj_index = d_VJ.generate(random_generator);
    std::size_t v_index = vj_index / nb_J;
    std::size_t j_index = vj_index % nb_J;

    dna seqV = seq_V_CDR3[v_index];
    std::size_t lenV = seqV.size();
    dna seqJ = seq_J_CDR3[ j_index ];

    std::size_t lenJ = seqJ.size();
    std::size_t delV = d_delV[v_index].generate(random_generator);

    if (lenV < delV) // maybe remove
      continue;

    std::size_t delJ = d_delJ[j_index].generate(random_generator);

    if(lenJ < delJ) //same
      continue;

    std::size_t insVJ = d_insVJ.generate(random_generator);

    // Check if the sequence is in-frame
    if(functional and (lenV - delV + lenJ - delJ + insVJ)%3 != 0)
      continue;

    // Generate the insertions
    dna ins_seq_VJ = markov_VJ.generate(insVJ, random_generator);

    // Generate the sequence
    dna seq;
    seq.insert(seq.end(), seqV.begin(), seqV.end() - delV);
    seq.insert(seq.end(), ins_seq_VJ.begin(), ins_seq_VJ.end());
    seq.insert(seq.end(), seqJ.begin() + delJ, seqJ.end());

    // translate
    amino seq_aa = translate(seq);

    // Check for stop codon
    if(functional and std::find(seq_aa.begin(), seq_aa.end(), 0) != seq_aa.end())
      continue;

    // Removed to match olga / sonia behaviour
    // // Check for conserved extremities
    // if(functional and seq_aa[0] != 2)
    //   continue;
    // amino_acid last = seq_aa.back();
    // if(functional and last != 5 and last != 18 and last != 19)
    //   continue;

    // apply error rate, on the nucleotide sequence
    add_error(seq, error_rate, random_generator);

    return std::make_tuple(dna_to_string(seq),
			   amino_to_string(translate(seq)),
			   v_index,
			   j_index);
  }
}
