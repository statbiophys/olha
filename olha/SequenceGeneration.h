#ifndef SEQUENCE_GENERATION
#define SEQUENCE_GENERATION

#include "utils.h"

class SequenceGenerationVDJ{
    private:
        bool verbose;

        std::mt19937_64 random_generator;

        std::vector<dna> seq_V_CDR3;
        std::vector<dna> seq_J_CDR3;
        std::vector<dna> seq_D_CDR3;

        std::size_t nb_V;
        std::size_t nb_J;
        std::size_t nb_D;
        std::size_t nb_delV;
        std::size_t nb_delDl;
        std::size_t nb_delDr;
        std::size_t nb_delJ;
        std::size_t nb_insDJ;
        std::size_t nb_insVD;
        DiscreteDistribution d_V;
        DiscreteDistribution d_DJ;
        DiscreteDistribution d_insVD;
        DiscreteDistribution d_insDJ;
        std::vector<DiscreteDistribution> d_delV;
        std::vector<DiscreteDistribution> d_delJ;
        std::vector<DiscreteDistribution> d_delD3_delD5;
        MarkovDNA markov_VD;
        MarkovDNA markov_DJ;
        double error_rate;
        double thymic_Q; // useless now, but removing it may be a pain

    public:
        SequenceGenerationVDJ(const std::string& file_gen);
        void load_file(const std::string& file_gen);
        std::tuple<std::string, std::string, std::size_t, std::size_t> generate(bool functional);

};


class SequenceGenerationVJ{
    private:
        bool verbose;

        std::mt19937_64 random_generator;

        std::vector<dna> seq_V_CDR3;
        std::vector<dna> seq_J_CDR3;

        std::size_t nb_V;
        std::size_t nb_J;
        std::size_t nb_delV;
        std::size_t nb_delJ;
        std::size_t nb_insVJ;
        DiscreteDistribution d_V;
        DiscreteDistribution d_VJ;
        DiscreteDistribution d_insVJ;
        std::vector<DiscreteDistribution> d_delV;
        std::vector<DiscreteDistribution> d_delJ;
        MarkovDNA markov_VJ;
        double error_rate;
        double thymic_Q; // useless now, but removing it may be a pain

    public:
        SequenceGenerationVJ(const std::string& file_gen);
        void load_file(const std::string& file_gen);
        std::tuple<std::string, std::string, std::size_t, std::size_t> generate(bool functional);

};


#endif //SEQUENCE_GENERATION
