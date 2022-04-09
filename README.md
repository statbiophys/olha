# olha
A package to generate TCR/BCR sequences fast, based on [olga](https://github.com/statbiophys/OLGA). Use the same syntax as olga but is up to 60x faster and can optionally generate non-functional sequences and include point-mutation "sequencing" errors.

Written in C++, interface with python3 with pybind11.

## Installation

```sh
pip install olha
```


## Example

```{py}
import olga
import olga.sequence_generation
import olga.load_model
import olha


## olga model loading
params_file_name = f'{olga.__path__[0]}/default_models/human_T_beta/model_params.txt'
marginals_file_name = f'{olga.__path__[0]}/default_models/human_T_beta/model_marginals.txt'
V_anchor_pos_file =f'{olga.__path__[0]}/default_models/human_T_beta/V_gene_CDR3_anchors.csv'
J_anchor_pos_file = f'{olga.__path__[0]}/default_models/human_T_beta/J_gene_CDR3_anchors.csv'

genomic_data = olga.load_model.GenomicDataVDJ()
genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
generative_model = olga.load_model.GenerativeModelVDJ()
generative_model.load_and_process_igor_model(marginals_file_name)

# sequence generation
olha_gen = olha.SequenceGeneration(genomic_data, generative_model, error_rate=0.1)
olga_gen.gen_rnd_prod_CDR3()
# ('TGCGCCAGCAGCTCCATGGACGGCTCCGAAAAACTGTTTTTT', 'CASSSMDGSEKLFF', 49, 3)

olga_gen = olga.sequence_generation.SequenceGenerationVDJ(generative_model, genomic_data)
```

## Comparison
```
import timeit
olha_gen = olha.SequenceGeneration(genomic_data, generative_model, error_rate=0.1)
olga_gen = olga.sequence_generation.SequenceGenerationVDJ(generative_model, genomic_data)

timeit.timeit(olha_gen.gen_rnd_prod_CDR3) # 3.31 μs
timeit.timeit(olga_gen.gen_rnd_prod_CDR3) # 103 μs
```
