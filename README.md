# bayesian-genotype-calling-pipeline
1) tissue_parser.py
Output:
The list of samples of the specific tissue
Train sample list
Test sample list
Genotype pcs

2) analyze_egenes_output.py
Output:
Find eqtls which meet selection criteria

3) genotype_finder.py
Output:
Genotypes of train and test set (according to the lists mentioned above)

4) bayesian_genotype_caller.py
Complete the computation procedure
Output:
The accuracy report w.r.t individuals and variants
